"""Incremental Multi-Extension FITS reader for HTTP streams.

Why: a SPHEREx L2 cutout trims IMAGE/FLAGS/VARIANCE/ZODI to the requested size,
but the PSF extension (a cube of oversampled PSFs over the focal plane) and the
WCS-WAVE lookup table are copied whole. For a 0.1 deg cutout the pixel data is
~60 KB while the file is ~5-6 MB. FITS is strictly sequential
(header blocks + data blocks, each padded to 2880 bytes), so we can parse HDUs
on the fly and *close the connection* as soon as every wanted EXTNAME has been
received - typically 50-100x less traffic.
"""
from __future__ import annotations

from collections.abc import AsyncIterator, Iterable

BLOCK = 2880
CARD = 80


class FitsStreamError(RuntimeError):
    pass


class NotFitsError(FitsStreamError):
    """Upstream answered with something that is not FITS.

    IRSA under load replies with a short text body ("Server busy, retry later.")
    - this must be treated as a retryable error, never cached.
    """

    def __init__(self, head: bytes):
        self.head = head
        super().__init__(f"response is not FITS: {head[:80]!r}")


def _parse_header(raw: bytes) -> dict[str, str]:
    out: dict[str, str] = {}
    for i in range(0, len(raw), CARD):
        card = raw[i : i + CARD].decode("ascii", errors="replace")
        key = card[:8].strip()
        if key == "END":
            break
        if card[8:10] == "= ":
            val = card[10:].split("/", 1)[0].strip() if not card[10:].lstrip().startswith("'") else card[10:].strip()
            if val.startswith("'"):
                val = val[1:].split("'", 1)[0].strip()
            out[key] = val
    return out


def _data_size(h: dict[str, str]) -> int:
    naxis = int(h.get("NAXIS", 0))
    if naxis == 0:
        return 0
    bitpix = abs(int(h["BITPIX"]))
    gcount = int(h.get("GCOUNT", 1))
    pcount = int(h.get("PCOUNT", 0))
    n = 1
    for i in range(1, naxis + 1):
        n *= int(h[f"NAXIS{i}"])
    nbytes = bitpix // 8 * gcount * (pcount + n)
    return (nbytes + BLOCK - 1) // BLOCK * BLOCK


async def read_hdus(
    chunks: AsyncIterator[bytes],
    wanted: Iterable[str] | None = None,
    keep_unwanted: bool = False,
) -> tuple[bytes, list[str], bool]:
    """Consume a FITS byte stream until all `wanted` EXTNAMEs were read.

    Returns (fits_bytes, extnames_in_output, truncated). `fits_bytes` is a valid
    FITS file containing the primary HDU plus the wanted extensions (unwanted
    extensions are skipped, i.e. not copied, unless keep_unwanted).
    `truncated` is True when we stopped before the end of the remote file -
    the caller must then close the HTTP response.
    """
    want = {w.upper() for w in wanted} if wanted else None
    buf = bytearray()
    out = bytearray()
    names: list[str] = []
    it = chunks.__aiter__()
    eof = False

    async def need(n: int) -> bool:
        nonlocal eof
        while len(buf) < n and not eof:
            try:
                buf.extend(await it.__anext__())
            except StopAsyncIteration:
                eof = True
        return len(buf) >= n

    hdu_index = 0
    while True:
        # ---- header: read 2880-byte blocks until the END card ----
        header = bytearray()
        while True:
            if hdu_index == 0 and not header:
                await need(BLOCK)
                if not bytes(buf[:9]) == b"SIMPLE  =":
                    raise NotFitsError(bytes(buf[:200]))
            if not await need(BLOCK):
                if not header and hdu_index > 0:
                    return bytes(out), names, False  # clean EOF between HDUs
                raise FitsStreamError("stream ended inside a FITS header")
            block = bytes(buf[:BLOCK])
            del buf[:BLOCK]
            header.extend(block)
            if any(block[i : i + 8] == b"END     " for i in range(0, BLOCK, CARD)):
                break
        h = _parse_header(bytes(header))
        name = "PRIMARY" if hdu_index == 0 else h.get("EXTNAME", f"HDU{hdu_index}").upper()
        size = _data_size(h)
        keep = hdu_index == 0 or want is None or name in want or keep_unwanted
        if keep:
            out.extend(header)
        # ---- data ----
        remaining = size
        while remaining:
            if not await need(1):
                raise FitsStreamError(f"stream ended inside data of HDU {name}")
            take = min(remaining, len(buf))
            if keep:
                out.extend(buf[:take])
            del buf[:take]
            remaining -= take
        if keep:
            names.append(name)
        hdu_index += 1
        if want is not None and want.issubset(set(names)):
            done_early = not (eof and not buf)
            return bytes(out), names, done_early
