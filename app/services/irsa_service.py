"""NASA IRSA access: SIAv2/TAP (ObsCore) -> epochs -> Datalink -> streamed cutouts.

Facts verified against the live archive (DR1, Nov 2026):
* `spherex.obscore` holds one row per (exposure, detector). `obs_id` is shared
  by the 6 detectors of an exposure -> frame identity is (obs_id, detector).
* Adding `ORDER BY t_min ASC` to ADQL on `spherex.obscore` makes IRSA's query
  planner occasionally stall for >60s. Omitting `ORDER BY` (sorting in Python)
  and racing IRSA SIAv2 (`/SIA`) with `/TAP/sync?FORMAT=csv` returns in ~2-5s.
* Cutouts: `{fits_url}?center={ra},{dec}&size={deg}`. IMAGE/FLAGS/VARIANCE/ZODI
  are trimmed (~158 KB); EPSF (~3.9 MB) is copied whole -> we stream and
  always stop immediately after ZODI (`fits_stream.read_hdus`).
* IRSA answers HTTP 503 "Server busy, retry later." under load -> bounded
  concurrency + exponential backoff with jitter.
"""
from __future__ import annotations

import asyncio
import csv
import io
import json
import logging
import math
import random
import re
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import httpx
import numpy as np

from app.config import Settings
from app.core.async_utils import race_first_success
from app.core.keys import normalize_dec, normalize_ra
from app.core.spherex import detector_from_bandname
from app.services.cache_manager import FitsCache
from app.services.fits_stream import FitsStreamError, NotFitsError, read_hdus

log = logging.getLogger(__name__)

CUTOUT_HDUS = ("IMAGE", "FLAGS", "VARIANCE", "ZODI")
RETRYABLE_STATUS = {429, 500, 502, 503, 504}
IRSA_SIA_URL = "https://irsa.ipac.caltech.edu/SIA"


@dataclass
class Frame:
    obs_id: str
    detector: int
    collection: str
    t_min: float  # MJD
    t_max: float
    access_url: str
    s_ra: float | None = None
    s_dec: float | None = None
    fits_url: str | None = None

    @property
    def key(self) -> str:
        return f"{self.obs_id}/D{self.detector}"

    @property
    def t_mid(self) -> float:
        return 0.5 * (self.t_min + self.t_max)


@dataclass
class Epoch:
    index: int
    frames: list[Frame] = field(default_factory=list)

    @property
    def mjd_start(self) -> float:
        return min(f.t_min for f in self.frames)

    @property
    def mjd_end(self) -> float:
        return max(f.t_max for f in self.frames)

    @property
    def mjd_mid(self) -> float:
        return float(np.median([f.t_mid for f in self.frames]))


def collection_priority(name: str) -> tuple[int, int]:
    """spherex_dr2 > spherex_dr1 > spherex_qr3 > spherex_qr2 ..."""
    m = re.match(r"spherex_(dr|qr)(\d+)", name or "")
    if not m:
        return (0, 0)
    return (2 if m.group(1) == "dr" else 1, int(m.group(2)))


class UpstreamError(RuntimeError):
    pass


class IrsaService:
    def __init__(self, settings: Settings, http: httpx.AsyncClient, cache: FitsCache):
        self.s = settings
        self.http = http
        self.cache = cache
        self._sem = asyncio.Semaphore(settings.download_concurrency)
        self._tap_cache_dir = cache.root / "tap"
        self._tap_cache_dir.mkdir(parents=True, exist_ok=True)

    # ======================================================= TAP / SIAv2 ====
    def _tap_cache_path(self, ra: float, dec: float) -> Path:
        return self._tap_cache_dir / f"{normalize_ra(ra):.4f}_{normalize_dec(dec):+.4f}.json"

    async def _query_tap_csv(self, ra: float, dec: float) -> list[Frame]:
        """Fast async HTTP query to IRSA TAP/sync in CSV format (no ORDER BY!)."""
        adql = (
            "SELECT obs_id, obs_collection, energy_bandpassname, t_min, t_max, s_ra, s_dec, access_url "
            f"FROM {self.s.spherex_obscore_table} "
            f"WHERE CONTAINS(POINT('ICRS', {ra:.8f}, {dec:.8f}), s_region) = 1"
        )
        url = f"{self.s.irsa_tap_url.rstrip('/')}/sync"
        r = await self.http.post(
            url,
            data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql},
            timeout=25.0,
        )
        r.raise_for_status()
        text = r.text
        if text.lstrip().startswith("<"):
            raise UpstreamError(f"TAP returned XML error: {text[:200]}")
        frames: list[Frame] = []
        for row in csv.DictReader(io.StringIO(text)):
            det = detector_from_bandname(row.get("energy_bandpassname"))
            if det is None:
                continue
            frames.append(
                Frame(
                    obs_id=row["obs_id"].strip(),
                    detector=det,
                    collection=row["obs_collection"].strip(),
                    t_min=float(row["t_min"]),
                    t_max=float(row["t_max"]),
                    access_url=row["access_url"].strip(),
                    s_ra=float(row["s_ra"]) if row.get("s_ra") else None,
                    s_dec=float(row["s_dec"]) if row.get("s_dec") else None,
                )
            )
        return frames

    @staticmethod
    def _parse_sia_votable(content: bytes) -> list[Frame]:
        from astropy.io import votable

        tab = votable.parse_single_table(io.BytesIO(content)).to_table(use_names_over_ids=True)
        frames: list[Frame] = []
        for row in tab:
            det = detector_from_bandname(str(row["energy_bandpassname"]))
            if det is None:
                continue
            frames.append(
                Frame(
                    obs_id=str(row["obs_id"]),
                    detector=det,
                    collection=str(row["obs_collection"]),
                    t_min=float(row["t_min"]),
                    t_max=float(row["t_max"]),
                    access_url=str(row["access_url"]),
                    s_ra=float(row["s_ra"]) if not np.ma.is_masked(row["s_ra"]) else None,
                    s_dec=float(row["s_dec"]) if not np.ma.is_masked(row["s_dec"]) else None,
                )
            )
        return frames

    async def _query_sia(self, ra: float, dec: float) -> list[Frame]:
        """IRSA SIAv2 spatial cone search (often faster than TAP when TAP is busy)."""
        r = await self.http.get(
            IRSA_SIA_URL,
            params={
                "COLLECTION": "spherex_dr1",
                "POS": f"CIRCLE {ra:.8f} {dec:.8f} 0.0005",
                "RESPONSEFORMAT": "VOTABLE",
            },
            timeout=25.0,
        )
        r.raise_for_status()
        return await asyncio.to_thread(self._parse_sia_votable, r.content)

    async def query_frames(self, ra: float, dec: float, force_refresh: bool = False) -> list[Frame]:
        cpath = self._tap_cache_path(ra, dec)
        if not force_refresh and self.cache.get(cpath):
            try:
                raw = json.loads(cpath.read_text())
                log.info("IRSA ObsCore disk cache hit for (%.4f, %+.4f): %d frames", ra, dec, len(raw))
                return [Frame(**item) for item in raw]
            except Exception:
                pass

        t0 = time.perf_counter()
        last_err: Exception | None = None
        for attempt in range(self.s.http_retries):
            # Race TAP CSV and SIAv2 until the FIRST SUCCESS (sibling failures don't cancel winner)
            tasks = [
                asyncio.create_task(self._query_tap_csv(ra, dec), name="irsa-tap-csv"),
                asyncio.create_task(self._query_sia(ra, dec), name="irsa-sia2"),
            ]
            try:
                raw_frames, winner = await race_first_success(tasks)
                frames = self.deduplicate(raw_frames)
                log.info(
                    "IRSA (%s) returned %d unique frames in %.2fs",
                    winner,
                    len(frames),
                    time.perf_counter() - t0,
                )
                self.cache.write_atomic(
                    cpath,
                    json.dumps([asdict(f) for f in frames]).encode(),
                )
                return frames
            except asyncio.CancelledError:
                raise
            except Exception as e:
                last_err = e
            log.warning("IRSA query attempt %d failed (%r), retrying...", attempt + 1, last_err)
            await asyncio.sleep(1.0 * 2**attempt + random.random() * 0.5)

        # Stale-if-Error fallback: serve previously cached ObsCore metadata even if force_refresh=True
        if cpath.is_file():
            try:
                raw = json.loads(cpath.read_text())
                log.warning(
                    "IRSA upstream failed (%r); serving stale disk cache (%d frames) for (%.4f, %+.4f)",
                    last_err,
                    len(raw),
                    ra,
                    dec,
                )
                self.cache.touch(cpath)
                return [Frame(**item) for item in raw]
            except Exception as cache_err:
                log.debug("Failed reading stale TAP cache %s: %r", cpath, cache_err)

        raise UpstreamError(f"IRSA query failed: {last_err!r}")

    @staticmethod
    def deduplicate(frames: list[Frame]) -> list[Frame]:
        best: dict[str, Frame] = {}
        for f in frames:
            cur = best.get(f.key)
            if cur is None or collection_priority(f.collection) > collection_priority(cur.collection):
                best[f.key] = f
        return sorted(best.values(), key=lambda f: f.t_min)

    # ============================================================= epochs ===
    def group_epochs(self, frames: list[Frame]) -> list[Epoch]:
        epochs: list[Epoch] = []
        cur: list[Frame] = []
        for f in sorted(frames, key=lambda x: x.t_min):
            if cur and (
                f.t_min - cur[-1].t_min > self.s.epoch_gap_days
                or f.t_min - cur[0].t_min > self.s.epoch_max_span_days
            ):
                epochs.append(Epoch(len(epochs), cur))
                cur = []
            cur.append(f)
        if cur:
            epochs.append(Epoch(len(epochs), cur))
        return epochs[-self.s.max_epochs :] if len(epochs) > self.s.max_epochs else epochs

    @staticmethod
    def _dist2(f: Frame, ra: float, dec: float) -> float:
        if f.s_ra is None or f.s_dec is None:
            return 1e9
        dra = (f.s_ra - ra + 180) % 360 - 180
        return (dra * math.cos(math.radians(dec))) ** 2 + (f.s_dec - dec) ** 2

    def select_display_frames(
        self,
        epoch: Epoch,
        detectors: tuple[int, ...],
        ra: float,
        dec: float,
        max_per_epoch: int | None = None,
    ) -> list[Frame]:
        """Select up to `max_per_epoch` (default 1-2) best-centered exposures per epoch
        closest to `(ra, dec)` so 4 epochs stream at most 4-8 FITS cutouts concurrently."""
        limit = max_per_epoch if max_per_epoch is not None else self.s.max_frames_per_epoch
        per_det: list[Frame] = []
        for d in detectors:
            cands = sorted((f for f in epoch.frames if f.detector == d), key=lambda f: self._dist2(f, ra, dec))
            per_det.extend(cands[: self.s.max_frames_per_detector_per_epoch])
        per_det.sort(key=lambda f: self._dist2(f, ra, dec))
        return per_det[: max(1, limit)]

    def select_spectrum_frames(self, epochs: list[Epoch], already: set[str]) -> list[Frame]:
        if self.s.spectrum_max_frames <= 0:
            return []
        pools = []
        for ep in epochs:
            for d in range(1, 7):
                pools.append([f for f in ep.frames if f.detector == d and f.key not in already])
        out: list[Frame] = []
        while len(out) < self.s.spectrum_max_frames and any(pools):
            for p in pools:
                if p and len(out) < self.s.spectrum_max_frames:
                    out.append(p.pop(len(p) // 2))
        return out

    # =========================================================== datalink ===
    async def _get(self, url: str, **kw) -> httpx.Response:
        for attempt in range(self.s.http_retries + 1):
            try:
                r = await self.http.get(url, **kw)
                if r.status_code in RETRYABLE_STATUS:
                    raise httpx.HTTPStatusError("retryable", request=r.request, response=r)
                r.raise_for_status()
                return r
            except (httpx.TransportError, httpx.HTTPStatusError) as e:
                if isinstance(e, httpx.HTTPStatusError) and e.response.status_code not in RETRYABLE_STATUS:
                    raise
                if attempt == self.s.http_retries:
                    raise
                await asyncio.sleep(0.5 * 2**attempt + random.random() * 0.5)
        raise AssertionError("unreachable")

    @staticmethod
    def _parse_datalink_votable(content: bytes) -> str | None:
        # Fast regex path for IRSA Datalink XML (<TD>https://...fits</TD>...<TD>#this</TD>)
        for m in re.finditer(rb"<TR[^>]*>(.*?)</TR>", content, re.DOTALL | re.IGNORECASE):
            tr = m.group(1)
            if b"#this" in tr:
                u = re.search(rb"<TD[^>]*>\s*(https?://[^<\s]+)\s*</TD>", tr, re.IGNORECASE)
                if u:
                    return u.group(1).decode("utf-8", errors="ignore").strip()

        from astropy.io import votable

        tab = votable.parse_single_table(io.BytesIO(content)).to_table(use_names_over_ids=True)
        for row in tab:
            if str(row["semantics"]) == "#this" and str(row["access_url"]).strip():
                return str(row["access_url"]).strip()
        return None

    async def resolve_one(self, frame: Frame) -> str | None:
        if frame.fits_url:
            return frame.fits_url
        if "datalink" not in frame.access_url:
            frame.fits_url = frame.access_url
            return frame.fits_url

        async with self._sem:
            try:
                r = await self._get(frame.access_url)
            except Exception as e:
                log.warning("datalink failed for %s: %s", frame.key, e)
                return None
        resolved = await asyncio.to_thread(self._parse_datalink_votable, r.content)
        if resolved:
            frame.fits_url = resolved
            return frame.fits_url
        return None

    async def resolve_fits_urls(self, frames: list[Frame], known: dict[str, str]) -> list[Frame]:
        for f in frames:
            f.fits_url = f.fits_url or known.get(f.key)
        missing = [f for f in frames if not f.fits_url]
        if missing:
            t0 = time.perf_counter()
            await asyncio.gather(*(self.resolve_one(f) for f in missing))
            log.info("Resolved %d Datalink URLs in %.2fs", len(missing), time.perf_counter() - t0)
        return [f for f in frames if f.fits_url]

    # ============================================================ cutouts ===
    def wavecal_key(self, frame: Frame) -> str:
        return f"{frame.collection}_D{frame.detector}"

    async def fetch_cutout(self, frame: Frame, ra: float, dec: float, size_deg: float) -> Path | None:
        """Stream a cutout and stop right after ZODI (~158 KB instead of 4 MB)."""
        assert frame.fits_url
        url = f"{frame.fits_url}?center={ra:.7f},{dec:.7f}&size={size_deg:.5f}"
        path = self.cache.raw_path(url)
        if self.cache.get(path):
            return path

        data = await self._stream_hdus(url, CUTOUT_HDUS)
        if data is None:
            return None
        self.cache.write_atomic(path, data)
        return path

    async def _stream_hdus(self, url: str, wanted: tuple[str, ...]) -> bytes | None:
        max_attempts = min(self.s.http_retries + 1, 2)
        for attempt in range(max_attempts):
            async with self._sem:
                try:
                    async with self.http.stream("GET", url, timeout=12.0) as r:
                        if r.status_code == 404:
                            return None
                        if r.status_code in RETRYABLE_STATUS:
                            raise UpstreamError(f"HTTP {r.status_code}")
                        r.raise_for_status()
                        data, names, _ = await read_hdus(r.aiter_bytes(), wanted)
                    if not set(CUTOUT_HDUS).issubset(names):
                        raise FitsStreamError(f"missing HDUs, got {names}")
                    return data
                except (httpx.TransportError, UpstreamError, NotFitsError, FitsStreamError) as e:
                    if attempt == max_attempts - 1:
                        log.warning("cutout skipped after %d attempts: %s (%s)", attempt + 1, url, e)
                        return None
            await asyncio.sleep(0.35 + random.random() * 0.35)
        return None

    @staticmethod
    def _split_wavecal(data: bytes, wpath: Path) -> bytes:
        from astropy.io import fits

        with fits.open(io.BytesIO(data)) as hdul:
            if "WCS-WAVE" in hdul:
                FitsCache.write_atomic(wpath, _hdul_bytes(fits.HDUList([fits.PrimaryHDU(), hdul["WCS-WAVE"].copy()])))
            keep = fits.HDUList([h.copy() for h in hdul if h.name in ("PRIMARY",) + CUTOUT_HDUS])
            return _hdul_bytes(keep)

    async def fetch_many(self, frames: list[Frame], ra: float, dec: float, size_deg: float) -> list[tuple[Frame, Path]]:
        t0 = time.perf_counter()
        paths = await asyncio.gather(*(self.fetch_cutout(f, ra, dec, size_deg) for f in frames))
        ok = [(f, p) for f, p in zip(frames, paths) if p is not None]
        log.info("Cutouts ready: %d/%d in %.2fs", len(ok), len(frames), time.perf_counter() - t0)
        return ok


def _hdul_bytes(hdul) -> bytes:
    buf = io.BytesIO()
    hdul.writeto(buf)
    return buf.getvalue()
