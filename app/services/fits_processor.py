"""FITS processing core: load -> reproject -> co-add -> stretch -> PNG; spectra.

All functions are synchronous & CPU-bound; the pipeline calls them through
`asyncio.to_thread` / a ThreadPoolExecutor (numpy, reproject's C core and
Pillow's encoder release the GIL for the heavy parts).

Design decisions
----------------
* ONE common output grid per target: a north-up TAN projection centred on the
  requested (ra, dec) at the native SPHEREx scale (6.15"/px). Every epoch is
  reprojected onto it, so static stars land on identical pixels and the
  frontend can "blink" PNGs without any client-side registration. Using the
  first epoch as reference would inherit its rotation and SIP distortion.
* Reproject at NATIVE resolution (e.g. 59x59 px for 0.1 deg), upscale only when
  encoding the PNG. Reprojecting straight to 512x512 costs ~75x more and adds
  zero information.
* Bitmasks are reprojected with nearest-neighbour; bilinear interpolation of
  FLAGS would invent bit combinations that never existed.
* Matplotlib is used only for its colormap LUTs; images are encoded with
  Pillow directly (no Figure => no axes/padding, ~10x faster, thread-safe).
"""
from __future__ import annotations

import asyncio
import base64
import copy
import io
import logging
import math
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning
from astropy.visualization import AsinhStretch, LinearStretch, LogStretch, SqrtStretch, ZScaleInterval
from astropy.wcs import WCS, FITSFixedWarning
from PIL import Image

from app.core.spherex import (
    CHANNEL_CENTERS_UM,
    CHANNEL_EDGES_UM,
    CHANNELS_PER_DETECTOR,
    DETECTORS,
    INFO_FLAG_NAMES,
    N_CHANNELS,
    bad_flag_mask,
    channel_index,
    parse_flag_bits,
    source_flag_mask,
)
from app.schemas.target import DataRange, DiffAnomaly, SpectrumPoint

log = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=FITSFixedWarning)

_RESAMPLE = {"nearest": Image.NEAREST, "bilinear": Image.BILINEAR, "bicubic": Image.BICUBIC, "lanczos": Image.LANCZOS}


# ===================================================================== grid ==
@dataclass
class DisplayGrid:
    ra: float
    dec: float
    fov: float
    n: int
    pixscale_deg: float
    png_size: int
    wcs: WCS = field(repr=False)

    def wcs_copy(self) -> WCS:
        """astropy WCS objects are not safe to share between threads (wcsset
        mutates internal state) - give each worker its own copy."""
        return copy.deepcopy(self.wcs)

    def world_to_png(self, ra: float, dec: float) -> tuple[float, float] | None:
        x, y = self.wcs_copy().world_to_pixel_values(ra, dec)
        if not (np.isfinite(x) and np.isfinite(y)):
            return None
        k = self.png_size / self.n
        px, py = (x + 0.5) * k, (self.n - 0.5 - y) * k  # FITS rows go up, PNG rows go down
        if 0 <= px <= self.png_size and 0 <= py <= self.png_size:
            return float(px), float(py)
        return None

    def png_wcs(self) -> dict[str, float | str]:
        """TAN WCS of the PNG itself (1-based pixels, row 1 = TOP)."""
        scale = self.pixscale_deg * self.n / self.png_size
        return {
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRVAL1": self.ra,
            "CRVAL2": self.dec,
            "CRPIX1": (self.png_size + 1) / 2,
            "CRPIX2": (self.png_size + 1) / 2,
            "CDELT1": -scale,
            "CDELT2": -scale,  # negative: PNG y grows towards SOUTH
            "NAXIS1": self.png_size,
            "NAXIS2": self.png_size,
        }


def make_grid(ra: float, dec: float, fov: float, pixscale_arcsec: float, png_size: int) -> DisplayGrid:
    pix = pixscale_arcsec / 3600.0
    n = int(min(max(math.ceil(fov / pix), 16), 600))
    n += (n + 1) % 2  # odd => target sits exactly on the central pixel
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-pix, pix]  # east-left, north-up
    w.wcs.radesys = "ICRS"
    w.array_shape = (n, n)
    return DisplayGrid(ra, dec, fov, n, pix, png_size, w)


# =================================================================== frames ==
@dataclass
class FrameData:
    key: str
    detector: int
    mjd: float
    image: np.ndarray  # IMAGE - ZODI, MJy/sr, bad pixels = NaN
    variance: np.ndarray
    flags: np.ndarray
    wcs: WCS = field(repr=False)
    wave_wcs: WCS | None = field(default=None, repr=False)
    crpix2w: float | None = None
    bits: dict[str, int] = field(default_factory=dict)
    pix_area_sr: float = 8.9e-10

    @property
    def bad_mask(self) -> int:
        return bad_flag_mask(self.bits)


def load_frame(path: Path, key: str, detector: int, mjd: float, wavecal_path: Path | None = None) -> FrameData:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        with fits.open(path, memmap=False) as hdul:
            hdr = hdul["IMAGE"].header
            img = hdul["IMAGE"].data.astype(np.float32)
            zodi = hdul["ZODI"].data.astype(np.float32) if "ZODI" in hdul else 0.0
            var = hdul["VARIANCE"].data.astype(np.float32)
            flg = hdul["FLAGS"].data.astype(np.int64)
            bits = parse_flag_bits(hdul["FLAGS"].header)
            wcs = WCS(hdr).celestial
            crpix2w = float(hdr["CRPIX2W"]) if "CRPIX2W" in hdr else None
            wave = None
            if "WCS-WAVE" not in hdul and wavecal_path is not None and wavecal_path.exists():
                with fits.open(wavecal_path, memmap=False) as wc:
                    hdul.append(wc["WCS-WAVE"].copy())
            if "WCS-WAVE" in hdul:
                try:
                    whdr = hdr.copy()
                    for k in list(whdr.keys()):
                        if k.startswith(("A_", "B_", "AP_", "BP_")):
                            del whdr[k]
                    wave = WCS(whdr, fobj=hdul, key="W")
                except Exception as e:
                    log.debug("spectral WCS unavailable for %s: %s", key, e)
    # Zodiacal light differs between epochs (solar elongation) -> must go,
    # otherwise the whole field "breathes" when blinking.
    img = img - zodi
    bad = (flg & bad_flag_mask(bits)) != 0
    img[bad | ~np.isfinite(img)] = np.nan
    var[bad | ~np.isfinite(var) | (var <= 0)] = np.nan
    try:
        area = float(wcs.proj_plane_pixel_area() * (math.pi / 180) ** 2)  # deg^2 -> sr
    except Exception:
        area = 8.9e-10
    return FrameData(key, detector, mjd, img, var, flg, wcs, wave, crpix2w, bits, area)


# ================================================================== co-adds ==
@dataclass
class EpochCoadd:
    index: int
    image: np.ndarray
    variance: np.ndarray
    flags: np.ndarray  # bitwise OR of contributing frames
    coverage: np.ndarray  # number of frames contributing per pixel
    n_frames: int
    detectors: list[int]
    mjd_start: float
    mjd_end: float
    mjd_mid: float
    rgb_bands: tuple[np.ndarray | None, np.ndarray | None, np.ndarray | None] = field(
        default=(None, None, None), repr=False
    )  # (blue=D1+D2, green=D3+D4, red=D5+D6)

    @property
    def date_iso(self) -> str:
        return Time(self.mjd_mid, format="mjd").isot[:19] + "Z"

    @property
    def label(self) -> str:
        return Time(self.mjd_mid, format="mjd").isot[:7]


def _robust_background(img: np.ndarray) -> float:
    v = img[np.isfinite(img)]
    if v.size < 10:
        return 0.0
    med = np.median(v)
    mad = np.median(np.abs(v - med)) * 1.4826 or 1.0
    return float(np.median(v[np.abs(v - med) < 3 * mad]))


def reproject_frame(fd: FrameData, grid: DisplayGrid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from reproject import reproject_interp

    out_wcs = grid.wcs_copy()
    shape = (grid.n, grid.n)
    img = fd.image - _robust_background(fd.image)  # equalise residual sky offsets between frames/detectors
    stack = np.stack([img, fd.variance])  # one call, broadcast over the leading axis
    (img_r, var_r), _ = reproject_interp((stack, fd.wcs), out_wcs, shape_out=shape, order="bilinear")
    flg_r, fp = reproject_interp((fd.flags.astype(np.float64), fd.wcs), out_wcs, shape_out=shape, order="nearest-neighbor")
    flg_r = np.where(np.isfinite(flg_r) & (fp > 0), flg_r, 0).astype(np.int64)
    return img_r.astype(np.float32), var_r.astype(np.float32), flg_r


def coadd_epoch(index: int, frames: list[FrameData], grid: DisplayGrid) -> EpochCoadd:
    """Inverse-variance weighted mean of all frames of one epoch, plus 3 sub-band
    co-adds (Blue=D1+D2, Green=D3+D4, Red=D5+D6) for WiseView-style RGB rendering."""
    n = grid.n
    num = np.zeros((n, n), np.float64)
    wsum = np.zeros((n, n), np.float64)
    flags = np.zeros((n, n), np.int64)
    cov = np.zeros((n, n), np.int16)

    # Sub-band accumulators: 0=Blue(D1,D2), 1=Green(D3,D4), 2=Red(D5,D6)
    b_num = [np.zeros((n, n), np.float64) for _ in range(3)]
    b_wsum = [np.zeros((n, n), np.float64) for _ in range(3)]

    for fd in frames:
        img, var, flg = reproject_frame(fd, grid)
        ok = np.isfinite(img) & np.isfinite(var) & (var > 0)
        w = np.where(ok, 1.0 / np.where(ok, var, 1.0), 0.0)
        weighted = np.where(ok, img, 0.0) * w
        num += weighted
        wsum += w
        flags |= flg
        cov += ok

        grp = 0 if fd.detector in (1, 2) else (1 if fd.detector in (3, 4) else 2)
        b_num[grp] += weighted
        b_wsum[grp] += w

    with np.errstate(invalid="ignore", divide="ignore"):
        image = np.where(wsum > 0, num / wsum, np.nan).astype(np.float32)
        variance = np.where(wsum > 0, 1.0 / wsum, np.nan).astype(np.float32)
        rgb_list: list[np.ndarray | None] = []
        for g in range(3):
            if np.any(b_wsum[g] > 0):
                rgb_list.append(np.where(b_wsum[g] > 0, b_num[g] / b_wsum[g], np.nan).astype(np.float32))
            else:
                rgb_list.append(None)

    mjds = [f.mjd for f in frames]
    return EpochCoadd(
        index=index,
        image=image,
        variance=variance,
        flags=flags,
        coverage=cov,
        n_frames=len(frames),
        detectors=sorted({f.detector for f in frames}),
        mjd_start=min(mjds),
        mjd_end=max(mjds),
        mjd_mid=float(np.median(mjds)),
        rgb_bands=(rgb_list[0], rgb_list[1], rgb_list[2]),
    )


def deep_stack(coadds: list[EpochCoadd]) -> tuple[np.ndarray, np.ndarray]:
    num = sum(np.nan_to_num(c.image / c.variance) for c in coadds)
    wsum = sum(np.nan_to_num(1.0 / c.variance) for c in coadds)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(wsum > 0, num / wsum, np.nan), np.where(wsum > 0, 1.0 / wsum, np.nan)


# ================================================================ stretching ==
def compute_cuts(images: list[np.ndarray]) -> tuple[float, float, DataRange]:
    """ONE set of cuts for ALL epochs - per-epoch normalisation makes static
    stars flicker in the blink and can hide real variability (RW Aur, R Dra)."""
    vals = np.concatenate([im[np.isfinite(im)].ravel() for im in images]) if images else np.array([])
    if vals.size < 4:
        return 0.0, 1.0, DataRange(p0_5=0, p50=0, p99_5=1, min=0, max=1)
    vmin, vmax = ZScaleInterval(n_samples=2000, contrast=0.15).get_limits(vals)
    p = np.percentile(vals, [0.5, 50, 99.5])
    if not vmax > vmin:
        vmin, vmax = float(p[0]), float(p[2]) if p[2] > p[0] else float(p[0]) + 1.0
    rng = DataRange(p0_5=float(p[0]), p50=float(p[1]), p99_5=float(p[2]), min=float(vals.min()), max=float(vals.max()))
    return float(vmin), float(vmax), rng


def _stretch_fn(stretch_type: str, asinh_a: float = 0.1, log_a: float = 1000.0):
    return {
        "linear": LinearStretch(),
        "zscale": LinearStretch(),  # ZScale is the interval; the curve is linear
        "sqrt": SqrtStretch(),
        "log": LogStretch(a=log_a),
        "asinh": AsinhStretch(a=asinh_a),
    }[stretch_type]


def _to_png(arr_u8: np.ndarray, size: int, resample: str, mode: str) -> bytes:
    im = Image.fromarray(np.ascontiguousarray(np.flipud(arr_u8)), mode=mode)
    if size and im.size != (size, size):
        im = im.resize((size, size), _RESAMPLE.get(resample, Image.NEAREST))
    buf = io.BytesIO()
    im.save(buf, format="PNG", compress_level=1)
    return buf.getvalue()


def render_image(
    data: np.ndarray, vmin: float, vmax: float, stretch_type: str = "asinh", size: int = 512,
    resample: str = "nearest", asinh_a: float = 0.1, log_a: float = 1000.0,
) -> bytes:
    x = (data.astype(np.float64) - vmin) / (vmax - vmin if vmax > vmin else 1.0)
    x = np.clip(np.nan_to_num(x, nan=0.0), 0.0, 1.0)
    y = _stretch_fn(stretch_type, asinh_a, log_a)(x, clip=True)
    return _to_png((y * 255 + 0.5).astype(np.uint8), size, resample, "L")


_LUT_CACHE: dict[str, np.ndarray] = {}


def _lut(name: str) -> np.ndarray:
    if name not in _LUT_CACHE:
        from matplotlib import colormaps

        _LUT_CACHE[name] = (colormaps[name](np.linspace(0, 1, 256))[:, :3] * 255).astype(np.uint8)
    return _LUT_CACHE[name]


def render_variance(variance: np.ndarray, size: int = 512, resample: str = "nearest") -> bytes:
    sigma = np.sqrt(variance)
    v = sigma[np.isfinite(sigma)]
    lo, hi = (np.percentile(v, [1, 99]) if v.size else (0, 1))
    x = np.clip(np.nan_to_num((sigma - lo) / (hi - lo if hi > lo else 1.0), nan=0.0), 0, 1)
    rgb = _lut("magma")[(x * 255).astype(np.uint8)]
    return _to_png(rgb, size, resample, "RGB")


def render_flags(flags: np.ndarray, coverage: np.ndarray, bits: dict[str, int], size: int = 512) -> bytes:
    """Transparent RGBA overlay: red = defect, orange = artifact (ghost/halo/
    crosstalk...), green = known SOURCE, grey = no data."""
    bad = (flags & bad_flag_mask(bits)) != 0
    art_mask = 0
    for name, b in bits.items():
        if name in INFO_FLAG_NAMES and name not in ("SOURCE", "FULLSAMPLE"):
            art_mask |= 1 << b
    art = (flags & art_mask) != 0
    src = (flags & source_flag_mask(bits)) != 0
    rgba = np.zeros(flags.shape + (4,), np.uint8)
    rgba[src] = (40, 220, 120, 90)
    rgba[art] = (255, 160, 0, 130)
    rgba[bad] = (255, 40, 40, 200)
    rgba[coverage == 0] = (90, 90, 90, 160)
    return _to_png(rgba, size, "nearest", "RGBA")


def png_data_url(png: bytes) -> str:
    return "data:image/png;base64," + base64.b64encode(png).decode()


def _normalize_band_snr(arr: np.ndarray | None) -> np.ndarray | None:
    """Convert a sub-band co-add into background-subtracted sigma units (`>= 0`)."""
    if arr is None:
        return None
    v = arr[np.isfinite(arr)]
    if v.size < 8:
        return None
    med = float(np.median(v))
    mad = float(1.4826 * np.median(np.abs(v - med)))
    clipped = v[np.abs(v - med) <= 3.0 * max(mad, 1e-9)]
    bg = float(np.median(clipped)) if clipped.size else med
    sig = float(1.4826 * np.median(np.abs(clipped - bg))) if clipped.size > 4 else mad
    sig = max(sig, 1e-6)
    return np.clip(np.nan_to_num((arr.astype(np.float64) - bg) / sig, nan=0.0), 0.0, None)


def render_rgb_composite(
    c: EpochCoadd,
    size: int = 512,
    resample: str = "nearest",
    asinh_q: float = 6.0,
    max_snr: float = 25.0,
) -> bytes:
    """Render WiseView-style false-color RGB PNG from the epoch's 3 IR sub-bands:
      * Blue  = D1 + D2 (0.75 - 1.62 um, stellar photosphere)
      * Green = D3 + D4 (1.63 - 3.82 um, H2O/CH4 molecular region)
      * Red   = D5 + D6 (3.83 - 5.00 um, 4.5 um thermal peak of brown dwarfs)
    Uses Lupton color-preserving arcsinh stretch so bright stars retain chromaticity
    and cold brown dwarfs / Planet X candidates glow bright orange-red.
    """
    b_raw, g_raw, r_raw = c.rgb_bands
    b = _normalize_band_snr(b_raw)
    g = _normalize_band_snr(g_raw)
    r = _normalize_band_snr(r_raw)
    tot = _normalize_band_snr(c.image)
    if tot is None:
        tot = np.zeros_like(c.image, dtype=np.float64)

    # Graceful fallback when only a subset of detector bands was requested/available
    present = sum(ch is not None for ch in (b, g, r))
    if present == 0:
        b = g = r = tot
    elif present == 1:
        if r is not None:
            # Red-only request -> warm WiseView W2-like amber palette
            g = 0.55 * r
            b = 0.18 * r
        elif b is not None:
            g = 0.70 * b
            r = 0.25 * b
        else:
            assert g is not None
            r = 0.65 * g
            b = 0.65 * g
    elif present == 2:
        if g is None and b is not None and r is not None:
            g = 0.5 * (b + r)
        elif b is None and g is not None and r is not None:
            b = 0.55 * g
        elif r is None and b is not None and g is not None:
            r = 0.55 * g

    assert b is not None and g is not None and r is not None
    # Lupton et al. (2004) intensity-preserving arcsinh stretch
    intensity = (r + g + b) / 3.0
    norm_factor = math.asinh(asinh_q)
    with np.errstate(invalid="ignore", divide="ignore"):
        stretched_i = np.arcsinh(asinh_q * (intensity / max_snr)) / norm_factor
        scale = np.where(intensity > 1e-9, stretched_i / intensity, 0.0)

    rgb_float = np.stack([r * scale, g * scale, b * scale], axis=-1)
    # Desaturate cores that exceed 1.0 in any channel while preserving hue
    max_ch = np.max(rgb_float, axis=-1, keepdims=True)
    rgb_float = np.where(max_ch > 1.0, rgb_float / np.maximum(max_ch, 1.0), rgb_float)
    rgb_u8 = (np.clip(rgb_float, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)
    return _to_png(rgb_u8, size, resample, "RGB")


# ============================================================ co-add FITS IO ==
def coadd_to_fits_bytes(c: EpochCoadd, grid: DisplayGrid) -> bytes:
    hdr = grid.wcs.to_header()
    hdr["BUNIT"] = "MJy/sr"
    hdr["MJD-MID"] = c.mjd_mid
    hdr["MJD-BEG"] = c.mjd_start
    hdr["MJD-END"] = c.mjd_end
    hdr["NFRAMES"] = c.n_frames
    hdr["EPOCH"] = c.index
    hdul = fits.HDUList(
        [
            fits.PrimaryHDU(c.image.astype(np.float32), header=hdr),
            fits.ImageHDU(c.variance.astype(np.float32), header=hdr, name="VARIANCE"),
        ]
    )
    buf = io.BytesIO()
    hdul.writeto(buf)
    return buf.getvalue()


def read_coadd(path_or_bytes: Path | bytes) -> np.ndarray:
    src = io.BytesIO(path_or_bytes) if isinstance(path_or_bytes, (bytes, bytearray)) else path_or_bytes
    with fits.open(src, memmap=False) as h:
        return h[0].data.astype(np.float32)


def read_coadd_with_var(path_or_bytes: Path | bytes) -> tuple[np.ndarray, np.ndarray, float]:
    """Read `(image, variance, mjd_mid)` from a cached epoch co-add FITS."""
    src = io.BytesIO(path_or_bytes) if isinstance(path_or_bytes, (bytes, bytearray)) else path_or_bytes
    with fits.open(src, memmap=False) as h:
        img = h[0].data.astype(np.float32)
        mjd_mid = float(h[0].header.get("MJD-MID", 61000.0))
        if "VARIANCE" in h and h["VARIANCE"].data is not None:
            var = h["VARIANCE"].data.astype(np.float32)
        else:
            v = img[np.isfinite(img)]
            mad = float(1.4826 * np.median(np.abs(v - np.median(v)))) if v.size else 0.05
            var = np.full_like(img, max(mad**2, 1e-6))
        return img, var, mjd_mid


# ============================================== difference imaging & dipoles ==
def compute_difference_map(
    img_a: np.ndarray,
    var_a: np.ndarray,
    img_b: np.ndarray,
    var_b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute `(diff_img, snr_map)` for `Epoch_B - Epoch_A` on the aligned WCS grid."""
    bg_a = _robust_background(img_a)
    bg_b = _robust_background(img_b)
    diff = (img_b.astype(np.float64) - bg_b) - (img_a.astype(np.float64) - bg_a)

    finite = np.isfinite(diff)
    if finite.sum() >= 10:
        med = float(np.median(diff[finite]))
        mad = float(1.4826 * np.median(np.abs(diff[finite] - med)))
        diff = diff - med
    else:
        mad = 0.02

    # Empirical MAD floor prevents tiny PSF-wing residuals on ultra-bright stars
    # from dominating the formal variance
    var_comb = np.where(np.isfinite(var_a) & (var_a > 0), var_a, mad**2) + np.where(
        np.isfinite(var_b) & (var_b > 0), var_b, mad**2
    )
    sigma = np.sqrt(np.maximum(var_comb, (0.75 * mad) ** 2 + 1e-12))
    snr_map = np.where(finite, diff / sigma, 0.0)
    return diff.astype(np.float32), snr_map.astype(np.float32)


def render_difference_png(
    snr_map: np.ndarray,
    size: int = 512,
    resample: str = "nearest",
    clip_sigma: float = 8.0,
) -> bytes:
    """Render a bipolar diverging difference image:
      * Blue/Cyan  (S < 0): source was at this pixel in Epoch A (or dimmed)
      * Dark Grey  (S ~ 0): stationary background & static stars cancel out
      * Red/Orange (S > 0): source moved to this pixel in Epoch B (or brightened)
    """
    s = np.nan_to_num(snr_map.astype(np.float64), nan=0.0)
    # Symmetric arcsinh compression around 0
    norm = math.asinh(clip_sigma / 2.0)
    u = np.clip(np.arcsinh(s / 2.0) / (norm if norm > 0 else 1.0), -1.0, 1.0)
    idx = ((u * 0.5 + 0.5) * 255.0 + 0.5).astype(np.uint8)
    rgb = _lut("coolwarm")[idx]
    return _to_png(rgb, size, resample, "RGB")


def detect_diff_anomalies(
    diff_img: np.ndarray,
    snr_map: np.ndarray,
    grid: DisplayGrid,
    mjd_a: float,
    mjd_b: float,
    snr_threshold: float = 3.8,
    max_anomalies: int = 12,
    coadds: list[EpochCoadd] | None = None,
) -> list[DiffAnomaly]:
    """Automatically detect moving dipoles (proper-motion objects) and variable
    stars (monopole brightenings / dimmings) on the epoch difference map."""
    ny, nx = snr_map.shape
    if ny < 7 or nx < 7:
        return []

    # 3x3 local extrema detection (excluding 2-pixel border)
    pos_peaks: list[tuple[int, int, float]] = []
    neg_peaks: list[tuple[int, int, float]] = []
    for y in range(2, ny - 2):
        for x in range(2, nx - 2):
            val = float(snr_map[y, x])
            if not np.isfinite(val):
                continue
            patch = snr_map[y - 1 : y + 2, x - 1 : x + 2]
            if val >= snr_threshold and val >= float(np.nanmax(patch)):
                pos_peaks.append((x, y, val))
            elif val <= -snr_threshold and val <= float(np.nanmin(patch)):
                neg_peaks.append((x, y, val))

    pos_peaks.sort(key=lambda t: t[2], reverse=True)
    neg_peaks.sort(key=lambda t: t[2])  # most negative first

    pix_area_sr = (grid.pixscale_deg * math.pi / 180.0) ** 2
    to_ujy = pix_area_sr * 1e12
    dt_days = max(abs(mjd_b - mjd_a), 1.0)
    dt_yr = dt_days / 365.25
    wcs = grid.wcs_copy()

    def _px_to_world_and_png(px: int, py: int) -> tuple[float, float, float, float]:
        ra_v, dec_v = wcs.pixel_to_world_values(px, py)
        k = grid.png_size / grid.n
        png_x = (px + 0.5) * k
        png_y = (grid.n - 0.5 - py) * k
        return float(ra_v), float(dec_v), round(float(png_x), 2), round(float(png_y), 2)

    used_pos: set[int] = set()
    used_neg: set[int] = set()
    anomalies: list[DiffAnomaly] = []

    # 1. Pair positive & negative peaks within [1.0, 10.0] native pixels -> moving dipoles!
    for ip, (xp, yp, sp) in enumerate(pos_peaks):
        best_in = None
        best_dist = 999.0
        fp_box = float(np.nansum(diff_img[yp - 1 : yp + 2, xp - 1 : xp + 2])) * to_ujy
        for ineg, (xn, yn, sn) in enumerate(neg_peaks):
            if ineg in used_neg:
                continue
            dist_px = math.hypot(xp - xn, yp - yn)
            if not (1.0 <= dist_px <= 10.0):
                continue
            ratio = sp / abs(sn)
            if not (0.25 <= ratio <= 4.0):
                continue
            if dist_px < best_dist:
                best_dist = dist_px
                best_in = (ineg, xn, yn, sn)
        if best_in is not None:
            ineg, xn, yn, sn = best_in
            used_pos.add(ip)
            used_neg.add(ineg)
            ra_to, dec_to, png_x_to, png_y_to = _px_to_world_and_png(xp, yp)
            ra_fr, dec_fr, png_x_fr, png_y_fr = _px_to_world_and_png(xn, yn)
            shift_arcsec = best_dist * (grid.pixscale_deg * 3600.0)
            pm_yr = shift_arcsec / dt_yr if dt_yr > 0.02 else None
            comb_snr = math.hypot(sp, sn)
            anomalies.append(
                DiffAnomaly(
                    anomaly_type="moving_dipole",
                    snr=round(comb_snr, 2),
                    ra=round(ra_to, 6),
                    dec=round(dec_to, 6),
                    x=png_x_to,
                    y=png_y_to,
                    ra_from=round(ra_fr, 6),
                    dec_from=round(dec_fr, 6),
                    x_from=png_x_fr,
                    y_from=png_y_fr,
                    shift_arcsec=round(shift_arcsec, 2),
                    pm_arcsec_yr=round(pm_yr, 2) if pm_yr is not None else None,
                    delta_flux_ujy=round(fp_box, 2),
                )
            )

    # 2. Remaining unpaired strong extrema -> variable stars / transients
    for ip, (xp, yp, sp) in enumerate(pos_peaks):
        if ip in used_pos or sp < snr_threshold + 0.5:
            continue
        ra_v, dec_v, png_x, png_y = _px_to_world_and_png(xp, yp)
        df = float(np.nansum(diff_img[yp - 1 : yp + 2, xp - 1 : xp + 2])) * to_ujy
        anomalies.append(
            DiffAnomaly(
                anomaly_type="variable_brightening",
                snr=round(sp, 2),
                ra=round(ra_v, 6),
                dec=round(dec_v, 6),
                x=png_x,
                y=png_y,
                delta_flux_ujy=round(df, 2),
            )
        )

    for ineg, (xn, yn, sn) in enumerate(neg_peaks):
        if ineg in used_neg or abs(sn) < snr_threshold + 0.5:
            continue
        ra_v, dec_v, png_x, png_y = _px_to_world_and_png(xn, yn)
        df = float(np.nansum(diff_img[yn - 1 : yn + 2, xn - 1 : xn + 2])) * to_ujy
        anomalies.append(
            DiffAnomaly(
                anomaly_type="variable_dimming",
                snr=round(abs(sn), 2),
                ra=round(ra_v, 6),
                dec=round(dec_v, 6),
                x=png_x,
                y=png_y,
                delta_flux_ujy=round(df, 2),
            )
        )

    anomalies.sort(key=lambda a: a.snr, reverse=True)
    result = anomalies[:max_anomalies]

    if coadds is not None and len(coadds) >= 3:
        from app.services.astrometry_fitter import fit_motion_model, measure_anomaly_centroids

        for anom in result:
            if anom.anomaly_type != "moving_dipole":
                continue
            centroids = measure_anomaly_centroids(coadds, grid, anom, mjd_a=mjd_a, mjd_b=mjd_b)
            model, plx = fit_motion_model(centroids)
            anom.motion_model = model
            anom.parallax_arcsec = plx

    return result


# ================================================================ detection ==
def aperture_snr(image: np.ndarray, variance: np.ndarray, grid: DisplayGrid, ra: float, dec: float, search_px: int = 1) -> float | None:
    """3x3 box S/N at the target (best within +-search_px to tolerate slight
    motion / centroid offsets), background from a 4-7 px annulus."""
    x, y = grid.wcs_copy().world_to_pixel_values(ra, dec)
    xi, yi = int(round(float(x))), int(round(float(y)))
    n = image.shape[0]
    yy, xx = np.mgrid[0:n, 0:n]
    best = None
    for dy in range(-search_px, search_px + 1):
        for dx in range(-search_px, search_px + 1):
            cx, cy = xi + dx, yi + dy
            if not (1 <= cx < n - 1 and 1 <= cy < n - 1):
                continue
            r = np.hypot(xx - cx, yy - cy)
            ann = image[(r >= 4) & (r <= 7) & np.isfinite(image)]
            ap = image[cy - 1 : cy + 2, cx - 1 : cx + 2]
            apv = variance[cy - 1 : cy + 2, cx - 1 : cx + 2]
            good = np.isfinite(ap) & np.isfinite(apv) & (apv > 0)
            n_good = int(good.sum())
            if ann.size < 8 or not good[1, 1] or n_good < 7:
                continue
            bg = float(np.median(ann))
            scale = 9.0 / n_good
            sig = math.sqrt(float(apv[good].sum()) * scale**2 + 9 * (1.2533 * float(np.std(ann))) ** 2 * 9 / ann.size)
            snr = float((ap[good] - bg).sum()) * scale / sig if sig > 0 else None
            if snr is not None and (best is None or snr > best):
                best = snr
    return best


# ================================================================= spectrum ==
@dataclass
class SpectralSample:
    detector: int
    wavelength_um: float
    bandwidth_um: float
    flux_ujy: float
    error_ujy: float
    mjd: float


def extract_sample(fd: FrameData, ra: float, dec: float) -> SpectralSample | None:
    """One LVF sample: forced 3x3 aperture photometry at (ra, dec).

    The wavelength comes from the per-pixel spectral WCS (WCS-WAVE lookup
    table, alternate WCS key 'W') evaluated at the *sub-pixel* source position.
    Background: median of a 3-5 px annulus (ZODI already removed).
    No aperture correction: fluxes are ~10-25% low but relative shape between
    channels of a detector is preserved (upgrade path: PSF fit with EPSF).
    """
    if fd.wave_wcs is None and fd.crpix2w is None:
        return None
    x, y = fd.wcs.world_to_pixel_values(ra, dec)
    if not (np.isfinite(x) and np.isfinite(y)):
        return None
    ny, nx = fd.image.shape
    xi, yi = int(round(float(x))), int(round(float(y)))
    if not (5 <= xi < nx - 5 and 5 <= yi < ny - 5):
        return None
    ap = fd.image[yi - 1 : yi + 2, xi - 1 : xi + 2]
    apv = fd.variance[yi - 1 : yi + 2, xi - 1 : xi + 2]
    good = np.isfinite(ap) & np.isfinite(apv) & (apv > 0)
    n_good = int(good.sum())
    if not good[1, 1] or n_good < 7:
        return None  # central pixel must be clean; up to 2 outer pixels may be masked
    yy, xx = np.mgrid[yi - 5 : yi + 6, xi - 5 : xi + 6]
    r = np.hypot(xx - x, yy - y)
    sub = fd.image[yi - 5 : yi + 6, xi - 5 : xi + 6]
    ann = sub[(r >= 3) & (r <= 5) & np.isfinite(sub)]
    if ann.size < 10:
        return None
    bg = float(np.median(ann))
    scale = 9.0 / n_good
    to_ujy = fd.pix_area_sr * 1e12  # (MJy/sr) * sr -> MJy -> uJy
    flux = float((ap[good] - bg).sum()) * scale * to_ujy
    err = math.sqrt(float(apv[good].sum()) * scale**2 + 9 * 9 * (1.2533 * float(np.std(ann))) ** 2 / ann.size) * to_ujy
    wl_um = bw_um = None
    if fd.wave_wcs is not None:
        try:
            wl, bw = fd.wave_wcs.pixel_to_world(float(x), float(y))
            wl_um, bw_um = float(wl.to_value("um")), float(bw.to_value("um"))
        except Exception:
            wl_um = None
    if wl_um is None and fd.crpix2w is not None and fd.detector in DETECTORS:
        # Analytic LVF dispersion along the 2048-px H2RG detector Y axis
        det = DETECTORS[fd.detector]
        det_y = min(max(float(y) + 1.0 - fd.crpix2w, 1.0), 2040.0)
        frac = (det_y - 1.0) / 2039.0
        wl_um = det.lambda_max_um * (det.lambda_min_um / det.lambda_max_um) ** frac
        bw_um = wl_um / det.resolving_power
    if wl_um is None or bw_um is None or not (np.isfinite(wl_um) and 0.6 < wl_um < 5.2):
        return None
    return SpectralSample(fd.detector, wl_um, bw_um, flux, err, fd.mjd)


def bin_spectrum(samples: list[SpectralSample]) -> list[SpectrumPoint]:
    """Bin LVF samples into the 102 canonical channels (inverse-variance mean,
    5-sigma MAD clipping against cosmic rays / unflagged artifacts)."""
    groups: dict[int, list[SpectralSample]] = {}
    for s in samples:
        ch = channel_index(s.wavelength_um)
        if ch is not None and s.error_ujy > 0:
            groups.setdefault(ch, []).append(s)
    points = []
    for ch in range(N_CHANNELS):
        lo, hi = CHANNEL_EDGES_UM[ch]
        g = groups.get(ch, [])
        flux = err = None
        if g:
            f = np.array([s.flux_ujy for s in g])
            e = np.array([s.error_ujy for s in g])
            if f.size >= 3:
                med = np.median(f)
                mad = 1.4826 * np.median(np.abs(f - med))
                keep = np.abs(f - med) <= 5 * max(mad, 1e-12) + 3 * e
                f, e = f[keep], e[keep]
            if f.size:
                w = 1.0 / e**2
                flux = float((f * w).sum() / w.sum())
                err = float(1.0 / math.sqrt(w.sum()))
        points.append(
            SpectrumPoint(
                channel=ch,
                detector=ch // CHANNELS_PER_DETECTOR + 1,
                wavelength_um=round(float(CHANNEL_CENTERS_UM[ch]), 4),
                bandwidth_um=round(float(hi - lo), 4),
                flux=None if flux is None else round(flux, 3),
                error=None if err is None else round(err, 3),
                n_samples=len(g),
            )
        )
    return points


# ========================================== non-blocking async wrappers ==
async def reproject_frame_async(fd: FrameData, grid: DisplayGrid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run `reproject_interp` in a worker thread without blocking the event loop."""
    return await asyncio.to_thread(reproject_frame, fd, grid)


async def coadd_epoch_async(index: int, frames: list[FrameData], grid: DisplayGrid) -> EpochCoadd:
    """Run WCS reprojection (`reproject_interp`) and inverse-variance co-adding via `asyncio.to_thread`."""
    return await asyncio.to_thread(coadd_epoch, index, frames, grid)


async def render_image_async(
    data: np.ndarray,
    vmin: float,
    vmax: float,
    stretch_type: str = "asinh",
    size: int = 512,
    resample: str = "nearest",
    asinh_a: float = 0.1,
    log_a: float = 1000.0,
) -> bytes:
    """Render stretched grayscale PNG in a worker thread via `asyncio.to_thread`."""
    return await asyncio.to_thread(
        render_image, data, vmin, vmax, stretch_type, size, resample, asinh_a, log_a
    )


async def render_rgb_composite_async(
    c: EpochCoadd,
    size: int = 512,
    resample: str = "nearest",
    asinh_q: float = 6.0,
    max_snr: float = 25.0,
) -> bytes:
    """Render false-color RGB PNG in a worker thread via `asyncio.to_thread`."""
    return await asyncio.to_thread(render_rgb_composite, c, size, resample, asinh_q, max_snr)

