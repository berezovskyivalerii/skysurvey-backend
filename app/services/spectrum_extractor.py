"""SPHEREx 102-channel LVF spectrophotometry & astrophysical summary extractor.

Provides:
* `extract_spectrum(fits_files, ra, dec, wavecal_dir=...) -> dict`
  Reads a collection of SPHEREx Level 2 Multi-Extension FITS cutouts, converts
  sky coordinates (ra, dec) to sub-pixel detector coordinates (x, y), resolves
  the LVF central wavelength & bandwidth at (x, y) via `WCS-WAVE` (with automatic
  analytic fallback to the detector's `CRPIX2W` LVF gradient), performs ZODI- and
  annulus-background-subtracted 3x3 aperture photometry, evaluates `FLAGS`,
  aggregates duplicate channel hits by inverse-variance weighting, and computes
  summary diagnostics (`ir_color_index`, `peak_wavelength_um`, `mean_snr`,
  `spectral_hint`).
* `find_covering_fits(raw_dir, ra, dec) -> list[Path]`
  Scans locally cached FITS cutouts in `/tmp/fits_cache/raw/` and returns all
  files whose spatial WCS footprint covers `(ra, dec)`.
"""
from __future__ import annotations

import logging
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS, FITSFixedWarning

from app.core.spherex import (
    CHANNEL_CENTERS_UM,
    CHANNEL_EDGES_UM,
    CHANNELS_PER_DETECTOR,
    DETECTORS,
    N_CHANNELS,
    bad_flag_mask,
    channel_index,
    detector_from_bandname,
    parse_flag_bits,
)
from app.schemas.spectrum import (
    DEFAULT_ABSORPTION_BANDS,
    ChartPoint,
    SpectralHint,
    SpectrumResponse,
    SpectrumSummary,
    TargetCoords,
)

log = logging.getLogger(__name__)
warnings.filterwarnings("ignore", category=FITSFixedWarning)

WARN_FLAG_NAMES: frozenset[str] = frozenset({"DICHROIC", "NONLINEAR", "GHOST", "GHOST_FPA", "GHOST_EXT", "CROSSTALK"})


@dataclass
class RawSpectralPoint:
    channel_id: int  # 1..102
    detector: int  # 1..6
    wavelength_um: float
    bandwidth_um: float
    flux_ujy: float
    flux_error_ujy: float
    is_flagged: bool


def _warn_flag_mask(bits: dict[str, int]) -> int:
    mask = 0
    for name, b in bits.items():
        if name in WARN_FLAG_NAMES:
            mask |= 1 << b
    return mask


def _resolve_detector(hdr: fits.Header, primary_hdr: fits.Header) -> int:
    for h in (hdr, primary_hdr):
        det = h.get("DETECTOR")
        if isinstance(det, (int, np.integer)) and 1 <= int(det) <= 6:
            return int(det)
        if isinstance(det, str):
            d = detector_from_bandname(det)
            if d:
                return d
            if det.strip().isdigit() and 1 <= int(det.strip()) <= 6:
                return int(det.strip())
        band = h.get("BAND") or h.get("FILTER")
        if band:
            d = detector_from_bandname(str(band))
            if d:
                return d
    return 1


def _resolve_wavelength(
    hdul: fits.HDUList,
    hdr: fits.Header,
    detector: int,
    x: float,
    y: float,
    wavecal_dir: Path | None,
) -> tuple[float, float] | None:
    """Resolve (wavelength_um, bandwidth_um) at pixel (x, y) with 3-tier fallback."""
    # Tier 1: WCS-WAVE lookup table (embedded or pinned in wavecal_dir)
    if "WCS-WAVE" not in hdul and wavecal_dir is not None and wavecal_dir.is_dir():
        for candidate in (
            wavecal_dir / f"spherex_dr1_D{detector}.fits",
            wavecal_dir / f"D{detector}.fits",
        ):
            if candidate.is_file():
                try:
                    with fits.open(candidate, memmap=False) as wc:
                        if "WCS-WAVE" in wc:
                            hdul.append(wc["WCS-WAVE"].copy())
                            break
                except Exception:
                    pass

    if "WCS-WAVE" in hdul:
        try:
            whdr = hdr.copy()
            for k in list(whdr.keys()):
                if k.startswith(("A_", "B_", "AP_", "BP_")):
                    del whdr[k]
            wave_wcs = WCS(whdr, fobj=hdul, key="W")
            wl, bw = wave_wcs.pixel_to_world(float(x), float(y))
            wl_um = float(wl.to_value("um"))
            bw_um = float(bw.to_value("um"))
            if np.isfinite(wl_um) and 0.65 <= wl_um <= 5.20:
                if not (np.isfinite(bw_um) and bw_um > 0):
                    bw_um = wl_um / DETECTORS.get(detector, DETECTORS[1]).resolving_power
                return wl_um, bw_um
        except Exception as e:
            log.debug("WCS-WAVE lookup failed, using analytic LVF fallback: %s", e)

    # Tier 2: Analytic LVF gradient along detector Y using CRPIX2W
    det_cfg = DETECTORS.get(detector, DETECTORS[1])
    if "CRPIX2W" in hdr:
        try:
            crpix2w = float(hdr["CRPIX2W"])
            det_y = min(max(float(y) + 1.0 - crpix2w, 1.0), 2040.0)
            frac = (det_y - 1.0) / 2039.0
            wl_um = det_cfg.lambda_max_um * (det_cfg.lambda_min_um / det_cfg.lambda_max_um) ** frac
            bw_um = wl_um / det_cfg.resolving_power
            return float(wl_um), float(bw_um)
        except Exception:
            pass

    # Tier 3: Explicit header wavelength keywords or detector midpoint
    for kw in ("WAVELENG", "EFF_WAVE", "LAMBDA", "CENWAVE"):
        if kw in hdr:
            try:
                wl_um = float(hdr[kw])
                if 0.65 <= wl_um <= 5.20:
                    return wl_um, wl_um / det_cfg.resolving_power
            except Exception:
                pass
    return None


def _header_near_target(hdr: fits.Header, ra: float, dec: float, max_deg: float = 4.0) -> bool:
    """Fast check whether a FITS header's tangent point (CRVAL1, CRVAL2) is within
    `max_deg` of `(ra, dec)` before invoking WCS/SIP Newton-Raphson inversion."""
    try:
        crval1 = float(hdr["CRVAL1"])
        crval2 = float(hdr["CRVAL2"])
    except (KeyError, TypeError, ValueError):
        return True
    if abs(crval2 - dec) > max_deg:
        return False
    dra = abs((crval1 - ra + 180.0) % 360.0 - 180.0) * math.cos(math.radians(dec))
    return dra <= max_deg


def _photometer_single_fits(
    path: Path,
    ra: float,
    dec: float,
    wavecal_dir: Path | None,
) -> RawSpectralPoint | None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with fits.open(path, memmap=False) as hdul:
            img_hdu = hdul["IMAGE"] if "IMAGE" in hdul else (hdul[1] if len(hdul) > 1 else hdul[0])
            if img_hdu.data is None:
                return None
            hdr = img_hdu.header
            if not _header_near_target(hdr, ra, dec):
                return None
            raw_img = img_hdu.data.astype(np.float64)
            ny, nx = raw_img.shape

            wcs = WCS(hdr).celestial
            x_arr, y_arr = wcs.world_to_pixel_values(ra, dec)
            x, y = float(x_arr), float(y_arr)
            if not (np.isfinite(x) and np.isfinite(y)):
                return None
            xi, yi = int(round(x)), int(round(y))
            if not (3 <= xi < nx - 3 and 3 <= yi < ny - 3):
                return None

            detector = _resolve_detector(hdr, hdul[0].header)
            wl_bw = _resolve_wavelength(hdul, hdr, detector, x, y, wavecal_dir)
            if wl_bw is None:
                return None
            wl_um, bw_um = wl_bw

            # Subtract ZODI model if present
            if "ZODI" in hdul and hdul["ZODI"].data is not None and hdul["ZODI"].data.shape == raw_img.shape:
                raw_img = raw_img - hdul["ZODI"].data.astype(np.float64)

            if "VARIANCE" in hdul and hdul["VARIANCE"].data is not None and hdul["VARIANCE"].data.shape == raw_img.shape:
                raw_var = hdul["VARIANCE"].data.astype(np.float64)
            else:
                finite_vals = raw_img[np.isfinite(raw_img)]
                mad_sig = 1.4826 * float(np.median(np.abs(finite_vals - np.median(finite_vals)))) if finite_vals.size else 0.05
                raw_var = np.full_like(raw_img, max(mad_sig**2, 1e-6))

            if "FLAGS" in hdul and hdul["FLAGS"].data is not None and hdul["FLAGS"].data.shape == raw_img.shape:
                flg = hdul["FLAGS"].data.astype(np.int64)
                bits = parse_flag_bits(hdul["FLAGS"].header)
            else:
                flg = np.zeros_like(raw_img, dtype=np.int64)
                bits = parse_flag_bits({})

            try:
                pix_area_sr = float(wcs.proj_plane_pixel_area() * (math.pi / 180.0) ** 2)
            except Exception:
                pix_area_sr = 8.91e-10

    bad_mask = bad_flag_mask(bits)
    warn_mask = _warn_flag_mask(bits)

    # 3x3 aperture around (xi, yi)
    ap_img = raw_img[yi - 1 : yi + 2, xi - 1 : xi + 2]
    ap_var = raw_var[yi - 1 : yi + 2, xi - 1 : xi + 2]
    ap_flg = flg[yi - 1 : yi + 2, xi - 1 : xi + 2]

    finite_mask = np.isfinite(ap_img) & np.isfinite(ap_var) & (ap_var > 0)
    if not finite_mask.any():
        return None

    clean_mask = finite_mask & ((ap_flg & bad_mask) == 0)
    center_bad = bool((ap_flg[1, 1] & bad_mask) != 0) or not bool(finite_mask[1, 1])
    any_warn = bool(((ap_flg & warn_mask) != 0).any())
    n_clean = int(clean_mask.sum())

    if n_clean < 4:
        # More than half the 3x3 aperture is destroyed by defects/transients -> drop exposure
        # so uncalibrated cosmic-ray spikes never blow up the frontend chart Y-axis.
        return None

    # Flag the point if the central pixel is defective, >2 aperture pixels are bad, or warning bits are set
    is_flagged = center_bad or (n_clean < 7) or any_warn
    use_mask = clean_mask
    n_used = n_clean

    # Local background from a 3..6 px annulus (excluding bad-flagged pixels)
    r_max = min(6, xi, yi, nx - 1 - xi, ny - 1 - yi)
    yy, xx = np.mgrid[yi - r_max : yi + r_max + 1, xi - r_max : xi + r_max + 1]
    rr = np.hypot(xx - x, yy - y)
    sub_img = raw_img[yi - r_max : yi + r_max + 1, xi - r_max : xi + r_max + 1]
    sub_flg = flg[yi - r_max : yi + r_max + 1, xi - r_max : xi + r_max + 1]
    ann = sub_img[(rr >= 2.8) & (rr <= max(r_max, 3.5)) & np.isfinite(sub_img) & ((sub_flg & bad_mask) == 0)]
    if ann.size >= 6:
        med = float(np.median(ann))
        mad = 1.4826 * float(np.median(np.abs(ann - med)))
        clipped = ann[np.abs(ann - med) <= 3.0 * max(mad, 1e-9)]
        bg = float(np.median(clipped)) if clipped.size else med
        bg_std = float(np.std(clipped)) if clipped.size > 1 else mad
        n_ann = max(int(clipped.size), 1)
    else:
        bg = 0.0
        bg_std = 0.0
        n_ann = 1

    scale = 9.0 / n_used
    to_ujy = pix_area_sr * 1e12  # MJy/sr * sr -> MJy -> uJy
    net_sb = float((ap_img[use_mask] - bg).sum()) * scale
    var_sb = float(ap_var[use_mask].sum()) * (scale**2) + 81.0 * (1.2533 * bg_std) ** 2 / n_ann
    flux_ujy = net_sb * to_ujy
    err_ujy = math.sqrt(max(var_sb, 1e-18)) * to_ujy
    if not (np.isfinite(flux_ujy) and np.isfinite(err_ujy) and err_ujy > 0):
        return None

    ch0 = channel_index(wl_um)
    if ch0 is None:
        # Clamp to nearest edge channel if slightly outside [0.75, 5.00]
        ch0 = 0 if wl_um < CHANNEL_CENTERS_UM[0] else N_CHANNELS - 1

    return RawSpectralPoint(
        channel_id=ch0 + 1,
        detector=detector,
        wavelength_um=wl_um,
        bandwidth_um=bw_um,
        flux_ujy=flux_ujy,
        flux_error_ujy=err_ujy,
        is_flagged=is_flagged,
    )


def _classify_spectrum(
    ir_color_index: float | None,
    mean_snr: float,
    max_snr: float,
    peak_wl: float | None,
    chart_data: list[ChartPoint],
) -> SpectralHint:
    if not chart_data or (max_snr < 2.5 and mean_snr < 1.2):
        return "no_significant_detection"
    if ir_color_index is not None:
        if ir_color_index >= 5.0 and (peak_wl or 0.0) >= 3.6:
            return "ultra_cool_y_dwarf_or_planet_x_candidate"
        if ir_color_index >= 2.0:
            return "cold_brown_dwarf_candidate"
        if ir_color_index >= 1.0:
            return "red_stellar_or_dusty_source"
        return "blue_main_sequence_star"
    # Only a subset of bands covered: check if red bands (3-5 um) show strong detection
    red_pts = [p for p in chart_data if p.wavelength_um >= 3.0 and not p.is_flagged and p.snr >= 3.0]
    if red_pts and (peak_wl or 0.0) >= 3.6:
        return "cold_brown_dwarf_candidate"
    return "insufficient_wavelength_coverage"


def extract_spectrum(
    fits_files: list[Path],
    ra: float,
    dec: float,
    wavecal_dir: Path | None = None,
) -> dict[str, Any]:
    """Extract, aggregate, and analyze the 102-channel SPHEREx spectrum at `(ra, dec)`.

    Returns a dictionary matching `SpectrumResponse` (`target`, `spectrum_summary`,
    `chart_data`, `absorption_bands_reference`). Never raises on missing `WCS-WAVE`
    or partial channel coverage.
    """
    if wavecal_dir is None:
        default_wc = Path("/tmp/fits_cache/wavecal")
        if default_wc.is_dir():
            wavecal_dir = default_wc

    raw_points: list[RawSpectralPoint] = []
    n_exposures_used = 0
    for fpath in fits_files:
        p = Path(fpath)
        if not p.is_file():
            continue
        try:
            pt = _photometer_single_fits(p, ra, dec, wavecal_dir)
            if pt is not None:
                raw_points.append(pt)
                n_exposures_used += 1
        except Exception as e:
            log.debug("Skipping FITS file %s in extract_spectrum: %s", p.name, e)

    # Group by channel_id (1..102)
    by_channel: dict[int, list[RawSpectralPoint]] = {}
    for pt in raw_points:
        by_channel.setdefault(pt.channel_id, []).append(pt)

    chart_data: list[ChartPoint] = []
    for ch_id in sorted(by_channel.keys()):
        samples = by_channel[ch_id]
        clean_samples = [s for s in samples if not s.is_flagged]
        if clean_samples:
            active = clean_samples
            flagged_channel = False
        else:
            active = samples
            flagged_channel = True

        fluxes = np.array([s.flux_ujy for s in active], dtype=np.float64)
        errs = np.array([s.flux_error_ujy for s in active], dtype=np.float64)
        wls = np.array([s.wavelength_um for s in active], dtype=np.float64)
        bws = np.array([s.bandwidth_um for s in active], dtype=np.float64)

        # 5-sigma MAD outlier clipping when >= 3 exposures hit the same channel
        if fluxes.size >= 3:
            med = float(np.median(fluxes))
            mad = 1.4826 * float(np.median(np.abs(fluxes - med)))
            keep = np.abs(fluxes - med) <= 5.0 * max(mad, 1e-9) + 3.0 * errs
            if keep.any():
                fluxes, errs, wls, bws = fluxes[keep], errs[keep], wls[keep], bws[keep]

        weights = 1.0 / np.maximum(errs**2, 1e-18)
        w_sum = float(weights.sum())
        flux_val = float((fluxes * weights).sum() / w_sum)
        err_val = float(1.0 / math.sqrt(w_sum))
        wl_val = float((wls * weights).sum() / w_sum)
        bw_val = float((bws * weights).sum() / w_sum)
        snr_val = flux_val / err_val if err_val > 0 else 0.0

        det = (ch_id - 1) // CHANNELS_PER_DETECTOR + 1
        chart_data.append(
            ChartPoint(
                channel_id=ch_id,
                wavelength_um=round(wl_val, 4),
                bandwidth_um=round(bw_val, 4),
                band=f"Band {det}",
                detector=det,
                flux=round(flux_val, 3),
                flux_error=round(err_val, 3),
                snr=round(snr_val, 2),
                is_flagged=flagged_channel,
                n_samples=len(samples),
            )
        )

    # Strictly sort by ascending wavelength_um
    chart_data.sort(key=lambda c: c.wavelength_um)

    valid_pts = [p for p in chart_data if not p.is_flagged]
    flagged_pts = [p for p in chart_data if p.is_flagged]
    ref_pts = valid_pts if valid_pts else chart_data

    if ref_pts:
        snrs = [p.snr for p in ref_pts]
        mean_snr = round(float(np.mean(snrs)), 2)
        max_snr = round(float(np.max(snrs)), 2)
        peak_pt = max(ref_pts, key=lambda p: p.flux)
        peak_wl = peak_pt.wavelength_um
        peak_flux = peak_pt.flux
    else:
        mean_snr = 0.0
        max_snr = 0.0
        peak_wl = None
        peak_flux = None

    # Infrared color index: F(3.0-5.0 um) / F(0.75-1.50 um)
    short_pts = [p for p in ref_pts if 0.75 <= p.wavelength_um <= 1.55 and p.flux_error > 0]
    long_pts = [p for p in ref_pts if 3.00 <= p.wavelength_um <= 5.05 and p.flux_error > 0]
    ir_color: float | None = None
    if short_pts and long_pts:
        w_s = np.array([1.0 / (p.flux_error**2) for p in short_pts])
        f_s = float(np.sum([p.flux * w for p, w in zip(short_pts, w_s)]) / np.sum(w_s))
        sig_s = float(1.0 / math.sqrt(float(np.sum(w_s))))
        w_l = np.array([1.0 / (p.flux_error**2) for p in long_pts])
        f_l = float(np.sum([p.flux * w for p, w in zip(long_pts, w_l)]) / np.sum(w_l))
        # Regularize denominator at the 1-sigma noise floor so cold brown dwarfs
        # with near-zero flux at 1 um don't divide by zero or negative noise
        denom = max(f_s, sig_s, 1.0)
        if f_l > 0:
            ir_color = round(f_l / denom, 3)

    hint = _classify_spectrum(ir_color, mean_snr, max_snr, peak_wl, chart_data)

    resp = SpectrumResponse(
        target=TargetCoords(ra=round(ra, 6), dec=round(dec, 6)),
        spectrum_summary=SpectrumSummary(
            total_channels=N_CHANNELS,
            valid_points=len(valid_pts),
            flagged_points=len(flagged_pts),
            n_exposures_used=n_exposures_used,
            mean_snr=mean_snr,
            max_snr=max_snr,
            peak_wavelength_um=peak_wl,
            peak_flux_ujy=peak_flux,
            ir_color_index=ir_color,
            spectral_hint=hint,
        ),
        chart_data=chart_data,
        absorption_bands_reference=list(DEFAULT_ABSORPTION_BANDS),
    )
    return resp.model_dump(mode="json")


def find_covering_fits(raw_dir: Path, ra: float, dec: float) -> list[Path]:
    """Scan local FITS cache directory and return all files whose WCS covers `(ra, dec)`."""
    if not raw_dir.is_dir():
        return []
    matched: list[Path] = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for p in sorted(raw_dir.glob("*.fits")):
            try:
                with fits.open(p, memmap=False) as hdul:
                    hdu = hdul["IMAGE"] if "IMAGE" in hdul else (hdul[1] if len(hdul) > 1 else hdul[0])
                    if hdu.data is None:
                        continue
                    if not _header_near_target(hdu.header, ra, dec):
                        continue
                    ny, nx = hdu.data.shape
                    wcs = WCS(hdu.header).celestial
                    x, y = wcs.world_to_pixel_values(ra, dec)
                    if np.isfinite(x) and np.isfinite(y) and 3 <= float(x) < nx - 3 and 3 <= float(y) < ny - 3:
                        matched.append(p)
            except Exception:
                continue
    return matched
