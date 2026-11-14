"""Multi-epoch astrometric fitter: linear proper motion vs annual parallax.

Distinguishes nearby brown dwarfs / high-PM stars (linear trajectory across 3+
epochs) from distant Solar System bodies like Planet X (300-500 AU, whose
apparent half-year displacement reverses sign over a 1-year baseline due to
Earth's 2 AU orbital parallax).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from astropy.coordinates import get_sun
from astropy.time import Time

from app.schemas.target import DiffAnomaly, MotionModel

if TYPE_CHECKING:
    from app.services.fits_processor import DisplayGrid


@dataclass
class EpochCentroid:
    epoch_index: int
    mjd: float
    ra: float
    dec: float
    x_px: float
    y_px: float
    sigma_arcsec: float


def _local_background(patch: np.ndarray) -> tuple[float, float]:
    vals = patch[np.isfinite(patch)]
    if vals.size < 4:
        return 0.0, 0.02
    med = float(np.median(vals))
    mad = float(1.4826 * np.median(np.abs(vals - med)))
    clipped = vals[np.abs(vals - med) <= 3.0 * max(mad, 1e-9)]
    bg = float(np.median(clipped)) if clipped.size else med
    sig = float(1.4826 * np.median(np.abs(clipped - bg))) if clipped.size > 2 else max(mad, 1e-4)
    return bg, max(sig, 1e-5)


def measure_subpixel_centroid(
    image: np.ndarray,
    variance: np.ndarray,
    grid: DisplayGrid,
    approx_x: float,
    approx_y: float,
    box_radius: int = 3,
) -> tuple[float, float, float, float, float] | None:
    """Measure sub-pixel flux-weighted centroid around `(approx_x, approx_y)` in native grid pixels.

    Returns `(ra, dec, x_sub, y_sub, sigma_arcsec)` or `None` if no significant source is present.
    """
    ny, nx = image.shape
    cx0 = int(round(approx_x))
    cy0 = int(round(approx_y))
    if not (1 <= cx0 < nx - 1 and 1 <= cy0 < ny - 1):
        return None

    x0 = max(1, cx0 - box_radius)
    x1 = min(nx - 1, cx0 + box_radius + 1)
    y0 = max(1, cy0 - box_radius)
    y1 = min(ny - 1, cy0 + box_radius + 1)

    search_patch = image[y0:y1, x0:x1]
    if not np.isfinite(search_patch).any():
        return None

    # Use broader annulus/region background around search box if available
    bx0 = max(0, cx0 - box_radius - 3)
    bx1 = min(nx, cx0 + box_radius + 4)
    by0 = max(0, cy0 - box_radius - 3)
    by1 = min(ny, cy0 + box_radius + 4)
    bg, bg_sig = _local_background(image[by0:by1, bx0:bx1])

    net_patch = search_patch - bg
    masked = np.where(np.isfinite(net_patch), net_patch, -np.inf)
    iy_rel, ix_rel = np.unravel_index(int(np.argmax(masked)), masked.shape)
    if not np.isfinite(masked[iy_rel, ix_rel]):
        return None

    px_peak = x0 + int(ix_rel)
    py_peak = y0 + int(iy_rel)

    win_img = image[py_peak - 1 : py_peak + 2, px_peak - 1 : px_peak + 2]
    win_var = variance[py_peak - 1 : py_peak + 2, px_peak - 1 : px_peak + 2]
    good = np.isfinite(win_img)
    if not good.any():
        return None

    net_win = np.where(good, win_img - bg, 0.0)
    weights = np.maximum(net_win, 0.0)
    w_sum = float(weights.sum())
    if w_sum <= 0.0:
        return None

    var_vals = np.where(good & np.isfinite(win_var) & (win_var > 0), win_var, bg_sig**2)
    snr = float(net_win.sum()) / math.sqrt(max(float(var_vals.sum()), 1e-12))
    if snr < 2.5:
        return None

    yy, xx = np.mgrid[py_peak - 1 : py_peak + 2, px_peak - 1 : px_peak + 2]
    x_sub = float((xx * weights).sum() / w_sum)
    y_sub = float((yy * weights).sum() / w_sum)

    ra_v, dec_v = grid.wcs_copy().pixel_to_world_values(x_sub, y_sub)
    if not (np.isfinite(ra_v) and np.isfinite(dec_v)):
        return None

    pixscale_arcsec = grid.pixscale_deg * 3600.0
    sigma_arcsec = min(max((pixscale_arcsec * 0.8) / max(snr, 2.0), 0.05), 5.0)
    return float(ra_v) % 360.0, float(dec_v), x_sub, y_sub, float(sigma_arcsec)


def measure_anomaly_centroids(
    coadds: list[Any],
    grid: DisplayGrid,
    anomaly: DiffAnomaly,
    mjd_a: float | None = None,
    mjd_b: float | None = None,
) -> list[EpochCentroid]:
    """Measure sub-pixel centroids of a moving anomaly across all available `coadds`."""
    if len(coadds) < 2:
        return []

    ordered = sorted(coadds, key=lambda c: float(c.mjd_mid))
    wcs = grid.wcs_copy()
    xb, yb = wcs.world_to_pixel_values(anomaly.ra, anomaly.dec)
    if anomaly.ra_from is not None and anomaly.dec_from is not None:
        xa, ya = wcs.world_to_pixel_values(anomaly.ra_from, anomaly.dec_from)
    else:
        xa, ya = xb, yb

    if not (np.isfinite(xa) and np.isfinite(ya) and np.isfinite(xb) and np.isfinite(yb)):
        return []

    ref_mjd_a = float(mjd_a) if mjd_a is not None else float(ordered[0].mjd_mid)
    ref_mjd_b = float(mjd_b) if mjd_b is not None else float(ordered[-1].mjd_mid)
    dt_span = ref_mjd_b - ref_mjd_a
    shift_px = math.hypot(float(xb) - float(xa), float(yb) - float(ya))
    search_r = min(max(4, int(math.ceil(shift_px)) + 2), 6)
    mid_x = 0.5 * (float(xa) + float(xb))
    mid_y = 0.5 * (float(ya) + float(yb))

    centroids: list[EpochCentroid] = []
    for c in ordered:
        frac = ((float(c.mjd_mid) - ref_mjd_a) / dt_span) if abs(dt_span) > 1e-3 else 0.5
        approx_x = float(xa) + frac * (float(xb) - float(xa))
        approx_y = float(ya) + frac * (float(yb) - float(ya))

        res = measure_subpixel_centroid(c.image, c.variance, grid, approx_x, approx_y, box_radius=search_r)
        if res is None and (frac < -0.05 or frac > 1.05):
            # Returning parallax trajectory may revisit the Epoch A / midpoint region
            res = measure_subpixel_centroid(c.image, c.variance, grid, mid_x, mid_y, box_radius=search_r)
        if res is None:
            continue
        ra_i, dec_i, x_sub, y_sub, sig = res
        centroids.append(
            EpochCentroid(
                epoch_index=int(c.index),
                mjd=float(c.mjd_mid),
                ra=ra_i,
                dec=dec_i,
                x_px=x_sub,
                y_px=y_sub,
                sigma_arcsec=sig,
            )
        )
    return centroids


def _parallax_factors(ra_deg: float, dec_deg: float, mjds: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute annual parallax projection factors `(P_alpha*, P_delta)` at `mjds`."""
    ra_rad = math.radians(ra_deg)
    dec_rad = math.radians(dec_deg)
    e_alpha = np.array([-math.sin(ra_rad), math.cos(ra_rad), 0.0], dtype=np.float64)
    e_delta = np.array(
        [
            -math.cos(ra_rad) * math.sin(dec_rad),
            -math.sin(ra_rad) * math.sin(dec_rad),
            math.cos(dec_rad),
        ],
        dtype=np.float64,
    )

    sun = get_sun(Time(mjds, format="mjd"))
    sun_cart = sun.cartesian.xyz.to_value("AU")  # (3, N)
    p_alpha = np.einsum("j,jn->n", e_alpha, sun_cart)
    p_delta = np.einsum("j,jn->n", e_delta, sun_cart)
    return p_alpha, p_delta


def _fit_1d_linear(t: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[float, float, float]:
    """Weighted linear fit `y = c0 + slope * t`, returning `(c0, slope, chi2)`."""
    w_sum = float(w.sum())
    t_bar = float((w * t).sum() / w_sum)
    y_bar = float((w * y).sum() / w_sum)
    dt = t - t_bar
    denom = float((w * dt * dt).sum())
    slope = float((w * dt * (y - y_bar)).sum() / denom) if denom > 1e-12 else 0.0
    c0 = y_bar - slope * t_bar
    resid = y - (c0 + slope * t)
    chi2 = float((w * resid * resid).sum())
    return c0, slope, chi2


def fit_motion_model(centroids: list[EpochCentroid]) -> tuple[MotionModel, float | None]:
    """Compare linear proper-motion vs annual-parallax models across `>= 3` epochs."""
    if len(centroids) < 3:
        return "ambiguous", None

    ordered = sorted(centroids, key=lambda c: c.mjd)
    ra0, dec0 = ordered[0].ra, ordered[0].dec
    cosd = max(math.cos(math.radians(dec0)), 1e-6)

    mjds = np.array([c.mjd for c in ordered], dtype=np.float64)
    t0 = float(np.mean(mjds))
    dt_yr = (mjds - t0) / 365.25
    if float(np.ptp(dt_yr)) < 0.05:
        return "ambiguous", None

    x = np.array(
        [((c.ra - ra0 + 180.0) % 360.0 - 180.0) * cosd * 3600.0 for c in ordered],
        dtype=np.float64,
    )
    y = np.array([(c.dec - dec0) * 3600.0 for c in ordered], dtype=np.float64)
    sig = np.array([max(c.sigma_arcsec, 0.05) for c in ordered], dtype=np.float64)
    w = 1.0 / (sig**2)

    # Reject completely static sources with no significant excursion
    total_span = math.hypot(float(np.ptp(x)), float(np.ptp(y)))
    if total_span < 1.5 * float(np.mean(sig)):
        return "ambiguous", None

    # Model 1: Linear proper motion (4 parameters: x0, mu_ra, y0, mu_dec)
    _, _, chi2_x = _fit_1d_linear(dt_yr, x, w)
    _, _, chi2_y = _fit_1d_linear(dt_yr, y, w)
    chi2_lin = chi2_x + chi2_y

    # Model 2: Annual parallax (3 parameters: x0, y0, parallax >= 0)
    p_alpha, p_delta = _parallax_factors(ra0, dec0, mjds)
    w_sum = float(w.sum())
    pa_bar = float((w * p_alpha).sum() / w_sum)
    pd_bar = float((w * p_delta).sum() / w_sum)
    x_bar = float((w * x).sum() / w_sum)
    y_bar = float((w * y).sum() / w_sum)

    dpa = p_alpha - pa_bar
    dpd = p_delta - pd_bar
    dx = x - x_bar
    dy = y - y_bar

    plx_den = float((w * (dpa * dpa + dpd * dpd)).sum())
    if plx_den <= 1e-9:
        return "linear_proper_motion", None

    plx_est = max(0.0, float((w * (dpa * dx + dpd * dy)).sum() / plx_den))
    resid_px = dx - plx_est * dpa
    resid_py = dy - plx_est * dpd
    chi2_plx = float((w * (resid_px**2 + resid_py**2)).sum())

    n_meas = 2 * len(ordered)
    bic_lin = chi2_lin + 4.0 * math.log(n_meas)
    bic_plx = chi2_plx + 3.0 * math.log(n_meas)

    # For N >= 4 epochs, also test combined slow orbital drift + annual parallax (5 parameters)
    if len(ordered) >= 4:
        # Build design matrix A * [x0, mu_x, y0, mu_y, plx] = [x, y]
        n_ep = len(ordered)
        a_mat = np.zeros((2 * n_ep, 5), dtype=np.float64)
        b_vec = np.concatenate([x, y])
        w_sqrt = np.sqrt(np.concatenate([w, w]))
        a_mat[:n_ep, 0] = 1.0
        a_mat[:n_ep, 1] = dt_yr
        a_mat[:n_ep, 4] = p_alpha
        a_mat[n_ep:, 2] = 1.0
        a_mat[n_ep:, 3] = dt_yr
        a_mat[n_ep:, 4] = p_delta

        params, _, _, _ = np.linalg.lstsq(a_mat * w_sqrt[:, None], b_vec * w_sqrt, rcond=None)
        plx_comb = max(0.0, float(params[4]))
        if plx_comb > 0.0:
            resid_comb = (b_vec - a_mat @ params) * w_sqrt
            chi2_comb = float((resid_comb**2).sum())
            bic_comb = chi2_comb + 5.0 * math.log(n_meas)
            if bic_comb < bic_plx and plx_comb > 2.0 * float(np.mean(sig)):
                bic_plx = bic_comb
                plx_est = plx_comb

    if abs(bic_lin - bic_plx) < 2.0:
        return "ambiguous", None
    if bic_plx < bic_lin and plx_est > 0.0:
        return "parallax_dominated", round(float(plx_est), 3)
    return "linear_proper_motion", None
