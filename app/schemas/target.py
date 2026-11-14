"""Pydantic schemas = the API contract shared with the frontend team.

Conventions
-----------
* Coordinates: ICRS degrees. Wavelengths: micrometres. Image units: MJy/sr.
* Spectrum flux: microjansky (uJy).
* Gaia source_id is a STRING (int64 > 2**53 loses precision in JS numbers).
* Pixel coordinates (x, y) refer to the delivered PNG: origin top-left,
  x -> right (East is LEFT, North is UP, standard sky orientation).
"""
from __future__ import annotations

from datetime import datetime
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

from app.schemas.spectrum import SpectrumResponse

TargetStatus = Literal["processing", "partial", "ready", "failed", "no_data"]
Classification = Literal["known_star", "candidate", "no_detection", "unknown"]
StretchType = Literal["linear", "sqrt", "log", "asinh", "zscale"]
BandName = Literal["all", "blue", "red", "D1", "D2", "D3", "D4", "D5", "D6"]
LayerName = Literal["image", "variance", "flags"]


class _Model(BaseModel):
    model_config = ConfigDict(extra="ignore")


# ----------------------------------------------------------------- display ---
class StretchInfo(_Model):
    stretch_type: StretchType = "asinh"
    vmin: float = Field(description="Lower cut applied to ALL epochs (MJy/sr)")
    vmax: float = Field(description="Upper cut applied to ALL epochs (MJy/sr)")
    unit: str = "MJy/sr"


class DataRange(_Model):
    """Robust statistics of the epoch stack - use to scale restretch sliders."""

    p0_5: float
    p50: float
    p99_5: float
    min: float
    max: float


class DisplayInfo(_Model):
    width: int = Field(512, description="PNG width (px)")
    height: int = Field(512, description="PNG height (px)")
    grid_px: int = Field(description="Native resampled grid size (SPHEREx pixels across FOV)")
    pixel_scale_arcsec: float = Field(description="Arcsec per PNG pixel")
    wcs: dict[str, float | str] = Field(
        description="FITS-style TAN WCS of the PNG (1-based CRPIX, y flipped => CDELT2<0 in PNG rows)"
    )
    stretch: StretchInfo
    data_range: DataRange


# ------------------------------------------------------------------ epochs ---
class LayerUrls(_Model):
    image: str = Field(description="Science image, ZScale/Asinh stretched, 8-bit grayscale PNG")
    variance: str = Field(description="sqrt(VARIANCE) noise map, magma colormap PNG")
    flags: str = Field(description="RGBA overlay PNG: transparent=clean, red=bad pixel, green=SOURCE")
    rgb: str | None = Field(
        None,
        description="WiseView-style false-color IR RGB PNG (Blue=D1+D2 0.75-1.6um, Green=D3+D4 1.6-3.8um, Red=D5+D6 3.8-5.0um)",
    )


class EpochOut(_Model):
    index: int
    label: str = Field(description="Human label, e.g. '2025-10'")
    date_obs: datetime = Field(description="Mid-epoch time (UTC)")
    mjd_start: float
    mjd_end: float
    mjd_mid: float
    n_frames: int = Field(description="Number of SPHEREx exposures co-added into this epoch")
    detectors: list[int]
    layers: LayerUrls
    coadd_fits_url: str | None = Field(None, description="Small float32 FITS of the co-add (for power users)")


# -------------------------------------------------------------- crossmatch ---
class GaiaMatch(_Model):
    source_id: str
    ra: float = Field(description="Gaia DR3 position at J2016.0")
    dec: float
    ra_obs: float = Field(description="Position propagated to the SPHEREx epoch using proper motion")
    dec_obs: float
    separation_arcsec: float = Field(description="Distance of propagated position from the query point")
    phot_g_mean_mag: float | None = None
    bp_rp: float | None = None
    parallax_mas: float | None = None
    parallax_error_mas: float | None = None
    pmra_mas_yr: float | None = None
    pmdec_mas_yr: float | None = None
    pm_total_mas_yr: float | None = None
    distance_pc: float | None = None
    ruwe: float | None = None


class OverlaySource(_Model):
    source_id: str
    ra: float = Field(description="Propagated to mean SPHEREx epoch")
    dec: float
    g_mag: float | None = None
    pm_total_mas_yr: float | None = None
    x: float = Field(description="PNG pixel x (0..width)")
    y: float = Field(description="PNG pixel y (0..height, top-left origin)")


class CrossMatch(_Model):
    classification: Classification
    reason: str
    catalog: str = "Gaia DR3"
    match: GaiaMatch | None = None
    spherex_snr: float | None = Field(None, description="IR detection S/N at target in the stacked co-add")
    high_proper_motion: bool = False
    overlay: list[OverlaySource] = []


# ---------------------------------------------------------------- spectrum ---
class SpectrumPoint(_Model):
    channel: int = Field(description="0..101")
    detector: int = Field(description="1..6")
    wavelength_um: float
    bandwidth_um: float
    flux: float | None = Field(None, description="uJy, null when no clean sample in this channel")
    error: float | None = None
    n_samples: int = 0


class Spectrum(_Model):
    flux_unit: str = "uJy"
    method: str = "3x3 aperture, ZODI- and annulus-background subtracted, inverse-variance binned"
    complete: bool = Field(description="False while only display frames were used (status=partial)")
    n_frames_used: int
    n_channels_filled: int
    points: list[SpectrumPoint]


# ------------------------------------------------------ difference imaging ---
AnomalyType = Literal["moving_dipole", "variable_brightening", "variable_dimming"]


class DiffAnomaly(_Model):
    anomaly_type: AnomalyType = Field(
        description="'moving_dipole' (proper-motion shift between epochs) or 'variable_brightening'/'variable_dimming'"
    )
    snr: float = Field(description="Combined significance (sigma) of the dipole or monopole residual")
    ra: float = Field(description="ICRS Right Ascension (deg) at Epoch B (current position)")
    dec: float = Field(description="ICRS Declination (deg) at Epoch B (current position)")
    x: float = Field(description="PNG pixel X (0..width) at Epoch B")
    y: float = Field(description="PNG pixel Y (0..height, top-left origin) at Epoch B")
    ra_from: float | None = Field(None, description="ICRS RA (deg) at Epoch A (origin of moving dipole)")
    dec_from: float | None = Field(None, description="ICRS Dec (deg) at Epoch A (origin of moving dipole)")
    x_from: float | None = Field(None, description="PNG pixel X at Epoch A (draw motion arrow from x_from,y_from -> x,y)")
    y_from: float | None = Field(None, description="PNG pixel Y at Epoch A")
    shift_arcsec: float | None = Field(None, description="Angular displacement between Epoch A and Epoch B (arcsec)")
    pm_arcsec_yr: float | None = Field(None, description="Estimated annualized proper motion (arcsec/yr)")
    delta_flux_ujy: float = Field(description="Peak aperture flux change between Epoch A and Epoch B (uJy)")


class DifferenceResult(_Model):
    epoch_a: int = Field(description="Index of reference (early) epoch")
    epoch_b: int = Field(description="Index of comparison (later) epoch")
    mjd_a: float
    mjd_b: float
    baseline_days: float = Field(description="Time baseline (MJD_B - MJD_A) in days")
    diff_image_url: str = Field(
        description="Bipolar difference PNG (Blue = Epoch A position / dimmed, Black = static background, Red/Gold = Epoch B position / brightened)"
    )
    deep_stack_url: str | None = Field(
        None,
        description="Inverse-variance deep stack PNG across all epochs for maximum faint-source sensitivity",
    )
    max_dipole_snr: float = Field(0.0, description="Highest moving-dipole S/N detected in the field")
    anomalies: list[DiffAnomaly] = Field(
        default_factory=list,
        description="Sorted list of detected moving dipoles and variable sources in the FOV",
    )


# ------------------------------------------------------------- responses -----
class Progress(_Model):
    stage: str
    percent: int = Field(ge=0, le=100)


class CacheInfo(_Model):
    hit: bool
    stale: bool = False
    created_at: datetime | None = None
    updated_at: datetime | None = None


class TargetResponse(_Model):
    target_key: str
    status: TargetStatus
    ra: float
    dec: float
    fov: float
    band: BandName = "all"
    progress: Progress | None = None
    poll_after_ms: int | None = Field(None, description="Present while status=processing|partial")
    cache: CacheInfo
    pipeline_version: str
    n_observations: int = 0
    display: DisplayInfo | None = None
    epochs: list[EpochOut] = []
    crossmatch: CrossMatch | None = None
    spectrum: Spectrum | None = None
    spectrum_data: SpectrumResponse | None = Field(
        None,
        description="Ready-to-plot 102-channel spectral chart data, summary metrics, and molecular absorption bands",
    )
    difference: DifferenceResult | None = Field(
        None,
        description="Epoch difference imaging (Epoch_last - Epoch_0), deep stack URL, and auto-detected moving dipoles / variable stars",
    )
    error: str | None = None


class RestretchRequest(_Model):
    model_config = ConfigDict(
        extra="ignore",
        json_schema_extra={
            "examples": [
                {
                    "ra": 269.4485,
                    "dec": 4.7394,
                    "fov": 0.1,
                    "band": "all",
                    "stretch_type": "asinh",
                    "min_cut": -0.02,
                    "max_cut": 0.35,
                    "asinh_a": 0.08,
                    "log_a": 1000.0,
                    "epochs": None,
                }
            ]
        },
    )

    ra: float = Field(269.4485, ge=0, lt=360)
    dec: float = Field(4.7394, ge=-90, le=90)
    fov: float = Field(0.1, gt=0, le=1.0)
    band: BandName = "all"
    stretch_type: StretchType = "asinh"
    min_cut: float | None = Field(None, description="MJy/sr. Omit (or null) to use automatic ZScale lower limit")
    max_cut: float | None = Field(None, description="MJy/sr. Omit (or null) to use automatic ZScale upper limit")
    asinh_a: float = Field(0.1, gt=0, le=1, description="AsinhStretch softening (smaller => stronger)")
    log_a: float = Field(1000.0, gt=1, description="LogStretch exponent")
    epochs: list[int] | None = Field(None, description="Subset of epoch indexes; null = all epochs")

    @field_validator("max_cut")
    @classmethod
    def _cuts_ordered(cls, v, info):
        lo = info.data.get("min_cut")
        if v is not None and lo is not None:
            # If Swagger UI sent placeholder 0 and 0, treat both as None (auto ZScale)
            if lo == 0 and v == 0:
                info.data["min_cut"] = None
                return None
            if v <= lo:
                raise ValueError("max_cut must be > min_cut")
        return v


class RestretchEpoch(_Model):
    index: int
    image: str = Field(description="data:image/png;base64,...")


class RestretchResponse(_Model):
    target_key: str
    stretch: StretchInfo
    width: int
    height: int
    epochs: list[RestretchEpoch]
    source: Literal["local_cache", "storage", "rebuilt"]
    elapsed_ms: float


class ErrorResponse(_Model):
    detail: str
