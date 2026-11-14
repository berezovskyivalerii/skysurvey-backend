"""Pydantic v2 schemas for the 102-channel SPHEREx spectral chart contract
(ready for direct rendering in Recharts / Chart.js)."""
from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

SpectralHint = Literal[
    "ultra_cool_y_dwarf_or_planet_x_candidate",
    "cold_brown_dwarf_candidate",
    "red_stellar_or_dusty_source",
    "blue_main_sequence_star",
    "no_significant_detection",
    "insufficient_wavelength_coverage",
]


class _Base(BaseModel):
    model_config = ConfigDict(extra="ignore")


class TargetCoords(_Base):
    ra: float = Field(description="Right Ascension (ICRS degrees)")
    dec: float = Field(description="Declination (ICRS degrees)")


class SpectralModelPoint(_Base):
    wavelength_um: float = Field(description="Wavelength in micrometres")
    model_flux: float = Field(description="Fitted model flux density in uJy")


class SpectrumSummary(_Base):
    total_channels: int = Field(102, description="Total SPHEREx LVF channels across Bands 1-6 (0.75-5.0 um)")
    valid_points: int = Field(description="Number of unflagged spectral channels with valid measurements")
    flagged_points: int = Field(0, description="Number of channels where only flagged pixels were available")
    n_exposures_used: int = Field(0, description="Number of individual SPHEREx FITS cutouts photometered")
    mean_snr: float = Field(description="Mean signal-to-noise ratio across valid channels")
    max_snr: float = Field(0.0, description="Maximum channel S/N")
    peak_wavelength_um: float | None = Field(None, description="Wavelength (um) of peak measured flux")
    peak_flux_ujy: float | None = Field(None, description="Peak flux density (uJy)")
    ir_color_index: float | None = Field(
        None,
        description="Long-to-short IR flux ratio F(3.0-5.0 um) / F(0.75-1.5 um). High (>2.0) indicates a cold brown dwarf / Planet X candidate.",
    )
    spectral_hint: SpectralHint = Field(
        description="Automated astrophysical classification hint derived from IR slope and S/N"
    )
    flux_unit: str = Field("uJy", description="Flux density unit (microjansky)")
    estimated_temp_k: float | None = Field(
        None,
        description="Estimated effective temperature T_eff (K) in [200, 6000]",
    )
    ch4_absorption_depth: float | None = Field(
        None,
        description="Relative CH4 absorption depth around 3.2-3.4 um (0.0..1.0)",
    )
    h2o_absorption_depth: float | None = Field(
        None,
        description="Relative H2O absorption depth around 2.55-2.85 um (0.0..1.0)",
    )
    model_curve: list[SpectralModelPoint] = Field(
        default_factory=list,
        description="Smooth fitted Planck + molecular absorption curve",
    )


class ChartPoint(_Base):
    channel_id: int = Field(ge=1, le=102, description="1-based SPHEREx channel index (1..102)")
    wavelength_um: float = Field(description="Central wavelength in micrometres (0.75..5.0)")
    bandwidth_um: float = Field(description="Effective LVF filter bandwidth (FWHM) in micrometres")
    band: str = Field(description="SPHEREx detector band label ('Band 1' .. 'Band 6')")
    detector: int = Field(ge=1, le=6, description="SPHEREx H2RG detector array number (1..6)")
    flux: float = Field(description="Background- & ZODI-subtracted aperture flux in microjansky (uJy)")
    flux_error: float = Field(description="1-sigma flux uncertainty in microjansky (uJy)")
    snr: float = Field(description="Signal-to-noise ratio (flux / flux_error)")
    is_flagged: bool = Field(
        description="True if central aperture pixels carried detector/optical quality flags (TRANSIENT, OVERFLOW, DICHROIC, etc.)"
    )
    n_samples: int = Field(1, description="Number of exposures co-added into this spectral channel")


class AbsorptionBandRef(_Base):
    molecule: str = Field(description="Chemical formula (H2O, CH4, CO, CO2)")
    range_um: list[float] = Field(min_length=2, max_length=2, description="[lambda_min_um, lambda_max_um]")
    label: str = Field(default="", description="Human-readable band label for chart tooltip")


DEFAULT_ABSORPTION_BANDS: list[AbsorptionBandRef] = [
    AbsorptionBandRef(molecule="H2O", range_um=[1.10, 1.20], label="H2O (J-band)"),
    AbsorptionBandRef(molecule="H2O", range_um=[1.35, 1.45], label="H2O (1.4 um)"),
    AbsorptionBandRef(molecule="CH4", range_um=[1.60, 1.70], label="CH4 (H-band T/Y dwarf)"),
    AbsorptionBandRef(molecule="H2O", range_um=[1.80, 1.95], label="H2O (1.9 um)"),
    AbsorptionBandRef(molecule="CH4", range_um=[2.15, 2.35], label="CH4 (K-band)"),
    AbsorptionBandRef(molecule="H2O", range_um=[2.55, 2.85], label="H2O (2.7 um fundamental)"),
    AbsorptionBandRef(molecule="CH4", range_um=[3.20, 3.40], label="CH4 (3.3 um fundamental)"),
    AbsorptionBandRef(molecule="CO2", range_um=[4.20, 4.35], label="CO2 (4.26 um)"),
    AbsorptionBandRef(molecule="CO", range_um=[4.55, 4.80], label="CO (4.67 um fundamental)"),
]


class SpectrumResponse(_Base):
    target: TargetCoords
    spectrum_summary: SpectrumSummary
    chart_data: list[ChartPoint]
    absorption_bands_reference: list[AbsorptionBandRef] = Field(
        default_factory=lambda: list(DEFAULT_ABSORPTION_BANDS)
    )
