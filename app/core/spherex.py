"""SPHEREx instrument constants.

SPHEREx has 6 H2RG detector arrays behind Linear Variable Filters (LVF). The
passband of an LVF changes continuously along one detector axis, so *every
pixel row has its own central wavelength*. Each array is split into 17
spectral channels -> 6 x 17 = 102 "bands". A single L2 image therefore holds
ONE wavelength sample per sky position; a full spectrum of a point is assembled
from many exposures that put the point on different parts of the LVF.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

CHANNELS_PER_DETECTOR = 17
N_DETECTORS = 6
N_CHANNELS = CHANNELS_PER_DETECTOR * N_DETECTORS  # 102


@dataclass(frozen=True)
class Detector:
    index: int  # 1..6
    lambda_min_um: float
    lambda_max_um: float
    resolving_power: float


# Nominal band edges (SPHEREx Explanatory Supplement, rounded).
DETECTORS: dict[int, Detector] = {
    1: Detector(1, 0.75, 1.09, 39),
    2: Detector(2, 1.10, 1.62, 41),
    3: Detector(3, 1.63, 2.41, 41),
    4: Detector(4, 2.42, 3.82, 35),
    5: Detector(5, 3.83, 4.41, 112),
    6: Detector(6, 4.42, 5.00, 130),
}

BAND_GROUPS: dict[str, tuple[int, ...]] = {
    "all": (1, 2, 3, 4, 5, 6),
    "blue": (1, 2, 3),  # 0.75-2.4 um (HgCdTe 2.5um cutoff arrays)
    "red": (4, 5, 6),  # 2.4-5.0 um: brown dwarfs / "W2-like", CH4 3.3um, CO 4.6um
    **{f"D{i}": (i,) for i in range(1, 7)},
}


def channel_edges_um() -> np.ndarray:
    """(102, 2) array of channel [lo, hi] edges.

    Within a detector the resolving power is ~constant, so channels are
    log-uniformly spaced between the band edges.
    """
    edges = []
    for d in DETECTORS.values():
        e = np.geomspace(d.lambda_min_um, d.lambda_max_um, CHANNELS_PER_DETECTOR + 1)
        edges.extend(zip(e[:-1], e[1:]))
    return np.asarray(edges)


CHANNEL_EDGES_UM = channel_edges_um()
CHANNEL_CENTERS_UM = np.sqrt(CHANNEL_EDGES_UM[:, 0] * CHANNEL_EDGES_UM[:, 1])


def channel_index(wavelength_um: float) -> int | None:
    """Map a wavelength to a global channel index 0..101 (None if outside)."""
    hit = np.nonzero((CHANNEL_EDGES_UM[:, 0] <= wavelength_um) & (wavelength_um < CHANNEL_EDGES_UM[:, 1]))[0]
    return int(hit[0]) if hit.size else None


def detector_from_bandname(name: str | None) -> int | None:
    """'SPHEREx-D3' -> 3."""
    if not name:
        return None
    name = str(name).strip()
    for i in range(1, 7):
        if name.endswith(f"D{i}"):
            return i
    return None


# ---------------------------------------------------------------------------
# FLAGS bitmask. The L2 FLAGS header self-documents every bit as
#   MSKN0021 = 'SOURCE'   MSKM0021 = 2097152   MSKD0021 = 'Pixel mapped to ...'
# (verified on DR1 cutouts). DEFAULT_FLAG_BITS is only a fallback.
# IMPORTANT: SOURCE (bit 21) means "pixel belongs to a known source" - it is
# informational, NOT a defect. Masking it would delete every star.
# ---------------------------------------------------------------------------
DEFAULT_FLAG_BITS: dict[str, int] = {
    "TRANSIENT": 0,
    "OVERFLOW": 1,
    "SUR_ERROR": 2,
    "NONFUNC": 6,
    "DICHROIC": 7,
    "MISSING_DATA": 9,
    "HOT": 10,
    "COLD": 11,
    "FULLSAMPLE": 12,
    "PHANMISS": 14,
    "NONLINEAR": 15,
    "PERSIST": 17,
    "OUTLIER": 19,
    "CROSSTALK": 20,
    "SOURCE": 21,
    "GHOST": 22,
    "GHOST_FPA": 23,
    "GHOST_EXT": 24,
    "STREAK": 25,
    "BLOOM": 26,
    "SNOWBALL": 27,
    "HALO": 28,
}

# Pixels with these bits are excluded from co-adds and photometry.
BAD_FLAG_NAMES: frozenset[str] = frozenset(
    {"TRANSIENT", "OVERFLOW", "SUR_ERROR", "NONFUNC", "MISSING_DATA", "HOT", "COLD", "PERSIST", "OUTLIER",
     "STREAK", "SNOWBALL"}
)
# Shown in the FLAGS overlay but kept in the data (they describe, not destroy).
INFO_FLAG_NAMES: frozenset[str] = frozenset(
    {"SOURCE", "NONLINEAR", "DICHROIC", "FULLSAMPLE", "PHANMISS", "CROSSTALK", "GHOST", "GHOST_FPA", "GHOST_EXT",
     "BLOOM", "HALO"}
)


def parse_flag_bits(header) -> dict[str, int]:
    """Read MSKNnnnn bit definitions from a FLAGS header (fallback to defaults)."""
    bits = {}
    for k, v in header.items():
        if k.startswith("MSKN") and k[4:].isdigit():
            bits[str(v).strip().upper()] = int(k[4:])
    return bits or dict(DEFAULT_FLAG_BITS)


def bad_flag_mask(bits: dict[str, int]) -> int:
    m = 0
    for name, bit in bits.items():
        if name in BAD_FLAG_NAMES:
            m |= 1 << bit
    return m


def source_flag_mask(bits: dict[str, int]) -> int:
    return (1 << bits["SOURCE"]) if "SOURCE" in bits else 0
