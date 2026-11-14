"""Canonical cache keys for targets.

Coordinates arrive as arbitrary floats from the UI (e.g. 269.4520751 vs
269.45207). We snap them to a 0.0001 deg (~0.36") grid so tiny float jitter maps
onto the same cache row and the same Storage folder. A tolerance search
(find_target RPC) additionally catches near-duplicates across the grid boundary.
"""
from __future__ import annotations

import re

COORD_DECIMALS = 4
FOV_DECIMALS = 3


def normalize_ra(ra: float) -> float:
    return round(ra % 360.0, COORD_DECIMALS) % 360.0


def normalize_dec(dec: float) -> float:
    if not -90.0 <= dec <= 90.0:
        raise ValueError("dec must be within [-90, 90]")
    return round(dec, COORD_DECIMALS)


def normalize_fov(fov: float) -> float:
    return round(fov, FOV_DECIMALS)


def target_key(ra: float, dec: float, fov: float, band: str = "all") -> str:
    """'269.4521_+4.6934_0.100_all' - safe for Storage paths and URLs."""
    ra_n, dec_n, fov_n = normalize_ra(ra), normalize_dec(dec), normalize_fov(fov)
    key = f"{ra_n:.{COORD_DECIMALS}f}_{dec_n:+.{COORD_DECIMALS}f}_{fov_n:.{FOV_DECIMALS}f}_{band}"
    assert re.fullmatch(r"[0-9A-Za-z_.+-]+", key)
    return key


def parse_sexagesimal(ra: str, dec: str) -> tuple[float, float]:
    """Accept '08h55m10.831s' / '-7d14m42.53s' or '08:55:10.8' style input."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame="icrs")
    return float(c.ra.deg), float(c.dec.deg)
