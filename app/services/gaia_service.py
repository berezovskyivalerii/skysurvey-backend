"""Gaia DR3 cross-match and target classification (pure async via httpx).

Key subtlety: Gaia DR3 positions are at epoch J2016.0 while SPHEREx observes in
2025+. Barnard's Star moves 10.4"/yr -> ~100" offset, i.e. 16 SPHEREx pixels.
A naive 5" cone search would MISS it and wrongly call it a 'candidate'.
So we:
  1. query ESA Gaia TAP (`/tap-server/tap/sync`, FORMAT=json) directly via our
     shared `httpx.AsyncClient` with a strict timeout (avoiding `astroquery.gaia`
     import-time synchronous blocking).
  2. propagate every source to the SPHEREx epoch (linear in pmra*/pmdec).
  3. match within r_match (~1.6 SPHEREx px) of the requested position.
"""
from __future__ import annotations

import asyncio
import logging
import math
import time
from dataclasses import dataclass
from datetime import datetime, timezone

import httpx
import numpy as np

from app.config import Settings
from app.schemas.target import CrossMatch, GaiaMatch, OverlaySource

log = logging.getLogger(__name__)

GAIA_EPOCH_MJD = 57388.5  # J2016.0
MAS_PER_DEG = 3.6e6
HIGH_PM_MAS_YR = 200.0
ESA_GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"


@dataclass
class GaiaResult:
    match_candidates: list[dict]  # sources that may coincide with target
    field: list[dict]  # brightest sources in the FOV (for the overlay)
    error: str | None = None


def _mjd_now() -> float:
    return datetime.now(timezone.utc).timestamp() / 86400.0 + 40587.0


def propagate(src: dict, mjd: float) -> tuple[float, float]:
    """Linear proper-motion propagation (parallax term < 1" -> ignored)."""
    pmra, pmdec = src.get("pmra"), src.get("pmdec")
    if pmra is None or pmdec is None or not np.isfinite(pmra) or not np.isfinite(pmdec):
        return src["ra"], src["dec"]
    dt = (mjd - GAIA_EPOCH_MJD) / 365.25
    cosd = max(math.cos(math.radians(src["dec"])), 1e-6)
    ra = (src["ra"] + pmra * dt / MAS_PER_DEG / cosd) % 360.0  # pmra is mu_alpha* (already * cos dec)
    dec = src["dec"] + pmdec * dt / MAS_PER_DEG
    return ra, dec


def ang_sep_arcsec(ra1, dec1, ra2, dec2) -> float:
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    s = math.sin((d2 - d1) / 2) ** 2 + math.cos(d1) * math.cos(d2) * math.sin((r2 - r1) / 2) ** 2
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(s)))) * 3600.0


def _clean(v):
    if v is None or (np.ma.is_masked(v)):
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if np.isfinite(f) else None


class GaiaService:
    COLS = "source_id, ra, dec, parallax, parallax_error, pmra, pmdec, pm, phot_g_mean_mag, bp_rp, ruwe"

    def __init__(self, settings: Settings, http: httpx.AsyncClient | None = None):
        self.s = settings
        self.http = http

    async def _query_esa_json(self, adql: str) -> list[dict]:
        params = {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "json",
            "QUERY": " ".join(adql.split()),
        }
        if self.http is not None:
            r = await self.http.post(ESA_GAIA_TAP_SYNC, data=params, timeout=12.0)
        else:
            async with httpx.AsyncClient(timeout=12.0) as client:
                r = await client.post(ESA_GAIA_TAP_SYNC, data=params)
        r.raise_for_status()
        payload = r.json()
        cols = [m["name"].lower() for m in payload.get("metadata", [])]
        rows = []
        for raw_row in payload.get("data", []):
            d = {}
            for col, val in zip(cols, raw_row):
                d[col] = str(val) if col == "source_id" and val is not None else _clean(val)
            rows.append(d)
        return rows

    def _run_vizier(self, ra: float, dec: float, radius_deg: float, limit: int) -> list[dict]:
        """Fallback when the ESA archive is slow/down: CDS VizieR I/355/gaiadr3."""
        from astropy import units as u
        from astropy.coordinates import SkyCoord
        from astroquery.vizier import Vizier

        v = Vizier(
            columns=["Source", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx", "pmRA", "pmDE", "PM", "Gmag", "BP-RP", "RUWE"],
            row_limit=limit,
        )
        res = v.query_region(SkyCoord(ra, dec, unit="deg"), radius=radius_deg * u.deg, catalog="I/355/gaiadr3")
        if not res:
            return []
        out = []
        for r in res[0]:
            out.append(
                {
                    "source_id": str(r["Source"]),
                    "ra": _clean(r["RA_ICRS"]),
                    "dec": _clean(r["DE_ICRS"]),
                    "parallax": _clean(r["Plx"]),
                    "parallax_error": _clean(r["e_Plx"]),
                    "pmra": _clean(r["pmRA"]),
                    "pmdec": _clean(r["pmDE"]),
                    "pm": _clean(r["PM"]),
                    "phot_g_mean_mag": _clean(r["Gmag"]),
                    "bp_rp": _clean(r["BP-RP"]),
                    "ruwe": _clean(r["RUWE"]),
                }
            )
        return out

    async def query(self, ra: float, dec: float, fov: float) -> GaiaResult:
        t0 = time.perf_counter()
        r_pad = self.s.gaia_search_pad_arcsec / 3600.0
        r_field = max(fov * math.sqrt(2) / 2, r_pad)
        adql = f"""
            SELECT TOP {self.s.gaia_overlay_max_sources} {self.COLS}
            FROM {self.s.gaia_table}
            WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.8f}, {dec:.8f}, {r_field:.8f}))
            ORDER BY phot_g_mean_mag ASC
        """
        try:
            rows = await self._query_esa_json(adql)
            log.info("Gaia DR3 returned %d sources in %.2fs", len(rows), time.perf_counter() - t0)
            return GaiaResult(rows, rows)
        except Exception as e:
            log.warning("ESA Gaia TAP failed (%s); falling back to VizieR", e)
            try:
                rows = await asyncio.wait_for(
                    asyncio.to_thread(self._run_vizier, ra, dec, r_field, self.s.gaia_overlay_max_sources),
                    timeout=10.0,
                )
                return GaiaResult(rows, rows)
            except Exception as e2:
                log.warning("Gaia VizieR fallback failed: %s", e2)
                return GaiaResult([], [], error=f"Gaia unavailable: {e2}")

    # ------------------------------------------------------- crossmatch ----
    def best_match(self, res: GaiaResult, ra: float, dec: float, mjd: float) -> GaiaMatch | None:
        best, best_sep = None, None
        for s in res.match_candidates:
            if s.get("ra") is None:
                continue
            ra_o, dec_o = propagate(s, mjd)
            sep = ang_sep_arcsec(ra_o, dec_o, ra, dec)
            if sep <= self.s.gaia_match_radius_arcsec and (best_sep is None or sep < best_sep):
                best, best_sep = (s, ra_o, dec_o), sep
        if best is None:
            return None
        s, ra_o, dec_o = best
        plx = s.get("parallax")
        pm_tot = s.get("pm") or (
            math.hypot(s["pmra"], s["pmdec"]) if s.get("pmra") is not None and s.get("pmdec") is not None else None
        )
        return GaiaMatch(
            source_id=s["source_id"],
            ra=s["ra"],
            dec=s["dec"],
            ra_obs=ra_o,
            dec_obs=dec_o,
            separation_arcsec=round(best_sep, 3),
            phot_g_mean_mag=s.get("phot_g_mean_mag"),
            bp_rp=s.get("bp_rp"),
            parallax_mas=plx,
            parallax_error_mas=s.get("parallax_error"),
            pmra_mas_yr=s.get("pmra"),
            pmdec_mas_yr=s.get("pmdec"),
            pm_total_mas_yr=pm_tot,
            distance_pc=(1000.0 / plx) if plx and plx > 0 else None,
            ruwe=s.get("ruwe"),
        )

    def overlay(self, res: GaiaResult, mjd: float, world_to_png) -> list[OverlaySource]:
        out = []
        for s in res.field:
            if s.get("ra") is None:
                continue
            ra_o, dec_o = propagate(s, mjd)
            xy = world_to_png(ra_o, dec_o)
            if xy is None:
                continue
            out.append(
                OverlaySource(
                    source_id=s["source_id"],
                    ra=round(ra_o, 7),
                    dec=round(dec_o, 7),
                    g_mag=s.get("phot_g_mean_mag"),
                    pm_total_mas_yr=s.get("pm"),
                    x=round(xy[0], 2),
                    y=round(xy[1], 2),
                )
            )
        return out

    def classify(self, match: GaiaMatch | None, snr: float | None, gaia_error: str | None) -> CrossMatch:
        hp = bool(match and match.pm_total_mas_yr and match.pm_total_mas_yr > HIGH_PM_MAS_YR)
        if match is not None:
            return CrossMatch(
                classification="known_star",
                reason=f"Gaia DR3 {match.source_id} within {match.separation_arcsec:.1f}\" after proper-motion propagation",
                match=match,
                spherex_snr=snr,
                high_proper_motion=hp,
            )
        if gaia_error:
            return CrossMatch(classification="unknown", reason=gaia_error, spherex_snr=snr)
        if snr is None:
            return CrossMatch(classification="unknown", reason="No SPHEREx coverage / photometry at this position")
        if snr >= self.s.detection_snr:
            return CrossMatch(
                classification="candidate",
                reason=(
                    f"IR source (S/N={snr:.1f}) with no Gaia DR3 counterpart. Could be a brown dwarf / cold "
                    "object - or a galaxy/QSO: check motion between epochs and the 3-5um spectrum."
                ),
                spherex_snr=snr,
            )
        return CrossMatch(classification="no_detection", reason=f"No Gaia source and SPHEREx S/N={snr:.1f}", spherex_snr=snr)
