"""CatWISE2020 / AllWISE infrared cross-validation and 10-year proper-motion baseline.

Ultracool Y/T brown dwarfs (T_eff < 400 K, e.g. WISE 0855) and Planet X candidates
are invisible in optical Gaia DR3, but many were imaged by NASA WISE / NEOWISE
in W1 (3.4 um) and W2 (4.6 um) around 2010-2020. Matching CatWISE2020 (~2015.4)
against SPHEREx (2025+) yields a ~10-year proper-motion baseline and W1-W2 color
diagnostics.
"""
from __future__ import annotations

import asyncio
import csv
import io
import logging
import math
import time
from dataclasses import dataclass

import httpx
import numpy as np

from app.config import Settings
from app.schemas.target import WiseSource
from app.services.gaia_service import ang_sep_arcsec

log = logging.getLogger(__name__)

CATWISE_EPOCH_MJD = 57170.0  # ~2015.4
ALLWISE_EPOCH_MJD = 55400.0  # ~2010.5
MAS_PER_DEG = 3.6e6


@dataclass
class WiseResult:
    sources: list[dict]
    catalog: str = "CatWISE2020"
    epoch_mjd: float = CATWISE_EPOCH_MJD
    error: str | None = None


def _clean(val: object) -> float | None:
    if val is None or val == "" or np.ma.is_masked(val):
        return None
    try:
        f = float(val)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def propagate_wise(src: dict, mjd: float, epoch_mjd: float = CATWISE_EPOCH_MJD) -> tuple[float, float]:
    """Linearly propagate WISE (ra, dec) to SPHEREx observation `mjd`."""
    ra0 = float(src["ra"])
    dec0 = float(src["dec"])
    pmra = src.get("pmra")
    pmdec = src.get("pmdec")
    if pmra is None or pmdec is None or not (math.isfinite(pmra) and math.isfinite(pmdec)):
        return ra0, dec0

    scale = 1000.0 if src.get("pm_unit") == "arcsec/yr" else 1.0
    pmra_mas = float(pmra) * scale
    pmdec_mas = float(pmdec) * scale

    dt_yr = (mjd - epoch_mjd) / 365.25
    cosd = max(math.cos(math.radians(dec0)), 1e-6)
    ra_obs = (ra0 + pmra_mas * dt_yr / MAS_PER_DEG / cosd) % 360.0
    dec_obs = max(-90.0, min(90.0, dec0 + pmdec_mas * dt_yr / MAS_PER_DEG))
    return ra_obs, dec_obs


class WiseService:
    def __init__(self, settings: Settings, http: httpx.AsyncClient | None = None):
        self.s = settings
        self.http = http

    async def _tap_csv(self, adql: str, tap_url: str | None = None) -> list[dict[str, str]]:
        base_url = (tap_url or self.s.irsa_tap_url).rstrip("/")
        url = f"{base_url}/sync"
        params = {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "csv",
            "QUERY": " ".join(adql.split()),
        }
        if self.http is not None:
            r = await self.http.post(url, data=params, timeout=self.s.wise_timeout_s)
        else:
            async with httpx.AsyncClient(timeout=self.s.wise_timeout_s) as client:
                r = await client.post(url, data=params)
        r.raise_for_status()
        text = r.text
        if text.lstrip().startswith("<"):
            raise RuntimeError(f"TAP returned XML error: {text[:200]}")
        return list(csv.DictReader(io.StringIO(text)))

    @staticmethod
    def _parse_catwise_rows(raw_rows: list[dict[str, str]]) -> list[dict]:
        sources: list[dict] = []
        for row in raw_rows:
            sid = (row.get("source_id") or row.get("objID") or "").strip()
            ra = _clean(row.get("ra") or row.get("RAdeg"))
            dec = _clean(row.get("dec") or row.get("DEdeg"))
            if not sid or ra is None or dec is None:
                continue
            pmra = _clean(row.get("pmra") or row.get("pmRA"))
            pmdec = _clean(row.get("pmdec") or row.get("pmDE"))
            # IRSA catwise_2020 and VizieR II/365/catwise store proper motion in arcsec/yr; convert to mas/yr
            if pmra is not None and abs(pmra) < 25.0:
                pmra *= 1000.0
            if pmdec is not None and abs(pmdec) < 25.0:
                pmdec *= 1000.0
            sources.append(
                {
                    "source_id": sid,
                    "ra": ra,
                    "dec": dec,
                    "pmra": pmra,
                    "pmdec": pmdec,
                    "w1mpro": _clean(row.get("w1mpro") or row.get("w1mpro_pm") or row.get("W1mproPM")),
                    "w2mpro": _clean(row.get("w2mpro") or row.get("w2mpro_pm") or row.get("W2mproPM")),
                    "w1snr": _clean(row.get("w1snr") or row.get("w1snr_pm") or row.get("snr1_pm") or row.get("snrW1pm")),
                    "w2snr": _clean(row.get("w2snr") or row.get("w2snr_pm") or row.get("snr2_pm") or row.get("snrW2pm")),
                }
            )
        return sources

    @staticmethod
    def _parse_allwise_rows(raw_rows: list[dict[str, str]]) -> list[dict]:
        sources: list[dict] = []
        for row in raw_rows:
            sid = (row.get("source_id") or row.get("designation") or "").strip()
            ra = _clean(row.get("ra"))
            dec = _clean(row.get("dec"))
            if not sid or ra is None or dec is None:
                continue
            sources.append(
                {
                    "source_id": sid,
                    "ra": ra,
                    "dec": dec,
                    "pmra": _clean(row.get("pmra")),
                    "pmdec": _clean(row.get("pmdec")),
                    "w1mpro": _clean(row.get("w1mpro")),
                    "w2mpro": _clean(row.get("w2mpro")),
                    "w1snr": _clean(row.get("w1snr")),
                    "w2snr": _clean(row.get("w2snr")),
                }
            )
        return sources

    async def query(self, ra: float, dec: float) -> WiseResult:
        t0 = time.perf_counter()
        radius = self.s.wise_search_radius_deg
        catwise_adql = f"""
            SELECT source_id, ra, dec, pmra, pmdec,
                   w1mpro_pm AS w1mpro, w2mpro_pm AS w2mpro,
                   w1snr_pm AS w1snr, w2snr_pm AS w2snr
            FROM {self.s.wise_table}
            WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.8f}, {dec:.8f}, {radius:.6f})) = 1
        """
        try:
            raw_rows = await self._tap_csv(catwise_adql)
            sources = self._parse_catwise_rows(raw_rows)
            log.info("CatWISE2020 returned %d sources in %.2fs", len(sources), time.perf_counter() - t0)
            return WiseResult(sources=sources, catalog="CatWISE2020", epoch_mjd=CATWISE_EPOCH_MJD)
        except asyncio.CancelledError:
            raise
        except Exception as exc:
            log.warning("CatWISE2020 query failed (%r); falling back to %s", exc, self.s.allwise_table)
            allwise_adql = f"""
                SELECT designation AS source_id, ra, dec, w1mpro, w2mpro, w1snr, w2snr
                FROM {self.s.allwise_table}
                WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.8f}, {dec:.8f}, {radius:.6f})) = 1
            """
            try:
                raw_rows = await self._tap_csv(allwise_adql)
                sources = self._parse_allwise_rows(raw_rows)
                log.info("AllWISE fallback returned %d sources in %.2fs", len(sources), time.perf_counter() - t0)
                return WiseResult(sources=sources, catalog="AllWISE", epoch_mjd=ALLWISE_EPOCH_MJD)
            except asyncio.CancelledError:
                raise
            except Exception as exc2:
                log.warning("AllWISE fallback failed (%r); trying VizieR CatWISE2020 mirror", exc2)
                vizier_adql = f"""
                    SELECT objID AS source_id, RAdeg AS ra, DEdeg AS dec,
                           pmRA AS pmra, pmDE AS pmdec,
                           W1mproPM AS w1mpro, W2mproPM AS w2mpro,
                           snrW1pm AS w1snr, snrW2pm AS w2snr
                    FROM "II/365/catwise"
                    WHERE CONTAINS(POINT('ICRS', RAdeg, DEdeg), CIRCLE('ICRS', {ra:.8f}, {dec:.8f}, {radius:.6f})) = 1
                """
                try:
                    raw_rows = await self._tap_csv(
                        vizier_adql, tap_url="https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
                    )
                    sources = self._parse_catwise_rows(raw_rows)
                    log.info("VizieR CatWISE2020 mirror returned %d sources in %.2fs", len(sources), time.perf_counter() - t0)
                    return WiseResult(sources=sources, catalog="CatWISE2020", epoch_mjd=CATWISE_EPOCH_MJD)
                except asyncio.CancelledError:
                    raise
                except Exception as exc3:
                    log.warning("VizieR CatWISE2020 fallback also failed: %r", exc3)
                    return WiseResult(sources=[], catalog="CatWISE2020", epoch_mjd=CATWISE_EPOCH_MJD, error=repr(exc2))

    def best_match(
        self,
        res: WiseResult,
        ra: float,
        dec: float,
        mjd: float,
        spherex_ra: float | None = None,
        spherex_dec: float | None = None,
    ) -> WiseSource | None:
        ref_ra = spherex_ra if spherex_ra is not None else ra
        ref_dec = spherex_dec if spherex_dec is not None else dec

        best: tuple[dict, float, float, float] | None = None
        for s in res.sources:
            if s.get("ra") is None or s.get("dec") is None:
                continue
            ra_o, dec_o = propagate_wise(s, mjd, res.epoch_mjd)
            sep_prop = ang_sep_arcsec(ra_o, dec_o, ref_ra, ref_dec)
            sep_cat = ang_sep_arcsec(float(s["ra"]), float(s["dec"]), ra, dec)
            sep = min(sep_prop, sep_cat)
            if sep <= self.s.wise_match_radius_arcsec and (best is None or sep < best[3]):
                best = (s, ra_o, dec_o, sep)

        if best is None:
            return None

        s, ra_o, dec_o, sep = best
        w1 = s.get("w1mpro")
        w2 = s.get("w2mpro")
        w1_w2_color = round(w1 - w2, 3) if (w1 is not None and w2 is not None) else None
        is_ultracool = bool(w1_w2_color is not None and w1_w2_color >= 1.5)

        pmra = s.get("pmra")
        pmdec = s.get("pmdec")
        if pmra is not None and s.get("pm_unit") == "arcsec/yr":
            pmra = float(pmra) * 1000.0
        if pmdec is not None and s.get("pm_unit") == "arcsec/yr":
            pmdec = float(pmdec) * 1000.0

        pos_ra = spherex_ra if spherex_ra is not None else ra_o
        pos_dec = spherex_dec if spherex_dec is not None else dec_o
        dt_yr = max(abs(mjd - res.epoch_mjd) / 365.25, 0.1)
        baseline_shift_arcsec = ang_sep_arcsec(float(s["ra"]), float(s["dec"]), pos_ra, pos_dec)
        if baseline_shift_arcsec < 0.05 and pmra is not None and pmdec is not None:
            baseline_pm = round(math.hypot(float(pmra), float(pmdec)) / 1000.0, 3)
        else:
            baseline_pm = round(baseline_shift_arcsec / dt_yr, 3)

        return WiseSource(
            source_id=str(s["source_id"]),
            catalog=res.catalog,
            ra=round(float(s["ra"]), 7),
            dec=round(float(s["dec"]), 7),
            ra_obs=round(ra_o, 7),
            dec_obs=round(dec_o, 7),
            separation_arcsec=round(sep, 3),
            w1mpro=round(float(w1), 3) if w1 is not None else None,
            w2mpro=round(float(w2), 3) if w2 is not None else None,
            w1_w2_color=w1_w2_color,
            w1snr=round(float(s["w1snr"]), 2) if s.get("w1snr") is not None else None,
            w2snr=round(float(s["w2snr"]), 2) if s.get("w2snr") is not None else None,
            pmra_mas_yr=round(float(pmra), 2) if pmra is not None else None,
            pmdec_mas_yr=round(float(pmdec), 2) if pmdec is not None else None,
            baseline_pm_arcsec_yr=baseline_pm,
            is_ultracool_color=is_ultracool,
        )
