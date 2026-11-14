"""IMCCE SkyBoT Solar System minor-planet & comet ephemeris filter.

Main-belt asteroids and comets shine brightly via thermal emission in SPHEREx
3.8-5.0 um channels without optical Gaia DR3 counterparts (because they move
rapidly across the sky). Querying SkyBoT at each epoch's `mjd_mid` lets the UI
label known Solar System bodies and separate them from Planet X candidates.
"""
from __future__ import annotations

import asyncio
import json
import logging
import math
import re
import time
from collections.abc import Callable
from typing import Any

import httpx

from app.config import Settings
from app.schemas.target import MinorPlanetOverlay
from app.services.gaia_service import ang_sep_arcsec

log = logging.getLogger(__name__)


def _to_float(val: Any) -> float | None:
    if val is None:
        return None
    s = str(val).strip()
    if not s or s in ("-", "null", "None", "nan", "NaN"):
        return None
    try:
        f = float(s)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def _parse_coord(val: Any, is_ra: bool) -> float | None:
    """Parse either decimal degrees or sexagesimal ('HH MM SS.s' / 'DD MM SS.s')."""
    if val is None:
        return None
    if isinstance(val, (int, float)):
        f = float(val)
        return (f % 360.0) if is_ra else f

    text = str(val).strip()
    if not text:
        return None
    try:
        f = float(text)
        if math.isfinite(f):
            return (f % 360.0) if is_ra else f
    except ValueError:
        pass

    sign = -1.0 if text.startswith("-") else 1.0
    cleaned = re.sub(r"[hdms:\s]+", " ", text.lstrip("+-")).strip()
    parts = cleaned.split()
    if not parts:
        return None
    try:
        nums = [float(p) for p in parts[:3]]
    except ValueError:
        return None

    deg = nums[0]
    if len(nums) > 1:
        deg += nums[1] / 60.0
    if len(nums) > 2:
        deg += nums[2] / 3600.0
    deg *= sign
    if is_ra and "d" not in text.lower():
        deg = (deg * 15.0) % 360.0
    return deg


def _normalize_col_name(col: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", col.strip().lower()).strip("_")


def _extract_rows_from_json(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [r for r in payload if isinstance(r, dict)]
    if not isinstance(payload, dict):
        return []

    for key in ("bodies", "objects", "results"):
        val = payload.get(key)
        if isinstance(val, list):
            return [r for r in val if isinstance(r, dict)]

    data = payload.get("data")
    if isinstance(data, list):
        if not data:
            return []
        if isinstance(data[0], dict):
            return data
        fields = payload.get("fields") or payload.get("columns") or payload.get("col")
        if isinstance(fields, list) and fields:
            norm_fields = [_normalize_col_name(str(f.get("name", f) if isinstance(f, dict) else f)) for f in fields]
            out: list[dict[str, Any]] = []
            for row in data:
                if isinstance(row, (list, tuple)):
                    out.append({norm_fields[i]: row[i] for i in range(min(len(norm_fields), len(row)))})
            return out
        default_cols = ["num", "name", "ra", "dec", "class", "vmag", "errpos", "d", "dra", "ddec"]
        out = []
        for row in data:
            if isinstance(row, (list, tuple)):
                out.append({default_cols[i]: row[i] for i in range(min(len(default_cols), len(row)))})
        return out
    return []


def _extract_rows_from_pipe_text(text: str) -> list[dict[str, Any]]:
    default_cols = ["num", "name", "ra", "dec", "class", "vmag", "errpos", "d", "dra", "ddec"]
    cols = default_cols
    rows: list[dict[str, Any]] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            header = line.lstrip("#").strip()
            if "|" in header and ("name" in header.lower() or "ra" in header.lower()):
                cols = [_normalize_col_name(c) for c in header.split("|")]
            continue
        if "|" not in line:
            continue
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 4:
            continue
        rows.append({cols[i]: parts[i] for i in range(min(len(cols), len(parts)))})
    return rows


def _pick(row: dict[str, Any], *keys: str) -> Any:
    for k in keys:
        if k in row and row[k] not in (None, ""):
            return row[k]
    return None


class SolarSystemService:
    def __init__(self, settings: Settings, http: httpx.AsyncClient | None = None):
        self.s = settings
        self.http = http

    async def _fetch_skybot(self, params: dict[str, str]) -> str:
        if self.http is not None:
            r = await self.http.get(self.s.skybot_url, params=params, timeout=self.s.skybot_timeout_s)
        else:
            async with httpx.AsyncClient(timeout=self.s.skybot_timeout_s) as client:
                r = await client.get(self.s.skybot_url, params=params)
        r.raise_for_status()
        return r.text

    def _parse_response(
        self,
        text: str,
        center_ra: float,
        center_dec: float,
        world_to_png: Callable[[float, float], tuple[float, float] | None] | None,
    ) -> list[MinorPlanetOverlay]:
        stripped = text.strip()
        if not stripped:
            return []

        rows: list[dict[str, Any]] = []
        if stripped.startswith("{") or stripped.startswith("["):
            try:
                rows = _extract_rows_from_json(json.loads(stripped))
            except Exception:
                rows = []
        if not rows and "|" in stripped:
            rows = _extract_rows_from_pipe_text(stripped)

        overlays: list[MinorPlanetOverlay] = []
        for row in rows:
            norm_row = {_normalize_col_name(str(k)): v for k, v in row.items()}
            name_raw = _pick(norm_row, "name", "designation", "object", "id")
            ra_raw = _pick(norm_row, "ra", "ra_deg", "raj2000", "ra_icrs")
            dec_raw = _pick(norm_row, "dec", "dec_deg", "dej2000", "dec_icrs")
            if not name_raw or ra_raw is None or dec_raw is None:
                continue

            mp_ra = _parse_coord(ra_raw, is_ra=True)
            mp_dec = _parse_coord(dec_raw, is_ra=False)
            if mp_ra is None or mp_dec is None or not (-90.0 <= mp_dec <= 90.0):
                continue

            if world_to_png is not None:
                xy = world_to_png(mp_ra, mp_dec)
                if xy is None:
                    continue
                px, py = round(float(xy[0]), 2), round(float(xy[1]), 2)
            else:
                px, py = 0.0, 0.0

            num_raw = _pick(norm_row, "num", "number")
            num_str = str(num_raw).strip() if num_raw not in (None, "", "-", "0") else None
            obj_type = str(_pick(norm_row, "class", "type", "object_type") or "Asteroid").strip()
            vmag = _to_float(_pick(norm_row, "vmag", "mv", "v_mag", "mag"))
            dra = _to_float(_pick(norm_row, "dra", "d_ra", "mu_ra", "dra_cosdec"))
            ddec = _to_float(_pick(norm_row, "ddec", "d_dec", "mu_dec"))
            sep = ang_sep_arcsec(center_ra, center_dec, mp_ra, mp_dec)

            overlays.append(
                MinorPlanetOverlay(
                    number=num_str,
                    name=str(name_raw).strip(),
                    object_type=obj_type,
                    ra=round(mp_ra, 6),
                    dec=round(mp_dec, 6),
                    v_mag=round(vmag, 2) if vmag is not None else None,
                    dra_arcsec_hr=round(dra, 3) if dra is not None else None,
                    ddec_arcsec_hr=round(ddec, 3) if ddec is not None else None,
                    separation_arcsec=round(sep, 2),
                    x=px,
                    y=py,
                )
            )

        overlays.sort(key=lambda o: o.separation_arcsec)
        return overlays

    async def query_epoch(
        self,
        ra: float,
        dec: float,
        fov_deg: float,
        mjd: float,
        world_to_png: Callable[[float, float], tuple[float, float] | None] | None = None,
    ) -> list[MinorPlanetOverlay]:
        if not self.s.skybot_enabled:
            return []

        params = {
            "-ep": f"{mjd:.6f}",
            "-ra": f"{ra:.6f}",
            "-dec": f"{dec:.6f}",
            "-rd": f"{max(fov_deg * 0.75, 0.05):.4f}",
            "-mime": "json",
            "-output": "basic",
            "-loc": "500",
        }
        t0 = time.perf_counter()
        try:
            text = await self._fetch_skybot(params)
            overlays = self._parse_response(text, ra, dec, world_to_png)
            log.debug(
                "SkyBoT MJD %.3f returned %d minor planets in %.2fs",
                mjd,
                len(overlays),
                time.perf_counter() - t0,
            )
            return overlays
        except Exception as exc:
            log.warning("SkyBoT query failed at MJD %.3f: %s", mjd, exc)
            return []

    async def query_epochs(
        self,
        ra: float,
        dec: float,
        fov_deg: float,
        mjds: list[float],
        world_to_png: Callable[[float, float], tuple[float, float] | None] | None = None,
    ) -> list[list[MinorPlanetOverlay]]:
        if not mjds or not self.s.skybot_enabled:
            return [[] for _ in mjds]
        results = await asyncio.gather(
            *(self.query_epoch(ra, dec, fov_deg, mjd, world_to_png) for mjd in mjds)
        )
        return list(results)
