"""Supabase PostgreSQL + Storage service (async, with local-disk fallback).

Features
--------
* Uses `supabase.create_async_client` (`AsyncClient`) for non-blocking PostgREST
  queries (`targets`, `observations`, and the `find_target` spherical-tolerance RPC).
* Uploads PNG layers and small float32 co-add FITS files in parallel to the
  public Storage bucket (`skysurvey`) with `x-upsert: true` and `Cache-Control`
  headers, returning deterministic public CDN URLs.
* If `SUPABASE_URL` / `SUPABASE_KEY` are not set (e.g. local hackathon demo or
  CI tests), transparently falls back to an in-process store + local disk directory
  served at `/static/storage/...`.
"""
from __future__ import annotations

import asyncio
import json
import logging
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import httpx

from app.config import Settings
from app.core.keys import normalize_dec, normalize_fov, normalize_ra, target_key
from app.services.gaia_service import ang_sep_arcsec

log = logging.getLogger(__name__)


class SupabaseService:
    def __init__(self, settings: Settings, http: httpx.AsyncClient):
        self.s = settings
        self.http = http
        self._client = None
        self._client_lock = asyncio.Lock()
        self._upload_sem = asyncio.Semaphore(settings.upload_concurrency)

        # In-memory fallback + L1 Datalink URL cache (persisted on disk across reloads)
        self._mem_targets: dict[str, dict[str, Any]] = {}
        self._obs_url_cache: dict[str, str] = {}
        self.local_storage_dir: Path = settings.cache_dir / "storage"
        self.local_storage_dir.mkdir(parents=True, exist_ok=True)
        self._local_db_file: Path = settings.cache_dir / "local_db.json"
        self._load_local_db()

    def _load_local_db(self) -> None:
        if self._local_db_file.is_file():
            try:
                data = json.loads(self._local_db_file.read_text())
                # Only load local target cache when Supabase is NOT configured;
                # when Supabase is enabled, PostgreSQL is the source of truth.
                if not self.enabled:
                    self._mem_targets = data.get("targets", {})
                self._obs_url_cache = data.get("observations", {})
            except Exception:
                pass

    def _save_local_db(self) -> None:
        try:
            self._local_db_file.write_text(
                json.dumps({"targets": self._mem_targets, "observations": self._obs_url_cache}, default=str)
            )
        except Exception:
            pass

    @property
    def enabled(self) -> bool:
        return bool(self.s.supabase_url and self.s.supabase_api_key)

    async def _get_client(self):
        if not self.enabled:
            return None
        if self._client is None:
            async with self._client_lock:
                if self._client is None:
                    from supabase import create_async_client

                    self._client = await create_async_client(
                        self.s.supabase_url.rstrip("/"),
                        self.s.supabase_api_key,
                    )
        return self._client

    # ======================================================= targets table ===
    async def get_target(
        self, ra: float, dec: float, fov: float, band: str = "all"
    ) -> dict[str, Any] | None:
        """Look up target by exact canonical `target_key`, falling back to the
        spherical-tolerance RPC `find_target` within `coord_match_tol_arcsec`."""
        key = target_key(ra, dec, fov, band)
        # Fast L1 in-process hit (includes targets whose BackgroundTasks cloud sync is still in flight)
        if key in self._mem_targets:
            return dict(self._mem_targets[key])
        for row in self._mem_targets.values():
            if (
                row.get("band") == band
                and abs(float(row.get("fov", 0)) - fov) < 5e-4
                and ang_sep_arcsec(float(row["ra"]), float(row["dec"]), ra, dec) <= self.s.coord_match_tol_arcsec
            ):
                return dict(row)

        client = await self._get_client()
        if client is not None:
            try:
                res = (
                    await client.table(self.s.targets_table)
                    .select("*")
                    .eq("target_key", key)
                    .limit(1)
                    .execute()
                )
                if res and getattr(res, "data", None):
                    rows = res.data
                    if isinstance(rows, list) and rows:
                        self._mem_targets[key] = dict(rows[0])
                        return dict(rows[0])
                # Tolerance search via RPC
                rpc = await client.rpc(
                    "find_target",
                    {
                        "p_ra": normalize_ra(ra),
                        "p_dec": normalize_dec(dec),
                        "p_fov": normalize_fov(fov),
                        "p_band": band,
                        "p_tol_arcsec": self.s.coord_match_tol_arcsec,
                    },
                ).execute()
                if rpc and getattr(rpc, "data", None):
                    rows = rpc.data
                    if isinstance(rows, list) and rows:
                        self._mem_targets[rows[0].get("target_key", key)] = dict(rows[0])
                        return dict(rows[0])
                # Supabase query succeeded and returned 0 rows -> true cache miss
                return None
            except Exception as e:
                log.warning("Supabase get_target(%s) failed, using fallback: %s", key, e)

        return None

    async def get_target_by_coords(self, ra: float, dec: float) -> dict[str, Any] | None:
        """Find any cached target within `coord_match_tol_arcsec` of (ra, dec) that has `spectrum_data`."""
        ra_n, dec_n = normalize_ra(ra), normalize_dec(dec)
        best_row: dict[str, Any] | None = None
        best_sep = 1e9
        tol = max(self.s.coord_match_tol_arcsec, 3.0)
        for row in self._mem_targets.values():
            if not row.get("spectrum_data"):
                continue
            sep = ang_sep_arcsec(float(row["ra"]), float(row["dec"]), ra_n, dec_n)
            if sep <= tol and sep < best_sep:
                best_row, best_sep = row, sep
        if best_row is not None:
            return dict(best_row)

        client = await self._get_client()
        if client is not None:
            try:
                rpc = await client.rpc(
                    "find_target_by_coords",
                    {
                        "p_ra": ra_n,
                        "p_dec": dec_n,
                        "p_tol_arcsec": tol,
                    },
                ).execute()
                if rpc and getattr(rpc, "data", None):
                    rows = rpc.data
                    if isinstance(rows, list) and rows:
                        return dict(rows[0])
                return None
            except Exception as e:
                log.debug("Supabase find_target_by_coords failed: %s", e)

        return None

    def upsert_target_local(self, record: dict[str, Any]) -> dict[str, Any]:
        """Update in-memory and local-disk target cache immediately without waiting for Supabase network I/O."""
        now_iso = datetime.now(timezone.utc).isoformat()
        key = record["target_key"]
        existing = self._mem_targets.get(key, {})
        merged = {
            **existing,
            **record,
            "created_at": existing.get("created_at") or record.get("created_at") or now_iso,
            "updated_at": now_iso,
        }
        self._mem_targets[key] = merged
        self._save_local_db()
        return merged

    async def upsert_target(self, record: dict[str, Any]) -> dict[str, Any]:
        """Insert or update a row in `targets` keyed by `target_key`."""
        merged = self.upsert_target_local(record)
        key = record["target_key"]

        client = await self._get_client()
        if client is not None:
            payload = {k: v for k, v in merged.items() if k not in ("id", "cx", "cy", "cz")}
            try:
                res = (
                    await client.table(self.s.targets_table)
                    .upsert(payload, on_conflict="target_key")
                    .execute()
                )
                if res and getattr(res, "data", None):
                    row = res.data[0] if isinstance(res.data, list) else res.data
                    self._mem_targets[key] = dict(row)
                    return dict(row)
            except Exception as e:
                log.warning("Supabase upsert_target(%s) failed: %s", key, e)
        return merged

    # ================================================== observations cache ===
    async def get_cached_fits_urls(self, obs_keys: list[str]) -> dict[str, str]:
        """Return `{obs_key: fits_url}` for already-resolved Datalink URLs."""
        if not obs_keys:
            return {}
        found = {k: self._obs_url_cache[k] for k in obs_keys if k in self._obs_url_cache}
        missing = [k for k in obs_keys if k not in found]
        if not missing:
            return found

        client = await self._get_client()
        if client is not None:
            try:
                res = (
                    await client.table(self.s.observations_table)
                    .select("obs_key, fits_url")
                    .in_("obs_key", missing[:500])
                    .execute()
                )
                for row in getattr(res, "data", None) or []:
                    k, u = row["obs_key"], row["fits_url"]
                    self._obs_url_cache[k] = u
                    found[k] = u
            except Exception as e:
                log.debug("Supabase get_cached_fits_urls failed: %s", e)
        return found

    async def save_observations(self, frames: list[Any]) -> None:
        """Persist newly resolved `(obs_key -> fits_url)` mappings."""
        rows = []
        for f in frames:
            if not getattr(f, "fits_url", None):
                continue
            self._obs_url_cache[f.key] = f.fits_url
            rows.append(
                {
                    "obs_key": f.key,
                    "obs_id": f.obs_id,
                    "detector": f.detector,
                    "collection": f.collection,
                    "t_min": f.t_min,
                    "t_max": f.t_max,
                    "s_ra": f.s_ra,
                    "s_dec": f.s_dec,
                    "fits_url": f.fits_url,
                }
            )
        if not rows:
            return
        self._save_local_db()
        client = await self._get_client()
        if client is not None:
            try:
                await (
                    client.table(self.s.observations_table)
                    .upsert(rows, on_conflict="obs_key")
                    .execute()
                )
            except Exception as e:
                log.debug("Supabase save_observations failed: %s", e)

    # ===================================================== storage bucket ===
    def _local_public_url(self, clean: str) -> str:
        import os

        base = (self.s.public_base_url or os.getenv("RENDER_EXTERNAL_URL", "")).rstrip("/")
        return f"{base}/static/storage/{clean}?v={self.s.pipeline_version}"

    def local_url(self, storage_path: str) -> str:
        return self._local_public_url(storage_path.lstrip("/"))

    def public_url(self, storage_path: str) -> str:
        clean = storage_path.lstrip("/")
        if self.enabled:
            base = self.s.supabase_url.rstrip("/")
            return f"{base}/storage/v1/object/public/{self.s.supabase_bucket}/{clean}?v={self.s.pipeline_version}"
        return self._local_public_url(clean)

    async def save_local_many(
        self, items: list[tuple[str, bytes, str]]
    ) -> dict[str, str]:
        """Write all PNG layers and co-add FITS files locally to `/tmp/fits_cache/storage`
        via `asyncio.to_thread` and immediately return `/static/storage/...` URLs."""
        def _write_batch() -> dict[str, str]:
            out: dict[str, str] = {}
            for path, data, _ctype in items:
                clean = path.lstrip("/")
                dest = self.local_storage_dir / clean
                self._write_local(dest, data)
                out[path] = self._local_public_url(clean)
            return out

        return await asyncio.to_thread(_write_batch)

    async def sync_target_to_cloud(
        self,
        key: str,
        row: dict[str, Any],
        uploads: list[tuple[str, bytes, str]],
        frames: list[Any] | None = None,
    ) -> None:
        """Background task (`fastapi.BackgroundTasks`): upload PNGs/FITS to Supabase
        Storage and upsert `observations` + `targets` into Supabase PostgreSQL."""
        try:
            if frames:
                await self.save_observations(frames)
            if not self.enabled:
                return
            cloud_row = dict(row)
            if uploads:
                await self.upload_many(uploads)
                # Replace ephemeral /static/storage/ URLs with durable Supabase CDN URLs in DB
                local_prefix = self._local_public_url("").split("?")[0]
                cloud_prefix = f"{self.s.supabase_url.rstrip('/')}/storage/v1/object/public/{self.s.supabase_bucket}/"
                raw_json = json.dumps(cloud_row, default=str).replace(local_prefix, cloud_prefix)
                cloud_row = json.loads(raw_json)
            await self.upsert_target(cloud_row)
            log.info("Background Supabase sync completed for %s (%d files)", key, len(uploads))
        except Exception as e:
            log.warning("Background Supabase sync failed for %s: %s", key, e)

    async def upload_bytes(
        self, storage_path: str, data: bytes, content_type: str = "image/png"
    ) -> str:
        """Upload raw bytes to Supabase Storage (and mirror locally if disabled)."""
        clean = storage_path.lstrip("/")
        if not self.enabled:
            dest = self.local_storage_dir / clean
            await asyncio.to_thread(self._write_local, dest, data)
            return self.public_url(clean)

        url = f"{self.s.supabase_url.rstrip('/')}/storage/v1/object/{self.s.supabase_bucket}/{clean}"
        headers = {
            "apikey": self.s.supabase_api_key,
            "Authorization": f"Bearer {self.s.supabase_api_key}",
            "Content-Type": content_type,
            "x-upsert": "true",
            "Cache-Control": "public, max-age=31536000",
        }
        async with self._upload_sem:
            try:
                r = await self.http.post(url, content=data, headers=headers, timeout=30.0)
                r.raise_for_status()
                return self.public_url(clean)
            except Exception as e:
                log.warning("Supabase Storage upload failed for %s (%s), saving locally", clean, e)
                dest = self.local_storage_dir / clean
                await asyncio.to_thread(self._write_local, dest, data)
                return self._local_public_url(clean)

    async def upload_many(
        self, items: list[tuple[str, bytes, str]]
    ) -> dict[str, str]:
        """Concurrently upload `[(storage_path, bytes, content_type), ...]`
        and return `{storage_path: public_url}`."""
        urls = await asyncio.gather(
            *(self.upload_bytes(path, data, ctype) for path, data, ctype in items)
        )
        return {path: u for (path, _, _), u in zip(items, urls)}

    async def download_storage_bytes(self, url_or_path: str) -> bytes | None:
        """Download an object back from Storage (used when `/tmp/fits_cache` was
        evicted or the container restarted before `/restretch`)."""
        if "/static/storage/" in url_or_path:
            rel = url_or_path.split("/static/storage/", 1)[1].split("?", 1)[0]
            p = self.local_storage_dir / rel
            return p.read_bytes() if p.is_file() else None
        url = url_or_path if url_or_path.startswith("http") else self.public_url(url_or_path)
        try:
            r = await self.http.get(url, timeout=20.0)
            if r.status_code == 200:
                return r.content
        except Exception as e:
            log.warning("Storage download failed for %s: %s", url, e)
        return None

    @staticmethod
    def _write_local(dest: Path, data: bytes) -> None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
