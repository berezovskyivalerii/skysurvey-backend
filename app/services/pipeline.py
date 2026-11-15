"""End-to-end processing pipeline & restretch orchestrator."""
from __future__ import annotations

import asyncio
import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from app.config import Settings
from app.core.keys import normalize_dec, normalize_fov, normalize_ra, target_key
from app.core.spherex import BAND_GROUPS, DEFAULT_FLAG_BITS
from app.schemas.spectrum import SpectrumResponse
from app.schemas.target import (
    BandName,
    CacheInfo,
    DataRange,
    DifferenceResult,
    DisplayInfo,
    EpochOut,
    LayerUrls,
    Progress,
    RestretchEpoch,
    RestretchRequest,
    RestretchResponse,
    Spectrum,
    StretchInfo,
    TargetResponse,
)
from app.services import fits_processor as fp
from app.services.cache_manager import FitsCache
from app.services.gaia_service import GaiaResult, GaiaService
from app.services.irsa_service import Epoch, Frame, IrsaService
from app.services.spectrum_extractor import extract_spectrum, find_covering_fits
from app.services.solar_system_service import SolarSystemService
from app.services.supabase_service import SupabaseService
from app.services.wise_service import WiseResult, WiseService

log = logging.getLogger(__name__)


class PipelineService:
    def __init__(
        self,
        settings: Settings,
        cache: FitsCache,
        irsa: IrsaService,
        gaia: GaiaService,
        db: SupabaseService,
        wise: WiseService | None = None,
        solar_system: SolarSystemService | None = None,
    ):
        self.s = settings
        self.cache = cache
        self.irsa = irsa
        self.gaia = gaia
        self.db = db
        self.wise = wise or WiseService(settings, irsa.http)
        self.solar_system = solar_system or SolarSystemService(settings, irsa.http)
        self._inflight: dict[str, asyncio.Task[dict[str, Any]]] = {}
        self._pending_sync: dict[str, tuple[dict[str, Any], list[tuple[str, bytes, str]], list[Frame]]] = {}
        self._bg_tasks: set[asyncio.Task[Any]] = set()

    # ===================================================== public entrypoint ==
    async def get_target(
        self,
        ra: float,
        dec: float,
        fov: float,
        band: BandName = "all",
        sync: bool = True,
        force_refresh: bool = False,
        fast: bool = True,
        background_tasks: Any | None = None,
    ) -> TargetResponse:
        ra_n, dec_n, fov_n = normalize_ra(ra), normalize_dec(dec), normalize_fov(fov)
        key = target_key(ra_n, dec_n, fov_n, band)

        if not force_refresh:
            cached = await self.db.get_target(ra_n, dec_n, fov_n, band)
            if cached and cached.get("pipeline_version") == self.s.pipeline_version:
                status = cached.get("status")
                if status in ("ready", "partial", "no_data"):
                    stale = self._is_stale(cached)
                    if stale and key not in self._inflight:
                        self._spawn_pipeline(ra_n, dec_n, fov_n, band, key, fast=fast)
                    return self._row_to_response(cached, hit=True, stale=stale)
                if status == "processing" and not self._is_lock_expired(cached):
                    if not sync:
                        return self._row_to_response(cached, hit=True, stale=False)
                    task = self._inflight.get(key)
                    if task:
                        row = await asyncio.wait_for(asyncio.shield(task), timeout=self.s.sync_timeout_s)
                        self._schedule_cloud_sync(key, background_tasks)
                        return self._row_to_response(row, hit=False, stale=False)

        # Cache miss (or forced refresh)
        task = self._inflight.get(key)
        if task is None:
            task = self._spawn_pipeline(ra_n, dec_n, fov_n, band, key, fast=fast, force_refresh=force_refresh)

        if not sync:
            # Give the task a brief moment to write the initial 'processing' row
            await asyncio.sleep(0.05)
            cached = await self.db.get_target(ra_n, dec_n, fov_n, band)
            if cached:
                return self._row_to_response(cached, hit=False, stale=False)

        row = await asyncio.wait_for(asyncio.shield(task), timeout=self.s.sync_timeout_s)
        self._schedule_cloud_sync(key, background_tasks)
        return self._row_to_response(row, hit=False, stale=False)

    def _schedule_cloud_sync(self, key: str, background_tasks: Any | None) -> None:
        sync_item = self._pending_sync.pop(key, None)
        if sync_item is None:
            return
        sync_row, sync_uploads, sync_frames = sync_item
        if background_tasks is not None and hasattr(background_tasks, "add_task"):
            background_tasks.add_task(
                self.db.sync_target_to_cloud,
                key,
                sync_row,
                sync_uploads,
                sync_frames,
            )
        else:
            bg = asyncio.create_task(
                self.db.sync_target_to_cloud(key, sync_row, sync_uploads, sync_frames),
                name=f"cloud-sync:{key}",
            )
            self._bg_tasks.add(bg)
            bg.add_done_callback(self._bg_tasks.discard)

    # ========================================================== restretch ====
    async def restretch(self, req: RestretchRequest) -> RestretchResponse:
        t0 = time.perf_counter()
        ra_n, dec_n, fov_n = normalize_ra(req.ra), normalize_dec(req.dec), normalize_fov(req.fov)
        key = target_key(ra_n, dec_n, fov_n, req.band)

        cached = await self.db.get_target(ra_n, dec_n, fov_n, req.band)
        source: str = "local_cache"

        if not cached or cached.get("status") not in ("ready", "partial"):
            await self.get_target(ra_n, dec_n, fov_n, req.band, sync=True)
            cached = await self.db.get_target(ra_n, dec_n, fov_n, req.band)
            source = "rebuilt"

        epochs_meta = (cached or {}).get("epochs") or []
        wanted_idx = set(req.epochs) if req.epochs is not None else {int(e["index"]) for e in epochs_meta}

        images: list[tuple[int, Any]] = []
        for ep in epochs_meta:
            idx = int(ep["index"])
            if idx not in wanted_idx:
                continue
            local_p = self.cache.get(self.cache.coadd_path(key, idx))
            if local_p is not None:
                arr = await asyncio.to_thread(fp.read_coadd, local_p)
                images.append((idx, arr))
                continue
            # Fallback: fetch small float32 co-add FITS from Supabase Storage
            coadd_url = ep.get("coadd_fits_url")
            if coadd_url:
                raw = await self.db.download_storage_bytes(coadd_url)
                if raw:
                    self.cache.write_atomic(self.cache.coadd_path(key, idx), raw)
                    arr = await asyncio.to_thread(fp.read_coadd, raw)
                    images.append((idx, arr))
                    if source == "local_cache":
                        source = "storage"

        if not images and epochs_meta:
            # Co-adds were missing both locally and in Storage -> rebuild once
            await self.get_target(ra_n, dec_n, fov_n, req.band, sync=True, force_refresh=True)
            source = "rebuilt"
            for ep in epochs_meta:
                idx = int(ep["index"])
                if idx in wanted_idx:
                    local_p = self.cache.get(self.cache.coadd_path(key, idx))
                    if local_p:
                        images.append((idx, await asyncio.to_thread(fp.read_coadd, local_p)))

        vmin, vmax = req.min_cut, req.max_cut
        if vmin is None or vmax is None:
            zmin, zmax, _ = await asyncio.to_thread(fp.compute_cuts, [im for _, im in images])
            vmin = zmin if vmin is None else vmin
            vmax = zmax if vmax is None else vmax

        def _render_all() -> list[RestretchEpoch]:
            out: list[RestretchEpoch] = []
            for idx, im in images:
                png = fp.render_image(
                    im,
                    vmin=float(vmin),
                    vmax=float(vmax),
                    stretch_type=req.stretch_type,
                    size=self.s.png_size,
                    resample=self.s.png_resample,
                    asinh_a=req.asinh_a,
                    log_a=req.log_a,
                )
                out.append(RestretchEpoch(index=idx, image=fp.png_data_url(png)))
            return out

        rendered = await asyncio.to_thread(_render_all)
        elapsed_ms = round((time.perf_counter() - t0) * 1000.0, 2)
        return RestretchResponse(
            target_key=key,
            stretch=StretchInfo(stretch_type=req.stretch_type, vmin=round(float(vmin), 5), vmax=round(float(vmax), 5)),
            width=self.s.png_size,
            height=self.s.png_size,
            epochs=rendered,
            source=source,  # type: ignore[arg-type]
            elapsed_ms=elapsed_ms,
        )

    # ==================================================== internal pipeline ==
    def _spawn_pipeline(
        self,
        ra: float,
        dec: float,
        fov: float,
        band: BandName,
        key: str,
        fast: bool = True,
        force_refresh: bool = False,
    ) -> asyncio.Task[dict[str, Any]]:
        task = asyncio.create_task(
            self._run_pipeline(ra, dec, fov, band, key, fast=fast, force_refresh=force_refresh),
            name=f"pipeline:{key}",
        )
        self._inflight[key] = task
        task.add_done_callback(lambda _: self._inflight.pop(key, None))
        return task

    async def _safe_gaia_query(self, ra: float, dec: float, fov: float) -> GaiaResult:
        try:
            return await self.gaia.query(ra, dec, fov)
        except asyncio.CancelledError:
            raise
        except Exception as exc:
            log.warning("Isolated Gaia error for (%.4f, %+.4f): %r", ra, dec, exc)
            return GaiaResult([], [], error=f"Gaia unavailable: {exc!r}")

    async def _safe_wise_query(self, ra: float, dec: float) -> WiseResult:
        try:
            return await self.wise.query(ra, dec)
        except asyncio.CancelledError:
            raise
        except Exception as exc:
            log.warning("Isolated WISE error for (%.4f, %+.4f): %r", ra, dec, exc)
            return WiseResult(sources=[], error=repr(exc))

    async def _run_pipeline(
        self,
        ra: float,
        dec: float,
        fov: float,
        band: BandName,
        key: str,
        fast: bool = True,
        force_refresh: bool = False,
    ) -> dict[str, Any]:
        t0 = time.perf_counter()
        prev_cached: dict[str, Any] | None = None
        if force_refresh:
            try:
                existing = await self.db.get_target(ra, dec, fov, band)
                if existing and existing.get("status") in ("ready", "partial", "no_data"):
                    prev_cached = existing
            except Exception:
                pass

        base_row: dict[str, Any] = {
            "target_key": key,
            "ra": ra,
            "dec": dec,
            "fov": fov,
            "band": band,
            "status": "processing",
            "classification": "unknown",
            "pipeline_version": self.s.pipeline_version,
            "n_observations": 0,
            "progress": {"stage": "Querying IRSA TAP, Gaia DR3 & CatWISE2020", "percent": 10},
            "epochs": [],
            "error": None,
        }
        self.db.upsert_target_local(base_row)

        gaia_task = asyncio.create_task(self._safe_gaia_query(ra, dec, fov), name=f"gaia:{key}")
        wise_task = asyncio.create_task(self._safe_wise_query(ra, dec), name=f"wise:{key}")
        skybot_task: asyncio.Task[Any] | None = None

        try:
            # Step 3: Start IRSA TAP/SIA immediately; Gaia & CatWISE run concurrently in background
            # so slow catalog queries never delay Datalink resolution or FITS cutout downloads.
            try:
                irsa_frames = self.irsa.query_frames(ra, dec, force_refresh=force_refresh)
            except TypeError:
                # Backwards compatibility for test fakes with signature (ra, dec)
                irsa_frames = self.irsa.query_frames(ra, dec)

            frames = await irsa_frames
            if not frames:
                gaia_res, wise_res = await asyncio.gather(gaia_task, wise_task)

                def _build_no_data() -> dict[str, Any]:
                    grid = fp.make_grid(ra, dec, fov, self.s.spherex_pixel_scale_arcsec, self.s.png_size)
                    mjd_ref = 61000.0
                    gmatch = self.gaia.best_match(gaia_res, ra, dec, mjd_ref)
                    xmatch = self.gaia.classify(gmatch, None, gaia_res.error)
                    xmatch.wise_match = self.wise.best_match(wise_res, ra, dec, mjd_ref)
                    xmatch.overlay = self.gaia.overlay(gaia_res, mjd_ref, grid.world_to_png)
                    return {
                        **base_row,
                        "status": "no_data",
                        "classification": xmatch.classification,
                        "progress": {"stage": "No SPHEREx coverage", "percent": 100},
                        "crossmatch": xmatch.model_dump(mode="json"),
                    }

                no_data_row = await asyncio.to_thread(_build_no_data)
                saved = self.db.upsert_target_local(no_data_row)
                self._pending_sync[key] = (saved, [], [])
                return saved

            epochs = self.irsa.group_epochs(frames)
            # On fast Cold Start, cap to at most 3 epochs x 1 best-centered frame (= 2-3 cutouts total)
            # so NASA IRSA IBE never throttles or queues per-IP FITS slicing requests.
            active_epochs = epochs[-3:] if (fast and len(epochs) > 3) else epochs
            detectors = BAND_GROUPS.get(band, BAND_GROUPS["all"])
            display_by_epoch: list[tuple[Epoch, list[Frame]]] = []
            display_frames: list[Frame] = []
            for ep in active_epochs:
                sel = self.irsa.select_display_frames(
                    ep,
                    detectors,
                    ra,
                    dec,
                    max_per_epoch=1 if fast else None,
                )
                if sel:
                    display_by_epoch.append((ep, sel))
                    display_frames.extend(sel)

            # Launch SkyBoT minor-planet ephemeris queries concurrently with FITS downloads
            early_grid = fp.make_grid(ra, dec, fov, self.s.spherex_pixel_scale_arcsec, self.s.png_size)
            early_mjds = [ep.mjd_mid for ep, _ in display_by_epoch]
            skybot_task = asyncio.create_task(
                self.solar_system.query_epochs(ra, dec, fov, early_mjds, early_grid.world_to_png),
                name=f"skybot:{key}",
            )

            self.db.upsert_target_local(
                {
                    **base_row,
                    "n_observations": len(frames),
                    "progress": {"stage": f"Resolving Datalink & downloading {len(display_frames)} cutouts", "percent": 35},
                }
            )

            # Resolve Datalink URLs & download cutouts (padded by 20% so rotated detector edges cover corners)
            known_urls = await self.db.get_cached_fits_urls([f.key for f in display_frames])
            cutout_deg = min(fov * 1.25, 0.5)
            if (
                hasattr(self.irsa, "resolve_and_fetch_many")
                and getattr(self.irsa.fetch_many, "__func__", None) is IrsaService.fetch_many
            ):
                downloaded = await self.irsa.resolve_and_fetch_many(
                    display_frames, known_urls, ra, dec, cutout_deg
                )
            else:
                await self.irsa.resolve_fits_urls(display_frames, known_urls)
                downloaded = await self.irsa.fetch_many(display_frames, ra, dec, cutout_deg)
            path_by_key: dict[str, Path] = {f.key: p for f, p in downloaded}

            self.db.upsert_target_local(
                {
                    **base_row,
                    "n_observations": len(frames),
                    "progress": {"stage": "WCS reprojecting epochs, rendering PNG layers & spectrum", "percent": 65},
                }
            )

            # Collect Gaia & CatWISE results (already running since t=0) before CPU cross-match
            gaia_res, wise_res = await asyncio.gather(gaia_task, wise_task)

            # CPU step in thread pool: load frames, WCS reproject, co-add, render PNGs, extract fast spectrum
            cpu_res = await asyncio.to_thread(
                self._process_display_cpu,
                ra,
                dec,
                fov,
                key,
                display_by_epoch,
                path_by_key,
                gaia_res,
                wise_res,
                fast,
            )

            # Save PNGs + co-add FITS locally to /tmp/fits_cache/storage (/static/storage/...)
            # and await SkyBoT (already running since before FITS downloads)
            url_map, all_minor_planets = await asyncio.gather(
                self.db.save_local_many(cpu_res["uploads"]),
                skybot_task,
            )
            mp_by_index = {
                ep.index: all_minor_planets[i]
                for i, (ep, _) in enumerate(display_by_epoch)
                if i < len(all_minor_planets)
            }
            minor_planets_by_epoch = [
                mp_by_index.get(ep_info["coadd"].index, []) for ep_info in cpu_res["epochs_info"]
            ]

            epochs_out: list[EpochOut] = []
            for idx_ep, ep_info in enumerate(cpu_res["epochs_info"]):
                c: fp.EpochCoadd = ep_info["coadd"]
                epochs_out.append(
                    EpochOut(
                        index=c.index,
                        label=c.label,
                        date_obs=datetime.fromisoformat(c.date_iso),
                        mjd_start=round(c.mjd_start, 5),
                        mjd_end=round(c.mjd_end, 5),
                        mjd_mid=round(c.mjd_mid, 5),
                        n_frames=c.n_frames,
                        detectors=c.detectors,
                        layers=LayerUrls(
                            image=url_map[ep_info["image_path"]],
                            variance=url_map[ep_info["variance_path"]],
                            flags=url_map[ep_info["flags_path"]],
                            rgb=url_map.get(ep_info["rgb_path"]),
                        ),
                        coadd_fits_url=url_map[ep_info["coadd_fits_path"]],
                        known_minor_planets=minor_planets_by_epoch[idx_ep],
                    )
                )

            diff_res: DifferenceResult | None = None
            if cpu_res.get("diff_meta") is not None:
                dm = cpu_res["diff_meta"]
                diff_res = DifferenceResult(
                    epoch_a=dm["epoch_a"],
                    epoch_b=dm["epoch_b"],
                    mjd_a=dm["mjd_a"],
                    mjd_b=dm["mjd_b"],
                    baseline_days=dm["baseline_days"],
                    diff_image_url=url_map[dm["diff_path"]],
                    deep_stack_url=url_map.get(dm["deep_path"]),
                    max_dipole_snr=dm["max_dipole_snr"],
                    anomalies=dm["anomalies"],
                )

            used_keys = set(path_by_key.keys())
            extra_spec_frames = [] if fast else self.irsa.select_spectrum_frames(epochs, used_keys)
            has_phase2 = len(extra_spec_frames) > 0

            spectrum_obj = Spectrum(
                complete=not has_phase2,
                n_frames_used=cpu_res["n_loaded"],
                n_channels_filled=sum(1 for p in cpu_res["spectrum_points"] if p.flux is not None),
                points=cpu_res["spectrum_points"],
            )
            xmatch: Any = cpu_res["crossmatch"]
            display_info: DisplayInfo = cpu_res["display"]

            row = self.db.upsert_target_local(
                {
                    "target_key": key,
                    "ra": ra,
                    "dec": dec,
                    "fov": fov,
                    "band": band,
                    "status": "partial" if has_phase2 else "ready",
                    "classification": xmatch.classification,
                    "pipeline_version": self.s.pipeline_version,
                    "n_observations": len(frames),
                    "progress": (
                        {"stage": "Refining 102-channel LVF spectrum in background", "percent": 85}
                        if has_phase2
                        else {"stage": "Complete", "percent": 100}
                    ),
                    "display": display_info.model_dump(mode="json"),
                    "epochs": [e.model_dump(mode="json") for e in epochs_out],
                    "crossmatch": xmatch.model_dump(mode="json"),
                    "spectrum": spectrum_obj.model_dump(mode="json"),
                    "spectrum_data": cpu_res["spectrum_data"],
                    "difference": diff_res.model_dump(mode="json") if diff_res else None,
                    "error": None,
                }
            )
            self._pending_sync[key] = (row, cpu_res["uploads"], display_frames)
            log.info(
                "Target %s phase-1 ready in %.2fs (%d epochs, %d frames, status=%s)",
                key,
                time.perf_counter() - t0,
                len(epochs_out),
                cpu_res["n_loaded"],
                row["status"],
            )

            if has_phase2:
                bg = asyncio.create_task(
                    self._complete_spectrum_bg(
                        key=key,
                        ra=cpu_res["phot_ra"],
                        dec=cpu_res["phot_dec"],
                        initial_samples=cpu_res["samples"],
                        initial_frames_count=cpu_res["n_loaded"],
                        extra_frames=extra_spec_frames,
                        row_snapshot=row,
                    ),
                    name=f"spectrum:{key}",
                )
                self._bg_tasks.add(bg)
                bg.add_done_callback(self._bg_tasks.discard)

            return row

        except Exception as e:
            if prev_cached is not None:
                log.warning("Pipeline refresh failed for %s (%r); serving stale cached target row", key, e)
                return self.db.upsert_target_local(prev_cached)
            log.exception("Pipeline failed for %s", key)
            failed_row = {
                **base_row,
                "status": "failed",
                "progress": {"stage": "Failed", "percent": 100},
                "error": str(e),
            }
            return self.db.upsert_target_local(failed_row)
        finally:
            for aux_task in (gaia_task, wise_task, skybot_task):
                if aux_task is not None and not aux_task.done():
                    aux_task.cancel()

    def _process_display_cpu(
        self,
        ra: float,
        dec: float,
        fov: float,
        key: str,
        display_by_epoch: list[tuple[Epoch, list[Frame]]],
        path_by_key: dict[str, Path],
        gaia_res: GaiaResult,
        wise_res: WiseResult | None = None,
        fast: bool = True,
    ) -> dict[str, Any]:
        grid = fp.make_grid(ra, dec, fov, self.s.spherex_pixel_scale_arcsec, self.s.png_size)
        coadds: list[fp.EpochCoadd] = []
        all_loaded: list[fp.FrameData] = []
        bits_merged: dict[str, int] = dict(DEFAULT_FLAG_BITS)

        for ep, flist in display_by_epoch:
            loaded_ep: list[fp.FrameData] = []
            for f in flist:
                p = path_by_key.get(f.key)
                if p is None:
                    continue
                wpath = self.cache.wavecal_path(self.irsa.wavecal_key(f))
                try:
                    fd = fp.load_frame(p, f.key, f.detector, f.t_mid, wpath)
                    loaded_ep.append(fd)
                    all_loaded.append(fd)
                    bits_merged.update(fd.bits)
                except Exception as e:
                    log.warning("Skipping corrupt cutout %s: %s", f.key, e)
            if loaded_ep:
                coadds.append(fp.coadd_epoch(ep.index, loaded_ep, grid))

        if not coadds:
            raise RuntimeError("No usable SPHEREx FITS cutouts could be processed for this target")

        vmin, vmax, data_range = fp.compute_cuts([c.image for c in coadds])

        uploads: list[tuple[str, bytes, str]] = []
        epochs_info: list[dict[str, Any]] = []
        for c in coadds:
            date_tag = c.date_iso[:10]
            base = f"targets/{key}/epoch_{c.index:02d}_{date_tag}"
            img_png = fp.render_image(c.image, vmin, vmax, "asinh", self.s.png_size, self.s.png_resample)
            var_png = fp.render_variance(c.variance, self.s.png_size, self.s.png_resample)
            flg_png = fp.render_flags(c.flags, c.coverage, bits_merged, self.s.png_size)
            rgb_png = fp.render_rgb_composite(c, self.s.png_size, self.s.png_resample)
            coadd_bytes = fp.coadd_to_fits_bytes(c, grid)

            # Persist co-add locally for <50ms /restretch and /diff
            self.cache.write_atomic(self.cache.coadd_path(key, c.index), coadd_bytes)

            p_img = f"{base}_image.png"
            p_var = f"{base}_variance.png"
            p_flg = f"{base}_flags.png"
            p_rgb = f"{base}_rgb.png"
            p_fit = f"{base}_coadd.fits"
            uploads.extend(
                [
                    (p_img, img_png, "image/png"),
                    (p_var, var_png, "image/png"),
                    (p_flg, flg_png, "image/png"),
                    (p_rgb, rgb_png, "image/png"),
                    (p_fit, coadd_bytes, "application/fits"),
                ]
            )
            epochs_info.append(
                {
                    "coadd": c,
                    "image_path": p_img,
                    "variance_path": p_var,
                    "flags_path": p_flg,
                    "rgb_path": p_rgb,
                    "coadd_fits_path": p_fit,
                }
            )

        mean_mjd = float(sum(c.mjd_mid for c in coadds) / len(coadds))
        gmatch = self.gaia.best_match(gaia_res, ra, dec, mean_mjd)
        # If Gaia matched a proper-motion star, center forced photometry at its
        # propagated SPHEREx-epoch coordinate!
        phot_ra = gmatch.ra_obs if gmatch else ra
        phot_dec = gmatch.dec_obs if gmatch else dec

        stack_img, stack_var = fp.deep_stack(coadds)
        snr = fp.aperture_snr(stack_img, stack_var, grid, phot_ra, phot_dec, search_px=1)
        xmatch = self.gaia.classify(gmatch, round(snr, 2) if snr is not None else None, gaia_res.error)
        wmatch = (
            self.wise.best_match(wise_res, ra, dec, mean_mjd, spherex_ra=phot_ra, spherex_dec=phot_dec)
            if wise_res is not None
            else None
        )
        xmatch.wise_match = wmatch
        if xmatch.classification == "candidate" and wmatch is not None:
            extras = [f"Matched {wmatch.catalog} {wmatch.source_id} (sep={wmatch.separation_arcsec:.1f}\")"]
            if wmatch.w1_w2_color is not None:
                extras.append(f"W1-W2={wmatch.w1_w2_color:.2f} mag")
            if wmatch.baseline_pm_arcsec_yr is not None:
                extras.append(f"10yr baseline PM={wmatch.baseline_pm_arcsec_yr:.2f}\"/yr")
            xmatch.reason = f"{xmatch.reason} | {'; '.join(extras)}"
        xmatch.overlay = self.gaia.overlay(gaia_res, mean_mjd, grid.world_to_png)

        # Difference imaging (Epoch_last - Epoch_0) + Deep Stack + Dipole/Variable detection
        c_a, c_b = coadds[0], coadds[-1]
        diff_img, snr_map = fp.compute_difference_map(c_a.image, c_a.variance, c_b.image, c_b.variance)
        diff_png = fp.render_difference_png(snr_map, self.s.png_size, self.s.png_resample)
        deep_png = fp.render_image(stack_img, vmin, vmax, "asinh", self.s.png_size, self.s.png_resample)
        anomalies = (
            fp.detect_diff_anomalies(diff_img, snr_map, grid, c_a.mjd_mid, c_b.mjd_mid, coadds=coadds)
            if len(coadds) >= 2
            else []
        )
        dipoles = [a.snr for a in anomalies if a.anomaly_type == "moving_dipole"]
        p_diff = f"targets/{key}/diff_e{c_a.index:02d}_e{c_b.index:02d}.png"
        p_deep = f"targets/{key}/deep_stack.png"
        uploads.extend(
            [
                (p_diff, diff_png, "image/png"),
                (p_deep, deep_png, "image/png"),
            ]
        )
        diff_meta = {
            "epoch_a": c_a.index,
            "epoch_b": c_b.index,
            "mjd_a": round(c_a.mjd_mid, 5),
            "mjd_b": round(c_b.mjd_mid, 5),
            "baseline_days": round(c_b.mjd_mid - c_a.mjd_mid, 2),
            "diff_path": p_diff,
            "deep_path": p_deep,
            "max_dipole_snr": round(max(dipoles), 2) if dipoles else 0.0,
            "anomalies": anomalies,
        }

        samples = [s for fd in all_loaded if (s := fp.extract_sample(fd, phot_ra, phot_dec)) is not None]
        points = fp.bin_spectrum(samples)

        # Fast spectrum on Cold Start: reuse already-loaded `samples` in memory without re-opening FITS files;
        # full directory scan / extra detector downloads are reserved for /target/spectrum or fast=False
        if fast:
            spectrum_data = extract_spectrum(
                [],
                phot_ra,
                phot_dec,
                wavecal_dir=self.cache.root / "wavecal",
                preloaded_samples=samples,
            )
        else:
            fits_set = {p for p in path_by_key.values() if p is not None}
            for cov_p in find_covering_fits(self.cache.root / "raw", phot_ra, phot_dec):
                fits_set.add(cov_p)
            spectrum_data = extract_spectrum(
                sorted(fits_set),
                phot_ra,
                phot_dec,
                wavecal_dir=self.cache.root / "wavecal",
            )
        # Preserve the requested query coordinates in the header block
        spectrum_data["target"] = {"ra": round(ra, 6), "dec": round(dec, 6)}

        display = DisplayInfo(
            width=self.s.png_size,
            height=self.s.png_size,
            grid_px=grid.n,
            pixel_scale_arcsec=round(grid.pixscale_deg * 3600.0 * grid.n / self.s.png_size, 4),
            wcs=grid.png_wcs(),
            stretch=StretchInfo(stretch_type="asinh", vmin=round(vmin, 5), vmax=round(vmax, 5)),
            data_range=DataRange(
                p0_5=round(data_range.p0_5, 5),
                p50=round(data_range.p50, 5),
                p99_5=round(data_range.p99_5, 5),
                min=round(data_range.min, 5),
                max=round(data_range.max, 5),
            ),
        )
        return {
            "uploads": uploads,
            "epochs_info": epochs_info,
            "grid": grid,
            "crossmatch": xmatch,
            "display": display,
            "samples": samples,
            "spectrum_points": points,
            "spectrum_data": spectrum_data,
            "diff_meta": diff_meta,
            "n_loaded": len(all_loaded),
            "phot_ra": phot_ra,
            "phot_dec": phot_dec,
        }

    async def _complete_spectrum_bg(
        self,
        key: str,
        ra: float,
        dec: float,
        initial_samples: list[fp.SpectralSample],
        initial_frames_count: int,
        extra_frames: list[Frame],
        row_snapshot: dict[str, Any],
    ) -> None:
        """Phase 2: download tiny 60" cutouts of the remaining exposures to fill
        all 102 LVF spectral channels, then mark the target `ready`."""
        try:
            known = await self.db.get_cached_fits_urls([f.key for f in extra_frames])
            await self.irsa.resolve_fits_urls(extra_frames, known)
            await self.db.save_observations(extra_frames)

            size_deg = self.s.spectrum_cutout_arcsec / 3600.0
            downloaded = await self.irsa.fetch_many(extra_frames, ra, dec, size_deg)

            def _extract_extra() -> tuple[list[fp.SpectralSample], int]:
                new_samples = list(initial_samples)
                n_extra = 0
                for f, p in downloaded:
                    wpath = self.cache.wavecal_path(self.irsa.wavecal_key(f))
                    try:
                        fd = fp.load_frame(p, f.key, f.detector, f.t_mid, wpath)
                        n_extra += 1
                        s = fp.extract_sample(fd, ra, dec)
                        if s is not None:
                            new_samples.append(s)
                    except Exception:
                        continue
                return new_samples, n_extra

            all_samples, n_extra = await asyncio.to_thread(_extract_extra)
            points = await asyncio.to_thread(fp.bin_spectrum, all_samples)
            cov_files = await asyncio.to_thread(find_covering_fits, self.cache.root / "raw", ra, dec)
            spec_data = await asyncio.to_thread(
                extract_spectrum, cov_files, ra, dec, self.cache.root / "wavecal"
            )
            spec_data["target"] = {"ra": round(float(row_snapshot["ra"]), 6), "dec": round(float(row_snapshot["dec"]), 6)}
            total_frames = initial_frames_count + n_extra
            spec = Spectrum(
                complete=True,
                n_frames_used=total_frames,
                n_channels_filled=sum(1 for p in points if p.flux is not None),
                points=points,
            )
            await self.db.upsert_target(
                {
                    **row_snapshot,
                    "status": "ready",
                    "progress": {"stage": "Complete", "percent": 100},
                    "spectrum": spec.model_dump(mode="json"),
                    "spectrum_data": spec_data,
                }
            )
            log.info(
                "Target %s phase-2 spectrum complete (%d/102 channels from %d frames)",
                key,
                spec.n_channels_filled,
                total_frames,
            )
        except Exception as e:
            log.warning("Phase-2 spectrum completion failed for %s: %s", key, e)
            await self.db.upsert_target(
                {
                    **row_snapshot,
                    "status": "ready",
                    "progress": {"stage": "Complete (partial spectrum)", "percent": 100},
                }
            )

    # ============================================================== helpers ==
    def _is_stale(self, row: dict[str, Any]) -> bool:
        ts = row.get("updated_at") or row.get("created_at")
        if not ts:
            return False
        dt = datetime.fromisoformat(str(ts).replace("Z", "+00:00")) if isinstance(ts, str) else ts
        age_days = (datetime.now(timezone.utc) - dt).total_seconds() / 86400.0
        return age_days > self.s.cache_max_age_days

    def _is_lock_expired(self, row: dict[str, Any]) -> bool:
        ts = row.get("updated_at") or row.get("created_at")
        if not ts:
            return True
        dt = datetime.fromisoformat(str(ts).replace("Z", "+00:00")) if isinstance(ts, str) else ts
        return (datetime.now(timezone.utc) - dt).total_seconds() > self.s.processing_lock_s

    def _row_to_response(self, row: dict[str, Any], hit: bool, stale: bool) -> TargetResponse:
        status = row.get("status", "processing")
        poll_ms = 1500 if status == "processing" else (3000 if status == "partial" else None)
        return TargetResponse(
            target_key=row["target_key"],
            status=status,
            ra=float(row["ra"]),
            dec=float(row["dec"]),
            fov=float(row["fov"]),
            band=row.get("band", "all"),
            progress=Progress.model_validate(row["progress"]) if row.get("progress") else None,
            poll_after_ms=poll_ms,
            cache=CacheInfo(
                hit=hit,
                stale=stale,
                created_at=row.get("created_at"),
                updated_at=row.get("updated_at"),
            ),
            pipeline_version=row.get("pipeline_version", self.s.pipeline_version),
            n_observations=int(row.get("n_observations", 0)),
            display=DisplayInfo.model_validate(row["display"]) if row.get("display") else None,
            epochs=[EpochOut.model_validate(e) for e in (row.get("epochs") or [])],
            crossmatch=row.get("crossmatch"),
            spectrum=Spectrum.model_validate(row["spectrum"]) if row.get("spectrum") else None,
            spectrum_data=SpectrumResponse.model_validate(row["spectrum_data"]) if row.get("spectrum_data") else None,
            difference=DifferenceResult.model_validate(row["difference"]) if row.get("difference") else None,
            error=row.get("error"),
        )

    # ============================================= standalone difference ===
    async def get_difference(
        self,
        ra: float,
        dec: float,
        fov: float = 0.1,
        band: BandName = "all",
        epoch_a: int = 0,
        epoch_b: int = -1,
    ) -> DifferenceResult:
        """Compute on-demand difference image (`Epoch_B - Epoch_A`) and detect
        moving dipoles & variable stars between any two epochs."""
        ra_n, dec_n, fov_n = normalize_ra(ra), normalize_dec(dec), normalize_fov(fov)
        key = target_key(ra_n, dec_n, fov_n, band)

        cached = await self.db.get_target(ra_n, dec_n, fov_n, band)
        if not cached or cached.get("status") not in ("ready", "partial"):
            await self.get_target(ra_n, dec_n, fov_n, band, sync=True)
            cached = await self.db.get_target(ra_n, dec_n, fov_n, band)

        epochs_meta = (cached or {}).get("epochs") or []
        if not epochs_meta:
            raise RuntimeError("No epochs available for difference imaging at these coordinates")

        idx_list = [int(e["index"]) for e in epochs_meta]
        idx_a = epoch_a if epoch_a in idx_list else idx_list[0]
        idx_b = (idx_list[-1] if epoch_b < 0 else epoch_b) if (epoch_b < 0 or epoch_b in idx_list) else idx_list[-1]

        # Fast path: if default (first, last) pair is already cached on the target row, return it
        if cached and cached.get("difference"):
            d_cached = DifferenceResult.model_validate(cached["difference"])
            if d_cached.epoch_a == idx_a and d_cached.epoch_b == idx_b:
                return d_cached

        async def _load_ep(ep_idx: int) -> tuple[ Any, Any, float]:
            local_p = self.cache.get(self.cache.coadd_path(key, ep_idx))
            if local_p is not None:
                return await asyncio.to_thread(fp.read_coadd_with_var, local_p)
            ep_obj = next((e for e in epochs_meta if int(e["index"]) == ep_idx), None)
            if ep_obj and ep_obj.get("coadd_fits_url"):
                raw = await self.db.download_storage_bytes(ep_obj["coadd_fits_url"])
                if raw:
                    self.cache.write_atomic(self.cache.coadd_path(key, ep_idx), raw)
                    return await asyncio.to_thread(fp.read_coadd_with_var, raw)
            raise RuntimeError(f"Co-add FITS for epoch {ep_idx} not found")

        (img_a, var_a, mjd_a), (img_b, var_b, mjd_b) = await asyncio.gather(
            _load_ep(idx_a), _load_ep(idx_b)
        )
        all_ep_coadds: list[fp.EpochCoadd] | None = None
        if len(idx_list) >= 3 and idx_a != idx_b:
            loaded_map = {idx_a: (img_a, var_a, mjd_a), idx_b: (img_b, var_b, mjd_b)}
            other_idxs = [i for i in idx_list if i not in loaded_map]
            other_data = await asyncio.gather(*(_load_ep(i) for i in other_idxs))
            for i, d in zip(other_idxs, other_data):
                loaded_map[i] = d
            all_ep_coadds = [
                fp.EpochCoadd(
                    index=i,
                    image=loaded_map[i][0],
                    variance=loaded_map[i][1],
                    flags=loaded_map[i][0].astype("int64") * 0,
                    coverage=loaded_map[i][0].astype("int16") * 0,
                    n_frames=1,
                    detectors=[],
                    mjd_start=loaded_map[i][2],
                    mjd_end=loaded_map[i][2],
                    mjd_mid=loaded_map[i][2],
                )
                for i in idx_list
            ]

        def _compute_cpu() -> DifferenceResult:
            grid = fp.make_grid(ra_n, dec_n, fov_n, self.s.spherex_pixel_scale_arcsec, self.s.png_size)
            diff_img, snr_map = fp.compute_difference_map(img_a, var_a, img_b, var_b)
            diff_png = fp.render_difference_png(snr_map, self.s.png_size, self.s.png_resample)
            anomalies = (
                fp.detect_diff_anomalies(diff_img, snr_map, grid, mjd_a, mjd_b, coadds=all_ep_coadds)
                if idx_a != idx_b
                else []
            )
            dipoles = [a.snr for a in anomalies if a.anomaly_type == "moving_dipole"]
            deep_url = (cached.get("difference") or {}).get("deep_stack_url") if cached else None
            return DifferenceResult(
                epoch_a=idx_a,
                epoch_b=idx_b,
                mjd_a=round(mjd_a, 5),
                mjd_b=round(mjd_b, 5),
                baseline_days=round(mjd_b - mjd_a, 2),
                diff_image_url=fp.png_data_url(diff_png),
                deep_stack_url=deep_url,
                max_dipole_snr=round(max(dipoles), 2) if dipoles else 0.0,
                anomalies=anomalies,
            )

        return await asyncio.to_thread(_compute_cpu)

    # ================================================ standalone spectrum ===
    async def get_spectrum(
        self,
        ra: float,
        dec: float,
        fov: float = 0.1,
        force_refresh: bool = False,
    ) -> SpectrumResponse:
        """Fast standalone endpoint logic for `GET /api/v1/target/spectrum`:
        1. Check Supabase / local DB cache (`spectrum_data` JSONB).
        2. Check `/tmp/fits_cache/raw/` for already-downloaded FITS cutouts covering `(ra, dec)`.
        3. Fallback: query IRSA and stream small cutouts around `(ra, dec)`, then extract.
        Never raises 500 on missing WCS-WAVE or partial channel coverage.
        """
        ra_n, dec_n = normalize_ra(ra), normalize_dec(dec)

        # If GET /api/v1/target is ALREADY in flight for these coordinates (e.g. frontend fired
        # /target and /target/spectrum concurrently in Promise.all / useEffect), coalesce onto
        # the in-flight pipeline task instead of launching a competing 24-frame IRSA download!
        if not force_refresh:
            coord_prefix = f"{ra_n:.4f}_{dec_n:+.4f}_"
            inflight_task = next(
                (t for k, t in self._inflight.items() if k.startswith(coord_prefix)),
                None,
            )
            if inflight_task is not None:
                row = await asyncio.wait_for(asyncio.shield(inflight_task), timeout=self.s.sync_timeout_s)
                if row and row.get("spectrum_data"):
                    return SpectrumResponse.model_validate(row["spectrum_data"])

            cached = await self.db.get_target_by_coords(ra_n, dec_n)
            if cached and cached.get("spectrum_data"):
                return SpectrumResponse.model_validate(cached["spectrum_data"])

        # Step 2: Check local FITS files in /tmp/fits_cache/raw/
        local_files = await asyncio.to_thread(find_covering_fits, self.cache.root / "raw", ra_n, dec_n)

        # Step 3: If no local files (or force_refresh), download cutouts from IRSA
        if not local_files or force_refresh:
            try:
                frames = await self.irsa.query_frames(ra_n, dec_n)
                if frames:
                    epochs = self.irsa.group_epochs(frames)
                    selected: list[Frame] = []
                    for ep in epochs[-2:]:
                        selected.extend(
                            self.irsa.select_display_frames(ep, BAND_GROUPS["all"], ra_n, dec_n, max_per_epoch=2)
                        )
                    known = await self.db.get_cached_fits_urls([f.key for f in selected])
                    cutout_deg = max(self.s.spectrum_cutout_arcsec / 3600.0, min(fov * 1.1, 0.15))
                    if (
                        hasattr(self.irsa, "resolve_and_fetch_many")
                        and getattr(self.irsa.fetch_many, "__func__", None) is IrsaService.fetch_many
                    ):
                        downloaded = await self.irsa.resolve_and_fetch_many(
                            selected, known, ra_n, dec_n, cutout_deg
                        )
                    else:
                        await self.irsa.resolve_fits_urls(selected, known)
                        downloaded = await self.irsa.fetch_many(selected, ra_n, dec_n, cutout_deg)
                    local_files = [p for _, p in downloaded]
                    # Also include any additional covering files present in raw/
                    extra_cov = await asyncio.to_thread(find_covering_fits, self.cache.root / "raw", ra_n, dec_n)
                    local_files = sorted(set(local_files).union(extra_cov))
            except Exception as e:
                log.warning("IRSA fallback inside get_spectrum(%.4f, %+.4f) failed: %s", ra_n, dec_n, e)

        spec_dict = await asyncio.to_thread(
            extract_spectrum,
            local_files,
            ra_n,
            dec_n,
            self.cache.root / "wavecal",
        )
        # Persist into `targets` if a row for (ra, dec) already exists
        existing = await self.db.get_target(ra_n, dec_n, normalize_fov(fov), "all")
        if existing:
            await self.db.upsert_target({**existing, "spectrum_data": spec_dict})

        return SpectrumResponse.model_validate(spec_dict)
