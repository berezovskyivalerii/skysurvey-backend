"""Warmup script to pre-compute and cache all PRESET_TARGETS into Supabase / local storage."""

from __future__ import annotations

import argparse
import asyncio
import logging
import sys
import time
from pathlib import Path
from typing import Any

import httpx

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from app.api.v1.targets import PRESET_TARGETS
from app.config import Settings, get_settings
from app.services.cache_manager import FitsCache
from app.services.gaia_service import GaiaService
from app.services.irsa_service import IrsaService
from app.services.pipeline import PipelineService
from app.services.solar_system_service import SolarSystemService
from app.services.supabase_service import SupabaseService
from app.services.wise_service import WiseService

logger = logging.getLogger("warmup_presets")


def _select_presets(target_id: str | None = None) -> list[dict[str, Any]]:
    if not target_id:
        return list(PRESET_TARGETS)
    selected = [p for p in PRESET_TARGETS if p["id"] == target_id]
    if not selected:
        available = ", ".join(p["id"] for p in PRESET_TARGETS)
        raise ValueError(f"Unknown preset target id {target_id!r}. Available: {available}")
    return selected


def _build_pipeline(cfg: Settings, http: httpx.AsyncClient) -> PipelineService:
    cache = FitsCache(
        root=cfg.cache_dir,
        ttl_s=cfg.cache_ttl_s,
        max_bytes=cfg.cache_max_bytes,
    )
    irsa = IrsaService(cfg, http, cache)
    gaia = GaiaService(cfg, http)
    wise = WiseService(cfg, http)
    solar_system = SolarSystemService(cfg, http)
    db = SupabaseService(cfg, http)
    return PipelineService(
        cfg,
        cache,
        irsa,
        gaia,
        db,
        wise=wise,
        solar_system=solar_system,
    )


async def warmup_all(
    force: bool = False,
    target_id: str | None = None,
    dry_run: bool = False,
    pipeline: PipelineService | None = None,
    settings: Settings | None = None,
) -> list[dict[str, Any]]:
    """Run synchronous full pipeline warmup for selected preset targets."""
    presets = _select_presets(target_id)
    results: list[dict[str, Any]] = []

    if dry_run:
        for preset in presets:
            results.append(
                {
                    "id": preset["id"],
                    "name": preset["name"],
                    "target_key": None,
                    "status": "dry_run",
                    "epochs": 0,
                    "classification": None,
                    "cache_hit": False,
                    "elapsed_s": 0.0,
                    "error": None,
                    "ok": True,
                }
            )
        return results

    async def _run_with_pipeline(active_pipeline: PipelineService) -> None:
        for preset in presets:
            t0 = time.perf_counter()
            logger.info(
                "Warming up %s (%s) ra=%.5f dec=%.5f fov=%.3f band=%s force=%s",
                preset["id"],
                preset["name"],
                preset["ra"],
                preset["dec"],
                preset["fov"],
                preset["band"],
                force,
            )
            try:
                resp = await active_pipeline.get_target(
                    ra=preset["ra"],
                    dec=preset["dec"],
                    fov=preset["fov"],
                    band=preset["band"],
                    sync=True,
                    force_refresh=force,
                )
                elapsed = round(time.perf_counter() - t0, 2)
                ok = resp.status != "error"
                results.append(
                    {
                        "id": preset["id"],
                        "name": preset["name"],
                        "target_key": resp.target_key,
                        "status": resp.status,
                        "epochs": len(resp.epochs),
                        "classification": resp.crossmatch.classification if resp.crossmatch else None,
                        "cache_hit": resp.cache.hit,
                        "elapsed_s": elapsed,
                        "error": resp.error,
                        "ok": ok,
                    }
                )
            except Exception as exc:
                elapsed = round(time.perf_counter() - t0, 2)
                logger.exception("Warmup failed for %s", preset["id"])
                results.append(
                    {
                        "id": preset["id"],
                        "name": preset["name"],
                        "target_key": None,
                        "status": "error",
                        "epochs": 0,
                        "classification": None,
                        "cache_hit": False,
                        "elapsed_s": elapsed,
                        "error": str(exc),
                        "ok": False,
                    }
                )

    if pipeline is not None:
        await _run_with_pipeline(pipeline)
    else:
        cfg = settings or get_settings()
        timeout = httpx.Timeout(connect=15.0, read=60.0, write=30.0, pool=10.0)
        async with httpx.AsyncClient(timeout=timeout, follow_redirects=True) as http:
            owned_pipeline = _build_pipeline(cfg, http)
            await _run_with_pipeline(owned_pipeline)

    return results


def _format_summary(results: list[dict[str, Any]]) -> str:
    header = f"{'ID':<16} {'STATUS':<10} {'EPOCHS':>6} {'HIT':<5} {'TIME(s)':>7} {'CLASS':<14} {'ERROR'}"
    lines = [header, "-" * len(header)]
    for item in results:
        lines.append(
            f"{item['id']:<16} "
            f"{str(item['status']):<10} "
            f"{int(item['epochs']):>6} "
            f"{('yes' if item['cache_hit'] else 'no'):<5} "
            f"{float(item['elapsed_s']):>7.2f} "
            f"{str(item['classification'] or '-'):<14} "
            f"{item['error'] or ''}"
        )
    return "\n".join(lines)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pre-compute and cache SPHEREx preset targets into Supabase/local cache."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Reprocess targets even if already cached (sets force_refresh=True).",
    )
    parser.add_argument(
        "--target",
        type=str,
        default=None,
        help="Warm up only a single preset id (e.g. wise_0855, barnards_star).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List selected presets without calling IRSA/Gaia/WISE/Supabase.",
    )
    return parser.parse_args(argv)


async def _async_main(args: argparse.Namespace) -> int:
    results = await warmup_all(
        force=args.force,
        target_id=args.target,
        dry_run=args.dry_run,
    )
    print(_format_summary(results))
    return 0 if all(item["ok"] for item in results) else 1


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    args = parse_args(argv)
    try:
        return asyncio.run(_async_main(args))
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
