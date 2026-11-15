"""FastAPI application entrypoint with shared connection pool and disk cache janitor."""
from __future__ import annotations

import logging
import os
from contextlib import asynccontextmanager

import httpx
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from app.api.v1.targets import router as targets_router
from app.config import Settings, get_settings
from app.services.cache_manager import FitsCache
from app.services.gaia_service import GaiaService
from app.services.irsa_service import IrsaService
from app.services.pipeline import PipelineService
from app.services.solar_system_service import SolarSystemService
from app.services.supabase_service import SupabaseService
from app.services.wise_service import WiseService

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")


def create_app(settings: Settings | None = None) -> FastAPI:
    cfg = settings or get_settings()
    level = getattr(logging, cfg.log_level.upper(), logging.INFO)
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
        force=True,
    )
    logging.getLogger("app").setLevel(level)
    logging.getLogger("httpx").setLevel(logging.WARNING)

    cache = FitsCache(
        root=cfg.cache_dir,
        ttl_s=cfg.cache_ttl_s,
        max_bytes=cfg.cache_max_bytes,
    )

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        limits = httpx.Limits(
            max_connections=cfg.http_max_connections,
            max_keepalive_connections=cfg.http_max_connections // 2,
        )
        timeout = httpx.Timeout(cfg.http_timeout_s, connect=10.0)
        transport = httpx.AsyncHTTPTransport(local_address="0.0.0.0", limits=limits)
        async with httpx.AsyncClient(
            transport=transport,
            limits=limits,
            timeout=timeout,
            follow_redirects=True,
        ) as http:
            irsa = IrsaService(cfg, http, cache)
            gaia = GaiaService(cfg, http)
            wise = WiseService(cfg, http)
            solar_system = SolarSystemService(cfg, http)
            db = SupabaseService(cfg, http)
            pipeline = PipelineService(cfg, cache, irsa, gaia, db, wise=wise, solar_system=solar_system)

            app.state.settings = cfg
            app.state.cache = cache
            app.state.irsa = irsa
            app.state.gaia = gaia
            app.state.wise = wise
            app.state.solar_system = solar_system
            app.state.db = db
            app.state.pipeline = pipeline

            cache.start(cfg.cache_sweep_interval_s)
            try:
                yield
            finally:
                await cache.stop()

    app = FastAPI(
        title=cfg.app_name,
        version=cfg.pipeline_version,
        description=(
            "High-performance Astro-informatics backend for NASA SPHEREx time-domain sky visualization "
            "(WCS epoch alignment, 3-layer PNG rendering, Gaia DR3 proper-motion cross-validation, "
            "and 102-channel LVF spectrophotometry)."
        ),
        lifespan=lifespan,
    )

    app.add_middleware(
        CORSMiddleware,
        allow_origins=cfg.cors_origins,
        allow_credentials=False,
        allow_methods=["*"],
        allow_headers=["*"],
        expose_headers=["*"],
    )

    local_storage = cfg.cache_dir / "storage"
    local_storage.mkdir(parents=True, exist_ok=True)
    app.mount("/static/storage", StaticFiles(directory=str(local_storage)), name="static-storage")

    app.include_router(targets_router)

    @app.get("/", tags=["system"])
    async def root():
        return {
            "service": cfg.app_name,
            "status": "ok",
            "version": cfg.pipeline_version,
            "docs": "/docs",
            "health": "/health",
            "presets": "/api/v1/presets",
            "target": "/api/v1/target?ra=269.4485&dec=4.7394&fov=0.1",
        }

    @app.get("/health", tags=["system"])
    async def health():
        return {
            "status": "ok",
            "pipeline_version": cfg.pipeline_version,
            "supabase_enabled": bool(cfg.supabase_url and cfg.supabase_api_key),
            "cache_usage_mb": round(cache.usage_bytes() / 1e6, 2),
        }

    return app


app = create_app()
