"""Application settings (12-factor: everything comes from env / .env)."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8", extra="ignore")

    # --- service ---
    app_name: str = "SkySurvey SPHEREx API"
    pipeline_version: str = "1.0.0"  # bump -> every cached target is treated as stale
    public_base_url: str = ""  # optional base URL prefix for /static/storage
    cors_origins: list[str] = ["*"]
    log_level: str = "INFO"

    # --- Supabase ---
    supabase_url: str = ""
    supabase_key: str = ""  # sb_publishable_... (needs RLS policies, see sql/001_schema.sql)
    supabase_secret_key: str = ""  # optional sb_secret_..., preferred server-side if set
    supabase_bucket: str = "skysurvey"
    targets_table: str = "targets"
    observations_table: str = "observations"

    # --- IRSA / SPHEREx ---
    irsa_tap_url: str = "https://irsa.ipac.caltech.edu/TAP"
    spherex_obscore_table: str = "spherex.obscore"
    spherex_pixel_scale_arcsec: float = 6.15
    epoch_gap_days: float = 45.0  # a gap larger than this starts a new epoch
    epoch_max_span_days: float = 150.0  # protects ecliptic-pole (deep field) continuous coverage
    max_frames_per_detector_per_epoch: int = 1  # display coadd: 1 best-centered frame per detector/epoch
    max_frames_per_epoch: int = 2  # cap total display frames per epoch to 1-2 closest to (ra, dec)
    max_epochs: int = 4  # at most 4 epochs * 2 frames = 4-8 concurrent FITS streams on cold start
    spectrum_max_frames: int = 0  # 0 = finish in Phase 1 (ready immediately); >0 = extra background frames
    spectrum_cutout_arcsec: float = 60.0
    tap_max_rows: int = 5000

    # --- HTTP ---
    http_timeout_s: float = 60.0
    http_max_connections: int = 32
    download_concurrency: int = 8  # allows 4-8 display cutouts to stream in a single parallel wave
    upload_concurrency: int = 8
    http_retries: int = 2

    # --- Gaia ---
    gaia_table: str = "gaiadr3.gaia_source"
    gaia_match_radius_arcsec: float = 10.0  # ~1.6 SPHEREx px, AFTER propagating proper motion
    gaia_search_pad_arcsec: float = 180.0  # covers 10 yr of motion for >15"/yr stars
    gaia_overlay_max_sources: int = 500
    gaia_timeout_s: float = 4.5

    # --- WISE (CatWISE2020 / AllWISE) ---
    wise_table: str = "catwise_2020"
    allwise_table: str = "allwise_p3as_psd"
    wise_search_radius_deg: float = 0.02
    wise_match_radius_arcsec: float = 10.0
    wise_timeout_s: float = 4.5

    # --- Solar System (IMCCE SkyBoT) ---
    skybot_url: str = "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"
    skybot_timeout_s: float = 2.0
    skybot_enabled: bool = True

    # --- rendering ---
    png_size: int = 512
    png_resample: str = "nearest"  # nearest keeps SPHEREx pixels crisp (WiseView style)
    detection_snr: float = 5.0

    # --- local cache ---
    cache_dir: Path = Path("/tmp/fits_cache")
    cache_ttl_s: int = 12 * 3600  # 12h TTL (still bounded by 500 MB LRU cap)
    cache_max_bytes: int = 500 * 1024 * 1024
    cache_sweep_interval_s: int = 60

    # --- DB cache policy ---
    coord_match_tol_arcsec: float = 2.0
    cache_max_age_days: float = 7.0  # stale-while-revalidate (new SPHEREx releases)
    processing_lock_s: int = 600  # a "processing" row older than this is considered dead

    # --- DEV ---
    sync_timeout_s: float = Field(default=120.0, description="max wait for ?sync=true")

    @property
    def supabase_api_key(self) -> str:
        return self.supabase_secret_key or self.supabase_key


@lru_cache
def get_settings() -> Settings:
    return Settings()
