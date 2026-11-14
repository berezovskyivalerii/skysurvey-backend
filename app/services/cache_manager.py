"""Local disk cache under /tmp/fits_cache with TTL + LRU size cap.

Layout
------
/tmp/fits_cache/
    raw/<sha1>.fits            streamed SPHEREx cutouts (IMAGE/FLAGS/VARIANCE/ZODI)
    wavecal/D<n>.fits          WCS-WAVE lookup table per detector (pinned, never evicted)
    targets/<target_key>/      epoch co-adds (float32 FITS, used by /restretch)

"Last access" == mtime: every read goes through `touch()`, so the sweeper is
both a TTL (age > ttl) and an LRU (oldest-first until under the byte cap).
On Render/Railway the disk is ephemeral and per-instance: anything here can
vanish, so callers must treat a miss as normal (Supabase Storage is the
durable fallback for co-adds).
"""
from __future__ import annotations

import asyncio
import hashlib
import logging
import os
import tempfile
import time
from pathlib import Path

log = logging.getLogger(__name__)

PINNED_DIRS = ("wavecal",)


class FitsCache:
    def __init__(self, root: Path, ttl_s: int, max_bytes: int):
        self.root = Path(root)
        self.ttl_s = ttl_s
        self.max_bytes = max_bytes
        for sub in ("raw", "wavecal", "targets"):
            (self.root / sub).mkdir(parents=True, exist_ok=True)
        self._task: asyncio.Task | None = None

    # ------------------------------------------------------------ paths ----
    @staticmethod
    def url_hash(url: str) -> str:
        return hashlib.sha1(url.encode()).hexdigest()

    def raw_path(self, url: str) -> Path:
        return self.root / "raw" / f"{self.url_hash(url)}.fits"

    def wavecal_path(self, key: str) -> Path:
        """WCS-WAVE table per (collection, detector), e.g. 'spherex_dr1_D6'."""
        return self.root / "wavecal" / f"{key}.fits"

    def target_dir(self, key: str) -> Path:
        d = self.root / "targets" / key
        d.mkdir(parents=True, exist_ok=True)
        return d

    def coadd_path(self, key: str, epoch_index: int) -> Path:
        return self.target_dir(key) / f"epoch_{epoch_index:02d}_coadd.fits"

    # ---------------------------------------------------------- access -----
    def get(self, path: Path) -> Path | None:
        if path.is_file():
            self.touch(path)
            return path
        return None

    @staticmethod
    def touch(path: Path) -> None:
        try:
            os.utime(path, None)
        except FileNotFoundError:
            pass

    @staticmethod
    def write_atomic(path: Path, data: bytes) -> Path:
        """Write via temp file + rename so readers never see half a FITS file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".part")
        try:
            with os.fdopen(fd, "wb") as f:
                f.write(data)
            os.replace(tmp, path)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise
        return path

    # ---------------------------------------------------------- sweeping ---
    def sweep(self) -> dict[str, int]:
        now = time.time()
        files: list[tuple[float, int, Path]] = []
        removed = freed = 0
        for p in self.root.rglob("*"):
            if not p.is_file() or p.relative_to(self.root).parts[0] in PINNED_DIRS:
                continue
            try:
                st = p.stat()
            except FileNotFoundError:
                continue
            if p.suffix == ".part" and now - st.st_mtime > 300:
                p.unlink(missing_ok=True)
                continue
            if now - st.st_mtime > self.ttl_s:
                p.unlink(missing_ok=True)
                removed += 1
                freed += st.st_size
            else:
                files.append((st.st_mtime, st.st_size, p))
        total = sum(s for _, s, _ in files)
        if total > self.max_bytes:
            for _, size, p in sorted(files):  # oldest first = LRU
                p.unlink(missing_ok=True)
                removed += 1
                freed += size
                total -= size
                if total <= self.max_bytes * 0.9:  # hysteresis
                    break
        for d in sorted((self.root / "targets").glob("*"), reverse=True):
            if d.is_dir() and not any(d.iterdir()):
                d.rmdir()
        if removed:
            log.info("cache sweep: removed %d files, freed %.1f MB, now %.1f MB", removed, freed / 1e6, total / 1e6)
        return {"removed": removed, "freed_bytes": freed, "total_bytes": total}

    def usage_bytes(self) -> int:
        return sum(p.stat().st_size for p in self.root.rglob("*") if p.is_file())

    async def _loop(self, interval_s: int) -> None:
        while True:
            try:
                await asyncio.to_thread(self.sweep)
            except Exception:  # never let the janitor die
                log.exception("cache sweep failed")
            await asyncio.sleep(interval_s)

    def start(self, interval_s: int) -> None:
        if self._task is None:
            self._task = asyncio.create_task(self._loop(interval_s), name="fits-cache-sweeper")

    async def stop(self) -> None:
        if self._task:
            self._task.cancel()
            try:
                await self._task
            except asyncio.CancelledError:
                pass
            self._task = None
