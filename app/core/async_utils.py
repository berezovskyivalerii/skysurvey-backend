"""Async concurrency helpers (race-to-first-success without cancelling siblings on error)."""
from __future__ import annotations

import asyncio
import logging
import time
from typing import TypeVar

log = logging.getLogger(__name__)

T = TypeVar("T")


async def race_first_success(
    tasks: list[asyncio.Task[T]],
    timeout: float | None = None,
) -> tuple[T, str]:
    """Wait for the FIRST task that returns without raising an exception.

    Unlike a naive `FIRST_COMPLETED` loop, if a task fails (e.g. HTTP 502),
    the error is logged at DEBUG level and the remaining `pending` tasks keep
    running until one succeeds, all fail, or `timeout` expires.
    """
    if not tasks:
        raise ValueError("race_first_success requires at least one task")

    pending_set: set[asyncio.Task[T]] = set(tasks)
    last_err: BaseException | None = None
    deadline = (time.perf_counter() + timeout) if timeout is not None else None

    try:
        while pending_set:
            wait_timeout: float | None = None
            if deadline is not None:
                wait_timeout = deadline - time.perf_counter()
                if wait_timeout <= 0:
                    raise TimeoutError(f"Race timed out after {timeout:.1f}s")

            done, pending_set = await asyncio.wait(
                pending_set,
                timeout=wait_timeout,
                return_when=asyncio.FIRST_COMPLETED,
            )
            if not done:
                raise TimeoutError(f"Race timed out after {timeout:.1f}s")

            for task in done:
                if task.cancelled():
                    continue
                exc = task.exception()
                if exc is not None:
                    last_err = exc
                    log.debug("Race branch %s failed (%r); continuing remaining branches", task.get_name(), exc)
                    continue
                result = task.result()
                for rem in pending_set:
                    rem.cancel()
                return result, task.get_name()

        if last_err is not None:
            raise last_err
        raise RuntimeError("All race tasks were cancelled")
    finally:
        for rem in pending_set:
            if not rem.done():
                rem.cancel()
        if pending_set:
            await asyncio.gather(*pending_set, return_exceptions=True)
