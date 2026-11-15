"""FastAPI v1 endpoints for `/api/v1/target`, `/api/v1/target/restretch`, and `/api/v1/presets`."""
from __future__ import annotations

from typing import Annotated

from fastapi import APIRouter, BackgroundTasks, Depends, HTTPException, Query, Request, Response

from app.schemas.spectrum import SpectrumResponse
from app.schemas.target import (
    BandName,
    DifferenceResult,
    ErrorResponse,
    RestretchRequest,
    RestretchResponse,
    TargetResponse,
)
from app.services.pipeline import PipelineService

router = APIRouter(prefix="/api/v1", tags=["targets"])


def get_pipeline(request: Request) -> PipelineService:
    return request.app.state.pipeline


PRESET_TARGETS = [
    {
        "id": "wise_0855",
        "name": "WISE J085510.83-071442.5",
        "category": "brown_dwarf_candidate",
        "ra": 133.795129,
        "dec": -7.245147,
        "fov": 0.1,
        "band": "red",
        "description": (
            "Coldest known brown dwarf (~250 K, Y4, d=2.28 pc, PM=8.1\"/yr). "
            "Completely invisible in Gaia DR3 optical -> classified as `candidate`, "
            "bright in SPHEREx D5-D6 (4-5 um CO/CH4 window)."
        ),
    },
    {
        "id": "barnards_star",
        "name": "Barnard's Star (2025.8 epoch)",
        "category": "high_proper_motion_star",
        "ra": 269.448500,
        "dec": 4.739400,
        "fov": 0.1,
        "band": "all",
        "description": (
            "Fastest-moving star in the sky (10.36\"/yr, M4V red dwarf at 1.83 pc). "
            "Shifts by ~0.84 SPHEREx pixels per 6-month epoch."
        ),
    },
    {
        "id": "eps_indi_b",
        "name": "Epsilon Indi Ba/Bb",
        "category": "brown_dwarf_binary",
        "ra": 331.068900,
        "dec": -56.782500,
        "fov": 0.1,
        "band": "all",
        "description": "Binary T1 + T6 brown dwarf pair at 3.64 pc (proper motion 4.7\"/yr).",
    },
    {
        "id": "rw_aurigae",
        "name": "RW Aurigae",
        "category": "variable_star",
        "ra": 76.956510,
        "dec": 30.401440,
        "fov": 0.1,
        "band": "all",
        "description": "Classical T Tauri variable star exhibiting deep circumstellar-disk occultations.",
    },
    {
        "id": "r_draconis",
        "name": "R Draconis",
        "category": "variable_star",
        "ra": 248.179540,
        "dec": 66.754230,
        "fov": 0.1,
        "band": "all",
        "description": "Mira-type pulsating AGB variable star (period 245.6 days) with strong H2O/CO bands.",
    },
]


@router.get(
    "/target",
    response_model=TargetResponse,
    responses={502: {"model": ErrorResponse}, 504: {"model": ErrorResponse}},
    summary="Fetch multi-epoch aligned PNG layers, Gaia DR3 cross-match, and 102-channel spectrum",
)
async def get_target_endpoint(
    response: Response,
    background_tasks: BackgroundTasks,
    pipeline: Annotated[PipelineService, Depends(get_pipeline)],
    ra: Annotated[float, Query(ge=0.0, lt=360.0, description="Right Ascension (ICRS degrees)", examples=[133.7951])],
    dec: Annotated[float, Query(ge=-90.0, le=90.0, description="Declination (ICRS degrees)", examples=[-7.2451])],
    fov: Annotated[float, Query(gt=0.01, le=0.5, description="Field of view in degrees (e.g. 0.05, 0.1, 0.2)")] = 0.1,
    band: Annotated[BandName, Query(description="Detector filter group (all, blue, red, D1..D6)")] = "all",
    sync: Annotated[
        bool,
        Query(description="Wait for Phase-1 completion before returning. Set false for instant 202 polling."),
    ] = True,
    fast: Annotated[
        bool,
        Query(description="Fast Cold Start: compute spectrum only from display cutouts without scanning/downloading extra exposures."),
    ] = True,
    force_refresh: Annotated[bool, Query(description="Bypass DB cache and re-run the pipeline")] = False,
) -> TargetResponse:
    try:
        res = await pipeline.get_target(
            ra=ra,
            dec=dec,
            fov=fov,
            band=band,
            sync=sync,
            force_refresh=force_refresh,
            fast=fast,
            background_tasks=background_tasks,
        )
    except TimeoutError as e:
        raise HTTPException(status_code=504, detail=f"Upstream archive timed out: {e}") from e
    except Exception as e:
        raise HTTPException(status_code=502, detail=str(e)) from e

    if res.status == "processing":
        response.status_code = 202
    elif res.status == "failed":
        raise HTTPException(status_code=502, detail=res.error or "Pipeline failed")
    return res


@router.get(
    "/target/spectrum",
    response_model=SpectrumResponse,
    responses={502: {"model": ErrorResponse}},
    summary="Fetch sorted 102-channel SPHEREx spectral chart data & molecular absorption bands",
)
async def get_target_spectrum_endpoint(
    pipeline: Annotated[PipelineService, Depends(get_pipeline)],
    ra: Annotated[float, Query(ge=0.0, lt=360.0, description="Right Ascension (ICRS degrees)", examples=[133.7951])],
    dec: Annotated[float, Query(ge=-90.0, le=90.0, description="Declination (ICRS degrees)", examples=[-7.2451])],
    fov: Annotated[float, Query(gt=0.01, le=0.5, description="Fallback cutout field of view in degrees")] = 0.1,
    force_refresh: Annotated[bool, Query(description="Bypass DB cache and recompute spectrum")] = False,
) -> SpectrumResponse:
    return await pipeline.get_spectrum(ra=ra, dec=dec, fov=fov, force_refresh=force_refresh)


@router.get(
    "/target/diff",
    response_model=DifferenceResult,
    responses={502: {"model": ErrorResponse}},
    summary="Compute difference image (Epoch_B - Epoch_A) and auto-detect moving dipoles & variable stars",
)
async def get_target_diff_endpoint(
    pipeline: Annotated[PipelineService, Depends(get_pipeline)],
    ra: Annotated[float, Query(ge=0.0, lt=360.0, description="Right Ascension (ICRS degrees)", examples=[133.7951])],
    dec: Annotated[float, Query(ge=-90.0, le=90.0, description="Declination (ICRS degrees)", examples=[-7.2451])],
    fov: Annotated[float, Query(gt=0.01, le=0.5, description="Field of view in degrees")] = 0.1,
    band: Annotated[BandName, Query(description="Detector filter group")] = "all",
    epoch_a: Annotated[int, Query(ge=0, description="Index of reference (early) epoch")] = 0,
    epoch_b: Annotated[int, Query(description="Index of comparison (later) epoch; -1 = latest epoch")] = -1,
) -> DifferenceResult:
    try:
        return await pipeline.get_difference(
            ra=ra,
            dec=dec,
            fov=fov,
            band=band,
            epoch_a=epoch_a,
            epoch_b=epoch_b,
        )
    except Exception as e:
        raise HTTPException(status_code=502, detail=str(e)) from e


@router.post(
    "/target/restretch",
    response_model=RestretchResponse,
    responses={502: {"model": ErrorResponse}},
    summary="Fast (<50ms) contrast restretch of cached epoch co-adds",
)
async def restretch_target_endpoint(
    req: RestretchRequest,
    pipeline: Annotated[PipelineService, Depends(get_pipeline)],
) -> RestretchResponse:
    try:
        return await pipeline.restretch(req)
    except Exception as e:
        raise HTTPException(status_code=502, detail=str(e)) from e


@router.get("/presets", summary="Official NASA Space Apps Challenge benchmark targets")
async def list_presets() -> dict[str, list[dict]]:
    return {"presets": PRESET_TARGETS}
