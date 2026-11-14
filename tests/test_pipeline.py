"""End-to-end tests for FITS streaming, WCS alignment, 3-layer PNG rendering,
102-channel LVF spectrum extraction, Gaia proper-motion propagation, LRU cache,
and the FastAPI `/api/v1/target` + `/api/v1/target/restretch` endpoints."""
from __future__ import annotations

import asyncio
import io
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from fastapi.testclient import TestClient

from app.config import Settings
from app.core.keys import target_key
from app.core.spherex import CHANNEL_CENTERS_UM, N_CHANNELS, parse_flag_bits
from app.main import create_app
from app.services import fits_processor as fp
from app.services.cache_manager import FitsCache
from app.services.fits_stream import NotFitsError, read_hdus
from app.services.gaia_service import GaiaResult, GaiaService, propagate
from app.services.irsa_service import Frame
from app.services.solar_system_service import SolarSystemService
from app.services.spectrum_extractor import extract_spectrum, fit_spectral_model
from app.services.wise_service import CATWISE_EPOCH_MJD, WiseResult, WiseService, propagate_wise


def _make_synthetic_spherex_me_fits(
    ra: float = 269.4521,
    dec: float = 4.6934,
    angle_deg: float = 15.0,
    detector: int = 6,
    mjd: float = 60902.8,
    star_offset_arcsec: tuple[float, float] = (0.0, 0.0),
    wavelength_um: float = 4.65,
    include_epsf: bool = True,
) -> bytes:
    """Create a realistic SPHEREx L2 Multi-Extension FITS file with rotated WCS,
    IMAGE, FLAGS (with MSKN* headers), VARIANCE, ZODI, optional bulky EPSF, and
    WCS-WAVE lookup table."""
    ny = nx = 41
    pix = 6.15 / 3600.0
    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [(nx + 1) / 2, (ny + 1) / 2]
    w.wcs.cdelt = [1.0, 1.0]
    w.wcs.pc = [[-pix * c, pix * s], [pix * s, pix * c]]
    w.wcs.radesys = "ICRS"
    hdr = w.to_header()
    hdr["BUNIT"] = "MJy / sr"
    hdr["DETECTOR"] = detector
    hdr["MJD-AVG"] = mjd

    # Alternate WCS 'W' (WAVE-TAB lookup table)
    hdr["WCSAXESW"] = 2
    hdr["CTYPE1W"] = "WAVE-TAB"
    hdr["CTYPE2W"] = "WAVE-TAB"
    hdr["CUNIT1W"] = "um"
    hdr["CUNIT2W"] = "um"
    hdr["CRPIX1W"] = 0.0
    hdr["CRPIX2W"] = 0.0
    hdr["CDELT1W"] = 1.0
    hdr["CDELT2W"] = 1.0
    hdr["CRVAL1W"] = 1.0
    hdr["CRVAL2W"] = 1.0
    hdr["PV1_3W"] = 1
    hdr["PV2_3W"] = 2
    hdr["PS1_0W"] = "WCS-WAVE"
    hdr["PS2_0W"] = "WCS-WAVE"
    hdr["PS1_1W"] = "VALUES"
    hdr["PS2_1W"] = "VALUES"
    hdr["PS1_2W"] = "X"
    hdr["PS2_2W"] = "Y"

    # Synthetic star at (ra + d_ra, dec + d_dec)
    star_ra = ra + star_offset_arcsec[0] / 3600.0
    star_dec = dec + star_offset_arcsec[1] / 3600.0
    sx, sy = w.world_to_pixel_values(star_ra, star_dec)
    yy, xx = np.mgrid[0:ny, 0:nx]
    psf = 8.0 * np.exp(-((xx - sx) ** 2 + (yy - sy) ** 2) / (2 * 1.1**2))
    zodi = np.full((ny, nx), 0.35, dtype=np.float32)
    rng = np.random.default_rng(42 + int(angle_deg))
    img = (zodi + psf + rng.normal(0.0, 0.03, (ny, nx))).astype(np.float32)
    var = np.full((ny, nx), 0.03**2, dtype=np.float32)
    flg = np.zeros((ny, nx), dtype=np.int32)
    # Mark star core as SOURCE (bit 21) and one corner pixel as TRANSIENT (bit 0)
    flg[int(round(float(sy))), int(round(float(sx)))] = 1 << 21
    flg[2, 2] = 1 << 0
    img[2, 2] = 999.0  # cosmic ray that must be masked out

    fhdr = w.to_header()
    fhdr["MSKN0000"] = "TRANSIENT"
    fhdr["MSKN0001"] = "OVERFLOW"
    fhdr["MSKN0021"] = "SOURCE"

    hdus = [
        fits.PrimaryHDU(),
        fits.ImageHDU(img, header=hdr, name="IMAGE"),
        fits.ImageHDU(flg, header=fhdr, name="FLAGS"),
        fits.ImageHDU(var, header=w.to_header(), name="VARIANCE"),
        fits.ImageHDU(zodi, header=w.to_header(), name="ZODI"),
    ]
    if include_epsf:
        # Simulate bulky EPSF extension after ZODI
        epsf_col = fits.Column(name="PSF", format="1024E", array=np.zeros((20, 1024), dtype=np.float32))
        hdus.append(fits.BinTableHDU.from_columns([epsf_col], name="EPSF"))

    # WCS-WAVE table: 2x2 grid covering detector coordinates with constant wavelength_um
    x_arr = np.array([[0, 100]], dtype=np.int32)
    y_arr = np.array([[0, 100]], dtype=np.int32)
    vals = np.zeros((1, 2, 2, 2), dtype=np.float32)
    vals[..., 0] = wavelength_um
    vals[..., 1] = wavelength_um / 41.0
    t = Table([x_arr, y_arr, vals], names=("X", "Y", "VALUES"))
    whdu = fits.table_to_hdu(t)
    whdu.name = "WCS-WAVE"
    hdus.append(whdu)

    buf = io.BytesIO()
    fits.HDUList(hdus).writeto(buf)
    return buf.getvalue()


def test_streaming_fits_early_stop():
    raw = _make_synthetic_spherex_me_fits(include_epsf=True)

    async def _go():
        async def _chunks():
            for i in range(0, len(raw), 4096):
                yield raw[i : i + 4096]

        return await read_hdus(_chunks(), wanted=("IMAGE", "FLAGS", "VARIANCE", "ZODI"))

    trimmed, names, truncated = asyncio.run(_go())
    assert truncated is True
    assert set(names) >= {"PRIMARY", "IMAGE", "FLAGS", "VARIANCE", "ZODI"}
    assert "EPSF" not in names
    assert len(trimmed) < len(raw)  # EPSF skipped!

    # Verify non-FITS "Server busy" raises NotFitsError
    async def _bad():
        async def _chunks():
            yield b"Server busy, retry later."

        await read_hdus(_chunks(), wanted=("IMAGE",))

    with pytest.raises(NotFitsError):
        asyncio.run(_bad())


def test_wcs_alignment_and_layers_and_spectrum(tmp_path: Path):
    ra, dec = 269.4521, 4.6934
    # Two frames rotated by 50 degrees relative to each other
    b1 = _make_synthetic_spherex_me_fits(ra, dec, angle_deg=10.0, detector=6, mjd=60900.0, wavelength_um=4.65)
    b2 = _make_synthetic_spherex_me_fits(ra, dec, angle_deg=60.0, detector=6, mjd=61085.0, wavelength_um=4.65)
    p1, p2 = tmp_path / "f1.fits", tmp_path / "f2.fits"
    p1.write_bytes(b1)
    p2.write_bytes(b2)

    fd1 = fp.load_frame(p1, "obs1/D6", 6, 60900.0)
    fd2 = fp.load_frame(p2, "obs2/D6", 6, 61085.0)
    # Cosmic ray at (2, 2) must be masked to NaN while SOURCE at centre stays finite
    assert np.isnan(fd1.image[2, 2])
    assert np.isfinite(fd1.image[20, 20])

    grid = fp.make_grid(ra, dec, fov=0.06, pixscale_arcsec=6.15, png_size=256)
    c1 = fp.coadd_epoch(0, [fd1], grid)
    c2 = fp.coadd_epoch(1, [fd2], grid)

    # Peak pixel of the stationary star must land at the exact same (y, x) despite 50 deg rotation!
    peak1 = np.unravel_index(np.nanargmax(c1.image), c1.image.shape)
    peak2 = np.unravel_index(np.nanargmax(c2.image), c2.image.shape)
    assert peak1 == peak2

    vmin, vmax, _ = fp.compute_cuts([c1.image, c2.image])
    img_png = fp.render_image(c1.image, vmin, vmax, "asinh", size=256)
    var_png = fp.render_variance(c1.variance, size=256)
    flg_png = fp.render_flags(c1.flags, c1.coverage, parse_flag_bits(fits.open(p1)["FLAGS"].header), size=256)
    assert img_png[:8] == b"\x89PNG\r\n\x1a\n"
    assert var_png[:8] == b"\x89PNG\r\n\x1a\n"
    assert flg_png[:8] == b"\x89PNG\r\n\x1a\n"

    # Spectrum extraction
    s1 = fp.extract_sample(fd1, ra, dec)
    assert s1 is not None
    assert s1.flux_ujy > 0
    pts = fp.bin_spectrum([s1])
    assert len(pts) == N_CHANNELS == len(CHANNEL_CENTERS_UM) == 102
    filled = [p for p in pts if p.flux is not None]
    assert len(filled) == 1
    assert abs(filled[0].wavelength_um - 4.65) < 0.05


def test_gaia_proper_motion_propagation():
    # Barnard's star motion (~ -801 mas/yr RA*, +10362 mas/yr Dec) over ~10 years
    src = {"source_id": "4472832130942575872", "ra": 269.454, "dec": 4.668, "pmra": -801.5, "pmdec": 10362.4, "pm": 10393.0}
    ra_obs, dec_obs = propagate(src, mjd=61041.0)  # ~2026.0 (dt = 10 yr from J2016.0)
    assert abs((dec_obs - src["dec"]) * 3600.0 - 103.6) < 1.0


def test_cache_lru_and_ttl(tmp_path: Path):
    cache = FitsCache(tmp_path, ttl_s=60, max_bytes=10_000)
    for i in range(5):
        p = cache.raw_path(f"https://example.org/{i}.fits")
        cache.write_atomic(p, b"X" * 3000)
    stats = cache.sweep()
    assert stats["removed"] >= 2
    assert stats["total_bytes"] <= 10_000


def test_api_target_and_restretch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    ra, dec, fov = 133.7951, -7.2451, 0.06
    cfg = Settings(
        cache_dir=tmp_path,
        png_size=128,
        max_frames_per_detector_per_epoch=1,
        supabase_url="",
        supabase_key="",
        supabase_anon_key="",
    )
    app = create_app(cfg)

    raw_e0 = _make_synthetic_spherex_me_fits(ra, dec, angle_deg=0.0, detector=6, mjd=60850.0, wavelength_um=4.55)
    raw_e1 = _make_synthetic_spherex_me_fits(
        ra, dec, angle_deg=35.0, detector=6, mjd=61040.0, star_offset_arcsec=(15.0, 0.0), wavelength_um=4.75
    )
    f0 = tmp_path / "e0.fits"
    f1 = tmp_path / "e1.fits"
    f0.write_bytes(raw_e0)
    f1.write_bytes(raw_e1)

    with TestClient(app) as client:
        pipeline = app.state.pipeline

        async def _fake_query_frames(r: float, d: float):
            return [
                Frame("2025W20_1A", 6, "spherex_dr1", 60850.0, 60850.01, "https://example/dl0", r, d, "https://example/f0.fits"),
                Frame("2025W48_1A", 6, "spherex_dr1", 61040.0, 61040.01, "https://example/dl1", r, d, "https://example/f1.fits"),
            ]

        async def _fake_fetch_many(frames, r, d, size_deg):
            return [(frames[0], f0), (frames[1], f1)][: len(frames)]

        async def _fake_gaia(r: float, d: float, fv: float):
            # Empty optical catalog -> WISE 0855 should be classified as `candidate`!
            return GaiaResult(match_candidates=[], field=[])

        async def _fake_wise(r: float, d: float):
            return WiseResult(
                sources=[
                    {
                        "source_id": "J085510.83-071442.5",
                        "ra": r + 0.005,
                        "dec": d,
                        "pmra": -1700.0,
                        "pmdec": 0.0,
                        "w1mpro": 17.8,
                        "w2mpro": 13.9,
                        "w1snr": 4.2,
                        "w2snr": 28.5,
                    }
                ]
            )

        async def _fake_skybot(r: float, d: float, fv: float, mjds: list[float], world_to_png=None):
            from app.schemas.target import MinorPlanetOverlay

            out = []
            for idx_m, _m in enumerate(mjds):
                if idx_m == 0:
                    xy = world_to_png(r, d) if world_to_png else (64.0, 64.0)
                    out.append(
                        [
                            MinorPlanetOverlay(
                                number="433",
                                name="Eros",
                                object_type="MBA",
                                ra=r,
                                dec=d,
                                v_mag=14.2,
                                separation_arcsec=0.0,
                                x=xy[0] if xy else 64.0,
                                y=xy[1] if xy else 64.0,
                            )
                        ]
                    )
                else:
                    out.append([])
            return out

        monkeypatch.setattr(pipeline.irsa, "query_frames", _fake_query_frames)
        monkeypatch.setattr(pipeline.irsa, "fetch_many", _fake_fetch_many)
        monkeypatch.setattr(pipeline.gaia, "query", _fake_gaia)
        monkeypatch.setattr(pipeline.wise, "query", _fake_wise)
        monkeypatch.setattr(pipeline.solar_system, "query_epochs", _fake_skybot)

        # 1. Cold call -> full pipeline runs
        r1 = client.get("/api/v1/target", params={"ra": ra, "dec": dec, "fov": fov})
        assert r1.status_code == 200
        body = r1.json()
        assert body["status"] == "ready"
        assert body["cache"]["hit"] is False
        assert len(body["epochs"]) == 2
        assert body["epochs"][0]["layers"]["rgb"] is not None
        assert len(body["epochs"][0]["known_minor_planets"]) == 1
        assert body["epochs"][0]["known_minor_planets"][0]["name"] == "Eros"
        assert body["epochs"][1]["known_minor_planets"] == []
        assert body["crossmatch"]["classification"] == "candidate"
        assert body["crossmatch"]["wise_match"] is not None
        assert body["crossmatch"]["wise_match"]["is_ultracool_color"] is True
        assert abs(body["crossmatch"]["wise_match"]["w1_w2_color"] - 3.9) < 0.01
        assert len(body["spectrum"]["points"]) == 102
        assert body["spectrum"]["n_channels_filled"] == 2
        assert body["spectrum_data"] is not None
        assert body["spectrum_data"]["spectrum_summary"]["total_channels"] == 102
        assert len(body["spectrum_data"]["chart_data"]) == 2
        assert len(body["spectrum_data"]["absorption_bands_reference"]) >= 3

        # Verify Difference Imaging & Moving Dipole auto-detection
        diff = body["difference"]
        assert diff is not None
        assert diff["epoch_a"] == 0 and diff["epoch_b"] == 1
        assert diff["baseline_days"] == 190.0
        assert diff["max_dipole_snr"] > 5.0
        dipoles = [a for a in diff["anomalies"] if a["anomaly_type"] == "moving_dipole"]
        assert len(dipoles) >= 1
        assert dipoles[0]["shift_arcsec"] is not None and 10.0 <= dipoles[0]["shift_arcsec"] <= 22.0
        assert dipoles[0]["pm_arcsec_yr"] is not None and dipoles[0]["pm_arcsec_yr"] > 15.0
        assert dipoles[0]["x_from"] is not None and dipoles[0]["y_from"] is not None

        # Verify standalone GET /api/v1/target/diff endpoint
        r_diff = client.get("/api/v1/target/diff", params={"ra": ra, "dec": dec, "fov": fov, "epoch_a": 0, "epoch_b": 1})
        assert r_diff.status_code == 200
        assert r_diff.json()["max_dipole_snr"] > 5.0

        # 2. Dedicated lightweight GET /api/v1/target/spectrum endpoint
        rs_spec = client.get("/api/v1/target/spectrum", params={"ra": ra, "dec": dec})
        assert rs_spec.status_code == 200
        spec_json = rs_spec.json()
        assert spec_json["target"] == {"ra": ra, "dec": dec}
        wls = [pt["wavelength_um"] for pt in spec_json["chart_data"]]
        assert wls == sorted(wls)
        assert spec_json["spectrum_summary"]["valid_points"] == 2

        # 3. Warm call -> instant Cache Hit (even with 0.3" coordinate jitter)
        r2 = client.get("/api/v1/target", params={"ra": ra + 0.00005, "dec": dec, "fov": fov})
        assert r2.status_code == 200
        assert r2.json()["cache"]["hit"] is True
        assert r2.json()["target_key"] == target_key(ra, dec, fov)

        # 3. Fast restretch from local `/tmp/fits_cache`
        r3 = client.post(
            "/api/v1/target/restretch",
            json={
                "ra": ra,
                "dec": dec,
                "fov": fov,
                "stretch_type": "log",
                "min_cut": -0.05,
                "max_cut": 4.0,
            },
        )
        assert r3.status_code == 200
        rs = r3.json()
        assert rs["source"] == "local_cache"
        assert len(rs["epochs"]) == 2
        assert rs["epochs"][0]["image"].startswith("data:image/png;base64,")


def test_extract_spectrum_fallback_without_wcs_wave(tmp_path: Path):
    """Verify extract_spectrum never crashes when WCS-WAVE HDU is missing and
    correctly classifies a red brown-dwarf-like SED."""
    ra, dec = 133.7951, -7.2451
    # Short-wave Band 1 frame (faint) and long-wave Band 6 frame (bright) without WCS-WAVE HDU
    b_short = _make_synthetic_spherex_me_fits(ra, dec, detector=1, wavelength_um=0.95, include_epsf=False)
    b_long = _make_synthetic_spherex_me_fits(ra, dec, detector=6, wavelength_um=4.65, include_epsf=False)

    p_short = tmp_path / "short_no_wave.fits"
    p_long = tmp_path / "long_no_wave.fits"
    with fits.open(io.BytesIO(b_short)) as h1:
        # Strip WCS-WAVE HDU and scale flux down so short-band is faint
        h1["IMAGE"].data = h1["ZODI"].data + (h1["IMAGE"].data - h1["ZODI"].data) * 0.1
        h1["IMAGE"].header["CRPIX2W"] = -800.0
        fits.HDUList([h for h in h1 if h.name != "WCS-WAVE"]).writeto(p_short)
    with fits.open(io.BytesIO(b_long)) as h6:
        h6["IMAGE"].header["CRPIX2W"] = -500.0
        fits.HDUList([h for h in h6 if h.name != "WCS-WAVE"]).writeto(p_long)

    out = extract_spectrum([p_short, p_long, tmp_path / "nonexistent.fits"], ra, dec, wavecal_dir=tmp_path)
    summary = out["spectrum_summary"]
    assert summary["total_channels"] == 102
    assert summary["valid_points"] == 2
    assert summary["ir_color_index"] is not None and summary["ir_color_index"] > 2.0
    assert summary["spectral_hint"] in ("cold_brown_dwarf_candidate", "ultra_cool_y_dwarf_or_planet_x_candidate")
    chart = out["chart_data"]
    assert len(chart) == 2
    assert chart[0]["wavelength_um"] < chart[1]["wavelength_um"]
    assert chart[0]["band"] == "Band 1"
    assert chart[1]["band"] == "Band 6"


def test_wise_crossmatch_and_propagation():
    """Verify CatWISE2020 proper-motion propagation, W1-W2 color, AllWISE fallback, and 10yr baseline PM."""
    import httpx

    cfg = Settings(supabase_url="", supabase_key="")
    wise = WiseService(cfg)

    # Simulate an ultracool Y dwarf observed by CatWISE2020 at MJD 57170 (~2015.4)
    # moving -8100 mas/yr in RA*cos(dec) over 10.5 years to SPHEREx MJD 61005.125
    target_ra = 133.7951
    target_dec = -7.2451
    dt_yr = (61005.125 - CATWISE_EPOCH_MJD) / 365.25
    cosd = np.cos(np.radians(target_dec))
    catwise_ra0 = (target_ra - (-8100.0) * dt_yr / 3.6e6 / cosd) % 360.0

    src = {
        "source_id": "WISEA_J0855",
        "ra": catwise_ra0,
        "dec": target_dec,
        "pmra": -8.1,  # arcsec/yr
        "pmdec": 0.0,
        "pm_unit": "arcsec/yr",
        "w1mpro": 17.9,
        "w2mpro": 13.8,
        "w1snr": 3.5,
        "w2snr": 32.0,
    }
    ra_prop, dec_prop = propagate_wise(src, mjd=61005.125)
    assert abs(ra_prop - target_ra) * 3600.0 < 0.1
    assert abs(dec_prop - target_dec) < 1e-6

    res = WiseResult(sources=[src], catalog="CatWISE2020", epoch_mjd=CATWISE_EPOCH_MJD)
    match = wise.best_match(res, target_ra, target_dec, mjd=61005.125, spherex_ra=target_ra, spherex_dec=target_dec)
    assert match is not None
    assert match.source_id == "WISEA_J0855"
    assert match.separation_arcsec < 0.2
    assert match.w1_w2_color == pytest.approx(4.1, abs=0.01)
    assert match.is_ultracool_color is True
    assert match.baseline_pm_arcsec_yr == pytest.approx(8.1, abs=0.1)

    # Test HTTP CSV parsing + automatic fallback from catwise_2020 to allwise_p3as_psd
    calls: list[str] = []

    def _handler(request: httpx.Request) -> httpx.Response:
        body = request.content.decode("utf-8")
        calls.append(body)
        if "catwise_2020" in body:
            return httpx.Response(503, text="CatWISE unavailable")
        csv_text = (
            "source_id,ra,dec,w1mpro,w2mpro,w1snr,w2snr\n"
            f"ALLWISE_001,{target_ra},{target_dec},15.2,13.4,12.0,25.0\n"
        )
        return httpx.Response(200, text=csv_text)

    async def _run_fallback():
        async with httpx.AsyncClient(transport=httpx.MockTransport(_handler)) as client:
            svc = WiseService(cfg, client)
            return await svc.query(target_ra, target_dec)

    fb_res = asyncio.run(_run_fallback())
    assert len(calls) == 2
    assert fb_res.catalog == "AllWISE"
    assert len(fb_res.sources) == 1
    assert fb_res.sources[0]["source_id"] == "ALLWISE_001"


def test_solar_system_skybot_overlay():
    """Verify IMCCE SkyBoT JSON + pipe-text parsing, sexagesimal conversion, and PNG coordinate projection."""
    import json

    import httpx

    ra, dec = 133.7951, -7.2451
    grid = fp.make_grid(ra, dec, fov=0.1, pixscale_arcsec=6.15, png_size=512)
    cfg = Settings(supabase_url="", supabase_key="", skybot_enabled=True)

    # JSON response with fields + sexagesimal RA (08 55 10.82 = 133.7951 deg) and Dec (-07 14 42.36 = -7.2451 deg)
    json_payload = json.dumps(
        {
            "fields": ["num", "name", "ra", "dec", "class", "vmag", "errpos", "d", "dra", "ddec"],
            "data": [
                ["1", "Ceres", "08 55 10.824", "-07 14 42.36", "MBA", 8.1, 0.1, 0.0, -25.4, 12.1],
                ["999", "FarAway", 140.0, 20.0, "MBA", 18.5, 0.2, 5000.0, 1.0, 1.0],
            ],
        }
    )
    pipe_payload = (
        "# Flag : 1\n"
        "# num | name | ra | dec | class | vmag | errpos | d | dra | ddec\n"
        f"4 | Vesta | {ra + 0.01} | {dec - 0.01} | MBA | 7.4 | 0.1 | 40.0 | -18.0 | 5.2\n"
    )

    def _handler(request: httpx.Request) -> httpx.Response:
        ep = float(request.url.params.get("-ep", "60800"))
        if ep < 60900:
            return httpx.Response(200, text=json_payload)
        return httpx.Response(200, text=pipe_payload)

    async def _run():
        async with httpx.AsyncClient(transport=httpx.MockTransport(_handler)) as client:
            svc = SolarSystemService(cfg, client)
            return await svc.query_epochs(ra, dec, 0.1, [60850.0, 61040.0], grid.world_to_png)

    res = asyncio.run(_run())
    assert len(res) == 2
    assert len(res[0]) == 1  # FarAway filtered out because outside PNG bounds
    assert res[0][0].name == "Ceres"
    assert res[0][0].number == "1"
    assert 250.0 <= res[0][0].x <= 262.0
    assert 250.0 <= res[0][0].y <= 262.0
    assert res[0][0].dra_arcsec_hr == pytest.approx(-25.4, abs=0.01)

    assert len(res[1]) == 1
    assert res[1][0].name == "Vesta"
    assert res[1][0].v_mag == pytest.approx(7.4, abs=0.01)


def test_fit_spectral_model_brown_dwarf_and_star():
    """Verify Planck Teff fitting and CH4/H2O absorption depth extraction."""
    from app.schemas.spectrum import ChartPoint

    def _make_sed(temp_k: float, scale: float, ch4_deep: float = 0.0, h2o_deep: float = 0.0) -> list[ChartPoint]:
        wls = np.linspace(0.8, 4.9, 45)
        pts: list[ChartPoint] = []
        raw = (wls**-3.0) / np.expm1(np.clip(14387.77 / (wls * temp_k), 1e-6, 700.0))
        raw = scale * raw / float(np.max(raw))
        for idx, (wl, fl) in enumerate(zip(wls, raw), start=1):
            if ch4_deep > 0.0 and 3.20 <= wl <= 3.40:
                fl *= 1.0 - ch4_deep
            if h2o_deep > 0.0 and 2.55 <= wl <= 2.85:
                fl *= 1.0 - h2o_deep
            det = min(6, max(1, int((wl - 0.75) / 0.71) + 1))
            pts.append(
                ChartPoint(
                    channel_id=min(102, idx * 2),
                    wavelength_um=float(wl),
                    bandwidth_um=0.04,
                    band=f"Band {det}",
                    detector=det,
                    flux=round(float(fl), 3),
                    flux_error=max(0.5, 0.05 * float(fl)),
                    snr=20.0,
                    is_flagged=False,
                    n_samples=2,
                )
            )
        return pts

    # 1. Cold 650 K T/Y brown dwarf with strong CH4 (3.3 um) and H2O (2.7 um) absorption
    bd_pts = _make_sed(temp_k=650.0, scale=250.0, ch4_deep=0.65, h2o_deep=0.50)
    bd_fit = fit_spectral_model(bd_pts)
    assert bd_fit["estimated_temp_k"] is not None
    assert 500.0 <= bd_fit["estimated_temp_k"] <= 800.0
    assert bd_fit["ch4_absorption_depth"] == pytest.approx(0.65, abs=0.12)
    assert bd_fit["h2o_absorption_depth"] == pytest.approx(0.50, abs=0.12)
    assert len(bd_fit["model_curve"]) == 50

    # 2. Hot 4500 K star without molecular absorption
    star_pts = _make_sed(temp_k=4500.0, scale=1000.0, ch4_deep=0.0, h2o_deep=0.0)
    star_fit = fit_spectral_model(star_pts)
    assert star_fit["estimated_temp_k"] is not None
    assert star_fit["estimated_temp_k"] >= 3500.0
    assert (star_fit["ch4_absorption_depth"] or 0.0) < 0.1
    assert (star_fit["h2o_absorption_depth"] or 0.0) < 0.1




