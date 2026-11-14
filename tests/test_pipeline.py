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
from app.services.spectrum_extractor import extract_spectrum


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

        monkeypatch.setattr(pipeline.irsa, "query_frames", _fake_query_frames)
        monkeypatch.setattr(pipeline.irsa, "fetch_many", _fake_fetch_many)
        monkeypatch.setattr(pipeline.gaia, "query", _fake_gaia)

        # 1. Cold call -> full pipeline runs
        r1 = client.get("/api/v1/target", params={"ra": ra, "dec": dec, "fov": fov})
        assert r1.status_code == 200
        body = r1.json()
        assert body["status"] == "ready"
        assert body["cache"]["hit"] is False
        assert len(body["epochs"]) == 2
        assert body["epochs"][0]["layers"]["rgb"] is not None
        assert body["crossmatch"]["classification"] == "candidate"
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

