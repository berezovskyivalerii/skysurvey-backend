"""Diagnostic script to inspect existing tables, columns, RPCs, and Storage buckets
in the configured Supabase project, and verify read/write permissions."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import httpx
from app.config import get_settings


def main() -> int:
    s = get_settings()
    if not s.supabase_url or not s.supabase_api_key:
        print("ERROR: SUPABASE_URL or SUPABASE_KEY is empty in .env")
        return 1

    base = s.supabase_url.rstrip("/")
    headers = {
        "apikey": s.supabase_api_key,
        "Authorization": f"Bearer {s.supabase_api_key}",
    }
    print(f"Checking Supabase project: {base}")

    with httpx.Client(timeout=15.0, headers=headers) as client:
        # 1. Probe tables directly via PostgREST (/rest/v1/<table>?limit=1)
        print("\n=== 1. Checking Required Tables ===")
        for tbl in ("targets", "observations"):
            rt = client.get(f"{base}/rest/v1/{tbl}", params={"select": "*", "limit": "1"})
            if rt.status_code == 200:
                rows = rt.json()
                info = f"{len(rows)} sample row(s)" if rows else "table ready (0 rows)"
                print(f"  [OK] public.{tbl}: HTTP 200 ({info})")
            else:
                print(f"  [MISSING] public.{tbl}: HTTP {rt.status_code} -> {rt.text[:200]}")

        # 2. Probe RPC functions
        print("\n=== 2. Checking RPC Functions ===")
        r_rpc = client.post(
            f"{base}/rest/v1/rpc/find_target",
            json={"p_ra": 133.7951, "p_dec": -7.2451, "p_fov": 0.1, "p_band": "all", "p_tol_arcsec": 2.0},
        )
        print(f"  [{'OK' if r_rpc.status_code == 200 else 'MISSING'}] rpc/find_target: HTTP {r_rpc.status_code}")

        # 3. Check Storage bucket & test upload to `skysurvey`
        # Note: sql/001_schema.sql restricts allowed_mime_types to
        # ['image/png', 'application/fits', 'application/octet-stream']
        print("\n=== 3. Checking Storage Bucket 'skysurvey' ===")
        tiny_png = (
            b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
            b"\x08\x02\x00\x00\x00\x90wS\xde\x00\x00\x00\x0cIDATx\x9cc\xf8\x0f\x00"
            b"\x00\x01\x01\x00\x05\x18\xd8N\x00\x00\x00\x00IEND\xaeB`\x82"
        )
        r_up = client.post(
            f"{base}/storage/v1/object/{s.supabase_bucket}/healthcheck.png",
            content=tiny_png,
            headers={**headers, "Content-Type": "image/png", "x-upsert": "true"},
        )
        if r_up.status_code == 200:
            pub_url = f"{base}/storage/v1/object/public/{s.supabase_bucket}/healthcheck.png"
            print(f"  [OK] Upload succeeded! Public URL: {pub_url}")
        else:
            print(f"  [MISSING] Upload failed: HTTP {r_up.status_code} -> {r_up.text[:200]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
