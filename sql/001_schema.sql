-- ============================================================================
-- SkySurvey SPHEREx Backend — Supabase PostgreSQL + Storage Schema
-- Run this script once in the Supabase SQL Editor.
-- ============================================================================

create extension if not exists pgcrypto;
create extension if not exists cube;
create extension if not exists earthdistance;

-- ----------------------------------------------------------------------------
-- 1. Main cache table: `targets`
-- Stores processed epoch layer URLs, Gaia DR3 cross-match, and 102-band spectra
-- keyed by canonical `target_key` (ra_dec_fov_band) and indexed on the sky
-- sphere for fast tolerance lookups.
-- ----------------------------------------------------------------------------
create table if not exists public.targets (
    id                uuid primary key default gen_random_uuid(),
    target_key        text not null unique,
    ra                double precision not null check (ra >= 0 and ra < 360),
    dec               double precision not null check (dec >= -90 and dec <= 90),
    fov               double precision not null check (fov > 0 and fov <= 2.0),
    band              text not null default 'all',

    -- Sky unit-vector components for fast index-assisted angular lookup
    cx                double precision generated always as (cos(radians(dec)) * cos(radians(ra))) stored,
    cy                double precision generated always as (cos(radians(dec)) * sin(radians(ra))) stored,
    cz                double precision generated always as (sin(radians(dec))) stored,

    status            text not null default 'processing'
                      check (status in ('processing', 'partial', 'ready', 'failed', 'no_data')),
    classification    text not null default 'unknown'
                      check (classification in ('known_star', 'candidate', 'no_detection', 'unknown')),
    pipeline_version  text not null default '1.0.0',
    n_observations    integer not null default 0,

    progress          jsonb,
    display           jsonb,
    epochs            jsonb not null default '[]'::jsonb,
    crossmatch        jsonb,
    spectrum          jsonb,
    spectrum_data     jsonb,
    difference        jsonb,
    error             text,

    created_at        timestamptz not null default now(),
    updated_at        timestamptz not null default now()
);

-- Safe migration if `targets` was already created before `spectrum_data` / `difference` were added
alter table public.targets
    add column if not exists spectrum_data jsonb,
    add column if not exists difference jsonb;

create index if not exists idx_targets_ra_dec_fov
    on public.targets (ra, dec, fov, band);

create index if not exists idx_targets_unit_vec
    on public.targets (cx, cy, cz);

create index if not exists idx_targets_classification
    on public.targets (classification)
    where status = 'ready';

-- Auto-maintain `updated_at`
create or replace function public.set_updated_at()
returns trigger
language plpgsql
as $$
begin
    new.updated_at = now();
    return new;
end;
$$;

drop trigger if exists trg_targets_updated_at on public.targets;
create trigger trg_targets_updated_at
before update on public.targets
for each row execute function public.set_updated_at();

-- ----------------------------------------------------------------------------
-- 2. Datalink URL cache table: `observations`
-- IRSA `spherex.obscore` returns a Datalink VOTable URL per (obs_id, detector)
-- which takes ~0.5-1.0s to resolve into the underlying physical FITS URL.
-- Caching `(obs_key -> fits_url)` eliminates Datalink round-trips across nearby
-- targets or different FOVs.
-- ----------------------------------------------------------------------------
create table if not exists public.observations (
    obs_key           text primary key,          -- e.g. '2025W33_2C_0321_1/D6'
    obs_id            text not null,
    detector          smallint not null check (detector between 1 and 6),
    collection        text not null,             -- e.g. 'spherex_dr1'
    t_min             double precision not null, -- MJD
    t_max             double precision not null,
    s_ra              double precision,
    s_dec             double precision,
    fits_url          text not null,
    created_at        timestamptz not null default now()
);

create index if not exists idx_observations_obs_id
    on public.observations (obs_id);

-- ----------------------------------------------------------------------------
-- 3. RPC: Spherical tolerance search (`find_target`)
-- Finds the closest cached target within `p_tol_arcsec` (handling RA=0/360
-- wrap-around cleanly via 3D unit vectors).
-- ----------------------------------------------------------------------------
create or replace function public.find_target(
    p_ra          double precision,
    p_dec         double precision,
    p_fov         double precision,
    p_band        text default 'all',
    p_tol_arcsec  double precision default 2.0
)
returns setof public.targets
language sql
stable
as $$
    with q as (
        select
            cos(radians(p_dec)) * cos(radians(p_ra)) as qx,
            cos(radians(p_dec)) * sin(radians(p_ra)) as qy,
            sin(radians(p_dec))                      as qz,
            2.0 * sin(radians(p_tol_arcsec / 3600.0) / 2.0) as chord_tol
    )
    select t.*
    from public.targets t, q
    where t.band = p_band
      and abs(t.fov - p_fov) < 0.0005
      and t.cx between q.qx - q.chord_tol and q.qx + q.chord_tol
      and t.cy between q.qy - q.chord_tol and q.qy + q.chord_tol
      and t.cz between q.qz - q.chord_tol and q.qz + q.chord_tol
      and sqrt(power(t.cx - q.qx, 2) + power(t.cy - q.qy, 2) + power(t.cz - q.qz, 2)) <= q.chord_tol
    order by sqrt(power(t.cx - q.qx, 2) + power(t.cy - q.qy, 2) + power(t.cz - q.qz, 2)) asc
    limit 1;
$$;

-- Coordinate-only lookup for GET /api/v1/target/spectrum?ra=...&dec=...
create or replace function public.find_target_by_coords(
    p_ra          double precision,
    p_dec         double precision,
    p_tol_arcsec  double precision default 3.0
)
returns setof public.targets
language sql
stable
as $$
    with q as (
        select
            cos(radians(p_dec)) * cos(radians(p_ra)) as qx,
            cos(radians(p_dec)) * sin(radians(p_ra)) as qy,
            sin(radians(p_dec))                      as qz,
            2.0 * sin(radians(p_tol_arcsec / 3600.0) / 2.0) as chord_tol
    )
    select t.*
    from public.targets t, q
    where t.spectrum_data is not null
      and t.cx between q.qx - q.chord_tol and q.qx + q.chord_tol
      and t.cy between q.qy - q.chord_tol and q.qy + q.chord_tol
      and t.cz between q.qz - q.chord_tol and q.qz + q.chord_tol
      and sqrt(power(t.cx - q.qx, 2) + power(t.cy - q.qy, 2) + power(t.cz - q.qz, 2)) <= q.chord_tol
    order by sqrt(power(t.cx - q.qx, 2) + power(t.cy - q.qy, 2) + power(t.cz - q.qz, 2)) asc,
             t.updated_at desc
    limit 1;
$$;

-- ----------------------------------------------------------------------------
-- 4. Row-Level Security (RLS) & Public Storage Bucket
-- Requirement: No user auth; backend connects using `SUPABASE_KEY`
-- (`sb_publishable_...`, role `anon`). We therefore enable RLS and grant
-- read/write policies on `targets`, `observations`, and the `skysurvey` bucket.
-- ----------------------------------------------------------------------------
alter table public.targets enable row level security;
alter table public.observations enable row level security;

drop policy if exists "targets_public_all" on public.targets;
create policy "targets_public_all"
    on public.targets
    for all
    to anon, authenticated, service_role
    using (true)
    with check (true);

drop policy if exists "observations_public_all" on public.observations;
create policy "observations_public_all"
    on public.observations
    for all
    to anon, authenticated, service_role
    using (true)
    with check (true);

insert into storage.buckets (id, name, public, file_size_limit, allowed_mime_types)
values (
    'skysurvey',
    'skysurvey',
    true,
    10485760, -- 10 MB per object
    array['image/png', 'application/fits', 'application/octet-stream']
)
on conflict (id) do update
set public = true;

drop policy if exists "skysurvey_storage_public_all" on storage.objects;
create policy "skysurvey_storage_public_all"
    on storage.objects
    for all
    to anon, authenticated, service_role
    using (bucket_id = 'skysurvey')
    with check (bucket_id = 'skysurvey');

-- Refresh PostgREST schema cache immediately
notify pgrst, 'reload schema';

