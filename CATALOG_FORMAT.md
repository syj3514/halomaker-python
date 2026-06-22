# Catalog Format

The pipeline produces up to two HDF5 outputs, **row-aligned and joined by `id`**:

1. `tree_bricks{iout:05d}.h5` — the **HaloMaker** halo/galaxy catalog (always).
2. `gas_bricks{iout:05d}.h5` — the **GasMaker** gas properties (optional
   post-processing step; see `README.md`).

Both share the same per-row ordering as `/catalog/halo`, so column `id` links a
halo's structural/stellar properties (output 1) to its gas properties (output 2).

> Unit conventions (read this first):
> - HaloMaker catalog: lengths/positions in **physical Mpc**, masses in
>   **10¹¹ M⊙**, velocities in **km/s**, temperatures in **K** (per the
>   HaloMaker convention).
> - GasMaker: gas/particle masses in **M⊙**, metallicity & chemistry are
>   **dimensionless mass fractions**, kinematics in **km/s**, gas angular
>   momentum in **10¹¹ M⊙·Mpc·km/s**.
> - ⚠️ **GasMaker `r200`/`r500` are in code units (box fraction, [0,1))** — not
>   Mpc. Multiply by `box_physical_mpc` (= `header.boxsize2 × aexp`) for physical
>   Mpc. (The catalog's `r`/`rvir`/`r50`/`r90`/`r*` are physical Mpc.)

---

## 1. HaloMaker catalog — `tree_bricks{iout:05d}.h5`

Groups: `/catalog/halo`, `/photometry/{CB07,BC03,FSPS}`, `/member`, `/header`,
`/input`.

### `/catalog/halo` (72 fields)

**Identity & hierarchy** — `id`, `timestep`, `nmem`, `ndm`, `nstar`, `nbsub`
(# subhalos), `hosthalo`, `hostsub`, `level`, `nextsub`.

**Position (Mpc) & velocity (km/s)** — `px,py,pz` (halo center),
`px*,py*,pz*` (stellar center), `vx,vy,vz`.

**Angular momentum** — `Lx,Ly,Lz` (total), `Lx*,Ly*,Lz*` (stars).

**Shape** — `sha,shb,shc` (principal-axis lengths).

**Mass (10¹¹ M⊙)** — `m` (total), `mdm` (dark matter), `m*` (stars).

**Radii (Mpc)** — `r` (extent), `r*` (stellar radius), `r50`, `r90` (stellar
half-/90%-mass radii), `rvir`.

**Stellar population** — `age`, `metal` (mass-weighted); star formation rate
over 100 Myr (`SFR`, `SFR_r50`, `SFR_r90`) and 10 Myr (`SFR10`, `SFR10_r50`,
`SFR10_r90`).

**Spin & dispersions** — `spin`; `sigma`, `sigma_dm`, `sigma*`; stellar
kinematics `vrot`, `sig3d`, `sigcyl` and their `_r50` / `_r90` apertures.

**Energy & virial** — `ek`, `ep`, `et`; `mvir`, `tvir` (K), `cvel`.

**Density profile (NFW)** — `rho_0`, `r_c`, `cNFW`, `cNFWerr`, `vmaxcir`,
`rmaxcir`, `inslope`, `inslopeerr` (inner-slope fit; note the
platform-sensitivity caveat in `WHATS_NEW.md`).

**Contamination** — `mcontam` (low-resolution / contaminant mass).

### `/photometry/{CB07,BC03,FSPS}`
Rest-frame stellar photometry, row-aligned with `/catalog/halo`: SDSS `umag`,
`gmag`, `rmag`, `imag`, `zmag` (AB); Johnson `Umag`, `Bmag`, `Vmag`, `Kmag`
(Vega; FSPS `Kmag` uses 2MASS Ks — see group metadata); plus r-band
luminosity-weighted `age_r`, `metal_r`, `r50_r`, `r90_r`. See `SSP_MODELS.md`.

### `/member`
Member particles: `index` and `count` (per-row offsets into the flat arrays) and
`pids` (particle IDs). When `dump_DMs` is enabled, `pos`, `vel`, `mass` are also
written in the same flat ordering.

### `/header`, `/input`
Run/cosmology metadata: box size (`boxsize2`), `aexp`, Hubble/cosmology
parameters, particle counts, and the input configuration.

---

## 2. GasMaker output — `gas_bricks{iout:05d}.h5`

Groups: `/header`, `/gas/{summary,processed,root_metrics}`.

### `/header` (attrs)
Schema and run metadata + aggregate progress: `schema_version`,
`source_catalog`, `run_mode`, `radius_field`, `root_padding`,
`cell_boundary_overlap`, `ownership`, `overlap_method`, `overlap_threshold`,
`overlap_depth`, `overlap_tolerance`, `requested_root_ids`,
`requested_root_count`, `completed_root_count`, `total_cells_read`,
`total_read_seconds`, `total_compute_seconds`. The `completed_root_count` and
`total_*` aggregates are recomputed from `/gas/root_metrics`.

### `/gas` (attrs + datasets)
Attrs: `row_alignment = /catalog/halo`, `completion_authority =
/gas/root_metrics`.

- **`/gas/summary`** — main results table, **row-aligned 1:1 with
  `/catalog/halo`** (full catalog length). Unprocessed rows are `NaN` (float) /
  `0` (int/bool); see §2.1.
- **`/gas/processed`** — catalog-length boolean mask; `True` where a row has been
  computed. (Multi-root runs also process a root's descendant/subhalo rows.)
- **`/gas/root_metrics`** — append-only per-root log and the **completion
  authority**: `root_id`, `descendant_count`, `cpu_count`, `cells_read`,
  `envelope_radius_code`, `read_seconds`, `compute_seconds`. A run is restartable
  — roots already present here are skipped on resume.

### 2.1 `/gas/summary` (58 fields)

**Identifiers** — `id` (this halo), `root_id` (its top-level root).

**Overlap / resolution diagnostics** (mostly on the `r*` aperture) — `ncells`,
`ncells_center`, `ncells_overlap`, `overlap_resolved` (bool),
`radius_to_min_cell`, `overlap_relative_bound`.

**Gas masses (M⊙)** — total `mgas`, cold (`T<10⁴ K`) `mcold`, dense (`T<10⁴ K` &
`n_H>5/cc`) `mdense`, each at six apertures: bare (= `r*`), `_r50`, `_r90`,
`_rvir`, `_r200`, `_r500`.

**Spherical overdensity** — `r200`, `r500` (⚠️ **code units, box fraction**),
`m200`, `m500` (M⊙, total enclosed). GasMaker uses threshold-crossing
interpolation (differs by design from RUR's nearest-shell selection — TASK-10
Class B diagnostic).

**Component masses (M⊙)** — `mdm_rvir/r200/r500` (dark matter),
`mstar_rvir/r200/r500` (stars).

**Metallicity & chemistry** (mass-weighted, dimensionless, on `r*`) — `metal_gas`
and per-element `H_gas`, `O_gas`, `Fe_gas`, `Mg_gas`, `C_gas`, `N_gas`, `Si_gas`,
`S_gas`, `D_gas` (filled when the corresponding cell fields exist).

**Gas kinematics (km/s)** — `vrot_gas`, `sig3d_gas`, `sigcyl_gas` at `r*`, plus
`_r50` / `_r90` apertures. (TASK-10 Class C diagnostic: axis/dispersion
definitions differ from RUR's scalar `vsig_gas`.)

**Gas angular momentum (10¹¹ M⊙·Mpc·km/s, on `r*`)** — `Lx_gas`, `Ly_gas`,
`Lz_gas`.

---

## 3. Conventions & notes

- **Aperture key:** `r*` = stellar radius · `r50`/`r90` = stellar half-/90%-mass
  radii · `rvir` = virial radius · `r200`/`r500` = spherical-overdensity radii.
- **Joining:** both files are row-aligned to `/catalog/halo`; match by `id`
  (order-independent). A row present in the catalog but unprocessed by GasMaker
  has `NaN`/`0` in `gas_bricks`.
- **Starless halos:** stellar-aperture GasMaker fields (bare `r*`, `_r50`,
  `_r90`, chemistry, kinematics) are left `NaN`; `_rvir` and SO quantities are
  still computed.
- **Validation:** GasMaker gas/particle masses, metallicity, and chemistry agree
  with the RUR reference at machine precision on a stratified sample (TASK-10
  Class A pass); `r200/r500` (Class B) and gas kinematics (Class C) are
  reported as documented definition differences. See `WHATS_NEW.md`.
