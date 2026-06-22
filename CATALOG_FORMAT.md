# Catalog Format

Field-by-field reference for the pipeline's HDF5 outputs, in table form. The
pipeline produces up to two files, **row-aligned and joined by `id`**: the
HaloMaker catalog (`tree_bricks*.h5`, always) and the optional GasMaker output
(`gas_bricks*.h5`). See `README.md` for how to run them and `WHATS_NEW.md` for
feature/validation context.

## 1. Output Files

| File | Producer | Required? | Main content | Row alignment / join key |
|---|---:|---:|---|---|
| `tree_bricks{iout:05d}.h5` | HaloMaker | Always | halo, galaxy, stellar population, photometry, member particles, run metadata | canonical `/catalog/halo`; join by `id` |
| `gas_bricks{iout:05d}.h5` | GasMaker | Optional post-processing | gas masses, gas chemistry, gas kinematics, SO component masses, processing diagnostics | same row order as `/catalog/halo`; join by `id` |

## 2. Global Unit Conventions

| Scope | Quantity type | Unit / convention | Notes |
|---|---|---|---|
| HaloMaker catalog | positions and radii | physical Mpc | applies to `px/py/pz`, `px*/py*/pz*`, `r`, `r*`, `r50`, `r90`, `rvir` |
| HaloMaker catalog | masses | `10^11 Msun` | HaloMaker convention |
| HaloMaker catalog | velocities | km/s | halo and stellar velocity fields |
| HaloMaker catalog | temperatures | K | e.g. `tvir` |
| GasMaker | gas / particle masses | `Msun` | not `10^11 Msun` |
| GasMaker | metallicity and chemistry | dimensionless mass fraction | mass-weighted where applicable |
| GasMaker | gas kinematics | km/s | `vrot_gas`, `sig3d_gas`, `sigcyl_gas` families |
| GasMaker | gas angular momentum | `10^11 Msun Mpc km/s` | `Lx_gas`, `Ly_gas`, `Lz_gas` |
| GasMaker | `r200`, `r500` | code units, box fraction `[0, 1)` | multiply by `box_physical_mpc = header.boxsize2 * aexp` to get physical Mpc |

## 3. HaloMaker File: `tree_bricks{iout:05d}.h5`

### 3.1 Top-Level Layout

| HDF5 path | Type | Row-aligned? | Purpose |
|---|---|---:|---|
| `/catalog/halo` | structured table | Yes | main halo/galaxy catalog |
| `/photometry/CB07` | structured table | Yes | rest-frame photometry from CB07 SSP |
| `/photometry/BC03` | structured table | Yes | rest-frame photometry from BC03 SSP |
| `/photometry/FSPS` | structured table | Yes | rest-frame photometry from FSPS SSP |
| `/member/index` | array | Yes | per-halo offset into flat member arrays |
| `/member/count` | array | Yes | per-halo member count |
| `/member/pids` | flat array | No, offset-based | member particle IDs |
| `/member/pos` | flat array, optional | No, offset-based | member particle positions when `dump_DMs` is enabled |
| `/member/vel` | flat array, optional | No, offset-based | member particle velocities when `dump_DMs` is enabled |
| `/member/mass` | flat array, optional | No, offset-based | member particle masses when `dump_DMs` is enabled |
| `/header` | attrs/group | N/A | run, box, cosmology, count metadata |
| `/input` | group | N/A | input configuration snapshot |

### 3.2 `/catalog/halo` Field Groups

The source document describes `/catalog/halo` as a 72-field structured table.

| Category | Fields | Unit / type | Meaning |
|---|---|---|---|
| Identity and hierarchy | `id`, `timestep`, `nmem`, `ndm`, `nstar`, `nbsub`, `hosthalo`, `hostsub`, `level`, `nextsub` | integer / counts | halo identity, member counts, host/subhalo relations, tree level |
| Position | `px`, `py`, `pz` | physical Mpc | halo center |
| Stellar center | `px*`, `py*`, `pz*` | physical Mpc | stellar center |
| Velocity | `vx`, `vy`, `vz` | km/s | halo velocity |
| Angular momentum | `Lx`, `Ly`, `Lz` | catalog convention | total angular momentum |
| Stellar angular momentum | `Lx*`, `Ly*`, `Lz*` | catalog convention | stellar angular momentum |
| Shape | `sha`, `shb`, `shc` | catalog convention | principal-axis lengths |
| Mass | `m`, `mdm`, `m*` | `10^11 Msun` | total, dark matter, and stellar mass |
| Radii | `r`, `r*`, `r50`, `r90`, `rvir` | physical Mpc | extent, stellar radius, stellar half/90%-mass radii, virial radius |
| Stellar population | `age`, `metal` | mass-weighted | stellar age and metallicity summaries |
| Star formation rate, 100 Myr | `SFR`, `SFR_r50`, `SFR_r90` | catalog convention | SFR in stellar apertures |
| Star formation rate, 10 Myr | `SFR10`, `SFR10_r50`, `SFR10_r90` | catalog convention | short-timescale SFR in stellar apertures |
| Spin and dispersions | `spin`, `sigma`, `sigma_dm`, `sigma*` | catalog convention | spin and velocity dispersion summaries |
| Stellar kinematics | `vrot`, `sig3d`, `sigcyl`, `vrot_r50`, `sig3d_r50`, `sigcyl_r50`, `vrot_r90`, `sig3d_r90`, `sigcyl_r90` | km/s-style kinematic quantities | stellar kinematics in `r*`, `r50`, `r90` apertures |
| Energy and virial | `ek`, `ep`, `et`, `mvir`, `tvir`, `cvel` | mixed; `tvir` in K | energy and virial diagnostics |
| NFW / density profile | `rho_0`, `r_c`, `cNFW`, `cNFWerr`, `vmaxcir`, `rmaxcir`, `inslope`, `inslopeerr` | catalog convention | profile fit and inner-slope diagnostics |
| Contamination | `mcontam` | catalog mass convention | low-resolution / contaminant mass |

Notes:

| Topic | Detail |
|---|---|
| `inslope` / `inslopeerr` | DM inner-density-slope fit; roundoff-degenerate shells are dropped from the fit so the value is platform-stable (residual within field-policy tolerance) |
| `r50` / `r90` | stellar half-/90%-mass radii in physical Mpc, not GasMaker code units |
| `level == 1` | top-level root halos for GasMaker root selection |

### 3.3 `/photometry/{CB07,BC03,FSPS}`

| Field group | Fields | Unit / convention | Notes |
|---|---|---|---|
| SDSS bands | `umag`, `gmag`, `rmag`, `imag`, `zmag` | AB magnitude | rest-frame stellar photometry |
| Johnson bands | `Umag`, `Bmag`, `Vmag`, `Kmag` | Vega magnitude | FSPS `Kmag` uses 2MASS Ks; see group metadata |
| r-band luminosity-weighted properties | `age_r`, `metal_r`, `r50_r`, `r90_r` | model/catalog convention | see `SSP_MODELS.md` |

## 4. GasMaker File: `gas_bricks{iout:05d}.h5`

### 4.1 Top-Level Layout

| HDF5 path | Type | Shape / length | Purpose |
|---|---|---:|---|
| `/header` | attrs/group | N/A | schema, run configuration, aggregate progress |
| `/gas` | group with attrs | N/A | GasMaker result namespace |
| `/gas/summary` | structured table | full `/catalog/halo` length | main gas result table; row-aligned 1:1 with HaloMaker catalog |
| `/gas/processed` | bool array | full `/catalog/halo` length | marks rows that have been computed |
| `/gas/root_metrics` | append-only structured table | one row per completed root | restart log and completion authority |

### 4.2 `/header` Attributes

| Attribute | Meaning | Notes |
|---|---|---|
| `schema_version` | GasMaker output schema version | compatibility gate on resume |
| `source_catalog` | path to source HaloMaker catalog | provenance |
| `run_mode` | GasMaker run mode | e.g. multi-root / single-root |
| `radius_field` | halo radius field used for root envelope | usually catalog radius field |
| `root_padding` | root read-envelope padding | GasMaker configuration |
| `cell_boundary_overlap` | whether boundary-overlap candidate cells are included | GasMaker configuration |
| `ownership` | ownership convention | current value documented as inclusive |
| `overlap_method` | overlap estimator | current method: adaptive radius fraction |
| `overlap_threshold` | adaptive overlap threshold | GasMaker configuration |
| `overlap_depth` | max adaptive subdivision depth | GasMaker configuration |
| `overlap_tolerance` | diagnostic tolerance for overlap bounds | GasMaker configuration |
| `requested_root_ids` | root IDs requested by the run | int array |
| `requested_root_count` | number of requested roots | aggregate |
| `completed_root_count` | number of completed roots | recomputed from `/gas/root_metrics` |
| `total_cells_read` | total cells read over completed roots | recomputed from `/gas/root_metrics` |
| `total_read_seconds` | total read time | recomputed from `/gas/root_metrics` |
| `total_compute_seconds` | total compute time | recomputed from `/gas/root_metrics` |

### 4.3 `/gas` Attributes and Datasets

| Item | Value / fields | Meaning |
|---|---|---|
| attr `row_alignment` | `/catalog/halo` | `summary` and `processed` use the same row order as the HaloMaker catalog |
| attr `completion_authority` | `/gas/root_metrics` | completed roots are determined from this append-only table |
| dataset `summary` | `SUMMARY_DTYPE` | main GasMaker results; unprocessed float fields are `NaN`, non-floats are `0` |
| dataset `processed` | boolean mask | `True` for rows written by GasMaker; multi-root processing includes descendants/subhalos |
| dataset `root_metrics` | `root_id`, `descendant_count`, `cpu_count`, `cells_read`, `envelope_radius_code`, `read_seconds`, `compute_seconds` | one row per completed root; roots present here are skipped on resume |

### 4.4 `/gas/summary` Field Groups

The source document describes `/gas/summary` as a 58-field structured table.

| Category | Fields | Unit / type | Aperture / convention |
|---|---|---|---|
| Identifiers | `id`, `root_id` | integer | `id` is this halo; `root_id` is its top-level root |
| Overlap diagnostics | `ncells`, `ncells_center`, `ncells_overlap`, `overlap_resolved`, `radius_to_min_cell`, `overlap_relative_bound` | counts, bool, floats | mostly on the stellar `r*` aperture (per-field meaning below) |

Overlap-diagnostic field meanings:

| Field | Meaning |
|---|---|
| `ncells` | candidate cells for the `r*` aperture |
| `ncells_center` | cells counted under the strict center-in rule |
| `ncells_overlap` | cells with a positive fractional-overlap weight |
| `overlap_resolved` | whether the overlap uncertainty bound is within tolerance |
| `radius_to_min_cell` | aperture radius relative to the smallest cell size |
| `overlap_relative_bound` | relative width between the lower/upper overlap-mass bounds |
| Total gas mass | `mgas`, `mgas_r50`, `mgas_r90`, `mgas_rvir`, `mgas_r200`, `mgas_r500` | `Msun` | bare field means `r*`; suffixes indicate aperture |
| Cold gas mass | `mcold`, `mcold_r50`, `mcold_r90`, `mcold_rvir`, `mcold_r200`, `mcold_r500` | `Msun` | cold definition: `T < 10^4 K` |
| Dense gas mass | `mdense`, `mdense_r50`, `mdense_r90`, `mdense_rvir`, `mdense_r200`, `mdense_r500` | `Msun` | dense definition: `T < 10^4 K` and `n_H > 5/cc` |
| Spherical overdensity radii | `r200`, `r500` | code units, box fraction `[0, 1)` | multiply by `box_physical_mpc` for physical Mpc |
| Spherical overdensity masses | `m200`, `m500` | `Msun` | total enclosed mass at SO radii |
| Dark matter component masses | `mdm_rvir`, `mdm_r200`, `mdm_r500` | `Msun` | component mass within aperture |
| Stellar component masses | `mstar_rvir`, `mstar_r200`, `mstar_r500` | `Msun` | component mass within aperture |
| Metallicity | `metal_gas` | dimensionless mass fraction | mass-weighted on `r*` |
| Chemistry | `H_gas`, `O_gas`, `Fe_gas`, `Mg_gas`, `C_gas`, `N_gas`, `Si_gas`, `S_gas`, `D_gas` | dimensionless mass fraction | mass-weighted on `r*`; filled when matching cell fields exist |
| Gas kinematics, `r*` | `vrot_gas`, `sig3d_gas`, `sigcyl_gas` | km/s | Class C diagnostic in TASK-10 |
| Gas kinematics, `r50` | `vrot_gas_r50`, `sig3d_gas_r50`, `sigcyl_gas_r50` | km/s | stellar `r50` aperture |
| Gas kinematics, `r90` | `vrot_gas_r90`, `sig3d_gas_r90`, `sigcyl_gas_r90` | km/s | stellar `r90` aperture |
| Gas angular momentum | `Lx_gas`, `Ly_gas`, `Lz_gas` | `10^11 Msun Mpc km/s` | on `r*` |

### 4.5 GasMaker Aperture Key

| Field pattern | Aperture |
|---|---|
| bare gas fields, e.g. `mgas`, `metal_gas`, `vrot_gas`, `Lx_gas` | stellar `r*` |
| `_r50` | stellar half-mass radius `r50` |
| `_r90` | stellar 90%-mass radius `r90` |
| `_rvir` | virial radius `rvir` |
| `_r200` | GasMaker spherical-overdensity `r200` |
| `_r500` | GasMaker spherical-overdensity `r500` |

## 5. Join, Restart, and Missing-Value Rules

| Rule | Detail |
|---|---|
| Primary join | use `id`; row order is also aligned to `/catalog/halo` |
| Order-independent matching | compare by `id` if row order is uncertain |
| GasMaker unprocessed rows | floats are `NaN`; integer/bool fields are `0` |
| Starless halos | stellar-aperture GasMaker fields (`r*`, `_r50`, `_r90`, chemistry, kinematics) remain `NaN`; `_rvir` and SO quantities are still computed |
| Restart authority | `/gas/root_metrics` decides completed roots |
| Resume behavior | roots already present in `/gas/root_metrics.root_id` are skipped |
| Aggregate progress | `/header.completed_root_count` and `/header.total_*` are recomputed from `/gas/root_metrics` |

## 6. Validation / Caveat Summary

| Topic | Status | Practical interpretation |
|---|---|---|
| GasMaker Class A fields | PASS on stratified TASK-10 sample | gas/particle masses, metallicity, and chemistry match RUR reference at machine precision under matched aperture definitions |
| Galaxy aperture scaling bug claim | rejected / acquitted | catalog `r50`/`r90` are physical Mpc; apples-to-apples comparison removes the apparent explosion |
| `r200` / `r500` | Class B definition difference | GasMaker uses threshold-crossing interpolation; RUR uses nearest-shell style selection |
| Gas kinematics | Class C diagnostic | axis/dispersion definitions differ from RUR scalar `vsig_gas` |
| Stellar chemistry correctness | follow-up candidate | not the main blocker for TASK-10 Class A GasMaker gate |

