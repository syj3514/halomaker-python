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
| All outputs | positions and radii | RAMSES code units, box fraction `[0, 1)` | positions are periodic; recover physical Mpc with `x_phys = x_code * box_physical_mpc` using the desired centered convention |
| HaloMaker catalog | masses | `Msun` | `10^11 Msun` legacy scaling has been removed from HDF5 output |
| HaloMaker catalog | velocities | km/s | halo and stellar velocity fields |
| HaloMaker catalog | angular momentum | `Msun Mpc km/s` | physical Mpc remains the length factor in angular momentum |
| HaloMaker catalog | energy | `Msun (km/s)^2` | kinetic, potential, total energies |
| HaloMaker catalog | density | `Msun/kpc^3` | `rho_0`; this is the only physical-volume exception to code-unit length output |
| HaloMaker catalog | stellar chemistry | dimensionless mass fraction | mass-weighted over the same full stellar-member aperture as `metal`; Ra3/no-chem inputs are `NaN` |
| HaloMaker catalog | temperatures | K | e.g. `tvir` |
| GasMaker | gas / particle masses | `Msun` | not `10^11 Msun` |
| GasMaker | metallicity and chemistry | dimensionless mass fraction | mass-weighted where applicable |
| GasMaker | gas kinematics | km/s | `vrot_gas`, `sig3d_gas`, `sigcyl_gas` families |
| GasMaker | gas angular momentum | `Msun Mpc km/s` | `Lx_gas`, `Ly_gas`, `Lz_gas` |
| Metadata | unit marker | `/header.units_version = halomaker_units_v2` | compound datasets also carry JSON `field_units` attrs |
| Metadata | box size | `box_comoving_mpc = boxsize2`, `box_physical_mpc = aexp * boxsize2` | legacy `boxsize2` and `aexp` are retained |

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
| `/member/pos` | flat array, optional | No, offset-based | member particle positions when `dump_members` is enabled (legacy `dump_DMs`) |
| `/member/vel` | flat array, optional | No, offset-based | member particle velocities when `dump_members` is enabled (legacy `dump_DMs`) |
| `/member/mass` | flat array, optional | No, offset-based | member particle masses when `dump_members` is enabled (legacy `dump_DMs`) |
| `/header` | attrs/group | N/A | run, box, cosmology, count metadata |
| `/input` | group | N/A | input configuration snapshot |

**Per-snapshot cosmology provenance** (`/header` and `/input` attrs): `H_f` is the
effective Hubble constant used for this snapshot; `info_H0` is the `H0` read from the
snapshot's `info_*.txt`; `H_f_source` records how `H_f` was set — `info_H0` (the
snapshot's authoritative H0, used when it differs materially from the config value,
`rtol>1e-6`), `config_equivalent_info` (config kept because it matches `info_H0`
within `rtol=1e-6`), or `config_fallback` (no `H0` in the header). `snapshot_dir`
records the source directory. When several snapshots are processed in one run and
their `numstep` labels collide, the output filename is disambiguated as
`tree_bricks{numstep}_{tag}.h5` (source-dir tag, then `_L{n}` line index); a single
snapshot keeps the plain `tree_bricks{numstep:05d}.h5` name. An optional 5th
`inputfiles` field sets an explicit per-line prefix that is always appended.

### 3.2 `/catalog/halo` Field Groups

The source document describes `/catalog/halo` as an 81-field structured table.

| Category | Fields | Unit / type | Meaning |
|---|---|---|---|
| Identity and hierarchy | `id`, `timestep`, `nmem`, `ndm`, `nstar`, `nbsub`, `hosthalo`, `hostsub`, `level`, `nextsub` | integer / counts | halo identity, member counts, host/subhalo relations, tree level |
| Position | `px`, `py`, `pz` | code unit `[0,1)` | halo center |
| Stellar center | `px*`, `py*`, `pz*` | code unit `[0,1)` | stellar center |
| Velocity | `vx`, `vy`, `vz` | km/s | halo velocity |
| Angular momentum | `Lx`, `Ly`, `Lz` | `Msun Mpc km/s` | total angular momentum |
| Stellar angular momentum | `Lx*`, `Ly*`, `Lz*` | `Msun Mpc km/s` | stellar angular momentum |
| Shape | `sha`, `shb`, `shc` | code unit `[0,1)` | principal-axis lengths, not axis ratios |
| Mass | `m`, `mdm`, `m*` | `Msun` | total, dark matter, and stellar mass |
| Radii | `r`, `r*`, `r50`, `r90`, `rvir` | code unit `[0,1)` | extent, stellar radius, stellar half/90%-mass radii, virial radius |
| Stellar population | `age`, `metal` | `Gyr`, mass fraction | stellar age and metallicity summaries |
| Stellar chemistry | `H_star`, `O_star`, `Fe_star`, `Mg_star`, `C_star`, `N_star`, `Si_star`, `S_star`, `D_star` | dimensionless mass fraction | mass-weighted over all stellar members, same aperture as `metal`; `NaN` when the snapshot has no stellar chemistry descriptor |
| Star formation rate, 100 Myr | `SFR`, `SFR_r50`, `SFR_r90` | `Msun/yr` | SFR in stellar apertures |
| Star formation rate, 10 Myr | `SFR10`, `SFR10_r50`, `SFR10_r90` | `Msun/yr` | short-timescale SFR in stellar apertures |
| Spin and dispersions | `spin`, `sigma`, `sigma_dm`, `sigma*` | dimensionless; km/s for dispersions | spin and velocity dispersion summaries |
| Stellar kinematics | `vrot`, `sig3d`, `sigcyl`, `vrot_r50`, `sig3d_r50`, `sigcyl_r50`, `vrot_r90`, `sig3d_r90`, `sigcyl_r90` | km/s | stellar kinematics in `r*`, `r50`, `r90` apertures |
| Energy and virial | `ek`, `ep`, `et`, `mvir`, `tvir`, `cvel` | `Msun (km/s)^2`; `mvir` Msun; `tvir` K; `cvel` km/s | energy and virial diagnostics |
| NFW / density profile | `rho_0`, `r_c`, `cNFW`, `cNFWerr`, `vmaxcir`, `rmaxcir`, `inslope`, `inslopeerr` | `rho_0` Msun/kpc^3; `r_c/rmaxcir` code unit; `vmaxcir` km/s; others dimensionless | profile fit and inner-slope diagnostics |
| Contamination | `mcontam` | `Msun` | low-resolution / contaminant mass |

Notes:

| Topic | Detail |
|---|---|
| `inslope` / `inslopeerr` | DM inner-density-slope fit; roundoff-degenerate shells are dropped from the fit so the value is platform-stable (residual within field-policy tolerance) |
| `r50` / `r90` | stellar half-/90%-mass radii in code units |
| `level == 1` | top-level root halos for GasMaker root selection |
| `*_star` chemistry | Ra4 snapshots with `chem_*` particle descriptors populate these fields; Ra3 or missing-descriptor snapshots keep them as `NaN` |

### 3.3 `/photometry/{CB07,BC03,FSPS}`

| Field group | Fields | Unit / convention | Notes |
|---|---|---|---|
| SDSS bands | `umag`, `gmag`, `rmag`, `imag`, `zmag` | AB magnitude | rest-frame stellar photometry |
| Johnson bands | `Umag`, `Bmag`, `Vmag`, `Kmag` | Vega magnitude | FSPS `Kmag` uses 2MASS Ks; see group metadata |
| r-band luminosity-weighted properties | `age_r`, `metal_r`, `r50_r`, `r90_r` | Gyr, mass fraction, code-unit radii | see `SSP_MODELS.md` |

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
| `units_version` | GasMaker output unit schema | current value: `halomaker_units_v2` |
| `source_catalog_units_version` | unit schema of the source catalog | used to interpret legacy vs current catalog position/radius fields |

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
| Gas angular momentum | `Lx_gas`, `Ly_gas`, `Lz_gas` | `Msun Mpc km/s` | on `r*` |

**Fixed physical definitions (not configurable):** the cold (`T < 10^4 K`) and
dense (`T < 10^4 K` and `n_H > 5/cc`) thresholds, the spherical-overdensity
contrasts (200× and 500× the critical density `ρ_crit`), and the 9-element
chemistry set (`H, O, Fe, Mg, C, N, Si, S, D`, shared with the catalog `*_star`
fields via `chem_species.py`) are fixed in source. Changing them is a
schema/feature change, not a run-time option.

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
| GasMaker Class A fields | PASS on stratified TASK-10 sample; A3 all-root extends it | gas/particle masses, metallicity, and chemistry match RUR reference at machine precision under matched aperture definitions |
| GasMaker gas masses at **low mass** (TASK-17) | characterized low-mass systematic (not a bug) | **High-mass / well-resolved halos: GasMaker ≈ RUR at machine-to-sub-% precision.** For **low-mass halos (npart ≲ 300)** the *diffuse* gas-cell quantities (`mgas*`, `mcold*`, gas metals) can differ from RUR by tens-of-% in the p99 tail: GasMaker weights boundary cells by **fractional volume overlap**, RUR uses **cell-center-in (binary) + a <8-cell radius relaxation**, and with only a few cells one boundary cell dominates. Particle-based (`mstar*`, `mdm*`) and *dense* gas match at machine precision regardless of mass. Median agreement stays ~0.3%. |
| Galaxy aperture scaling bug claim | rejected / acquitted | legacy catalog `r50`/`r90` were physical Mpc; in `halomaker_units_v2` they are stored as code-unit radii, so downstream readers must branch on `units_version` |
| `r200` / `r500` | Class B definition difference | GasMaker uses threshold-crossing interpolation; RUR uses nearest-shell style selection |
| Gas kinematics | Class C diagnostic | axis/dispersion definitions differ from RUR scalar `vsig_gas` |
| Stellar chemistry correctness | TASK-13 schema | `*_star` fields are catalog-level stellar mass-fraction averages when Ra4 particle descriptors provide per-element chemistry; otherwise `NaN` fallback applies |
