# HaloMaker Python

Python + Fortran implementation of HaloMaker / AdaptaHOP for RAMSES
snapshots. The Fortran extensions provide the memory-sensitive neighbor and
structure-tree routines. The pipeline runs in **full-box** mode; the legacy
zoom-in path has been removed.

This release is staged from the dev3 memory-optimized implementation. The
public filename remains canonical (`compute_adaptahop.f90`) so existing scripts
do not need development branch names.

Compared with the local reference implementation used during development,
dev3 reduced the peak RSS of the 39990 full-box test from about 104 GiB to
about 24 GiB while preserving catalog-level results. Runtime stayed close to
the reference path on the tested workloads. The experimental dev4 zoom-in
read-compact path is not included in this release staging copy.

## Requirements

- Linux
- Python 3.10 or newer
- `gfortran`
- OpenMP runtime (`libgomp`)
- RAMSES snapshot files

## Python Environment

Choose either `uv`, conda, or another Python environment manager. Python 3.10
or newer is supported. The examples below use Python 3.12, but you may replace
it with another supported version.

### uv

Create a project environment with your preferred Python version:

```bash
uv venv --python 3.12
source .venv/bin/activate
uv sync
```

`uv sync` creates a local `uv.lock` file. This lock is intentionally not
distributed so that users can resolve dependencies for their chosen Python
version and platform.

### conda

Create and activate an environment with your preferred Python version:

```bash
conda create -n halomaker-python python=3.12
conda activate halomaker-python
python -m pip install -e .
```

Alternatively, let conda resolve a supported Python version from the included
environment file:

```bash
conda env create -f environment.yml
conda activate halomaker-python
```

## Build

Compile both f2py extensions:

```bash
bash build.sh
```

This creates:

- `compute_adaptahop*.so`

Build and run from the same activated Python environment. If you do not want
to activate the environment, pass the interpreter explicitly:

```bash
PYTHON=.venv/bin/python bash build.sh
PYTHON=.venv/bin/python bash run.sh
```

## Configure

For production runs, create a dedicated ignored run directory and copy the
example inputs there:

```bash
mkdir -p runs/my_run
cd runs/my_run
cp ../../examples/input_HaloMaker.dat.example input_HaloMaker.dat
cp ../../examples/inputfiles_HaloMaker.dat.example inputfiles_HaloMaker.dat
```

Edit `inputfiles_HaloMaker.dat` so that each active line points to an existing
RAMSES snapshot. Each line is `'<snapshot dir>' <format> <worker count>
<output number> [prefix]`; the 5th `prefix` field is optional. Lines are
processed in order, so several snapshots (even from different simulations) can
be chained in one run. The pipeline runs in periodic full-box mode; the legacy
zoom-in mode was removed, so a stray `zoomin` key in `input_HaloMaker.dat` is
ignored.

For RAMSES snapshots, `lbox`, `omega_f`, `lambda_f`, and `H_f` are read
per-snapshot from the RAMSES header/`info_*.txt` during `read_data()`, so they
are optional in `input_HaloMaker.dat`. `H_f` honours the config value when it
matches the snapshot's `H0` to within `rtol=1e-6` and otherwise uses the
snapshot's authoritative `H0` (recorded as `H_f_source`/`info_H0` in the output
header) — this lets a single run chain simulations with different `H0`. If
config values are omitted, box-size- and cosmology-dependent quantities are
finalized only after the snapshot header is read.

When chained snapshots share an output number, output filenames are
auto-disambiguated (`tree_bricks{n}_{tag}.h5`) so none are overwritten; a single
snapshot keeps the plain `tree_bricks{n:05d}.h5` name.

Recommended print levels:

- `verbose = .false.`, `megaverbose = .false.`: compact production log.
- `verbose = .true.`, `megaverbose = .false.`: normal diagnostic log with
  major counts and timings.
- `megaverbose = .true.`: development log with detailed Fortran and memory
  tracking output.

## Run

```bash
bash ../../run.sh
```

Run directories should live under `runs/<name>/` so local configs, logs, and
large outputs do not accumulate in the release root. `runs/` is ignored by Git.
The main HDF5 catalog output is written in the current run directory as
`tree_bricks*.h5`.

### Parallelism (`nbPes`)

The third field of each line in `inputfiles_HaloMaker.dat` is `nbPes`, the
number of worker processes/threads used for both the Python stages (particle
reading, density, halo properties) and the Fortran OpenMP stages (neighbor
search, mean density, saddle connections, node tree). Larger `nbPes` speeds up
most phases, but total speedup **saturates** (≈12× at 32 cores on a 19M-particle
box): the `create nodes` step of the structure tree is effectively serial and does
not scale with `nbPes`, so it dominates wall time on very large boxes. This is an
algorithmic limit, not a tuning gap — the node tree is built by a sequential
density-percolation of the most massive structure (a data-dependency chain), and
on the 07206 box the single largest halo alone is ~69% of that step, capping any
same-output parallel speedup at ≈1.45×. Reducing it requires an output-changing
algorithm change, not more threads.

Setting `optimize_nodes = .true.` in `input_HaloMaker.dat` (default `.false.` =
classic path) enables such a path: large components are summarized through a
per-group **blocked prefix-moment index** (replacing the per-level rescan with
O(log n) lookups), while small (`< 1000` particles) and threshold-near components
fall back to the exact legacy computation, so the catalog stays **bit-for-bit on
the classic output** (verified `real_regression=0` on the frozen goldens). On a
39990 full box this cuts the `create nodes` step from ~1428 s to ~102 s (≈14×) at
**+4.8 MB** peak memory.

### Output units (breaking change: `halomaker_units_v2`)

Catalog and GasMaker outputs now use one consistent unit system: **masses in
`Msun`**, **positions / radii / shape axes in RAMSES code units `[0, 1)`**,
angular momentum in `Msun Mpc km/s`, energies in `Msun (km/s)^2`, `rho_0` in
`Msun/kpc^3`, velocities in `km/s`, ages in `Gyr`. The previous `10^11 Msun`
mass scaling and physical-Mpc positions/radii are **no longer** in the HDF5
output. Files are tagged `/header.units_version = "halomaker_units_v2"`, carry
`box_comoving_mpc` / `box_physical_mpc` attributes, and per-field `field_units`
JSON attributes. Recover physical Mpc with `x_phys = x_code * box_physical_mpc`.
Analysis scripts that assumed the old units must branch on `units_version`. See
**`docs/CATALOG_FORMAT.md`** for the full field-by-field unit table.

Each run computes intrinsic rest-frame stellar photometry with all bundled SSP
models and writes the aligned results under:

- `/photometry/CB07/data`
- `/photometry/BC03/data`
- `/photometry/FSPS/data`

The datasets contain SDSS `ugriz`, Johnson `UBV`, `Kmag`, and r-band
luminosity-weighted stellar age, metallicity, `r50`, and `r90`. FSPS uses the
2MASS Ks response for its `Kmag`; this distinction is recorded in the group
metadata. See `docs/SSP_MODELS.md` for model definitions and data provenance.

For Ra4 snapshots that store per-element stellar chemistry, the catalog also
includes mass-weighted per-element stellar abundances `H_star, O_star, Fe_star,
Mg_star, C_star, N_star, Si_star, S_star, D_star` (mass fraction); snapshots
without stellar chemistry leave these `NaN`. Set `dump_members` (which writes
per-member pos/vel/mass for all members; the old keys `dump_DMs`/`dump_stars`
still work) to also export the flat `/member` arrays. See `docs/CATALOG_FORMAT.md`.

If a run is interrupted with Ctrl-C, killed by a scheduler, or leaves Python
`forkserver` / `resource_tracker` processes behind, inspect and clean runtime
leftovers with:

```bash
bash clean_runtime.sh            # dry-run: list only, delete nothing
bash clean_runtime.sh --orphan   # reclaim provably-dead orphan shm only
bash clean_runtime.sh --force    # kill this repo's runtime procs, then reclaim
```

Shared-memory files are tagged with the owning run's `(pid, start-time)`, so the
helper classifies each as **DEAD** (owning process gone / PID recycled / zombie),
**LIVE** (process still running), or **UNKNOWN** (legacy/generic, e.g. `psm_*`).
The bare command is a true dry-run (deletes nothing). `--orphan` removes only
DEAD orphans (and their manifest), never touching LIVE or UNKNOWN files. `--force`
first terminates this repository's own runtime processes, re-classifies, then
removes the now-DEAD segments plus UNKNOWN/legacy files. A LIVE run owned by
another working tree or agent is never killed here, so it is never removed — the
"never touch a live run" invariant holds even under `--force`. (DEAD is judged on
the owning PID; a hard-killed run may briefly leave a fork child mapping the
segment, but unlinking only drops the name and leaves existing mappings valid.)

## GasMaker (gas post-processor)

GasMaker is a separate tool that reads an existing HaloMaker/GalaxyMaker catalog
plus the RAMSES AMR/hydro data and adds **gas properties** per halo and galaxy:
total / cold (T<10⁴ K) / dense gas mass, gas metallicity and per-element
chemistry, gas kinematics and angular momentum (within r\*, r50, r90, r_vir),
and spherical-overdensity quantities (r200/m200, r500/m500 plus enclosed
DM/star/gas masses). It is restartable.

```bash
# specific roots:
python GasMaker.py <catalog.h5> <ramses_repo> <iout> --root-ids 3,11,15
# or all top-level halos:
python GasMaker.py <catalog.h5> <ramses_repo> <iout> --roots all
```

GasMaker can also be run from HaloMaker-style config files. Put
`input_GasMaker.dat` (parameters) and `inputfiles_GasMaker.dat` (one job per
line) in the working directory, then run:

```bash
python GasMaker.py
```

See `examples/input_GasMaker.dat.example` and
`examples/inputfiles_GasMaker.dat.example`. `inputfiles_GasMaker.dat` follows
the same whitespace-splitting convention as HaloMaker's inputfiles parser:
leading `#`/`!` comments and blank lines are ignored, quotes around paths are
decorative, and paths with spaces are not supported. Multiple lines are
processed sequentially; duplicate output filenames are rejected before running.
For production use, keep GasMaker config files, logs, and outputs in a dedicated
`runs/<name>/` directory as well.

Long runs show run-wide progress: a startup banner, stage timings, a per-root
progress display, and an end summary (wall time, read/compute totals, slowest
root). Progress goes to **stderr** while stdout keeps only the final summary
lines, so `> result.log` stays machine-parseable and `> run.log 2>&1` captures
the full timeline. Control it with `--progress {auto,bar,plain,quiet}` (config
key `progress`; `auto` picks a tqdm bar on a terminal and plain timestamped
lines when redirected) and `--progress-every N` / `progress_every`.

The output defaults to `gas_bricks{iout:05d}.h5` (override with `--output`). It
is row-aligned with the catalog and joined by `id`; see **`docs/CATALOG_FORMAT.md`**
for the full field list of both outputs. A `gas_bricks` file is derived from one
specific catalog (recorded as `source_catalog` in its `/header`): the `_rvir`
apertures and SO fields inherit that catalog's `rvir`, so regenerate it whenever
the catalog changes (e.g. after the TASK-21 virial refinement). The fixed physical thresholds (cold
`T<10⁴ K`, dense `n_H>5/cc`, spherical overdensities `200`/`500 ρ_crit`) are
recorded in the output `/header` for provenance.

GasMaker writes its output one root at a time, which trips HDF5's file locking on
some network/parallel filesystems (Lustre/NFS) with `BlockingIOError [errno 11]`.
Since GasMaker is the sole writer, it sets `HDF5_USE_FILE_LOCKING=FALSE` by
default (safe); export that variable yourself to override.

Snapshot reading is **pluggable**. The default reader (`gasmaker/readers/rur.py`)
uses the `rur` package and is imported **lazily** — the GasMaker core does not
depend on `rur`, so it installs and imports without it. By default it uses an
installed `rur` (or `$RUR_PATH`); pass `--rur-path` for a checkout, or implement
the small `gasmaker.readers.base.CellReader` interface (incl. `read_boxes` and
`hydro_fields`) to read another simulation/format.

> Status: validated against the RUR reference at machine precision on a
> stratified NH2 sample (gas/particle masses, metallicity, chemistry — see
> `docs/WHATS_NEW.md`). `r200/r500` use threshold-crossing interpolation, which
> differs by design from RUR's nearest-shell selection.

## Files

- `HaloMaker.py`: command-line entry point
- `compute_halo_props.py`: HaloMaker workflow and halo properties
- `input_output.py`: RAMSES readers and HDF5 catalog writer
- `halo_defs.py`: shared runtime state and utilities
- `num_rec.py`: numerical helpers
- `compute_neiKDtree_mod.py`: Python-to-Fortran bridge
- `compute_adaptahop.f90`: optimized full-box AdaptaHOP extension
- `compute_adaptahop.pyf`: explicit f2py interface for portable builds
- `examples/hdf_output_example.py`: simple HDF5 catalog reader example
- `ssp_photometry.py`: compact SSP-table interpolation
- `halomaker_data/ssp_tables`: generated SSP runtime tables; ignored by Git
- `clean_runtime.sh`: dry-run / cleanup helper for interrupted runs
- `GasMaker.py`: gas post-processor command-line entry point
- `gasmaker/`: GasMaker package (`pipeline`, `catalog`, `geometry`, `overlap`)
- `gasmaker/readers/`: pluggable snapshot readers (`base` interface + `rur` adapter)

## Preparing SSP tables

The compact BC03, CB07, and FSPS tables are not redistributed with HaloMaker.
Provide the source model data before the first build. Local development copies,
when present, use the standard ignored locations `assets/ssp_originals/bc03`
and `assets/ssp_originals/cb07`. Explicit paths override those defaults:

- `BC03_PATH`: BC03 Chabrier/Padova 1994 source tarball or extracted directory
- `CB07_PATH`: CB07 source-table directory used by RUR
- `FSPS_PATH`: FSPS source/data installation (`SPS_HOME` is also accepted)

BC03 can be obtained from the Bruzual & Charlot 2003 original release page:
`https://www.bruzual.org/bc03/Original_version_2003/`

Install the optional FSPS generator dependency when generating FSPS:

```bash
uv sync --extra ssp-generation
BC03_PATH=/path/to/bc03 \
CB07_PATH=/path/to/cb07 \
FSPS_PATH=/path/to/fsps \
PYTHON=.venv/bin/python bash build.sh
```

`build.sh` calls `tools/prepare_ssp_tables.sh` to generate missing tables
under `halomaker_data/ssp_tables/`, then compiles the Fortran extensions. It
reuses existing table files on later builds. The generated files are ignored by
Git. To prepare or regenerate only the SSP tables:

```bash
PYTHON=.venv/bin/python bash tools/prepare_ssp_tables.sh
HALOMAKER_TABLES_ONLY=1 PYTHON=.venv/bin/python bash build.sh
```

To regenerate individual tables explicitly:

```bash
uv run python tools/generate_bc03_table.py \
    --bc03-path /path/to/bc03 --force
uv run python tools/generate_cb07_table.py \
    --cb07-path /path/to/cb07 --force
uv run python tools/generate_fsps_table.py \
    --fsps-path /path/to/fsps --force
```

The FSPS source/data installation remains separate because `python-fsps`
alone does not provide the full model data. The BC03 and CB07 source tables
also remain separate because they have their own upstream distribution terms.

## Release Checklist

Before publishing, choose a license and add a `LICENSE` file. Also add a small
redistributable RAMSES fixture and an automated smoke test if redistribution
rights permit it.
