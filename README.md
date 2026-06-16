# HaloMaker Python

Python + Fortran implementation of HaloMaker / AdaptaHOP for RAMSES
snapshots. The Fortran extensions provide the memory-sensitive neighbor and
structure-tree routines. Both full-box and zoom-in workflows are included.

This release is staged from the dev3 memory-optimized implementation. The
public filenames remain canonical (`compute_adaptahop.f90` and
`compute_adaptahop_zoomin.f90`) so existing scripts do not need development
branch names.

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
- `compute_adaptahop_zoomin*.so`

Build and run from the same activated Python environment. If you do not want
to activate the environment, pass the interpreter explicitly:

```bash
PYTHON=.venv/bin/python bash build.sh
PYTHON=.venv/bin/python bash run.sh
```

## Configure

Copy the example inputs into the repository root:

```bash
cp examples/input_HaloMaker.dat.example input_HaloMaker.dat
cp examples/inputfiles_HaloMaker.dat.example inputfiles_HaloMaker.dat
```

Edit `inputfiles_HaloMaker.dat` so that each active line points to an existing
RAMSES snapshot. Set `zoomin = .true.` in `input_HaloMaker.dat` for zoom-in
processing and `zoomin = .false.` for periodic full-box processing.

For RAMSES snapshots, `lbox`, `omega_f`, and `lambda_f` are optional. The code
reads the authoritative box size and snapshot cosmology from the RAMSES AMR
header during `read_data()`. If those values are omitted, box-size-dependent
and cosmology-dependent quantities are finalized only after the snapshot header
is read.

Recommended print levels:

- `verbose = .false.`, `megaverbose = .false.`: compact production log.
- `verbose = .true.`, `megaverbose = .false.`: normal diagnostic log with
  major counts and timings.
- `megaverbose = .true.`: development log with detailed Fortran and memory
  tracking output.

## Run

```bash
bash run.sh
```

The main HDF5 catalog output is written as `tree_bricks*.h5`.

Each run computes intrinsic rest-frame stellar photometry with all bundled SSP
models and writes the aligned results under:

- `/photometry/CB07/data`
- `/photometry/BC03/data`
- `/photometry/FSPS/data`

The datasets contain SDSS `ugriz`, Johnson `UBV`, `Kmag`, and r-band
luminosity-weighted stellar age, metallicity, `r50`, and `r90`. FSPS uses the
2MASS Ks response for its `Kmag`; this distinction is recorded in the group
metadata. See `SSP_MODELS.md` for model definitions and data provenance.

If a run is interrupted with Ctrl-C, killed by a scheduler, or leaves Python
`forkserver` / `resource_tracker` processes behind, inspect and clean runtime
leftovers with:

```bash
bash clean_runtime.sh
bash clean_runtime.sh --force
```

The first command is a dry run. The `--force` command terminates matching
HaloMaker runtime processes in this repository and removes matching shared
memory files owned by the current user.

## Files

- `HaloMaker.py`: command-line entry point
- `compute_halo_props.py`: HaloMaker workflow and halo properties
- `input_output.py`: RAMSES readers and HDF5 catalog writer
- `halo_defs.py`: shared runtime state and utilities
- `num_rec.py`: numerical helpers
- `compute_neiKDtree_mod.py`: Python-to-Fortran bridge
- `compute_adaptahop.f90`: optimized full-box AdaptaHOP extension
- `compute_adaptahop_zoomin.f90`: optimized zoom-in AdaptaHOP extension
- `compute_adaptahop*.pyf`: explicit f2py interfaces for portable builds
- `hdf_output_example.py`: simple HDF5 catalog reader example
- `ssp_photometry.py`: compact SSP-table interpolation
- `halomaker_data/ssp_tables`: generated SSP runtime tables; ignored by Git
- `clean_runtime.sh`: dry-run / cleanup helper for interrupted runs

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
