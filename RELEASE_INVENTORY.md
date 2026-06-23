# Release Inventory

This directory is a clean public-release staging area. It intentionally does
not include local simulation data, generated catalogs, binary extensions, or
development experiments.

This staging copy uses the dev3 memory-optimized AdaptaHOP implementation as
the canonical release source. Development branch names are intentionally not
part of the public filenames.

## 1. Existing Files Included

Runtime Python code:

- `HaloMaker.py`
- `compute_halo_props.py`
- `compute_neiKDtree_mod.py`
- `halo_defs.py`
- `input_output.py`
- `chem_species.py` (shared 9-element chemistry list for HaloMaker + GasMaker)
- `num_rec.py`
- `hdf_output_example.py`
- `ssp_photometry.py`

GasMaker (gas post-processor):

- `GasMaker.py` (entry point)
- `gasmaker/` package: `pipeline.py`, `catalog.py`, `geometry.py`, `overlap.py`
- `gasmaker/readers/`: `base.py` (reader interface), `rur.py` (optional rur
  adapter, lazily imported — the core does not depend on `rur`)

Final Fortran sources:

- `compute_adaptahop.f90`
- `compute_adaptahop.pyf`

The pipeline is now **full-box only**: the legacy zoom-in Fortran module
(`compute_adaptahop_zoomin.f90`/`.pyf`) was removed, and the fullbox interface
no longer carries the zoom-in argument. The release Fortran filename is
canonical, staged from `compute_adaptahop_dev3.f90` so existing Python imports
work without branch-specific names.

## 2. Existing Files Excluded

Development-only material:

- `AGENTS.md`
- `ADAPTAHOP_ZOOMIN_DEV_NOTES.md`
- `codex_eval/` logs, comparators, diagnostics, and intermediate runners
- intermediate `compute_adaptahop_*dev*.f90`, test, balanced, and debug files
- experimental dev4 zoom-in read-compact patch files
- `trashcan/`
- exploratory scripts such as `test.py`, `main.py`, `fortrun.sh`,
  `param_test.sh`, and `nc_test.sh`
- visualization notebooks, generated plots, and `draw/`

Machine-specific or generated material:

- `.venv/`, `__pycache__/`, `*.pyc`
- compiled `*.so`, `*.o`, and `*.mod`
- `tree_bricks*`, `struct_tree.dat`, `ncontam_halos*.dat`, and
  `resim_masses*.dat`
- local reference logs and HDF5 catalogs
- local experiment output directories
- active `input_HaloMaker.dat` and `inputfiles_HaloMaker.dat`, because they
  contain environment-specific settings and filesystem paths

Repository metadata not copied:

- the existing `f2py.sh`, because it compiles development files
- the existing `pyrun.sh`, because it is tied to the local `uv` workflow
- the existing `uv.lock`, because it captures development dependencies

## 3. New Files Added

- `README.md`: installation, build, configuration, and run instructions
- `pyproject.toml`: minimal runtime package metadata
- `requirements.txt`: minimal runtime dependency list
- `environment.yml`: optional conda environment definition
- `.gitignore`: generated-output exclusions
- `build.sh`: reproducible f2py build for the two final Fortran modules
- `run.sh`: minimal execution wrapper
- `examples/input_HaloMaker.dat.example`: sanitized parameter example
- `examples/inputfiles_HaloMaker.dat.example`: sanitized snapshot-list example

## Still Required Before Public Release

- Choose and add a `LICENSE` file. This cannot be inferred safely.
- Confirm authorship and citation text for the README.
- Add a small redistributable RAMSES fixture and smoke test if licensing and
  storage constraints allow it.
- Run the staged release against a small redistributable full-box dataset once
  such a fixture can be included.

`uv.lock` and `.python-version` are intentionally not distributed. Users may
generate them locally if they want to pin their selected interpreter and exact
dependency resolution.
