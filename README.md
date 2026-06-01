# HaloMaker Python

Python + Fortran implementation of HaloMaker / AdaptaHOP for RAMSES
snapshots. The Fortran extensions provide the memory-sensitive neighbor and
structure-tree routines. Both full-box and zoom-in workflows are included.

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

## Configure

Copy the example inputs into the repository root:

```bash
cp examples/input_HaloMaker.dat.example input_HaloMaker.dat
cp examples/inputfiles_HaloMaker.dat.example inputfiles_HaloMaker.dat
```

Edit `inputfiles_HaloMaker.dat` so that each active line points to an existing
RAMSES snapshot. Set `zoomin = .true.` in `input_HaloMaker.dat` for zoom-in
processing and `zoomin = .false.` for periodic full-box processing.

## Run

```bash
bash run.sh
```

The main HDF5 catalog output is written as `tree_bricks*.h5`.

## Files

- `HaloMaker.py`: command-line entry point
- `compute_halo_props.py`: HaloMaker workflow and halo properties
- `input_output.py`: RAMSES readers and HDF5 catalog writer
- `halo_defs.py`: shared runtime state and utilities
- `num_rec.py`: numerical helpers
- `compute_neiKDtree_mod.py`: Python-to-Fortran bridge
- `compute_adaptahop.f90`: optimized full-box AdaptaHOP extension
- `compute_adaptahop_zoomin.f90`: optimized zoom-in AdaptaHOP extension
- `hdf_output_example.py`: simple HDF5 catalog reader example

## Release Checklist

Before publishing, choose a license and add a `LICENSE` file. Also add a small
redistributable RAMSES fixture and an automated smoke test if redistribution
rights permit it.
