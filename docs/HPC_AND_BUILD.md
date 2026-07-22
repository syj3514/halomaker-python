# Building and Running HaloMaker on HPC Systems

This guide covers a first build, OpenMP verification, SSP-table choices, and
single-node PBS or Slurm execution. Commands are run from the repository root
unless a section says otherwise.

## Quick Start

1. Create and activate a Python 3.10+ environment; install the project
   dependencies.
2. Choose how to provide the three SSP tables, or disable photometry.
3. Build with GNU Fortran (default) or the separate Intel build script.
4. Verify the extension's OpenMP runtime before submitting a long job.
5. Put each run's input files in a dedicated run directory and make `nbPes`
   match the CPU allocation.

GNU build with existing SSP tables:

```bash
PYTHON=.venv/bin/python bash build.sh
```

Halo-only build without SSP tables:

```bash
HALOMAKER_SKIP_SSP=1 PYTHON=.venv/bin/python bash build.sh
```

The second form must be paired with `photometry = .false.` in the run's
`input_HaloMaker.dat`. Otherwise the build succeeds, but the run correctly
stops when it cannot load the tables.

## Choose a Compiler

### GNU Fortran (default)

`build.sh` supports **gfortran only**. It passes `-fopenmp` through both the
classic f2py and Python 3.12+ Meson paths, links `libgomp`, and rejects an
extension that lacks the GNU OpenMP symbol.

```bash
PYTHON=/path/to/python F90=/path/to/gfortran bash build.sh
```

Use the same Python environment to build and run. A login-node `python` can be
different from the interpreter in a batch job, so set `PYTHON` explicitly in
both places.

### Intel ifort or ifx

Do not set `F90=ifort` on `build.sh`. Meson/f2py combinations can misidentify
Intel Fortran, and `build.sh` intentionally remains gfortran-only. Use the
Meson-bypass script instead:

```bash
PYTHON=/path/to/python FC=ifx CC=icx bash tools/build_intel.sh
```

For an older Intel Classic installation:

```bash
PYTHON=/path/to/python FC=ifort CC=icx bash tools/build_intel.sh
```

The script generates the f2py wrappers, compiles the Fortran implementation
directly with `-qopenmp`, compiles the C wrapper with `icx`, and links
`libiomp5`. `ifort` is deprecated upstream; prefer `ifx` for a new toolchain.
If the Intel runtime is outside the default loader path, set `INTEL_LIBDIR` to
the directory containing `libiomp5`.

Useful optional variables are:

- `BRANCH`: source branch under `patch/` (`main` by default).
- `OUTPUT_DIR`: extension installation directory (`src/` by default).
- `BUILD_DIR`: retained work directory for debugging (otherwise temporary).
- `PYTHON`, `FC`, `CC`, `INTEL_LIBDIR`: toolchain paths.

## SSP Tables and Photometry

With photometry enabled (the default), HaloMaker needs all three compact files
at build time **and** runtime:

```text
halomaker_data/ssp_tables/bc03.npz
halomaker_data/ssp_tables/cb07.npz
halomaker_data/ssp_tables/fsps.npz
```

`build.sh` runs SSP preparation before compiling. If a required table and its
source are both absent, the build stops before the Fortran extension is built.
Choose one of the following workflows.

### Generate from original model data

Provide the model sources and let `build.sh` generate missing compact tables:

```bash
uv sync --extra ssp-generation
BC03_PATH=/path/to/bc03 \
CB07_PATH=/path/to/cb07 \
FSPS_PATH=/path/to/fsps \
PYTHON=.venv/bin/python bash build.sh
```

`SPS_HOME` is accepted in place of `FSPS_PATH`. Source information:

- BC03: <https://www.bruzual.org/bc03/Original_version_2003/>
- CB07: <https://www.bruzual.org/cb07/>
- FSPS: <https://github.com/cconroy20/fsps>

The model data have their own distribution terms and are not included in this
repository. See [SSP_MODELS.md](SSP_MODELS.md) for the adopted models.

### Reuse trusted prebuilt tables

If a collaborator or an existing installation provides compact tables, place
all three `.npz` files in `halomaker_data/ssp_tables/`. Verify them against the
provider's checksum manifest before building, for example:

```bash
mkdir -p halomaker_data/ssp_tables
cp /trusted/location/{bc03,cb07,fsps}.npz halomaker_data/ssp_tables/
sha256sum halomaker_data/ssp_tables/{bc03,cb07,fsps}.npz
PYTHON=.venv/bin/python bash build.sh
```

The repository does not ship a checksum manifest because it does not ship
these tables. Compare all three printed digests with the manifest supplied by
the same trusted provider as the files.

### Disable photometry

Users who only need halo properties can build and run without SSP data:

```bash
HALOMAKER_SKIP_SSP=1 PYTHON=.venv/bin/python bash build.sh
```

Add this to `input_HaloMaker.dat`:

```text
photometry = .false.
```

The default is `.true.`. With photometry off, the catalog records
`/input.photometry = false` and omits the entire `/photometry` group. Analysis
code must inspect `/input.photometry` before requiring photometry datasets.
The `/catalog` and `/member` science datasets retain their normal schema.

`HALOMAKER_SKIP_SSP=1` controls only the build-time preparation; the config
flag controls runtime behavior. `HALOMAKER_SKIP_SSP=1` and
`HALOMAKER_TABLES_ONLY=1` are mutually exclusive.

## Verify OpenMP

A successful import is not sufficient. If a build backend silently drops the
OpenMP compile flag, the code still runs and can still print
`[OMP] ... ncore=N`; that line echoes `nbPes`, not the observed thread count.

`build.sh` and `tools/build_intel.sh` run this check automatically. To repeat
it for the extension imported by a specific Python:

```bash
PYTHON=/path/to/python
EXT=$(PYTHONPATH=src "$PYTHON" -c \
  'import compute_adaptahop; print(compute_adaptahop.__file__)')
bash tools/check_openmp_extension.sh "$EXT"
```

The check requires either:

- `GOMP_parallel` for a GNU `libgomp` build; or
- `__kmpc_fork_call` and a resolved `libiomp5` for an Intel build.

For manual inspection:

```bash
nm -D "$EXT" | grep -E 'GOMP_parallel|__kmpc_fork_call'
ldd "$EXT" | grep -E 'libgomp|libiomp5'
```

Also verify a real compute-node run. During an OpenMP-heavy phase, the main
process should show more than one lightweight process/thread and CPU usage can
exceed 100%:

```bash
ps -o pid,nlwp,pcpu,etime,cmd -p <halo-maker-pid>
```

Check on a compute node as well as the login node: compiler runtimes, loader
paths, CPU affinity, and the selected Python can differ. A missing symbol is a
build failure; a present symbol but one-thread runtime points instead to the
allocation, affinity, environment, or workload phase.

## Threading Model

The third `inputfiles_HaloMaker.dat` field, `nbPes`, is HaloMaker's shared
concurrency knob:

- Python stages create up to `nbPes` multiprocessing workers.
- SciPy/KD-tree calls receive `workers=nbPes` where implemented.
- Fortran OpenMP stages call `omp_set_num_threads(nbPes)` or use
  `NUM_THREADS(nbPes)`.

These stages are mostly stage-wise rather than a nested `nbPes x nbPes` model.
`nbPes` is unrelated to RAMSES `ncpu`: `ncpu` is the number of RAMSES domain
files and is read from the snapshot header.

Set `nbPes` no larger than the physical/logical CPUs allocated to the batch
job. Prevent NumPy-linked libraries from adding another thread pool inside
Python workers:

```bash
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

It is reasonable to set `OMP_NUM_THREADS` to the same value as `nbPes` in the
job environment, but the Fortran routines explicitly select `nbPes`; the input
file remains the application-level source of truth.

## Inputfiles Reference

Each non-comment line of `inputfiles_HaloMaker.dat` has exactly four or five
whitespace-separated fields:

```text
'<snapshot_dir>' <format> <nbPes> <numstep> [suffix]
```

| Field | Meaning |
|---|---|
| `snapshot_dir` | RAMSES snapshot directory. Quotes are decorative; paths containing spaces are not supported. |
| `format` | Snapshot format understood by HaloMaker, such as `Ra3` or `Ra4`. |
| `nbPes` | Python-worker and Fortran-OpenMP concurrency for this snapshot. |
| `numstep` | Integer output number used in `tree_bricksNNNNN.h5`. |
| `suffix` | Optional output tag. `nh2` produces `tree_bricksNNNNN_nh2.h5`. |

Blank lines and lines beginning with `#` or `!` are ignored. Multiple snapshot
lines run sequentially. If output numbers collide and no suffix is supplied,
HaloMaker derives a disambiguating tag so an earlier catalog is not silently
overwritten.

## PBS Example

Scheduler syntax varies by site; adjust resource and queue directives locally.
The run directory must contain `input_HaloMaker.dat` and
`inputfiles_HaloMaker.dat`.

```bash
#!/usr/bin/env bash
#PBS -N halomaker
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -l walltime=06:00:00
#PBS -l place=excl
#PBS -j oe

set -euo pipefail

REPO=/path/to/halomaker-python-release
RUN_DIR=/path/to/runs/my_run
PYTHON=/path/to/environment/bin/python

cd "$PBS_O_WORKDIR"
cd "$RUN_DIR"

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=32

EXT=$(PYTHONPATH="$REPO/src" "$PYTHON" -c \
  'import compute_adaptahop; print(compute_adaptahop.__file__)')
bash "$REPO/tools/check_openmp_extension.sh" "$EXT"
PYTHON="$PYTHON" bash "$REPO/run.sh" HaloMaker.log
```

Submit with `qsub job.pbs`. Set the inputfiles `nbPes` field to `32` for this
allocation.

## Slurm Example

```bash
#!/usr/bin/env bash
#SBATCH --job-name=halomaker
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --exclusive
#SBATCH --output=halomaker-%j.log

set -euo pipefail

REPO=/path/to/halomaker-python-release
RUN_DIR=/path/to/runs/my_run
PYTHON=/path/to/environment/bin/python

cd "$RUN_DIR"

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

EXT=$(PYTHONPATH="$REPO/src" "$PYTHON" -c \
  'import compute_adaptahop; print(compute_adaptahop.__file__)')
bash "$REPO/tools/check_openmp_extension.sh" "$EXT"
PYTHON="$PYTHON" bash "$REPO/run.sh" HaloMaker.log
```

Submit with `sbatch job.slurm`. Make the inputfiles `nbPes` value equal to
`SLURM_CPUS_PER_TASK`.

## Filesystems, Timing, and Memory

- Request one full/exclusive node for reproducible benchmarks and large-memory
  runs. Concurrent jobs on the same node can distort runtime and memory data.
- Shared NFS/Lustre metadata and I/O can dominate snapshot reading. For clean
  compute timing, stage inputs and the run directory to scheduler-provided
  node-local scratch, run there, then copy the catalog and log back. Do not
  delete the scratch copy until the output checksum and destination are
  verified.
- Build in the target software environment. At minimum, import and run the
  OpenMP self-check inside the compute job. Rebuild on the compute node when
  login and compute nodes use different compiler/runtime stacks.
- `/usr/bin/time -v` reports useful process-tree timing and a maximum RSS, but
  it is not a sum of simultaneously resident Python workers and shared-memory
  mappings. Do not treat it as the job's aggregate peak.
- Prefer scheduler accounting for aggregate memory: PBS
  `resources_used.mem` from `qstat -xf <jobid>`, or the site's Slurm accounting
  fields (commonly `MaxRSS` via `sacct`). Confirm the local scheduler's exact
  aggregation semantics before comparing systems.

## Preflight Checklist

- The batch job uses the intended Python (`PYTHON=/absolute/path`).
- The extension imports under that Python on the compute node.
- `tools/check_openmp_extension.sh` passes.
- A short run shows `NLWP > 1` or CPU usage above 100% in an OpenMP phase.
- `nbPes` matches the scheduler CPU allocation.
- BLAS/NumExpr thread counts are pinned to one.
- Either all three SSP tables exist, or both `HALOMAKER_SKIP_SSP=1` was used at
  build time and `photometry = .false.` is set at runtime.
- The run directory and output destination are writable and have enough free
  space.
- The scheduler's aggregate-memory field, not only `/usr/bin/time -v`, is
  captured for memory profiling.
