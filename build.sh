#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
F90="${F90:-gfortran}"
BRANCH="${BRANCH:-main}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SSP_TABLE_DIR="$ROOT/halomaker_data/ssp_tables"
BC03_TABLE="$SSP_TABLE_DIR/bc03.npz"
CB07_TABLE="$SSP_TABLE_DIR/cb07.npz"
FSPS_TABLE="$SSP_TABLE_DIR/fsps.npz"
BC03_SOURCE="${BC03_PATH:-}"
CB07_SOURCE="${CB07_PATH:-}"
FSPS_SOURCE="${FSPS_PATH:-${SPS_HOME:-}}"

if [[ -f "$BC03_TABLE" ]]; then
    echo "Using existing BC03 table: $BC03_TABLE"
else
    if [[ -z "$BC03_SOURCE" ]]; then
        cat >&2 <<EOF
Missing BC03 table: $BC03_TABLE
Set BC03_PATH to a BC03 source tarball or extracted directory, then rerun build.sh.
The BC03 source is not redistributed with HaloMaker.
EOF
        exit 1
    fi
    if [[ ! -e "$BC03_SOURCE" ]]; then
        echo "BC03_PATH does not exist: $BC03_SOURCE" >&2
        exit 1
    fi

    echo "Generating BC03 table from: $BC03_SOURCE"
    "$PYTHON" "$ROOT/tools/generate_bc03_table.py" \
        --bc03-path "$BC03_SOURCE"
fi

if [[ -f "$CB07_TABLE" ]]; then
    echo "Using existing CB07 table: $CB07_TABLE"
else
    if [[ -z "$CB07_SOURCE" ]]; then
        cat >&2 <<EOF
Missing CB07 table: $CB07_TABLE
Set CB07_PATH to a CB07 source-table directory, then rerun build.sh.
The CB07 source tables are not redistributed with HaloMaker.
EOF
        exit 1
    fi
    if [[ ! -d "$CB07_SOURCE" ]]; then
        echo "CB07_PATH is not a directory: $CB07_SOURCE" >&2
        exit 1
    fi

    echo "Generating CB07 table from: $CB07_SOURCE"
    "$PYTHON" "$ROOT/tools/generate_cb07_table.py" \
        --cb07-path "$CB07_SOURCE"
fi

if [[ -f "$FSPS_TABLE" ]]; then
    echo "Using existing FSPS table: $FSPS_TABLE"
else
    if [[ -z "$FSPS_SOURCE" ]]; then
        cat >&2 <<EOF
Missing FSPS table: $FSPS_TABLE
Set FSPS_PATH to an FSPS source/data installation, then rerun build.sh.
Install the table generator first if needed:
  uv sync --extra ssp-generation
EOF
        exit 1
    fi
    if [[ ! -d "$FSPS_SOURCE" ]]; then
        echo "FSPS_PATH is not a directory: $FSPS_SOURCE" >&2
        exit 1
    fi
    if ! SPS_HOME="$FSPS_SOURCE" "$PYTHON" -c "import fsps" >/dev/null 2>&1; then
        cat >&2 <<EOF
python-fsps is required to generate $FSPS_TABLE.
Install it with:
  uv sync --extra ssp-generation
EOF
        exit 1
    fi

    echo "Generating FSPS table from: $FSPS_SOURCE"
    "$PYTHON" "$ROOT/tools/generate_fsps_table.py" \
        --fsps-path "$FSPS_SOURCE"
fi

if [[ "${HALOMAKER_TABLES_ONLY:-0}" == "1" ]]; then
    echo "HALOMAKER_TABLES_ONLY=1: skipping Fortran extension build"
    exit 0
fi

if [[ "$BRANCH" == "main" ]]; then
    BUILD_DIR="$ROOT"
else
    BUILD_DIR="$ROOT/patch/$BRANCH"
    if [[ ! -d "$BUILD_DIR" ]]; then
        echo "Unknown patch branch: $BRANCH" >&2
        exit 1
    fi
fi

cd "$BUILD_DIR"

for source in compute_adaptahop.f90 compute_adaptahop_zoomin.f90; do
    module="${source%.f90}"
    source_path="$BUILD_DIR/$source"
    interface_path="$BUILD_DIR/${module}.pyf"
    [[ -f "$source_path" ]] || source_path="$ROOT/$source"
    [[ -f "$interface_path" ]] || interface_path="$ROOT/${module}.pyf"

    echo "Compiling ${module} [branch=${BRANCH}]"
    "$PYTHON" -m numpy.f2py \
        -lgomp \
        --f90exec="$F90" \
        --f77exec="$F90" \
        --f90flags="-fopenmp -O3 -x f95-cpp-input" \
        -c "$interface_path" "$source_path"
done

rm -f ./*.o ./*.mod
