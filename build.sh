#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
F90="${F90:-gfortran}"
BRANCH="${BRANCH:-main}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FSPS_TABLE="$ROOT/halomaker_data/ssp_tables/fsps.npz"
FSPS_SOURCE="${FSPS_PATH:-${SPS_HOME:-}}"

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
