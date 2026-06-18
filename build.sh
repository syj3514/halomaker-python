#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
F90="${F90:-gfortran}"
BRANCH="${BRANCH:-main}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "$ROOT/tools/prepare_ssp_tables.sh"

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

for source in compute_adaptahop.f90; do
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
