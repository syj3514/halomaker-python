#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
F90="${F90:-gfortran}"
BRANCH="${BRANCH:-main}"
SCRIPT="$(readlink -f "${BASH_SOURCE[0]}")"
ROOT="$(cd "$(dirname "$SCRIPT")" && pwd)"

if ! python_command="$(command -v "$PYTHON")"; then
    echo "Python executable not found: $PYTHON" >&2
    exit 1
fi
if [[ "$python_command" == /* ]]; then
    PYTHON="$python_command"
else
    PYTHON="$PWD/$python_command"
fi

if [[ "${HALOMAKER_SKIP_SSP:-0}" == "1" && \
      "${HALOMAKER_TABLES_ONLY:-0}" == "1" ]]; then
    echo "HALOMAKER_SKIP_SSP=1 and HALOMAKER_TABLES_ONLY=1 are mutually exclusive." >&2
    exit 1
fi

if [[ "${HALOMAKER_SKIP_SSP:-0}" == "1" ]]; then
    echo "HALOMAKER_SKIP_SSP=1: skipping SSP table preparation"
else
    PYTHON="$PYTHON" bash "$ROOT/tools/prepare_ssp_tables.sh"
fi

if [[ "${HALOMAKER_TABLES_ONLY:-0}" == "1" ]]; then
    echo "HALOMAKER_TABLES_ONLY=1: skipping Fortran extension build"
    exit 0
fi

if ! fortran_command="$(command -v "$F90")"; then
    echo "Fortran compiler not found: $F90" >&2
    exit 1
fi
if [[ "$fortran_command" == /* ]]; then
    F90="$fortran_command"
else
    F90="$PWD/$fortran_command"
fi
fortran_version="$("$F90" --version 2>&1)"
if ! grep -qi 'GNU Fortran' <<<"$fortran_version"; then
    echo "build.sh supports gfortran only (got F90=$F90)." >&2
    echo "Use tools/build_intel.sh for ifort or ifx." >&2
    exit 1
fi

if [[ "$BRANCH" == "main" ]]; then
    BUILD_DIR="$ROOT/src"
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
    [[ -f "$source_path" ]] || source_path="$ROOT/fortran/$source"
    [[ -f "$interface_path" ]] || interface_path="$ROOT/fortran/${module}.pyf"

    echo "Compiling ${module} [branch=${BRANCH}]"
    openmp_fflags="-fopenmp -O3 -x f95-cpp-input"
    FC="$F90" \
    F77="$F90" \
    FFLAGS="${FFLAGS:+$FFLAGS }$openmp_fflags" \
    LDFLAGS="${LDFLAGS:+$LDFLAGS }-fopenmp" \
    "$PYTHON" -m numpy.f2py \
        -lgomp \
        --f90exec="$F90" \
        --f77exec="$F90" \
        --f90flags="$openmp_fflags" \
        -c "$interface_path" "$source_path"

    extension_suffix="$("$PYTHON" -c \
        'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"
    artifact="$BUILD_DIR/${module}${extension_suffix}"
    if [[ ! -f "$artifact" ]]; then
        echo "Build completed without the expected extension: $artifact" >&2
        exit 1
    fi
    "$ROOT/tools/check_openmp_extension.sh" "$artifact"
done

rm -f ./*.o ./*.mod
