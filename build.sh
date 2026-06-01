#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
F90="${F90:-gfortran}"

for source in compute_adaptahop.f90 compute_adaptahop_zoomin.f90; do
    module="${source%.f90}"
    echo "Compiling ${module}"
    "$PYTHON" -m numpy.f2py \
        -lgomp \
        --f90exec="$F90" \
        --f77exec="$F90" \
        --f90flags="-fopenmp -O3 -x f95-cpp-input" \
        -c "${module}.pyf" "$source"
done

rm -f ./*.o ./*.mod
