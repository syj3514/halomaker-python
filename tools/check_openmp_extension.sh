#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 /path/to/extension.so" >&2
    exit 2
fi

artifact="$1"
if [[ ! -f "$artifact" ]]; then
    echo "OpenMP check failed: extension not found: $artifact" >&2
    exit 1
fi

if ! command -v nm >/dev/null 2>&1; then
    echo "OpenMP check failed: nm is required to inspect $artifact" >&2
    exit 1
fi

symbols="$(nm -D "$artifact" 2>/dev/null || true)"
if grep -q 'GOMP_parallel' <<<"$symbols"; then
    echo "OpenMP check passed (GNU libgomp): $artifact"
    exit 0
fi

if grep -q '__kmpc_fork_call' <<<"$symbols"; then
    if ! command -v ldd >/dev/null 2>&1; then
        echo "OpenMP check failed: ldd is required for an Intel extension." >&2
        exit 1
    fi
    libraries="$(ldd "$artifact" 2>&1 || true)"
    if grep -q 'libiomp5' <<<"$libraries" && \
       ! grep -q 'libiomp5.*not found' <<<"$libraries"; then
        echo "OpenMP check passed (Intel libiomp5): $artifact"
        exit 0
    fi
    echo "OpenMP check failed: __kmpc_fork_call is present but libiomp5 is not resolved." >&2
    exit 1
fi

cat >&2 <<EOF
OpenMP check failed: no parallel runtime symbol found in:
  $artifact
Expected GOMP_parallel (gfortran) or __kmpc_fork_call with libiomp5 (ifort/ifx).
The extension may have been built without OpenMP; refusing a silent serial build.
EOF
exit 1
