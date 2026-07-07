#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
LOG="${1:-GasMaker.log}"
ROOT="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
ENTRYPOINT="$ROOT/GasMaker.py"

echo "Running GasMaker"
PYTHONPATH="$ROOT/src:$ROOT${PYTHONPATH:+:$PYTHONPATH}" \
    "$PYTHON" "$ENTRYPOINT" 2>&1 | tee "$LOG"
