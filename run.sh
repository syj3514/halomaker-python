#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
LOG="${1:-HaloMaker.log}"

"$PYTHON" HaloMaker.py 2>&1 | tee "$LOG"
