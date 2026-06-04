#!/usr/bin/env bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
LOG="${1:-HaloMaker.log}"
BRANCH="${BRANCH:-main}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ "$BRANCH" == "main" ]]; then
    PATCH_PATH=""
else
    PATCH_PATH="$ROOT/patch/$BRANCH"
    if [[ ! -d "$PATCH_PATH" ]]; then
        echo "Unknown patch branch: $BRANCH" >&2
        exit 1
    fi
fi

ENTRYPOINT="$ROOT/HaloMaker.py"
if [[ -n "$PATCH_PATH" && -f "$PATCH_PATH/HaloMaker.py" ]]; then
    ENTRYPOINT="$PATCH_PATH/HaloMaker.py"
fi

cd "$ROOT"
echo "Running HaloMaker [branch=${BRANCH}]"
PYTHONPATH="${PATCH_PATH:+$PATCH_PATH:}$ROOT${PYTHONPATH:+:$PYTHONPATH}" \
    "$PYTHON" "$ENTRYPOINT" 2>&1 | tee "$LOG"
