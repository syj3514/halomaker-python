#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${PYTHON:-python}"
SSP_TABLE_DIR="$ROOT/halomaker_data/ssp_tables"
BC03_TABLE="$SSP_TABLE_DIR/bc03.npz"
CB07_TABLE="$SSP_TABLE_DIR/cb07.npz"
FSPS_TABLE="$SSP_TABLE_DIR/fsps.npz"
BC03_SOURCE="${BC03_PATH:-$ROOT/assets/ssp_originals/bc03}"
CB07_SOURCE="${CB07_PATH:-$ROOT/assets/ssp_originals/cb07}"
FSPS_SOURCE="${FSPS_PATH:-${SPS_HOME:-}}"

missing_table_help() {
    local model="$1"
    local path="$2"
    PYTHONPATH="$ROOT/src:$ROOT${PYTHONPATH:+:$PYTHONPATH}" \
        "$PYTHON" -c \
        'from ssp_photometry import missing_table_message; import sys; print(missing_table_message(sys.argv[1], sys.argv[2]), end="")' \
        "$model" "$path" >&2
}

if [[ -f "$BC03_TABLE" ]]; then
    echo "Using existing BC03 table: $BC03_TABLE"
else
    if [[ ! -e "$BC03_SOURCE" ]]; then
        missing_table_help BC03 "$BC03_TABLE"
        exit 1
    fi

    echo "Generating BC03 table from: $BC03_SOURCE"
    "$PYTHON" "$ROOT/tools/generate_bc03_table.py" \
        --bc03-path "$BC03_SOURCE"
fi

if [[ -f "$CB07_TABLE" ]]; then
    echo "Using existing CB07 table: $CB07_TABLE"
else
    if [[ ! -d "$CB07_SOURCE" ]]; then
        missing_table_help CB07 "$CB07_TABLE"
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
        missing_table_help FSPS "$FSPS_TABLE"
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
