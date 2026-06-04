#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORCE=0
if [[ "${1:-}" == "--force" ]]; then
    FORCE=1
elif [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    cat <<'EOF'
Usage:
  bash clean_runtime.sh
  bash clean_runtime.sh --force

Without --force, this script only lists HaloMaker-related runtime leftovers.
With --force, it terminates matching processes and removes HaloMaker shared
memory files owned by the current user.
EOF
    exit 0
fi

echo "[clean] repository: $ROOT"
echo "[clean] mode: $([[ "$FORCE" == 1 ]] && echo force || echo dry-run)"

mapfile -t PIDS < <(
python - "$ROOT" <<'PY'
import os
import sys

root = os.path.realpath(sys.argv[1])
me = os.getpid()
matches = []

for name in os.listdir("/proc"):
    if not name.isdigit():
        continue
    pid = int(name)
    if pid == me:
        continue
    base = f"/proc/{pid}"
    try:
        with open(f"{base}/cmdline", "rb") as f:
            cmd = f.read().replace(b"\0", b" ").decode("utf-8", "replace").strip()
        cwd = os.path.realpath(os.readlink(f"{base}/cwd"))
    except (FileNotFoundError, ProcessLookupError, PermissionError, OSError):
        continue

    in_repo = cwd == root or cwd.startswith(root + os.sep)
    looks_halomaker = any(token in cmd for token in (
        "HaloMaker.py",
        "run.sh",
        "multiprocessing.resource_tracker",
        "multiprocessing.forkserver",
        "from multiprocessing",
    ))
    if in_repo and looks_halomaker:
        matches.append((pid, cwd, cmd))

for pid, cwd, cmd in sorted(matches):
    print(pid)
    print(f"  cwd={cwd}", file=sys.stderr)
    print(f"  cmd={cmd}", file=sys.stderr)
PY
)

if [[ "${#PIDS[@]}" -eq 0 ]]; then
    echo "[clean] no matching runtime processes"
else
    echo "[clean] matching runtime process IDs: ${PIDS[*]}"
    if [[ "$FORCE" == 1 ]]; then
        kill "${PIDS[@]}" 2>/dev/null || true
        sleep 1
        kill -9 "${PIDS[@]}" 2>/dev/null || true
        echo "[clean] terminated matching processes"
    else
        echo "[clean] dry-run only; use --force to terminate them"
    fi
fi

mapfile -t SHM_FILES < <(
    find /dev/shm -maxdepth 1 -user "$(id -un)" \
        \( -name 'HaloMaker_*' -o -name 'psm_*' \) -print 2>/dev/null | sort
)

if [[ "${#SHM_FILES[@]}" -eq 0 ]]; then
    echo "[clean] no matching shared-memory files"
else
    printf '[clean] matching shared-memory files:\n'
    printf '  %s\n' "${SHM_FILES[@]}"
    if [[ "$FORCE" == 1 ]]; then
        rm -f -- "${SHM_FILES[@]}"
        echo "[clean] removed matching shared-memory files"
    else
        echo "[clean] dry-run only; use --force to remove them"
    fi
fi
