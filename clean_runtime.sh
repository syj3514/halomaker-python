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

# Shared-memory sweep. Files created by this pipeline are named
#   HaloMaker_u<uid>_t<time>_p<pid>_s<starttime>_<array>
# plus a co-located  ..._MANIFEST  breadcrumb. The (pid, starttime) pair is unique per
# host across time, so we can tell a DEAD orphan (owning process gone / PID recycled)
# from a LIVE run (process still running with the same start-time). Dead orphans are
# auto-removed even in dry-run mode — they are provably safe. Live runs are NEVER
# touched automatically (only --force, which first kills the process). Legacy files
# with no _p<pid>_s<starttime>_ tag can't be liveness-judged, so they keep the old
# semantics: listed in dry-run, removed only with --force.
SHM_CLASSIFY="$(
python - <<'PY'
import os, re, glob
uid = os.getuid()
pat = re.compile(r'_p(\d+)_s(\d+)_')

def alive(pid, starttime):
    try:
        with open(f"/proc/{pid}/stat") as fh:
            data = fh.read()
        return data[data.rindex(")") + 2:].split()[19] == starttime
    except (OSError, IndexError, ValueError):
        return False

for path in sorted(glob.glob("/dev/shm/HaloMaker_*") + glob.glob("/dev/shm/psm_*")):
    try:
        if os.stat(path).st_uid != uid:
            continue
    except OSError:
        continue
    m = pat.search(os.path.basename(path))
    if not m:
        print("UNKNOWN\t" + path)
    else:
        print(("LIVE" if alive(m.group(1), m.group(2)) else "DEAD") + "\t" + path)
PY
)"

SHM_DEAD=(); SHM_LIVE=(); SHM_UNKNOWN=()
while IFS=$'\t' read -r cls path; do
    [[ -z "$path" ]] && continue
    case "$cls" in
        DEAD)    SHM_DEAD+=("$path") ;;
        LIVE)    SHM_LIVE+=("$path") ;;
        UNKNOWN) SHM_UNKNOWN+=("$path") ;;
    esac
done <<< "$SHM_CLASSIFY"

if [[ "${#SHM_LIVE[@]}" -gt 0 ]]; then
    echo "[clean] live-run shared-memory (owning process alive — NOT touched):"
    printf '  %s\n' "${SHM_LIVE[@]}"
    if [[ "$FORCE" == 1 ]]; then
        rm -f -- "${SHM_LIVE[@]}"
        echo "[clean] --force: removed live-run shared-memory (process was killed above)"
    fi
fi

if [[ "${#SHM_DEAD[@]}" -gt 0 ]]; then
    echo "[clean] dead-orphan shared-memory (owning process gone — auto-removing):"
    printf '  %s\n' "${SHM_DEAD[@]}"
    rm -f -- "${SHM_DEAD[@]}"
    echo "[clean] removed ${#SHM_DEAD[@]} dead-orphan shared-memory file(s)"
fi

if [[ "${#SHM_UNKNOWN[@]}" -gt 0 ]]; then
    echo "[clean] untagged/legacy shared-memory (liveness unknown):"
    printf '  %s\n' "${SHM_UNKNOWN[@]}"
    if [[ "$FORCE" == 1 ]]; then
        rm -f -- "${SHM_UNKNOWN[@]}"
        echo "[clean] --force: removed untagged/legacy shared-memory files"
    else
        echo "[clean] dry-run only; use --force to remove untagged/legacy files"
    fi
fi

if [[ "${#SHM_DEAD[@]}" -eq 0 && "${#SHM_LIVE[@]}" -eq 0 && "${#SHM_UNKNOWN[@]}" -eq 0 ]]; then
    echo "[clean] no matching shared-memory files"
fi
