#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORCE=0
ORPHAN=0
case "${1:-}" in
    --force)  FORCE=1 ;;
    --orphan) ORPHAN=1 ;;
    "")       ;;  # bare: true dry-run, deletes nothing
    -h|--help)
        cat <<'EOF'
Usage:
  bash clean_runtime.sh            # dry-run: list only, delete NOTHING
  bash clean_runtime.sh --orphan   # safe mode: reclaim provably-dead orphan shm only
  bash clean_runtime.sh --force    # aggressive: kill this repo's runtime procs, then reclaim

Shared-memory files are named  HaloMaker_u<uid>_t<time>_p<pid>_s<starttime>_<array>
(+ a co-located ..._MANIFEST breadcrumb). The (pid, starttime) pair is unique per host
across time, so a segment is classified as:
  DEAD    owning process gone / PID recycled (starttime mismatch) / zombie  -> orphan
  LIVE    owning process still running with the same start-time            -> a real run
  UNKNOWN legacy / generic (no _p<pid>_s<starttime>_ tag, e.g. psm_*)       -> can't judge

Deletion policy:
  dry-run  : nothing removed (lists what --orphan / --force would do)
  --orphan : DEAD orphans (+their manifest) removed; LIVE and UNKNOWN untouched
  --force  : kills this repository's runtime processes FIRST, re-classifies, then removes
             the now-DEAD segments plus UNKNOWN/legacy files. A LIVE run owned by another
             working tree / another agent is never killed here and is therefore never
             removed — the "never touch a live run" invariant holds under --force too.

Note: DEAD is judged on the owning (parent) PID. A hard-killed run may briefly leave a
fork child mapping the segment; unlinking only drops the name (existing mappings stay
valid, no crash), but "owner PID gone" is not strictly identical to "logical run fully
ended".
EOF
        exit 0 ;;
    *)
        echo "clean_runtime.sh: unknown option '$1' (see --help)" >&2
        exit 2 ;;
esac

# --force also reclaims dead orphans (after killing this repo's procs turns them dead).
RECLAIM_DEAD=0
[[ "$FORCE" == 1 || "$ORPHAN" == 1 ]] && RECLAIM_DEAD=1

MODE=dry-run
[[ "$ORPHAN" == 1 ]] && MODE=orphan
[[ "$FORCE" == 1 ]] && MODE=force
echo "[clean] repository: $ROOT"
echo "[clean] mode: $MODE"

# --- runtime processes (only --force kills, and only this repository's own runs) ---
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
    echo "[clean] matching runtime process IDs (this repository): ${PIDS[*]}"
    if [[ "$FORCE" == 1 ]]; then
        kill "${PIDS[@]}" 2>/dev/null || true
        sleep 1
        kill -9 "${PIDS[@]}" 2>/dev/null || true
        echo "[clean] terminated matching processes"
    else
        echo "[clean] not terminated; use --force to terminate them"
    fi
fi

# --- shared-memory sweep. Classified AFTER the (optional) kill above, so a run this
# repo just terminated is now judged DEAD and reclaimed, while a LIVE run owned by
# some other working tree / agent (which we never killed) stays LIVE and is preserved. ---
SHM_CLASSIFY="$(
python - <<'PY'
import os, re, glob
uid = os.getuid()
pat = re.compile(r'_p(\d+)_s(\d+)_')

def alive(pid, starttime):
    """Live iff /proc/<pid>/stat is readable, not a zombie, and its start-time
    (field 22) matches — so a recycled PID or a just-killed zombie reads as dead."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            data = fh.read()
        fields = data[data.rindex(")") + 2:].split()
        state = fields[0]        # field 3
        if state == "Z":         # zombie: owning process has exited
            return False
        return fields[19] == starttime  # field 22: kernel start-time
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

# LIVE: never removed by any mode. Reported so the operator knows a run is active.
if [[ "${#SHM_LIVE[@]}" -gt 0 ]]; then
    echo "[clean] live-run shared-memory (owning process alive — preserved, never removed):"
    printf '  %s\n' "${SHM_LIVE[@]}"
fi

# DEAD orphans: removed by --orphan and --force; only listed in dry-run.
if [[ "${#SHM_DEAD[@]}" -gt 0 ]]; then
    if [[ "$RECLAIM_DEAD" == 1 ]]; then
        echo "[clean] dead-orphan shared-memory (owning process gone — reclaiming):"
        printf '  %s\n' "${SHM_DEAD[@]}"
        rm -f -- "${SHM_DEAD[@]}"
        echo "[clean] reclaimed ${#SHM_DEAD[@]} dead-orphan shared-memory file(s)"
    else
        echo "[clean] dead-orphan shared-memory (would reclaim; use --orphan or --force):"
        printf '  %s\n' "${SHM_DEAD[@]}"
    fi
fi

# UNKNOWN/legacy: can't judge liveness, so only --force removes them.
if [[ "${#SHM_UNKNOWN[@]}" -gt 0 ]]; then
    echo "[clean] untagged/legacy shared-memory (liveness unknown):"
    printf '  %s\n' "${SHM_UNKNOWN[@]}"
    if [[ "$FORCE" == 1 ]]; then
        rm -f -- "${SHM_UNKNOWN[@]}"
        echo "[clean] --force: removed untagged/legacy shared-memory files"
    else
        echo "[clean] not removed; use --force to remove untagged/legacy files"
    fi
fi

if [[ "${#SHM_DEAD[@]}" -eq 0 && "${#SHM_LIVE[@]}" -eq 0 && "${#SHM_UNKNOWN[@]}" -eq 0 ]]; then
    echo "[clean] no matching shared-memory files"
fi
