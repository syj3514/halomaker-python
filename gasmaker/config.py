from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PARAM_FILE = "input_GasMaker.dat"
INPUTFILES_FILE = "inputfiles_GasMaker.dat"


SUPPORTED_KEYS = {
    "mode",
    "radius_field",
    "padding",
    "nthread",
    "overlap_depth",
    "overlap_tolerance",
    "overlap_threshold",
    "roots",
    "root_ids",
    "max_roots",
    "read_grav",
    "rur_path",
    "overwrite",
    "progress",
    "progress_every",
}


@dataclass(frozen=True)
class GasMakerJob:
    line_no: int
    catalog: str
    repo: str
    iout: int
    output: Path


def _strip_quotes(value: str) -> str:
    if len(value) >= 2 and value[0] in ("'", '"') and value[-1] == value[0]:
        return value[1:-1]
    return value


def _strip_value_comment(value: str) -> str:
    if "!" in value:
        value = value.split("!", 1)[0]
    return value.strip()


def _parse_bool(value: str, *, key: str) -> bool:
    lowered = value.strip().lower()
    if lowered in (".true.", "true", "t", "1", "yes", "y"):
        return True
    if lowered in (".false.", "false", "f", "0", "no", "n"):
        return False
    raise ValueError(f"{PARAM_FILE}: key {key!r} expects .true. or .false., got {value!r}")


def _parse_root_ids(value: str) -> list[int]:
    root_ids = []
    for item in value.split(","):
        item = item.strip()
        if item:
            root_ids.append(int(item))
    if not root_ids:
        raise ValueError(f"{PARAM_FILE}: key 'root_ids' expects at least one root id")
    return root_ids


def _coerce_value(key: str, value: str):
    if key in {"padding", "overlap_tolerance", "overlap_threshold"}:
        return float(value)
    if key in {"nthread", "overlap_depth", "max_roots", "progress_every"}:
        return int(value)
    if key == "progress":
        from .progress import MODES
        if value not in MODES:
            raise ValueError(
                f"{PARAM_FILE}: key 'progress' expects one of {MODES}, got {value!r}"
            )
        return value
    if key in {"read_grav", "overwrite"}:
        return _parse_bool(value, key=key)
    if key == "root_ids":
        return _parse_root_ids(value)
    if key == "rur_path" and value == "":
        return None
    return _strip_quotes(value)


def read_params(path: str | Path = PARAM_FILE) -> dict:
    params = {}
    with open(path, "r") as handle:
        for lineno, line in enumerate(handle, start=1):
            stripped = line.strip()
            if (not stripped) or stripped[0] in ("#", "!"):
                continue
            if "=" not in stripped:
                raise ValueError(
                    f"{path}: expected 'key = value' on line {lineno}: {stripped!r}"
                )
            key, value = map(str.strip, stripped.split("=", 1))
            if key not in SUPPORTED_KEYS:
                print(
                    f"Warning: {path}: ignoring unknown key {key!r} on line {lineno}",
                    file=sys.stderr,
                )
                continue
            params[key] = _coerce_value(key, _strip_value_comment(value))
    return params


def read_inputfiles(path: str | Path = INPUTFILES_FILE) -> list[GasMakerJob]:
    jobs = []
    with open(path, "r") as handle:
        for lineno, line in enumerate(handle, start=1):
            stripped = line.strip()
            if (not stripped) or stripped[0] in ("#", "!"):
                continue
            fields = stripped.split()
            if len(fields) not in (3, 4):
                raise ValueError(
                    f"{path}: expected 3 or 4 fields "
                    "('<catalog.h5>' '<ramses_repo>' <iout> [output]), got "
                    f"{len(fields)} on line {lineno}: {stripped!r}"
                )
            catalog, repo, iout_text = fields[:3]
            iout = int(iout_text)
            output = fields[3] if len(fields) == 4 else f"gas_bricks{iout:05d}.h5"
            jobs.append(
                GasMakerJob(
                    line_no=lineno,
                    catalog=_strip_quotes(catalog),
                    repo=_strip_quotes(repo),
                    iout=iout,
                    output=Path(_strip_quotes(output)),
                )
            )
    if not jobs:
        raise ValueError(
            f"{path}: reached end of file without a job line "
            "(expected: '<catalog.h5>' '<ramses_repo>' <iout> [output])"
        )
    return jobs


def validate_unique_outputs(jobs: Iterable[GasMakerJob]) -> None:
    seen = {}
    for job in jobs:
        key = str(job.output)
        if key in seen:
            raise ValueError(
                f"{INPUTFILES_FILE}: output collision for {key!r} on lines "
                f"{seen[key]} and {job.line_no}"
            )
        seen[key] = job.line_no
