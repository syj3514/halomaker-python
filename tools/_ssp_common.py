"""Shared helpers for the BC03 / CB07 / FSPS SSP compact-table generators.

These build-time generators (`generate_{bc03,cb07,fsps}_table.py`) all emit the
same 9-column magnitude layout and, for the text-table models (BC03/CB07), parse
the same `#log...`-headered source files. The identical pieces live here so the
per-model scripts only carry what genuinely differs (source globs and the
model-specific magnitude derivation in their own `_build_grid`).
"""

import re

import numpy as np


# Output magnitude columns, in order — shared by every SSP compact table.
FIELDS = (
    "umag", "gmag", "rmag", "imag", "zmag",
    "Umag", "Bmag", "Vmag", "Kmag",
)

METALLICITY_PATTERN = re.compile(r"\bZ\s*=\s*([0-9.eE+-]+)")


def read_table(path):
    """Read a `#log...`-headered SSP text table into a numpy structured array."""
    lines = path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("#log"))
    names = lines[header][1:].split()
    return np.genfromtxt(path, names=names, skip_header=header)


def read_metallicity(path):
    """Extract the `Z=...` metallicity value from an SSP table header."""
    with path.open() as handle:
        for line in handle:
            match = METALLICITY_PATTERN.search(line)
            if match:
                return float(match.group(1))
    raise ValueError(f"Could not find metallicity Z=... in {path}")


def discover_sources(source_dir, glob_pattern, source_pattern, model_label):
    """Find (tag, Z, ABmag_path, color_path) source pairs, sorted by tag.

    `glob_pattern` selects the ABmag files, `source_pattern` is a compiled regex
    with a `tag` group, and `model_label` names the model in error messages.
    """
    if not source_dir.is_dir():
        return []
    pairs = []
    for ab_path in source_dir.glob(glob_pattern):
        match = source_pattern.match(ab_path.name)
        if not match:
            continue
        tag = int(match.group("tag"))
        color_path = ab_path.with_suffix(".1color")
        if not color_path.is_file():
            raise SystemExit(
                f"Missing {model_label} color table for {ab_path}: {color_path}"
            )
        pairs.append((tag, read_metallicity(ab_path), ab_path, color_path))
    return sorted(pairs, key=lambda item: item[0])
