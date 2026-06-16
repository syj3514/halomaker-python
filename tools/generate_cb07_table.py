import argparse
import os
import re
from pathlib import Path

import numpy as np


FIELDS = (
    "umag", "gmag", "rmag", "imag", "zmag",
    "Umag", "Bmag", "Vmag", "Kmag",
)
SOURCE_PATTERN = re.compile(r"^cb2007_lr_BaSeL_m(?P<tag>\d+)_chab_ssp\.1ABmag$")
METALLICITY_PATTERN = re.compile(r"\bZ\s*=\s*([0-9.eE+-]+)")


def _read_table(path):
    lines = path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("#log"))
    names = lines[header][1:].split()
    return np.genfromtxt(path, names=names, skip_header=header)


def _read_metallicity(path):
    with path.open() as handle:
        for line in handle:
            match = METALLICITY_PATTERN.search(line)
            if match:
                return float(match.group(1))
    raise ValueError(f"Could not find metallicity Z=... in {path}")


def _discover_sources(source_dir):
    if not source_dir.is_dir():
        return []
    pairs = []
    for ab_path in source_dir.glob("cb2007_lr_BaSeL_m*_chab_ssp.1ABmag"):
        match = SOURCE_PATTERN.match(ab_path.name)
        if not match:
            continue
        tag = int(match.group("tag"))
        color_path = ab_path.with_suffix(".1color")
        if not color_path.is_file():
            raise SystemExit(f"Missing CB07 color table for {ab_path}: {color_path}")
        pairs.append((tag, _read_metallicity(ab_path), ab_path, color_path))
    return sorted(pairs, key=lambda item: item[0])


def _find_source_dir(source):
    source = source.expanduser().resolve()
    candidates = [source, source / "cb07"]
    for candidate in candidates:
        if _discover_sources(candidate):
            return candidate
    raise SystemExit(
        "CB07 source files were not found. Expected pairs named "
        "cb2007_lr_BaSeL_m*_chab_ssp.{1ABmag,1color} under "
        f"{source} or {source / 'cb07'}"
    )


def _build_grid(source_dir):
    grids = []
    age_grid = None
    metallicities = []
    sources = _discover_sources(source_dir)
    if not sources:
        raise SystemExit(f"No CB07 source table pairs found in {source_dir}")
    for tag, metallicity, ab_path, color_path in sources:
        ab = _read_table(ab_path)
        colors = _read_table(color_path)

        if len(ab) != len(colors):
            raise ValueError(f"Row count mismatch for {ab_path} and {color_path}")
        current_age = np.asarray(ab["logageyr"], dtype=np.float64)
        if age_grid is None:
            age_grid = current_age
        elif not np.array_equal(age_grid, current_age):
            raise ValueError(f"Age grid mismatch for CB07 metallicity tag {tag}")

        values = np.column_stack(
            (
                ab["u"],
                ab["g"],
                ab["r"],
                ab["i"],
                ab["z"],
                colors["Umag"],
                colors["Bmag"],
                colors["Vmag"],
                colors["Kmag"],
            )
        )
        metallicities.append(metallicity)
        grids.append(values)
    return (
        age_grid,
        np.asarray(metallicities, dtype=np.float64),
        np.asarray(grids, dtype=np.float64),
    )


def generate(source, output):
    source_dir = _find_source_dir(source)
    log_age, metallicities, magnitudes = _build_grid(source_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        log_age=log_age,
        metallicity=metallicities,
        field_names=np.asarray(FIELDS),
        magnitudes=magnitudes,
        meta_model=np.asarray("Charlot & Bruzual 2007"),
        meta_model_version=np.asarray("CB07 local table set"),
        meta_IMF=np.asarray("Chabrier"),
        meta_stellar_tracks=np.asarray("Padova 1994 + Marigo 2007"),
        meta_spectral_library=np.asarray("BaSeL"),
        meta_mass_normalization=np.asarray("1 Msun initially formed"),
        meta_K_filter=np.asarray("Johnson K"),
        meta_source=np.asarray("RUR CB07 table set"),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cb07-path",
        type=Path,
        help=(
            "CB07 source directory "
            "(priority: CLI > CB07_PATH > assets/ssp_originals/cb07)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="output compact table path",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing compact table",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    output = args.output or root / "halomaker_data" / "ssp_tables" / "cb07.npz"
    if output.exists() and not args.force:
        print(f"Using existing {output}")
        return

    source = args.cb07_path
    if source is None:
        value = os.environ.get("CB07_PATH")
        source = Path(value) if value else root / "assets" / "ssp_originals" / "cb07"
    if not source.expanduser().is_dir():
        raise SystemExit(f"CB07 source directory does not exist: {source}")

    generate(source, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
