import argparse
import os
from pathlib import Path

import numpy as np


FIELDS = (
    "umag", "gmag", "rmag", "imag", "zmag",
    "Umag", "Bmag", "Vmag", "Kmag",
)
METALLICITIES = (0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05, 0.1)
TAGS = (22, 32, 42, 52, 62, 72, 82)


def _read_table(path):
    lines = path.read_text().splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("#log"))
    names = lines[header][1:].split()
    return np.genfromtxt(path, names=names, skip_header=header)


def _required_files(source_dir):
    files = []
    for tag in TAGS:
        stem = source_dir / f"cb2007_lr_BaSeL_m{tag}_chab_ssp"
        files.append(stem.with_suffix(".1ABmag"))
        files.append(stem.with_suffix(".1color"))
    return files


def _find_source_dir(source):
    source = source.expanduser().resolve()
    candidates = [source, source / "cb07"]
    for candidate in candidates:
        missing = [path for path in _required_files(candidate) if not path.is_file()]
        if not missing:
            return candidate
    missing_text = "\n".join(str(path) for path in missing[:6])
    if len(missing) > 6:
        missing_text += f"\n... and {len(missing) - 6} more"
    raise SystemExit(
        "CB07 source files were not found. Expected 14 files named "
        "cb2007_lr_BaSeL_m{22,32,42,52,62,72,82}_chab_ssp.{1ABmag,1color}. "
        f"Missing examples:\n{missing_text}"
    )


def _build_grid(source_dir):
    grids = []
    age_grid = None
    for tag in TAGS:
        stem = source_dir / f"cb2007_lr_BaSeL_m{tag}_chab_ssp"
        ab = _read_table(stem.with_suffix(".1ABmag"))
        colors = _read_table(stem.with_suffix(".1color"))

        if len(ab) != len(colors):
            raise ValueError(f"Row count mismatch for {stem}")
        current_age = np.asarray(ab["logageyr"], dtype=np.float64)
        if age_grid is None:
            age_grid = current_age
        elif not np.array_equal(age_grid, current_age):
            raise ValueError(f"Age grid mismatch for {stem}")

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
        grids.append(values)
    return age_grid, np.asarray(grids, dtype=np.float64)


def generate(source, output):
    source_dir = _find_source_dir(source)
    log_age, magnitudes = _build_grid(source_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        log_age=log_age,
        metallicity=np.asarray(METALLICITIES, dtype=np.float64),
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
        help="CB07 source directory (also accepted via CB07_PATH)",
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
        source = Path(value) if value else None
    if source is None:
        raise SystemExit("Set CB07_PATH or pass --cb07-path with a CB07 source directory")
    if not source.expanduser().is_dir():
        raise SystemExit(f"CB07 source directory does not exist: {source}")

    generate(source, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
