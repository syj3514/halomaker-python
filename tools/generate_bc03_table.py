import argparse
import os
import re
import tarfile
import tempfile
from pathlib import Path

import numpy as np


FIELDS = (
    "umag", "gmag", "rmag", "imag", "zmag",
    "Umag", "Bmag", "Vmag", "Kmag",
)
SOURCE_PATTERN = re.compile(r"^bc2003_lr_m(?P<tag>\d+)_chab_ssp\.1ABmag$")
METALLICITY_PATTERN = re.compile(r"\bZ\s*=\s*([0-9.eE+-]+)")
RELATIVE_MODEL_DIR = Path("bc03") / "models" / "Padova1994" / "chabrier"


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


def _discover_sources(path):
    if not path.is_dir():
        return []
    pairs = []
    for ab_path in path.glob("bc2003_lr_m*_chab_ssp.1ABmag"):
        match = SOURCE_PATTERN.match(ab_path.name)
        if not match:
            continue
        tag = int(match.group("tag"))
        color_path = ab_path.with_suffix(".1color")
        if not color_path.is_file():
            raise SystemExit(f"Missing BC03 color table for {ab_path}: {color_path}")
        pairs.append((tag, _read_metallicity(ab_path), ab_path, color_path))
    return sorted(pairs, key=lambda item: item[0])


def _find_source_dir(source):
    source = source.expanduser().resolve()
    if source.is_file() and tarfile.is_tarfile(source):
        return source
    candidates = [
        source,
        source / RELATIVE_MODEL_DIR,
        source / "models" / "Padova1994" / "chabrier",
    ]
    for candidate in candidates:
        if _discover_sources(candidate):
            return candidate
    raise SystemExit(
        "BC03 source files were not found. Expected a directory containing "
        f"{RELATIVE_MODEL_DIR} or the chabrier model files themselves: {source}"
    )


def _model_dir_from_tar(tar_path, temp_dir):
    with tarfile.open(tar_path) as archive:
        try:
            archive.extractall(temp_dir, filter="data")
        except TypeError:
            archive.extractall(temp_dir)
    model_dir = Path(temp_dir) / RELATIVE_MODEL_DIR
    if not _discover_sources(model_dir):
        raise SystemExit(
            f"BC03 tarball does not contain expected files under {RELATIVE_MODEL_DIR}: {tar_path}"
        )
    return model_dir


def _build_grid(source_dir):
    grids = []
    age_grid = None
    metallicities = []
    sources = _discover_sources(source_dir)
    if not sources:
        raise SystemExit(f"No BC03 source table pairs found in {source_dir}")
    for tag, metallicity, ab_path, color_path in sources:
        ab = _read_table(ab_path)
        colors = _read_table(color_path)

        if len(ab) != len(colors):
            raise ValueError(f"Row count mismatch for {ab_path} and {color_path}")
        current_age = np.asarray(ab["logageyr"], dtype=np.float64)
        if age_grid is None:
            age_grid = current_age
        elif not np.array_equal(age_grid, current_age):
            raise ValueError(f"Age grid mismatch for BC03 metallicity tag {tag}")

        gmag = np.asarray(ab["g_AB"], dtype=np.float64)
        values = np.column_stack(
            (
                gmag + np.asarray(ab["ugAB"], dtype=np.float64),
                gmag,
                gmag - np.asarray(ab["grAB"], dtype=np.float64),
                gmag - np.asarray(ab["giAB"], dtype=np.float64),
                gmag - np.asarray(ab["gzAB"], dtype=np.float64),
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
    source = _find_source_dir(source)
    if source.is_file():
        with tempfile.TemporaryDirectory() as temp_dir:
            source_dir = _model_dir_from_tar(source, temp_dir)
            _write_table(source_dir, output)
    else:
        _write_table(source, output)


def _write_table(source_dir, output):
    log_age, metallicities, magnitudes = _build_grid(source_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        log_age=log_age,
        metallicity=metallicities,
        field_names=np.asarray(FIELDS),
        magnitudes=magnitudes,
        meta_model=np.asarray("Bruzual & Charlot 2003"),
        meta_model_version=np.asarray("Original BC03 release, 2003-08-15"),
        meta_IMF=np.asarray("Chabrier (2003), 0.1-100 Msun"),
        meta_stellar_tracks=np.asarray("Padova 1994 + S. Charlot 1997"),
        meta_spectral_library=np.asarray("STELIB/BaSeL low-resolution product"),
        meta_mass_normalization=np.asarray("1 Msun initially formed"),
        meta_K_filter=np.asarray("Johnson K"),
        meta_source=np.asarray(
            "https://www.bruzual.org/bc03/Original_version_2003/"
            "bc03.models.padova_1994_chabrier_imf.tar.gz"
        ),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bc03-path",
        type=Path,
        help=(
            "BC03 source tarball or extracted directory "
            "(priority: CLI > BC03_PATH > assets/ssp_originals/bc03)"
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
    output = args.output or root / "halomaker_data" / "ssp_tables" / "bc03.npz"
    if output.exists() and not args.force:
        print(f"Using existing {output}")
        return

    source = args.bc03_path
    if source is None:
        value = os.environ.get("BC03_PATH")
        source = Path(value) if value else root / "assets" / "ssp_originals" / "bc03"
    if source is None:
        raise SystemExit(
            "Set BC03_PATH or pass --bc03-path with a BC03 source tarball/directory"
        )
    if not source.expanduser().exists():
        raise SystemExit(f"BC03 source path does not exist: {source}")

    generate(source, output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
