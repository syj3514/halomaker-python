import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from _ssp_common import FIELDS


AB_BANDS = ("sdss_u", "sdss_g", "sdss_r", "sdss_i", "sdss_z")
VEGA_BANDS = ("u", "b", "v", "2mass_ks")
FSPS_COMMIT = "82a873508d500ca353bbb922459bf928498f7a72"


def generate_partial(system, output):
    import fsps

    use_vega = system == "vega"
    bands = VEGA_BANDS if use_vega else AB_BANDS
    population = fsps.StellarPopulation(
        compute_vega_mags=use_vega,
        zcontinuous=0,
        imf_type=1,
        sfh=0,
        dust1=0.0,
        dust2=0.0,
        add_dust_emission=False,
        add_agb_dust_model=False,
        add_neb_emission=False,
    )
    grids = [
        population.get_mags(zmet=zmet, tage=0.0, bands=bands)
        for zmet in range(1, len(population.zlegend) + 1)
    ]
    np.savez_compressed(
        output,
        log_age=np.asarray(population.ssp_ages, dtype=np.float64),
        metallicity=np.asarray(population.zlegend, dtype=np.float64),
        magnitudes=np.asarray(grids, dtype=np.float64),
        fsps_version=np.asarray(fsps.__version__),
        libraries=np.asarray([
            item.decode() if isinstance(item, bytes) else str(item)
            for item in population.libraries
        ]),
    )


def combine(ab_path, vega_path, output):
    with np.load(ab_path, allow_pickle=False) as ab:
        with np.load(vega_path, allow_pickle=False) as vega:
            if not np.array_equal(ab["log_age"], vega["log_age"]):
                raise ValueError("FSPS AB and Vega age grids differ")
            if not np.array_equal(ab["metallicity"], vega["metallicity"]):
                raise ValueError("FSPS AB and Vega metallicity grids differ")
            np.savez_compressed(
                output,
                log_age=ab["log_age"],
                metallicity=ab["metallicity"],
                field_names=np.asarray(FIELDS),
                magnitudes=np.concatenate(
                    (ab["magnitudes"], vega["magnitudes"]), axis=2
                ),
                meta_model=np.asarray(
                    "Flexible Stellar Population Synthesis"
                ),
                meta_model_version=np.asarray(
                    f"python-fsps {str(ab['fsps_version'])}"
                ),
                meta_FSPS_commit=np.asarray(FSPS_COMMIT),
                meta_IMF=np.asarray(
                    "Chabrier (2003), FSPS imf_type=1"
                ),
                meta_libraries=np.asarray(
                    ",".join(str(item) for item in ab["libraries"])
                ),
                meta_stellar_tracks=np.asarray("MIST"),
                meta_spectral_library=np.asarray("MILES"),
                meta_mass_normalization=np.asarray(
                    "1 Msun initially formed"
                ),
                meta_K_filter=np.asarray("2MASS Ks (stored as Kmag)"),
                meta_source=np.asarray(
                    "https://github.com/cconroy20/fsps"
                ),
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--partial", choices=("ab", "vega"))
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--fsps-path",
        type=Path,
        help="FSPS source/data installation (also accepted via FSPS_PATH or SPS_HOME)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing compact table",
    )
    args = parser.parse_args()
    if args.partial:
        generate_partial(args.partial, args.output)
        return

    root = Path(__file__).resolve().parents[1]
    output = root / "halomaker_data" / "ssp_tables" / "fsps.npz"
    if output.exists() and not args.force:
        print(f"Using existing {output}")
        return

    fsps_path = args.fsps_path
    if fsps_path is None:
        value = os.environ.get("FSPS_PATH") or os.environ.get("SPS_HOME")
        fsps_path = Path(value) if value else None
    if fsps_path is None:
        raise SystemExit(
            "Set FSPS_PATH or pass --fsps-path with an FSPS source/data installation"
        )
    fsps_path = fsps_path.expanduser().resolve()
    if not fsps_path.is_dir():
        raise SystemExit(f"FSPS source/data directory does not exist: {fsps_path}")
    os.environ["SPS_HOME"] = str(fsps_path)

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        for system in ("ab", "vega"):
            subprocess.run(
                [
                    sys.executable, __file__, "--partial", system,
                    "--output", str(temp / f"{system}.npz"),
                ],
                check=True,
            )
        combine(temp / "ab.npz", temp / "vega.npz", output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
