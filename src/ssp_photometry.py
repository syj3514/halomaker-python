from importlib.resources import files

import numpy as np


BANDS = (
    ("umag", "SDSS u", "AB"),
    ("gmag", "SDSS g", "AB"),
    ("rmag", "SDSS r", "AB"),
    ("imag", "SDSS i", "AB"),
    ("zmag", "SDSS z", "AB"),
    ("Umag", "Johnson U", "Vega"),
    ("Bmag", "Johnson B", "Vega"),
    ("Vmag", "Johnson V", "Vega"),
    ("Kmag", "K", "Vega"),
)
BAND_INDEX = {name: index for index, (name, _, _) in enumerate(BANDS)}
MODELS = ("CB07", "BC03", "FSPS")

SSP_SOURCES = {
    "BC03": {
        "environment": "BC03_PATH",
        "url": "https://www.bruzual.org/bc03/Original_version_2003/",
    },
    "CB07": {
        "environment": "CB07_PATH",
        "url": "https://www.bruzual.org/cb07/",
    },
    "FSPS": {
        "environment": "FSPS_PATH (or SPS_HOME)",
        "url": "https://github.com/cconroy20/fsps",
    },
}

_TABLES = {}


def missing_table_message(model, path):
    model = model.upper()
    source = SSP_SOURCES[model]
    filename = f"{model.lower()}.npz"
    return f"""Missing {model} SSP photometry table.
Expected compact table: {path}
Photometry is enabled and requires all three compact SSP tables.

Resolve this in one of these ways:
  1. Generate it from the original model data: set {source['environment']},
     then rerun build.sh. Source information: {source['url']}
  2. Copy a trusted prebuilt {filename} into halomaker_data/ssp_tables/
     and verify its checksum against the provider's manifest.
  3. For a halo-only run, set `photometry = .false.` in input_HaloMaker.dat.
     Build without SSP preparation using `HALOMAKER_SKIP_SSP=1 bash build.sh`.
"""


def _text(value):
    value = np.asarray(value)
    if value.ndim == 0:
        return str(value.item())
    return ",".join(str(item) for item in value.tolist())


def load_table(model):
    model = model.upper()
    if model not in MODELS:
        raise ValueError(f"Unknown SSP model: {model}")
    if model in _TABLES:
        return _TABLES[model]

    path = files("halomaker_data.ssp_tables").joinpath(
        f"{model.lower()}.npz"
    )
    if not path.is_file():
        raise FileNotFoundError(missing_table_message(model, path))
    try:
        source_context = np.load(path, allow_pickle=False)
    except FileNotFoundError:
        raise FileNotFoundError(missing_table_message(model, path)) from None
    with source_context as source:
        field_names = tuple(str(name) for name in source["field_names"])
        if field_names != tuple(name for name, _, _ in BANDS):
            raise ValueError(f"Unexpected photometry fields in {path}: {field_names}")
        metadata = {
            name.removeprefix("meta_"): _text(source[name])
            for name in source.files
            if name.startswith("meta_")
        }
        table = {
            "model": model,
            "log_age": np.asarray(source["log_age"], dtype=np.float64),
            "log_metal": np.log10(
                np.asarray(source["metallicity"], dtype=np.float64)
            ),
            "magnitudes": np.asarray(source["magnitudes"], dtype=np.float64),
            "metadata": metadata,
        }
    _TABLES[model] = table
    return table


def load_all_tables():
    return {model: load_table(model) for model in MODELS}


def model_metadata(model):
    table = load_table(model)
    return table["model"], table["metadata"]


def interpolation_coordinates(model, age_gyr, metal):
    table = load_table(model)
    age_grid = table["log_age"]
    metal_grid = table["log_metal"]
    log_age = np.log10(
        np.clip(age_gyr * 1.0e9, 10.0**age_grid[0], 10.0**age_grid[-1])
    )
    log_metal = np.log10(
        np.clip(metal, 10.0**metal_grid[0], 10.0**metal_grid[-1])
    )

    ia1 = np.clip(np.searchsorted(age_grid, log_age, side="right"), 1, len(age_grid) - 1)
    iz1 = np.clip(np.searchsorted(metal_grid, log_metal, side="right"), 1, len(metal_grid) - 1)
    ia0 = ia1 - 1
    iz0 = iz1 - 1

    age0 = age_grid[ia0]
    age1 = age_grid[ia1]
    metal0 = metal_grid[iz0]
    metal1 = metal_grid[iz1]
    wa = np.divide(
        log_age - age0, age1 - age0,
        out=np.zeros_like(log_age), where=age1 > age0,
    )
    wz = np.divide(
        log_metal - metal0, metal1 - metal0,
        out=np.zeros_like(log_metal), where=metal1 > metal0,
    )
    return ia0, ia1, iz0, iz1, wa, wz


def interpolate_magnitude(model, band, coordinates):
    ia0, ia1, iz0, iz1, wa, wz = coordinates
    iband = BAND_INDEX[band]
    grid = load_table(model)["magnitudes"][:, :, iband]
    mag0 = grid[iz0, ia0] * (1.0 - wa) + grid[iz0, ia1] * wa
    mag1 = grid[iz1, ia0] * (1.0 - wa) + grid[iz1, ia1] * wa
    return mag0 * (1.0 - wz) + mag1 * wz
