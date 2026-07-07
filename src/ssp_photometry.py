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

_TABLES = {}


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
    with np.load(path, allow_pickle=False) as source:
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
