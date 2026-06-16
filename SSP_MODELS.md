# SSP Photometry Models

HaloMaker computes all three SSP variants in one run. Each HDF5 photometry
dataset is row-aligned with `/catalog/halo`; the `id` field provides an
additional explicit check.

## Common Calculation

- Intrinsic rest-frame magnitudes
- No foreground or internal dust attenuation
- No nebular emission
- SDSS `ugriz` in AB magnitudes
- Johnson `UBV` in Vega magnitudes
- Stellar particle mass is treated as the SSP normalization mass
- Bilinear interpolation in log age and log metallicity, in magnitude space
- Values outside a model grid are clipped to the nearest grid boundary

For Ra4 snapshots with `initial_mass` in `part_file_descriptor.txt`,
photometry uses the stored initial stellar particle mass. For formats such as
Ra3 that do not provide it, current stellar mass is used as a fallback. The
choice is recorded as `mass_source=initial_mass` or
`mass_source=current_mass_fallback` in every HDF5 model group.

## CB07

- Model: Charlot & Bruzual 2007
- IMF: Chabrier
- Tracks: Padova 1994 + Marigo 2007
- Spectral library: BaSeL
- Metallicities: 7
- Ages: 220
- `Kmag`: Johnson K, Vega

The compact runtime table is generated at build time from the CB07 source
tables used by RUR. The standard local source location is
`assets/ssp_originals/cb07` when available; set `CB07_PATH` to override it.
The CB07 source tables and generated compact table are not redistributed with
HaloMaker.

## BC03

- Model: Bruzual & Charlot 2003
- Release: original 2003 model set
- IMF: Chabrier, 0.1-100 solar masses
- Tracks: Padova 1994 + S. Charlot 1997
- Metallicities: 6
- Ages: 220
- `Kmag`: Johnson K, Vega

Source:
`https://www.bruzual.org/bc03/Original_version_2003/`

The compact runtime table is generated at build time from the BC03
Chabrier/Padova 1994 source tarball or extracted directory. The standard local
source location is `assets/ssp_originals/bc03` when available; set `BC03_PATH`
to override it. The BC03 source tables and generated compact table are not
redistributed with HaloMaker.

## FSPS

- Model: Flexible Stellar Population Synthesis
- python-fsps: 0.4.7
- FSPS commit: `82a873508d500ca353bbb922459bf928498f7a72`
- IMF: Chabrier (`imf_type=1`)
- Tracks: MIST
- Spectral library: MILES
- Metallicities: 12
- Ages: 107
- `Kmag`: 2MASS Ks, Vega

Source: `https://github.com/cconroy20/fsps`

The FSPS compact table is generated locally at build time and is not
redistributed with HaloMaker. Set `FSPS_PATH` when running `build.sh` for the
first time. Once generated, FSPS itself is not needed at runtime.

## Redistribution Note

Generated BC03, CB07, and FSPS compact tables are intentionally excluded from
this repository because their upstream model data have separate distribution
terms. Before publishing or redistributing derived tables, confirm permission
from the relevant model authors.
