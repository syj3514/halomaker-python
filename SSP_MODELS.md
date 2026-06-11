# Bundled SSP Photometry Models

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

The current RAMSES Ra3 input does not provide initial stellar particle mass.
Consequently, photometry uses current stellar mass as a fallback. This is
recorded as `mass_source=current_mass_fallback` in every HDF5 model group.

## CB07

- Model: Charlot & Bruzual 2007
- IMF: Chabrier
- Tracks: Padova 1994 + Marigo 2007
- Spectral library: BaSeL
- Metallicities: 7
- Ages: 220
- `Kmag`: Johnson K, Vega

The compact runtime table was generated from the CB07 tables used by RUR.

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

FSPS is needed only to regenerate the compact table, not to run HaloMaker.

## Redistribution Note

FSPS source code is MIT licensed. BC03 and CB07 model tables have separate
upstream terms and citations. Before publishing or redistributing this
repository, the maintainer must confirm that distribution of the compact
derived BC03 and CB07 tables is permitted by the model authors.
