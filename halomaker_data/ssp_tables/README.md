# SSP photometry tables

This directory holds the SSP (simple stellar population) photometry tables used
by `ssp_photometry.py` for rest-frame magnitudes (CB07 / BC03 / FSPS).

The `.npz` tables are **not committed** (they are large and license-bound) — only
this folder and its `__init__.py` are tracked. Generate the tables before first
use:

```bash
FSPS_PATH=/path/to/fsps bash build.sh
```

See `SSP_MODELS.md` (repo root) for the full table-generation and FSPS setup
instructions.
