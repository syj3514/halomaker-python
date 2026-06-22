"""Reader interface for GasMaker.

GasMaker's core (pipeline/catalog/geometry/overlap) is reader-agnostic: it never
imports a specific simulation reader. Instead it consumes the `ReadResult`
contract below and calls the small `CellReader` interface. This keeps the heavy
RAMSES I/O dependency (e.g. the optional `rur` adapter in `readers/rur.py`) out
of the core import path — install/import GasMaker without `rur`, and only pull
a reader in when you actually read snapshot data.

To support another simulation/format, implement a class with the same three
members and return a `ReadResult`.
"""
from dataclasses import dataclass
from typing import Protocol, Sequence

import numpy as np


@dataclass
class ReadResult:
    """Data contract returned by a reader for a set of bounding boxes.

    - cells: structured ndarray with at least x, y, z, level (+ requested hydro
      fields such as rho, P, vx, vy, vz, metal, chemical fractions). Positions
      are in code units [0, 1).
    - dms / stars: structured ndarrays (x, y, z, m, cpu) or None.
    - cpus: int32 array of the CPU/domain files that were read.
    - units: dict of conversion factors and cosmology, with keys
      Msol, K, H/cc, kpc, km/s, box_mpc, H0, aexp, rho_crit_msol_pkpc3.
    """

    cells: np.ndarray
    dms: np.ndarray | None
    stars: np.ndarray | None
    cpus: np.ndarray
    units: dict


class CellReader(Protocol):
    """Minimal interface GasMaker needs from a snapshot reader."""

    @property
    def maximum_cell_half_diagonal(self) -> float:
        ...

    def read_boxes(
        self,
        boxes: Sequence,
        fields: Sequence[str],
        read_grav: bool = False,
        nthread: int = 1,
        read_particles: bool = False,
    ) -> ReadResult:
        ...

    def close(self) -> None:
        ...
