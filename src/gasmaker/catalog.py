from dataclasses import dataclass

import h5py
import numpy as np


@dataclass(frozen=True)
class HaloCatalog:
    halos: np.ndarray
    box_physical_mpc: float
    box_comoving_mpc: float
    aexp: float
    units_version: str

    @classmethod
    def read(cls, path):
        with h5py.File(path, "r") as hdf:
            halos = hdf["catalog/halo"][:]
            box_comoving_mpc = float(hdf["header"].attrs["boxsize2"])
            aexp = float(hdf["header"].attrs["aexp"])
            units_version = str(hdf["header"].attrs.get("units_version", "legacy"))
        return cls(
            halos=halos,
            box_physical_mpc=box_comoving_mpc * aexp,
            box_comoving_mpc=box_comoving_mpc,
            aexp=aexp,
            units_version=units_version,
        )

    def row_for_id(self, halo_id):
        rows = np.flatnonzero(self.halos["id"] == halo_id)
        if rows.size != 1:
            raise ValueError(f"Expected one halo id={halo_id}, found {rows.size}")
        return int(rows[0])

    def descendants(self, root_id):
        root_row = self.row_for_id(root_id)
        root = self.halos[root_row]
        if root["level"] != 1:
            raise ValueError(f"Halo id={root_id} is not level 1")

        selected = np.zeros(self.halos.size, dtype=bool)
        selected[root_row] = True
        changed = True
        while changed:
            parent_ids = self.halos["id"][selected]
            children = (
                np.isin(self.halos["hosthalo"], parent_ids)
                | np.isin(self.halos["hostsub"], parent_ids)
            )
            children[root_row] = True
            updated = selected | children
            changed = not np.array_equal(updated, selected)
            selected = updated
        return np.flatnonzero(selected)

    def root_rows(self):
        return np.flatnonzero(self.halos["level"] == 1)
