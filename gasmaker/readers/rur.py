import os
import sys

import numpy as np

from .base import ReadResult


class RurCellReader:
    def __init__(self, repo, iout, mode="nh2", rur_path=None):
        # Default: use an installed `rur`. Otherwise fall back to RUR_PATH, then
        # to an explicit rur_path. No local user path is hard-coded.
        if rur_path is None:
            rur_path = os.environ.get("RUR_PATH")
        if rur_path and rur_path not in sys.path:
            sys.path.insert(0, rur_path)
        from rur import uri

        self.snapshot = uri.RamsesSnapshot(
            repo, iout, mode=mode, verbose=0
        )

    @property
    def maximum_cell_half_diagonal(self):
        levelmin = int(self.snapshot.params["levelmin"])
        return np.sqrt(3.0) * 0.5 * np.exp2(-levelmin)

    def hydro_fields(self):
        """Return (all_fields, chem_elements) available in the hydro data.

        Part of the CellReader contract so the pipeline never reaches into a
        backend-specific snapshot object.
        """
        return self.snapshot.hydro_desc()

    def read_boxes(
        self,
        boxes,
        fields,
        read_grav=False,
        nthread=1,
        read_particles=False,
    ):
        snap = self.snapshot
        cpu_lists = []
        for box in boxes:
            snap.box = np.asarray(box)
            cpu_lists.append(snap.get_involved_cpu())
        cpus = np.unique(np.concatenate(cpu_lists)).astype(np.int32)
        if read_grav:
            missing = [
                int(cpu) for cpu in cpus
                if not os.path.exists(snap.get_path("grav", int(cpu)))
            ]
            if missing:
                raise FileNotFoundError(
                    f"RAMSES grav files are unavailable for CPUs {missing}"
                )
        requested = list(dict.fromkeys(
            ["x", "y", "z", "level", "cpu", *fields]
        ))
        cells = snap.get_cell(
            box=None,
            target_fields=requested,
            domain_slicing=True,
            exact_box=False,
            cpulist=cpus,
            read_grav=read_grav,
            python=True,
            nthread=nthread,
            use_cache=False,
        ).table.copy()
        dms = None
        stars = None
        part_fields = ["x", "y", "z", "m", "cpu"]
        if read_particles:
            dms = snap.get_part(
                box=None,
                pname="dm",
                target_fields=part_fields,
                domain_slicing=True,
                exact_box=False,
                cpulist=cpus,
                python=True,
                nthread=nthread,
                use_cache=False,
            ).table.copy()
            snap.clear()
            try:
                stars = snap.get_part(
                    box=None,
                    pname="star",
                    target_fields=part_fields,
                    domain_slicing=True,
                    exact_box=False,
                    cpulist=cpus,
                    python=True,
                    nthread=nthread,
                    use_cache=False,
                ).table.copy()
            except Exception:
                stars = None
        units = {
            "Msol": float(snap.unit["Msol"]),
            "K": float(snap.unit["K"]),
            "H/cc": float(snap.unit["H/cc"]),
            "kpc": float(snap.unit["kpc"]),
            "box_mpc": float(snap.params["boxsize_physical"]),
            "km/s": float(snap.unit["km/s"]),
            "H0": float(snap.params["H0"]),
            "aexp": float(snap.params["aexp"]),
        }
        h0 = units["H0"]
        aexp = units["aexp"]
        h02 = (h0 * 3.24078e-20) ** 2
        grav = 6.6743e-11
        rho_crit = 3.0 * h02 / 8.0 / np.pi / grav
        rho_crit *= 5.02785e-31 * (3.086e19) ** 3
        units["rho_crit_msol_pkpc3"] = rho_crit / (aexp ** 3)
        snap.clear()
        return ReadResult(cells=cells, dms=dms, stars=stars, cpus=cpus, units=units)

    def close(self):
        self.snapshot.clear()
