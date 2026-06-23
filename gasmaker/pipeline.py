from dataclasses import dataclass
from pathlib import Path
import json
import time

import h5py
import numpy as np

from .catalog import HaloCatalog
from .geometry import (
    halo_center_code,
    halo_radius_code,
    periodic_boxes,
    sphere_cell_mask,
    periodic_delta,
)
from .overlap import fractional_sphere_cell_overlap


SUMMARY_DTYPE = np.dtype([
    ("id", "i4"),
    ("root_id", "i4"),
    
    # Diagnostics (calculated for r* aperture, or defaults if starless)
    ("ncells", "i8"),
    ("ncells_center", "i8"),
    ("ncells_overlap", "i8"),
    ("overlap_resolved", "?"),
    ("radius_to_min_cell", "f8"),
    ("overlap_relative_bound", "f8"),
    
    # Gas masses (4 apertures)
    ("mgas", "f8"),
    ("mgas_r50", "f8"),
    ("mgas_r90", "f8"),
    ("mgas_rvir", "f8"),
    
    ("mcold", "f8"),
    ("mcold_r50", "f8"),
    ("mcold_r90", "f8"),
    ("mcold_rvir", "f8"),
    
    ("mdense", "f8"),
    ("mdense_r50", "f8"),
    ("mdense_r90", "f8"),
    ("mdense_rvir", "f8"),

    # Tier2 spherical overdensity and component masses
    ("r200", "f8"),
    ("m200", "f8"),
    ("r500", "f8"),
    ("m500", "f8"),

    ("mgas_r200", "f8"),
    ("mgas_r500", "f8"),
    ("mcold_r200", "f8"),
    ("mcold_r500", "f8"),
    ("mdense_r200", "f8"),
    ("mdense_r500", "f8"),

    ("mdm_rvir", "f8"),
    ("mdm_r200", "f8"),
    ("mdm_r500", "f8"),
    ("mstar_rvir", "f8"),
    ("mstar_r200", "f8"),
    ("mstar_r500", "f8"),
    
    # Metallicity (representative r* aperture)
    ("metal_gas", "f8"),
    
    # Kinematics (3 apertures: r*, r50, r90)
    ("vrot_gas", "f8"),
    ("sig3d_gas", "f8"),
    ("sigcyl_gas", "f8"),
    
    ("vrot_gas_r50", "f8"),
    ("sig3d_gas_r50", "f8"),
    ("sigcyl_gas_r50", "f8"),
    
    ("vrot_gas_r90", "f8"),
    ("sig3d_gas_r90", "f8"),
    ("sigcyl_gas_r90", "f8"),
    
    # Gas Angular Momentum (r* aperture)
    ("Lx_gas", "f8"),
    ("Ly_gas", "f8"),
    ("Lz_gas", "f8"),
    
    # Chemistry (r* aperture)
    ("H_gas", "f8"),
    ("O_gas", "f8"),
    ("Fe_gas", "f8"),
    ("Mg_gas", "f8"),
    ("C_gas", "f8"),
    ("N_gas", "f8"),
    ("Si_gas", "f8"),
    ("S_gas", "f8"),
    ("D_gas", "f8"),
])


ROOT_METRICS_DTYPE = np.dtype([
    ("root_id", "i4"),
    ("descendant_count", "i8"),
    ("cpu_count", "i8"),
    ("cells_read", "i8"),
    ("envelope_radius_code", "f8"),
    ("read_seconds", "f8"),
    ("compute_seconds", "f8"),
])

SCHEMA_VERSION = 3
UNITS_VERSION = "halomaker_units_v2"
GAS_FIELD_UNITS = {
    "id": "int",
    "root_id": "int",
    "ncells": "count",
    "ncells_center": "count",
    "ncells_overlap": "count",
    "overlap_resolved": "bool",
    "radius_to_min_cell": "dimensionless",
    "overlap_relative_bound": "dimensionless",
    "mgas": "Msun",
    "mgas_r50": "Msun",
    "mgas_r90": "Msun",
    "mgas_rvir": "Msun",
    "mcold": "Msun",
    "mcold_r50": "Msun",
    "mcold_r90": "Msun",
    "mcold_rvir": "Msun",
    "mdense": "Msun",
    "mdense_r50": "Msun",
    "mdense_r90": "Msun",
    "mdense_rvir": "Msun",
    "r200": "code_unit",
    "m200": "Msun",
    "r500": "code_unit",
    "m500": "Msun",
    "mgas_r200": "Msun",
    "mgas_r500": "Msun",
    "mcold_r200": "Msun",
    "mcold_r500": "Msun",
    "mdense_r200": "Msun",
    "mdense_r500": "Msun",
    "mdm_rvir": "Msun",
    "mdm_r200": "Msun",
    "mdm_r500": "Msun",
    "mstar_rvir": "Msun",
    "mstar_r200": "Msun",
    "mstar_r500": "Msun",
    "metal_gas": "mass_fraction",
    "vrot_gas": "km/s",
    "sig3d_gas": "km/s",
    "sigcyl_gas": "km/s",
    "vrot_gas_r50": "km/s",
    "sig3d_gas_r50": "km/s",
    "sigcyl_gas_r50": "km/s",
    "vrot_gas_r90": "km/s",
    "sig3d_gas_r90": "km/s",
    "sigcyl_gas_r90": "km/s",
    "Lx_gas": "Msun Mpc km/s",
    "Ly_gas": "Msun Mpc km/s",
    "Lz_gas": "Msun Mpc km/s",
    "H_gas": "mass_fraction",
    "O_gas": "mass_fraction",
    "Fe_gas": "mass_fraction",
    "Mg_gas": "mass_fraction",
    "C_gas": "mass_fraction",
    "N_gas": "mass_fraction",
    "Si_gas": "mass_fraction",
    "S_gas": "mass_fraction",
    "D_gas": "mass_fraction",
}


@dataclass
class RootMetrics:
    root_id: int
    descendant_count: int
    cpu_count: int
    cells_read: int
    envelope_radius_code: float
    read_seconds: float
    compute_seconds: float

    def as_record(self):
        record = np.zeros((), dtype=ROOT_METRICS_DTYPE)
        record["root_id"] = self.root_id
        record["descendant_count"] = self.descendant_count
        record["cpu_count"] = self.cpu_count
        record["cells_read"] = self.cells_read
        record["envelope_radius_code"] = self.envelope_radius_code
        record["read_seconds"] = self.read_seconds
        record["compute_seconds"] = self.compute_seconds
        return record


class GasMaker:
    def __init__(
        self,
        catalog_path,
        reader,
        radius_field="r",
        padding=1.1,
        include_boundary=True,
        overlap_depth=8,
        overlap_tolerance=0.05,
        overlap_chunk_size=2048,
        overlap_threshold=0.1,
    ):
        self.catalog_path = catalog_path
        self.catalog = HaloCatalog.read(catalog_path)
        self.reader = reader
        self.radius_field = radius_field
        self.padding = padding
        self.include_boundary = include_boundary
        self.overlap_depth = overlap_depth
        self.overlap_tolerance = overlap_tolerance
        self.overlap_chunk_size = overlap_chunk_size
        self.overlap_threshold = overlap_threshold

    def _fractional_overlap(self, cells, center, radius):
        count = cells.size
        if count <= self.overlap_chunk_size:
            return fractional_sphere_cell_overlap(
                cells, center, radius, max_depth=self.overlap_depth, threshold_ratio=self.overlap_threshold
            )
        fraction = np.zeros(count, dtype=np.float64)
        lower = np.zeros(count, dtype=np.float64)
        upper = np.zeros(count, dtype=np.float64)
        for start in range(0, count, self.overlap_chunk_size):
            stop = min(start + self.overlap_chunk_size, count)
            f_chunk, l_chunk, u_chunk = fractional_sphere_cell_overlap(
                cells[start:stop],
                center,
                radius,
                max_depth=self.overlap_depth,
                threshold_ratio=self.overlap_threshold,
            )
            fraction[start:stop] = f_chunk
            lower[start:stop] = l_chunk
            upper[start:stop] = u_chunk
        return fraction, lower, upper

    def _gas_kinematics_1d(self, mass, dx, dy, dz, dvx, dvy, dvz, Lvec):
        sig3d2 = np.average(dvx**2 + dvy**2 + dvz**2, weights=mass)
        sig3d = np.sqrt(max(sig3d2, 0.0))
        
        L_norm = np.linalg.norm(Lvec)
        if L_norm <= 0 or not np.isfinite(L_norm):
            return np.nan, sig3d, np.nan
            
        axis = Lvec / L_norm
        zpos = dx * axis[0] + dy * axis[1] + dz * axis[2]
        rperp_x = dx - zpos * axis[0]
        rperp_y = dy - zpos * axis[1]
        rperp_z = dz - zpos * axis[2]
        rperp = np.sqrt(rperp_x**2 + rperp_y**2 + rperp_z**2)
        
        valid = rperp > np.finfo(np.float64).eps
        if np.count_nonzero(valid) < 2:
            return np.nan, sig3d, np.nan
            
        inv_rperp = 1.0 / rperp[valid]
        erad_x = rperp_x[valid] * inv_rperp
        erad_y = rperp_y[valid] * inv_rperp
        erad_z = rperp_z[valid] * inv_rperp
        
        ephi_x = axis[1] * erad_z - axis[2] * erad_y
        ephi_y = axis[2] * erad_x - axis[0] * erad_z
        ephi_z = axis[0] * erad_y - axis[1] * erad_x
        
        vrad = dvx[valid] * erad_x + dvy[valid] * erad_y + dvz[valid] * erad_z
        vphi = dvx[valid] * ephi_x + dvy[valid] * ephi_y + dvz[valid] * ephi_z
        vzyl = dvx[valid] * axis[0] + dvy[valid] * axis[1] + dvz[valid] * axis[2]
        
        cyl_weight = mass[valid]
        vrot = np.average(vphi, weights=cyl_weight)
        mean_vrad = np.average(vrad, weights=cyl_weight)
        mean_vzyl = np.average(vzyl, weights=cyl_weight)
        
        sigcyl2 = np.average(
            (vrad - mean_vrad)**2
            + (vphi - vrot)**2
            + (vzyl - mean_vzyl)**2,
            weights=cyl_weight,
        )
        sigcyl = np.sqrt(max(sigcyl2, 0.0))
        
        return vrot, sig3d, sigcyl

    def _distance_to_center(self, table, center):
        if table is None or len(table) == 0:
            return np.empty(0, dtype=np.float64)
        dx = periodic_delta(table["x"], center[0])
        dy = periodic_delta(table["y"], center[1])
        dz = periodic_delta(table["z"], center[2])
        return np.sqrt(dx * dx + dy * dy + dz * dz)

    def _particle_mass_msol(self, particles, units):
        if particles is None or len(particles) == 0:
            return np.empty(0, dtype=np.float64)
        return particles["m"] / units["Msol"]

    def _full_cell_mass_msol(self, cells, units):
        if cells is None or len(cells) == 0:
            return np.empty(0, dtype=np.float64)
        volume = np.exp2(-3.0 * cells["level"])
        return cells["rho"] * volume / units["Msol"]

    def _component_mass_within(self, mass, distance, radius):
        if not (np.isfinite(radius) and radius > 0):
            return np.nan
        if mass.size == 0:
            return 0.0
        return np.sum(mass[distance < radius])

    def _spherical_overdensity(self, halo, cells, dms, stars, units):
        center = halo_center_code(
            halo, self.catalog.box_physical_mpc, self.catalog.units_version
        )
        halo_radius = halo_radius_code(
            halo, self.catalog.box_physical_mpc, self.radius_field,
            units_version=self.catalog.units_version,
        )

        cell_dist = self._distance_to_center(cells, center)
        dm_dist = self._distance_to_center(dms, center)
        star_dist = self._distance_to_center(stars, center)

        cell_inside = cell_dist <= halo_radius
        dm_inside = dm_dist <= halo_radius
        star_inside = star_dist <= halo_radius

        cell_dist = cell_dist[cell_inside]
        dm_dist = dm_dist[dm_inside]
        star_dist = star_dist[star_inside]
        cells_inside = cells[cell_inside]
        cell_mass = self._full_cell_mass_msol(cells_inside, units)
        if cells_inside.size:
            temperature = cells_inside["P"] / cells_inside["rho"] / units["K"]
            cell_cold = temperature < 1.0e4
            cell_dense = cell_cold & (cells_inside["rho"] / units["H/cc"] > 5.0)
        else:
            cell_cold = np.zeros(0, dtype=bool)
            cell_dense = np.zeros(0, dtype=bool)
        dm_mass = self._particle_mass_msol(dms[dm_inside] if dms is not None else None, units)
        star_mass = self._particle_mass_msol(
            stars[star_inside] if stars is not None else None, units
        )

        total_count = cell_dist.size + dm_dist.size + star_dist.size
        if total_count == 0:
            return {
                "r200": np.nan,
                "m200": np.nan,
                "r500": np.nan,
                "m500": np.nan,
                "cell_dist": cell_dist,
                "cell_mass": cell_mass,
                "cell_cold": cell_cold,
                "cell_dense": cell_dense,
                "dm_dist": dm_dist,
                "dm_mass": dm_mass,
                "star_dist": star_dist,
                "star_mass": star_mass,
            }

        distance = np.empty(total_count, dtype=np.float64)
        mass = np.empty(total_count, dtype=np.float64)
        cursor = 0
        count = cell_dist.size
        distance[cursor:cursor + count] = cell_dist
        mass[cursor:cursor + count] = cell_mass
        cursor += count
        count = star_dist.size
        distance[cursor:cursor + count] = star_dist
        mass[cursor:cursor + count] = star_mass
        cursor += count
        count = dm_dist.size
        distance[cursor:cursor + count] = dm_dist
        mass[cursor:cursor + count] = dm_mass

        valid = np.isfinite(distance) & np.isfinite(mass) & (distance > 0) & (mass > 0)
        if np.count_nonzero(valid) == 0:
            r200 = m200 = r500 = m500 = np.nan
        else:
            distance = distance[valid]
            mass = mass[valid]
            order = np.argsort(distance)
            distance = distance[order]
            mass = mass[order]
            distance_pkpc = distance / units["kpc"]
            cumulative_mass = np.cumsum(mass)
            density = cumulative_mass / (4.0 / 3.0 * np.pi * distance_pkpc ** 3)
            r200, m200 = self._overdensity_radius_mass(
                distance,
                distance_pkpc,
                cumulative_mass,
                density,
                200.0 * units["rho_crit_msol_pkpc3"],
            )
            r500, m500 = self._overdensity_radius_mass(
                distance,
                distance_pkpc,
                cumulative_mass,
                density,
                500.0 * units["rho_crit_msol_pkpc3"],
            )

        return {
            "r200": r200,
            "m200": m200,
            "r500": r500,
            "m500": m500,
            "cell_dist": cell_dist,
            "cell_mass": cell_mass,
            "cell_cold": cell_cold,
            "cell_dense": cell_dense,
            "dm_dist": dm_dist,
            "dm_mass": dm_mass,
            "star_dist": star_dist,
            "star_mass": star_mass,
        }

    def _overdensity_radius_mass(
        self, distance_code, distance_pkpc, cumulative_mass, density, target_density
    ):
        if distance_code.size == 0:
            return np.nan, np.nan
        above = density >= target_density
        if not np.any(above):
            return np.nan, np.nan
        index = int(np.flatnonzero(above)[-1])
        if index >= distance_code.size - 1:
            return np.nan, np.nan

        x0 = density[index]
        x1 = density[index + 1]
        if not (np.isfinite(x0) and np.isfinite(x1)) or x0 == x1:
            fraction = 0.0
        else:
            fraction = (target_density - x0) / (x1 - x0)
            fraction = float(np.clip(fraction, 0.0, 1.0))

        radius_code = distance_code[index] + fraction * (
            distance_code[index + 1] - distance_code[index]
        )
        mass = cumulative_mass[index] + fraction * (
            cumulative_mass[index + 1] - cumulative_mass[index]
        )
        return radius_code, mass

    def _apply_tier2(self, output, halo, cells, dms, stars, units):
        tier2 = self._spherical_overdensity(halo, cells, dms, stars, units)
        for field in ("r200", "m200", "r500", "m500"):
            output[field] = tier2[field]

        radius_rvir = halo_radius_code(
            halo, self.catalog.box_physical_mpc, "rvir",
            units_version=self.catalog.units_version,
        )
        aperture_radii = {
            "_rvir": radius_rvir,
            "_r200": tier2["r200"],
            "_r500": tier2["r500"],
        }
        dm_dist = tier2["dm_dist"]
        dm_mass = tier2["dm_mass"]
        star_dist = tier2["star_dist"]
        star_mass = tier2["star_mass"]
        cell_dist = tier2["cell_dist"]
        cell_mass = tier2["cell_mass"]
        cold = tier2["cell_cold"]
        dense = tier2["cell_dense"]

        for suffix, radius in aperture_radii.items():
            output[f"mdm{suffix}"] = self._component_mass_within(
                dm_mass, dm_dist, radius
            )
            output[f"mstar{suffix}"] = self._component_mass_within(
                star_mass, star_dist, radius
            )

        for suffix in ("_r200", "_r500"):
            radius = aperture_radii[suffix]
            gas_mask = cell_dist < radius if np.isfinite(radius) else np.zeros(
                cell_dist.size, dtype=bool
            )
            if np.isfinite(radius):
                output[f"mgas{suffix}"] = np.sum(cell_mass[gas_mask])
                output[f"mcold{suffix}"] = np.sum(cell_mass[gas_mask & cold])
                output[f"mdense{suffix}"] = np.sum(cell_mass[gas_mask & dense])

    def _summarize(self, halo, root_id, cells, dms, stars, units):
        box_physical_mpc = self.catalog.box_physical_mpc
        output = np.zeros((), dtype=SUMMARY_DTYPE)
        output["id"] = halo["id"]
        output["root_id"] = root_id
        
        for field in SUMMARY_DTYPE.names:
            if np.issubdtype(output.dtype[field], np.floating):
                output[field] = np.nan
                
        # 1. Halo-scale: rvir (always computed)
        center_halo = halo_center_code(
            halo, box_physical_mpc, self.catalog.units_version
        )
        radius_rvir = halo_radius_code(
            halo, box_physical_mpc, "rvir",
            units_version=self.catalog.units_version,
        )
        candidate_mask_rvir = sphere_cell_mask(cells, center_halo, radius_rvir, self.include_boundary)
        selected_rvir = cells[candidate_mask_rvir]
        
        if selected_rvir.size > 0:
            fraction_rvir, lower_rvir, upper_rvir = self._fractional_overlap(
                selected_rvir, center_halo, radius_rvir
            )
            volume_rvir = np.exp2(-3.0 * selected_rvir["level"])
            mass_code_rvir = selected_rvir["rho"] * volume_rvir * fraction_rvir
            mass_msol_rvir = mass_code_rvir / units["Msol"]
            
            output["mgas_rvir"] = np.sum(mass_msol_rvir)
            
            temp_rvir = selected_rvir["P"] / selected_rvir["rho"] / units["K"]
            rho_hcc_rvir = selected_rvir["rho"] / units["H/cc"]
            cold_mask_rvir = temp_rvir < 1e4
            dense_mask_rvir = cold_mask_rvir & (rho_hcc_rvir > 5.0)
            
            output["mcold_rvir"] = np.sum(mass_msol_rvir[cold_mask_rvir])
            output["mdense_rvir"] = np.sum(mass_msol_rvir[dense_mask_rvir])
        else:
            output["mgas_rvir"] = np.nan
            output["mcold_rvir"] = np.nan
            output["mdense_rvir"] = np.nan
            
        # 2. Galaxy-scale: r*, r50, r90 (only if has stars)
        has_stars = (halo["nstar"] > 0) and np.isfinite(halo["r*"]) and (halo["r*"] > 0)
        if has_stars:
            center_star = np.array([halo["px*"], halo["py*"], halo["pz*"]])
            if self.catalog.units_version == "halomaker_units_v2":
                center_star = np.mod(center_star, 1.0)
            else:
                center_star = np.mod(center_star / box_physical_mpc, 1.0)
            
            radius_rstar = halo_radius_code(
                halo, box_physical_mpc, "r*",
                units_version=self.catalog.units_version,
            )
            r50_val = float(halo["r50"])
            r90_val = float(halo["r90"])
            if self.catalog.units_version == "halomaker_units_v2":
                radius_r50 = r50_val if (np.isfinite(r50_val) and r50_val > 0) else 0.0
                radius_r90 = r90_val if (np.isfinite(r90_val) and r90_val > 0) else 0.0
            else:
                radius_r50 = r50_val / box_physical_mpc if (np.isfinite(r50_val) and r50_val > 0) else 0.0
                radius_r90 = r90_val / box_physical_mpc if (np.isfinite(r90_val) and r90_val > 0) else 0.0
            
            candidate_mask_rstar = sphere_cell_mask(cells, center_star, radius_rstar, self.include_boundary)
            selected_rstar = cells[candidate_mask_rstar]
            
            if selected_rstar.size > 0:
                fraction_rstar, lower_rstar, upper_rstar = self._fractional_overlap(
                    selected_rstar, center_star, radius_rstar
                )
                volume_rstar = np.exp2(-3.0 * selected_rstar["level"])
                mass_code_rstar = selected_rstar["rho"] * volume_rstar * fraction_rstar
                mass_msol_rstar = mass_code_rstar / units["Msol"]
                mass_sum_rstar = np.sum(mass_msol_rstar)
                
                # Diagnostics (on r*)
                output["ncells"] = selected_rstar.size
                center_mask_rstar = sphere_cell_mask(selected_rstar, center_star, radius_rstar, include_boundary=False)
                output["ncells_center"] = np.count_nonzero(center_mask_rstar)
                output["ncells_overlap"] = np.count_nonzero(fraction_rstar > 0)
                output["radius_to_min_cell"] = radius_rstar / np.min(np.exp2(-selected_rstar["level"]))
                
                mgas_overlap_lower_rstar = np.sum(selected_rstar["rho"] * volume_rstar * lower_rstar) / units["Msol"]
                mgas_overlap_upper_rstar = np.sum(selected_rstar["rho"] * volume_rstar * upper_rstar) / units["Msol"]
                if mass_sum_rstar > 0:
                    output["overlap_relative_bound"] = (mgas_overlap_upper_rstar - mgas_overlap_lower_rstar) / mass_sum_rstar
                    output["overlap_resolved"] = (output["overlap_relative_bound"] <= self.overlap_tolerance)
                else:
                    output["overlap_resolved"] = False
                    
                # Cold/dense masks for selected_rstar
                temp_rstar = selected_rstar["P"] / selected_rstar["rho"] / units["K"]
                rho_hcc_rstar = selected_rstar["rho"] / units["H/cc"]
                cold_mask_rstar = temp_rstar < 1e4
                dense_mask_rstar = cold_mask_rstar & (rho_hcc_rstar > 5.0)
                
                # Metallicity & chemistry (on r*)
                if mass_sum_rstar > 0:
                    if "metal" in selected_rstar.dtype.names:
                        output["metal_gas"] = np.average(selected_rstar["metal"], weights=mass_code_rstar)
                    
                    # Dynamically discovered chem element fractions
                    for elem in ["H", "O", "Fe", "Mg", "C", "N", "Si", "S", "D"]:
                        if elem in selected_rstar.dtype.names:
                            output[f"{elem}_gas"] = np.average(selected_rstar[elem], weights=mass_code_rstar)
                            
                # Kinematics relative to gas bulk velocity in r*
                if mass_sum_rstar > 0:
                    # Physical delta (Mpc)
                    dx = periodic_delta(selected_rstar["x"], center_star[0]) * box_physical_mpc
                    dy = periodic_delta(selected_rstar["y"], center_star[1]) * box_physical_mpc
                    dz = periodic_delta(selected_rstar["z"], center_star[2]) * box_physical_mpc
                    
                    # Velocities (km/s)
                    vx_kms = selected_rstar["vx"] / units["km/s"]
                    vy_kms = selected_rstar["vy"] / units["km/s"]
                    vz_kms = selected_rstar["vz"] / units["km/s"]
                    
                    bulk_vx = np.average(vx_kms, weights=mass_code_rstar)
                    bulk_vy = np.average(vy_kms, weights=mass_code_rstar)
                    bulk_vz = np.average(vz_kms, weights=mass_code_rstar)
                    
                    dvx = vx_kms - bulk_vx
                    dvy = vy_kms - bulk_vy
                    dvz = vz_kms - bulk_vz
                    
                    # Gas angular momentum
                    Lx_val = np.sum(mass_msol_rstar * (dy * dvz - dz * dvy))
                    Ly_val = np.sum(mass_msol_rstar * (dz * dvx - dx * dvz))
                    Lz_val = np.sum(mass_msol_rstar * (dx * dvy - dy * dvx))
                    Lvec = np.array([Lx_val, Ly_val, Lz_val])
                    
                    output["Lx_gas"] = Lx_val
                    output["Ly_gas"] = Ly_val
                    output["Lz_gas"] = Lz_val
                    
                    L_norm = np.linalg.norm(Lvec)
                else:
                    Lvec = np.zeros(3)
                    L_norm = 0.0
                    dx = dy = dz = dvx = dvy = dvz = None
                    
                # Loop over apertures to assign masses and compute kinematics
                for ap_name, radius in [("r*", radius_rstar), ("r50", radius_r50), ("r90", radius_r90)]:
                    suffix = "" if ap_name == "r*" else f"_{ap_name}"
                    if radius <= 0:
                        output[f"mgas{suffix}"] = 0.0
                        output[f"mcold{suffix}"] = 0.0
                        output[f"mdense{suffix}"] = 0.0
                        continue
                    fraction, lower, upper = self._fractional_overlap(
                        selected_rstar, center_star, radius
                    )
                    mass_code = selected_rstar["rho"] * volume_rstar * fraction
                    mass_msol = mass_code / units["Msol"]
                    mass_sum = np.sum(mass_msol)
                    
                    output[f"mgas{suffix}"] = mass_sum
                    output[f"mcold{suffix}"] = np.sum(mass_msol[cold_mask_rstar])
                    output[f"mdense{suffix}"] = np.sum(mass_msol[dense_mask_rstar])
                    
                    if mass_sum > 0 and np.count_nonzero(mass_code > 0) >= 2 and dvx is not None:
                        vrot, sig3d, sigcyl = self._gas_kinematics_1d(
                            mass_code, dx, dy, dz, dvx, dvy, dvz, Lvec
                        )
                        output[f"vrot_gas{suffix}"] = vrot
                        output[f"sig3d_gas{suffix}"] = sig3d
                        output[f"sigcyl_gas{suffix}"] = sigcyl
            else:
                # selected_rstar is empty
                output["ncells"] = 0
                output["ncells_center"] = 0
                output["ncells_overlap"] = 0
                output["overlap_resolved"] = False
        
        self._apply_tier2(output, halo, cells, dms, stars, units)
        return output

    def process_root(self, root_id, read_grav=False, nthread=1):
        rows = self.catalog.descendants(root_id)
        halos = self.catalog.halos[rows]
        root = self.catalog.halos[self.catalog.row_for_id(root_id)]
        center = halo_center_code(
            root, self.catalog.box_physical_mpc, self.catalog.units_version
        )
        
        required_radius = 0.0
        # TASK-07: keep the envelope bounded by the catalog halo max radius `r`.
        for halo in halos:
            child_center_halo = halo_center_code(
                halo, self.catalog.box_physical_mpc, self.catalog.units_version
            )
            dist_halo = np.linalg.norm(
                child_center_halo - center
                - np.rint(child_center_halo - center)
            )
            child_radius = halo_radius_code(
                halo, self.catalog.box_physical_mpc, self.radius_field,
                units_version=self.catalog.units_version,
            )
            required_radius = max(required_radius, dist_halo + child_radius)

            has_stars = (
                (halo["nstar"] > 0)
                and np.isfinite(halo["r*"])
                and (halo["r*"] > 0)
            )
            if has_stars:
                child_center_star = np.array([halo["px*"], halo["py*"], halo["pz*"]])
                if self.catalog.units_version == "halomaker_units_v2":
                    child_center_star = np.mod(child_center_star, 1.0)
                else:
                    child_center_star = np.mod(
                        child_center_star / self.catalog.box_physical_mpc, 1.0
                    )
                dist_star = np.linalg.norm(
                    child_center_star - center
                    - np.rint(child_center_star - center)
                )
                radius_rstar = halo_radius_code(
                    halo, self.catalog.box_physical_mpc, "r*",
                    units_version=self.catalog.units_version,
                )
                required_radius = max(required_radius, dist_star + radius_rstar)

        root_radius = halo_radius_code(
            root, self.catalog.box_physical_mpc, self.radius_field,
            units_version=self.catalog.units_version,
        )
        envelope_radius = max(
            root_radius * self.padding,
            required_radius + self.reader.maximum_cell_half_diagonal,
        )
        boxes = periodic_boxes(center, envelope_radius)
        
        # Discover hydro fields dynamically
        all_fields, chem_elements = self.reader.hydro_fields()
        base_fields = ["rho", "P", "metal", "vx", "vy", "vz"]
        fields = [f for f in base_fields + chem_elements if f in all_fields]

        started = time.perf_counter()
        read = self.reader.read_boxes(
            boxes,
            fields=fields,
            read_grav=read_grav,
            nthread=nthread,
            read_particles=True,
        )
        read_seconds = time.perf_counter() - started

        started = time.perf_counter()
        results = np.empty(halos.size, dtype=SUMMARY_DTYPE)
        for index, halo in enumerate(halos):
            results[index] = self._summarize(
                halo, root_id, read.cells, read.dms, read.stars, read.units
            )
        compute_seconds = time.perf_counter() - started
        metrics = RootMetrics(
            root_id=root_id,
            descendant_count=halos.size - 1,
            cpu_count=read.cpus.size,
            cells_read=read.cells.size,
            envelope_radius_code=envelope_radius,
            read_seconds=read_seconds,
            compute_seconds=compute_seconds,
        )
        return rows, results, metrics

    def _empty_summary(self):
        full = np.zeros(self.catalog.halos.size, dtype=SUMMARY_DTYPE)
        full["id"] = self.catalog.halos["id"]
        for field in SUMMARY_DTYPE.names[1:]:
            dtype = full.dtype[field]
            if np.issubdtype(dtype, np.floating):
                full[field] = np.nan
            else:
                full[field] = 0
        return full

    def _write_common_header(self, header, run_mode, requested_root_ids):
        header.attrs["schema_version"] = SCHEMA_VERSION
        header.attrs["units_version"] = UNITS_VERSION
        header.attrs["source_catalog_units_version"] = self.catalog.units_version
        header.attrs["source_catalog"] = str(self.catalog_path)
        header.attrs["radius_field"] = self.radius_field
        header.attrs["root_padding"] = self.padding
        header.attrs["cell_boundary_overlap"] = self.include_boundary
        header.attrs["overlap_method"] = "adaptive_r_fraction"
        header.attrs["overlap_threshold"] = self.overlap_threshold
        header.attrs["overlap_depth"] = self.overlap_depth
        header.attrs["overlap_tolerance"] = self.overlap_tolerance
        header.attrs["ownership"] = "inclusive"
        header.attrs["run_mode"] = run_mode
        header.attrs["requested_root_ids"] = np.asarray(
            requested_root_ids, dtype=np.int32
        )
        header.attrs["requested_root_count"] = len(requested_root_ids)
        header.attrs["completed_root_count"] = 0
        header.attrs["total_cells_read"] = 0
        header.attrs["total_read_seconds"] = 0.0
        header.attrs["total_compute_seconds"] = 0.0

    def _create_output(self, output_path, root_ids, run_mode):
        with h5py.File(output_path, "w") as hdf:
            header = hdf.create_group("header")
            self._write_common_header(header, run_mode, root_ids)
            gas = hdf.create_group("gas")
            gas.attrs["row_alignment"] = "/catalog/halo"
            gas.attrs["completion_authority"] = "/gas/root_metrics"
            gas.create_dataset("processed", data=np.zeros(self.catalog.halos.size, dtype=bool))
            summary_ds = gas.create_dataset(
                "summary", data=self._empty_summary(), compression="lzf"
            )
            summary_ds.attrs["field_units"] = json.dumps(GAS_FIELD_UNITS, sort_keys=True)
            summary_ds.attrs["units_version"] = UNITS_VERSION
            gas.create_dataset(
                "root_metrics",
                shape=(0,),
                maxshape=(None,),
                dtype=ROOT_METRICS_DTYPE,
                chunks=True,
            )

    def _recompute_header_aggregates(self, hdf):
        rm = hdf["gas/root_metrics"][:]
        header = hdf["header"]
        n = int(rm.shape[0])
        header.attrs["completed_root_count"] = n
        header.attrs["total_cells_read"] = int(np.sum(rm["cells_read"])) if n else 0
        header.attrs["total_read_seconds"] = float(np.sum(rm["read_seconds"])) if n else 0.0
        header.attrs["total_compute_seconds"] = (
            float(np.sum(rm["compute_seconds"])) if n else 0.0
        )

    def _require_compatible_output(self, hdf):
        if "header" not in hdf or "gas" not in hdf:
            raise ValueError("Output file is missing required header/gas groups")
        version = int(hdf["header"].attrs.get("schema_version", -1))
        if version != SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported GasMaker schema_version={version}; expected {SCHEMA_VERSION}"
            )
        gas = hdf["gas"]
        for name in ("summary", "processed", "root_metrics"):
            if name not in gas:
                raise ValueError(f"Output file is missing /gas/{name}")
        if gas["summary"].dtype != SUMMARY_DTYPE:
            raise ValueError("Output /gas/summary dtype does not match current SUMMARY_DTYPE")
        if gas["processed"].shape != (self.catalog.halos.size,):
            raise ValueError("Output /gas/processed shape does not match catalog rows")
        if gas["root_metrics"].dtype != ROOT_METRICS_DTYPE:
            raise ValueError("Output /gas/root_metrics dtype does not match current schema")

    def initialize_output(self, output_path, root_ids, overwrite=False, run_mode="multi_root"):
        output_path = Path(output_path)
        if output_path.exists() and overwrite:
            output_path.unlink()
        if not output_path.exists():
            output_path.parent.mkdir(parents=True, exist_ok=True)
            self._create_output(output_path, root_ids, run_mode)
        else:
            with h5py.File(output_path, "r+") as hdf:
                self._require_compatible_output(hdf)
                self._recompute_header_aggregates(hdf)
                hdf.flush()

    def completed_roots(self, output_path):
        with h5py.File(output_path, "r") as hdf:
            self._require_compatible_output(hdf)
            metrics = hdf["gas/root_metrics"][:]
        return set(int(root_id) for root_id in metrics["root_id"])

    def write_root_result(
        self,
        output_path,
        rows,
        results,
        metrics,
        simulate_crash_after_summary=False,
    ):
        with h5py.File(output_path, "r+") as hdf:
            self._require_compatible_output(hdf)
            gas = hdf["gas"]
            root_id = int(metrics.root_id)
            completed = set(int(item) for item in gas["root_metrics"]["root_id"])
            if root_id in completed:
                return False

            gas["summary"][rows] = results
            gas["processed"][rows] = True
            hdf.flush()
            if simulate_crash_after_summary:
                raise RuntimeError(
                    f"Simulated crash after summary write for root_id={root_id}"
                )

            root_metrics = gas["root_metrics"]
            size = root_metrics.shape[0]
            root_metrics.resize((size + 1,))
            root_metrics[size] = metrics.as_record()
            self._recompute_header_aggregates(hdf)
            hdf.flush()
        return True

    def process_roots(
        self,
        root_ids,
        output_path,
        overwrite=False,
        read_grav=False,
        nthread=1,
        stop_after_roots=None,
        simulate_crash_after_summary_root=None,
    ):
        root_ids = [int(root_id) for root_id in root_ids]
        self.initialize_output(output_path, root_ids, overwrite=overwrite)
        completed = self.completed_roots(output_path)
        processed_now = []
        skipped = []
        for root_id in root_ids:
            if root_id in completed:
                skipped.append(root_id)
                continue
            rows, results, metrics = self.process_root(
                root_id,
                read_grav=read_grav,
                nthread=nthread,
            )
            self.write_root_result(
                output_path,
                rows,
                results,
                metrics,
                simulate_crash_after_summary=(
                    simulate_crash_after_summary_root == root_id
                ),
            )
            completed.add(root_id)
            processed_now.append(root_id)
            if stop_after_roots is not None and len(processed_now) >= stop_after_roots:
                break
        return {
            "requested": root_ids,
            "processed": processed_now,
            "skipped": skipped,
            "remaining": [
                root_id
                for root_id in root_ids
                if root_id not in completed
            ],
        }

    def write_result(self, output_path, rows, results, metrics):
        self.initialize_output(
            output_path,
            [metrics.root_id],
            overwrite=True,
            run_mode="single_root",
        )
        self.write_root_result(output_path, rows, results, metrics)
