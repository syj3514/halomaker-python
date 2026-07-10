import numpy as np


def halo_center_code(halo, box_mpc, units_version="legacy"):
    center_mpc = np.array([halo["px"], halo["py"], halo["pz"]])
    if units_version == "halomaker_units_v2":
        return np.mod(center_mpc, 1.0)
    return np.mod(center_mpc / box_mpc, 1.0)


def halo_radius_code(halo, box_mpc, radius_field="r", padding=1.0, units_version="legacy",
                     allow_nonfinite=False):
    radius = float(halo[radius_field])
    if not np.isfinite(radius) or radius <= 0:
        # Non-virialized/degenerate halos carry rvir=NaN in the catalog
        # (HaloMaker det_vir_props_1b4 fix, dc71d71). For virial-radius apertures
        # propagate the NaN so downstream masks skip the aperture and the gas
        # fields stay NaN, instead of aborting the whole run.
        if allow_nonfinite:
            return float("nan")
        raise ValueError(
            f"Invalid {radius_field}={radius} for halo id={halo['id']}"
        )
    if units_version == "halomaker_units_v2":
        return radius * padding
    return radius * padding / box_mpc


def periodic_delta(values, center, boxlen=1.0):
    delta = values - center
    return delta - np.rint(delta / boxlen) * boxlen


def sphere_cell_mask(cells, center, radius, include_boundary=True):
    dx = periodic_delta(cells["x"], center[0])
    dy = periodic_delta(cells["y"], center[1])
    dz = periodic_delta(cells["z"], center[2])
    distance = np.sqrt(dx * dx + dy * dy + dz * dz)
    if include_boundary:
        half_diagonal = np.sqrt(3.0) * 0.5 * np.exp2(-cells["level"])
        return distance <= radius + half_diagonal
    return distance <= radius


class PeriodicSpatialIndex:
    """Periodic uniform-grid index returning source-order candidate rows."""

    def __init__(self, table, grid_size=256):
        self.table = table
        self.grid_size = int(grid_size)
        self.size = 0 if table is None else len(table)
        self._all = None
        self._keys = np.empty(0, dtype=np.uint32)
        self._starts = np.empty(0, dtype=np.int64)
        self._order = np.empty(0, dtype=np.int64)

        if self.size == 0:
            return

        self._all = np.arange(self.size, dtype=np.int64)
        g = self.grid_size

        key = np.floor(np.mod(table["x"], 1.0) * g).astype(np.uint32)
        np.minimum(key, g - 1, out=key)
        key *= np.uint32(g)

        component = np.floor(np.mod(table["y"], 1.0) * g).astype(np.uint32)
        np.minimum(component, g - 1, out=component)
        key += component
        key *= np.uint32(g)

        component = np.floor(np.mod(table["z"], 1.0) * g).astype(np.uint32)
        np.minimum(component, g - 1, out=component)
        key += component

        self._order = np.argsort(key, kind="stable")
        sorted_key = key[self._order]
        self._keys, self._starts = np.unique(sorted_key, return_index=True)

    def query(self, center, radius):
        """Return source-order row indices in a periodic AABB around center."""
        if self.size == 0:
            return self._order
        if not (np.isfinite(radius) and radius > 0):
            return np.empty(0, dtype=np.int64)
        if radius >= 0.5:
            return self._all

        g = self.grid_size
        span = int(np.ceil(radius * g)) + 1
        if 2 * span + 1 >= g:
            return self._all

        center = np.mod(np.asarray(center, dtype=np.float64), 1.0)
        cbin = np.floor(center * g).astype(np.int64)
        offsets = np.arange(-span, span + 1, dtype=np.int64)
        xs = np.mod(cbin[0] + offsets, g)
        ys = np.mod(cbin[1] + offsets, g)
        zs = np.mod(cbin[2] + offsets, g)

        keys = (
            (xs[:, None, None] * g + ys[None, :, None]) * g
            + zs[None, None, :]
        ).ravel()
        keys = np.unique(keys).astype(self._keys.dtype, copy=False)
        positions = np.searchsorted(self._keys, keys)
        in_range = positions < self._keys.size
        if not np.any(in_range):
            return np.empty(0, dtype=np.int64)
        keys = keys[in_range]
        positions = positions[in_range]
        matched = self._keys[positions] == keys
        if not np.any(matched):
            return np.empty(0, dtype=np.int64)

        positions = positions[matched]
        lengths = np.diff(np.append(self._starts, self._order.size))[positions]
        total = int(np.sum(lengths))
        if total == 0:
            return np.empty(0, dtype=np.int64)

        indices = np.empty(total, dtype=np.int64)
        cursor = 0
        for start, length in zip(self._starts[positions], lengths):
            stop = cursor + int(length)
            indices[cursor:stop] = self._order[start:start + int(length)]
            cursor = stop
        indices.sort()
        return indices


def periodic_boxes(center, radius):
    intervals = []
    for value in center:
        lower = value - radius
        upper = value + radius
        if lower < 0:
            intervals.append(((0.0, upper), (lower + 1.0, 1.0)))
        elif upper > 1:
            intervals.append(((lower, 1.0), (0.0, upper - 1.0)))
        else:
            intervals.append(((lower, upper),))

    boxes = []
    for xlim in intervals[0]:
        for ylim in intervals[1]:
            for zlim in intervals[2]:
                boxes.append(np.array((xlim, ylim, zlim), dtype=np.float64))
    return boxes
