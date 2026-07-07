import numpy as np


def halo_center_code(halo, box_mpc, units_version="legacy"):
    center_mpc = np.array([halo["px"], halo["py"], halo["pz"]])
    if units_version == "halomaker_units_v2":
        return np.mod(center_mpc, 1.0)
    return np.mod(center_mpc / box_mpc, 1.0)


def halo_radius_code(halo, box_mpc, radius_field="r", padding=1.0, units_version="legacy"):
    radius = float(halo[radius_field])
    if not np.isfinite(radius) or radius <= 0:
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
