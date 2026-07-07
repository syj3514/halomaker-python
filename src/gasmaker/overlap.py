import numpy as np

from .geometry import periodic_delta


OFFSETS = np.array(
    [
        (x, y, z)
        for x in (-1.0, 1.0)
        for y in (-1.0, 1.0)
        for z in (-1.0, 1.0)
    ],
    dtype=np.float64,
)


def _classify(offsets, half_size, radius_squared):
    absolute = np.abs(offsets)
    minimum = np.maximum(absolute - half_size[:, None], 0.0)
    maximum = absolute + half_size[:, None]
    fully_inside = np.sum(maximum * maximum, axis=1) <= radius_squared
    fully_outside = np.sum(minimum * minimum, axis=1) >= radius_squared
    return fully_inside, fully_outside


def fractional_sphere_cell_overlap(cells, center, radius, max_depth=4, threshold_ratio=0.1):
    """Estimate sphere/cube overlap with deterministic adaptive subdivision.

    Returns estimate, lower bound, and upper bound for every input cell.
    Bounds differ only for subcubes unresolved at ``max_depth`` or stopped early.
    """
    if max_depth < 0:
        raise ValueError("max_depth must be non-negative")
    count = cells.size
    estimate = np.zeros(count, dtype=np.float64)
    lower = np.zeros(count, dtype=np.float64)
    upper = np.zeros(count, dtype=np.float64)
    if count == 0:
        return estimate, lower, upper

    offsets = np.column_stack((
        periodic_delta(cells["x"], center[0]),
        periodic_delta(cells["y"], center[1]),
        periodic_delta(cells["z"], center[2]),
    ))
    half_size = 0.5 * np.exp2(-cells["level"].astype(np.float64))
    owners = np.arange(count, dtype=np.int64)
    weights = np.ones(count, dtype=np.float64)
    radius_squared = radius * radius

    for depth in range(max_depth + 1):
        inside, outside = _classify(
            offsets, half_size, radius_squared
        )
        if np.any(inside):
            np.add.at(estimate, owners[inside], weights[inside])
            np.add.at(lower, owners[inside], weights[inside])
            np.add.at(upper, owners[inside], weights[inside])

        unresolved = ~(inside | outside)
        if not np.any(unresolved):
            break
        offsets = offsets[unresolved]
        half_size = half_size[unresolved]
        owners = owners[unresolved]
        weights = weights[unresolved]

        if depth == max_depth:
            center_inside = (
                np.sum(offsets * offsets, axis=1) <= radius_squared
            )
            if np.any(center_inside):
                np.add.at(
                    estimate,
                    owners[center_inside],
                    weights[center_inside],
                )
            np.add.at(upper, owners, weights)
            break

        # Adaptive early stop
        if threshold_ratio is not None and threshold_ratio > 0:
            # Child's dx = child_half * 2.0 = half_size
            early_stop = half_size <= radius * threshold_ratio
            if np.any(early_stop):
                center_inside = (
                    np.sum(offsets[early_stop] * offsets[early_stop], axis=1) <= radius_squared
                )
                if np.any(center_inside):
                    np.add.at(
                        estimate,
                        owners[early_stop][center_inside],
                        weights[early_stop][center_inside],
                    )
                np.add.at(
                    upper,
                    owners[early_stop],
                    weights[early_stop],
                )
                keep = ~early_stop
                if not np.any(keep):
                    break
                offsets = offsets[keep]
                half_size = half_size[keep]
                owners = owners[keep]
                weights = weights[keep]

        child_half = half_size * 0.5
        offsets = (
            offsets[:, None, :]
            + OFFSETS[None, :, :] * child_half[:, None, None]
        ).reshape(-1, 3)
        half_size = np.repeat(child_half, 8)
        owners = np.repeat(owners, 8)
        weights = np.repeat(weights / 8.0, 8)

    return estimate, lower, upper

