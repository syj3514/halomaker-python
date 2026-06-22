from pathlib import Path

import numpy as np

from gasmaker.catalog import HaloCatalog
from gasmaker.geometry import (
    halo_center_code,
    periodic_delta,
    periodic_boxes,
    sphere_cell_mask,
)
from gasmaker.overlap import fractional_sphere_cell_overlap


CATALOG = (
    Path(__file__).parents[1]
    / "nh2_00060_stellar_chemistry_run"
    / "tree_bricks00060.h5"
)


def main():
    catalog = HaloCatalog.read(CATALOG)
    rows = catalog.descendants(3)
    halos = catalog.halos[rows]
    assert halos[0]["id"] == 3
    assert np.all(
        (halos["id"] == 3)
        | (halos["hosthalo"] == 3)
        | np.isin(halos["hostsub"], halos["id"])
    )
    assert np.count_nonzero(halos["level"] == 2) == 21
    assert np.count_nonzero(halos["level"] == 3) == 6
    assert halos.size == 28

    np.testing.assert_allclose(
        catalog.box_physical_mpc,
        catalog.box_comoving_mpc * catalog.aexp,
    )
    np.testing.assert_allclose(catalog.box_physical_mpc, 19.947893436430064)
    center = halo_center_code(halos[0], catalog.box_physical_mpc)
    recovered = (center - 0.5) * catalog.box_physical_mpc
    original = np.array(
        [halos[0]["px"], halos[0]["py"], halos[0]["pz"]]
    )
    np.testing.assert_allclose(
        periodic_delta(recovered, original, catalog.box_physical_mpc),
        0.0,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        periodic_delta(np.array([0.99, 0.01]), 0.0),
        [-0.01, 0.01],
    )
    split = periodic_boxes(np.array([0.5, 0.01, 0.99]), 0.02)
    assert len(split) == 4
    for box in split:
        assert np.all(box >= 0)
        assert np.all(box <= 1)
        assert np.all(box[:, 0] <= box[:, 1])

    cells = np.zeros(
        2,
        dtype=[("x", "f8"), ("y", "f8"), ("z", "f8"), ("level", "i4")],
    )
    cells["x"] = [0.99, 0.5]
    cells["level"] = 10
    mask = sphere_cell_mask(cells, np.zeros(3), 0.02, False)
    np.testing.assert_array_equal(mask, [True, False])

    unit_cell = np.zeros(
        1,
        dtype=[("x", "f8"), ("y", "f8"), ("z", "f8"), ("level", "i4")],
    )
    unit_cell["x"] = unit_cell["y"] = unit_cell["z"] = 0.5
    fraction, lower, upper = fractional_sphere_cell_overlap(
        unit_cell, np.full(3, 0.5), 0.25, max_depth=6
    )
    analytic = 4.0 / 3.0 * np.pi * 0.25 ** 3
    assert lower[0] <= analytic <= upper[0]
    assert abs(fraction[0] - analytic) < 5.0e-4

    inside, inside_lower, inside_upper = fractional_sphere_cell_overlap(
        unit_cell, np.full(3, 0.5), 1.0, max_depth=0
    )
    np.testing.assert_array_equal(inside, [1.0])
    np.testing.assert_array_equal(inside_lower, [1.0])
    np.testing.assert_array_equal(inside_upper, [1.0])
    print("GASMAKER_UNIT_TEST=PASS")
    print(f"ROOT3_ROWS={rows.size}")
    print(f"ROOT3_DESCENDANTS={rows.size - 1}")


if __name__ == "__main__":
    main()
