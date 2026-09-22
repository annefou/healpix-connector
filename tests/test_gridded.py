import numpy as np
import pytest

from healpix_connector.cells import assign_cells
from healpix_connector.sources.gridded import sample_cells_from_grid

RES = 0.25  # ERA5


def grid():
    lon = np.arange(-11, 5.01, RES)
    lat = np.arange(45, 33.99, -RES)
    return lon, lat


def test_takes_the_containing_box_value():
    lon, lat = grid()
    values = np.broadcast_to(lat[:, None], (lat.size, lon.size)).astype(float).copy()
    cells = assign_cells([-3.7, -3.7], [40.4, 43.0], 8)
    got = sample_cells_from_grid(values, lon, lat, cells, 8, RES)
    # The field equals latitude, so each cell gets its own latitude, to within half a box.
    assert abs(got.value[0] - 40.4) <= RES
    assert abs(got.value[1] - 43.0) <= RES


def test_support_is_the_source_box_and_is_declared_when_wider_than_the_cell():
    lon, lat = grid()
    values = np.ones((lat.size, lon.size))
    cells = assign_cells([-3.7], [40.4], 8)  # ~25 km cells, ERA5 box ~21-28 km
    got = sample_cells_from_grid(values, lon, lat, cells, 8, RES)
    assert got.support_m > 15_000
    assert "support is the source box" in got.method or "grid box containing" in got.method
    # At depth 10 (~6 km) the source box is far wider than the cell: it must say so.
    deep = sample_cells_from_grid(values, lon, lat, assign_cells([-3.7], [40.4], 10), 10, RES)
    assert "wider than the cell" in deep.method


def test_shape_mismatch_is_refused():
    lon, lat = grid()
    with pytest.raises(ValueError):
        sample_cells_from_grid(np.ones((3, 3)), lon, lat, assign_cells([-3.7], [40.4], 8), 8, RES)
