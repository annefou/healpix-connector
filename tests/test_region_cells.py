import numpy as np
import pytest

from healpix_connector.cells import assign_cells, footprint_cells
from healpix_connector.region import Region


def test_bbox_contains_and_ccw_wkt():
    r = Region.from_bbox(-4, 40, -3, 41)
    assert r.contains([-3.5, -5.0], [40.5, 40.5]).tolist() == [True, False]
    wkt = r.query_wkt()
    assert wkt.startswith("POLYGON((-4.000000 40.000000,-3.000000 40.000000")


def test_polygon_is_made_counter_clockwise_and_tested_exactly():
    cw = [(0, 0), (0, 2), (2, 2), (2, 0)]  # clockwise
    r = Region.from_polygon(cw)
    assert r.query_wkt().startswith("POLYGON((2.000000 0.000000")  # reversed to CCW
    # L-shape: the notch must be outside.
    L = Region.from_polygon([(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)])
    assert L.contains([0.5, 1.5, 1.5], [0.5, 0.5, 1.5]).tolist() == [True, True, False]


def test_cell_region_query_covers_cells_and_filters_exactly():
    cells = assign_cells([-3.7, -3.4], [40.4, 40.4], 8)
    r = Region.from_cells(cells, 8)
    x0, y0, x1, y1 = r.query_bounds()
    assert x0 < -3.7 < x1 and y0 < 40.4 < y1
    assert r.contains([-3.7, -3.4, 5.0], [40.4, 40.4, 40.4]).tolist() == [True, True, False]


def test_footprint_without_uncertainty_is_own_cell_flagged_unknown():
    cells, known = footprint_cells(-3.7, 40.4, None, 8)
    assert known is False
    assert cells.tolist() == assign_cells([-3.7], [40.4], 8).tolist()


def test_footprint_grows_with_uncertainty_and_contains_own_cell():
    own = int(assign_cells([-3.7], [40.4], 8)[0])
    small, _ = footprint_cells(-3.7, 40.4, 100.0, 8)
    large, known = footprint_cells(-3.7, 40.4, 30_000.0, 8)  # wider than a ~25 km cell
    assert known and own in small and own in large
    assert large.size > small.size >= 1
    assert large.size >= 4


def test_empty_cell_region_is_refused():
    with pytest.raises(ValueError):
        Region.from_cells(np.array([], dtype=np.uint64), 8)
