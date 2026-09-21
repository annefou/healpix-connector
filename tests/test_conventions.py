import numpy as np
import pytest

from healpix_connector.conventions import ELLIPSOID, cell_size_m


@pytest.mark.parametrize(
    "depth, km",
    [(7, 50.93), (8, 25.47), (10, 6.37), (12, 1.59), (13, 0.80)],
)
def test_cell_sizes_match_the_design_note(depth, km):
    assert cell_size_m(depth) / 1000 == pytest.approx(km, abs=0.01)


def test_depth_out_of_range_is_refused():
    with pytest.raises(ValueError):
        cell_size_m(30)


def test_healpix_geo_is_used_on_wgs84_nested():
    # The Canberra point used throughout the design note, at depth 8.
    from healpix_geo import nested

    cells = nested.lonlat_to_healpix(
        np.array([149.13]), np.array([-35.28]), np.uint8(8), ellipsoid=ELLIPSOID
    )
    assert cells.shape == (1,)
    # NESTED: the parent at depth 7 is the cell id shifted right by two bits.
    parent = nested.lonlat_to_healpix(
        np.array([149.13]), np.array([-35.28]), np.uint8(7), ellipsoid=ELLIPSOID
    )
    assert int(cells[0]) >> 2 == int(parent[0])
