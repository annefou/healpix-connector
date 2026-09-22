import json

import numpy as np

from healpix_connector.binning import bin_to_cells
from healpix_connector.sources import chelsa

RES = 1 / 120


def small_dataset():
    lon = np.arange(-3.9958, -3.0, RES)
    lat = np.arange(40.9958, 40.0, -RES)
    vals = np.full((lat.size, lon.size), 14.5)
    return chelsa.to_dataset(bin_to_cells(vals, lon, lat, 8, RES), 1)


def test_dataset_declares_grid4earth_dggs_convention():
    ds = small_dataset()
    assert ds.attrs["zarr_conventions"][0]["name"] == "dggs"
    assert ds.attrs["dggs"] == {
        "name": "healpix", "refinement_level": 8, "indexing_scheme": "nested",
        "ellipsoid": {"name": "wgs84", "semi_major_axis": 6378137.0, "inverse_flattening": 298.257223563},
        "spatial_dimension": "cells", "coordinate": "cell_ids", "compression": "none",
    }
    assert ds["crs"].attrs["grid_mapping_name"] == "healpix"
    assert ds["cell_ids"].attrs["standard_name"] == "healpix_index"
    assert ds["bio1"].attrs["grid_mapping"] == "crs"


def test_dataset_carries_layer3_fields():
    a = small_dataset().attrs
    assert a["support_depth"] == 8
    assert a["valid_time"] == "climatology"
    assert a["climatology_bounds"] == ["1981-01-01", "2010-12-31"]
    assert a["source_doi"] == "10.16904/envidat.228"
    assert a["license"] == "CC0-1.0"
    assert "binning" in a["resampling_method"]


def test_write_roundtrip_and_stac_item(tmp_path):
    import xarray as xr

    ds = small_dataset()
    z = chelsa.write(ds, tmp_path, (-4, 40, -3, 41))
    back = xr.open_zarr(z, consolidated=False)
    assert back.attrs["dggs"]["refinement_level"] == 8
    np.testing.assert_allclose(back["bio1"].values, 14.5)
    item = json.loads(next(tmp_path.glob("*.stac.json")).read_text())
    assert item["properties"]["start_datetime"].startswith("1981")
    assert item["properties"]["sci:doi"] == "10.16904/envidat.228"
    assert item["assets"]["data"]["href"].endswith(".zarr")


def test_all_model_variables_have_verified_units():
    # Units as documented in CHELSA's file specification, except bio4: its
    # values are 100x the standard deviation (checked against ERA5).
    assert chelsa.BIO[4][1] == "0.01 degC" and chelsa.BIO[12][1] == "kg m-2"
    assert "100 x standard deviation" in chelsa.BIO[4][2]
    assert set(chelsa.BIO) == {1, 4, 5, 6, 12, 15}
    ds = chelsa.to_dataset(bin_to_cells(
        np.full((60, 60), 800.0), np.arange(-3.9958, -3.5, RES), np.arange(40.9958, 40.5, -RES), 8, RES), 12)
    assert ds["bio12"].attrs["units"] == "kg m-2"
