"""GRID4EARTH DGGS-Zarr metadata, matching healpix-convert's output exactly.

Group attributes carry the zarr-conventions/dggs v1 declaration; a scalar
``crs`` variable carries the CF HEALPix grid mapping; the cell-id coordinate is
``cell_ids`` with the CF ``healpix_index`` standard name.
"""

from __future__ import annotations

DGGS_CONVENTION = {
    "uuid": "7b255807-140c-42ca-97f6-7a1cfecdbc38",
    "name": "dggs",
    "schema_url": "https://raw.githubusercontent.com/zarr-conventions/dggs/refs/tags/v1/schema.json",
    "spec_url": "https://github.com/zarr-conventions/dggs/blob/v1/README.md",
    "description": "Discrete Global Grid Systems convention for zarr",
}

WGS84 = {"name": "wgs84", "semi_major_axis": 6378137.0, "inverse_flattening": 298.257223563}


def dggs_attrs(depth: int) -> dict:
    """Group attributes declaring a NESTED WGS84 HEALPix grid at ``depth``."""
    return {
        "zarr_conventions": [DGGS_CONVENTION],
        "dggs": {
            "name": "healpix",
            "refinement_level": int(depth),
            "indexing_scheme": "nested",
            "ellipsoid": dict(WGS84),
            "spatial_dimension": "cells",
            "coordinate": "cell_ids",
            "compression": "none",
        },
    }


def cell_dataset(cell_ids, depth: int, variables: dict, attrs: dict, *, time_bounds=None):
    """A DGGS-Zarr dataset of ``variables`` over ``cell_ids`` (sorted on the way in).

    ``variables`` maps a name to ``(values, var_attrs)``; values are ``(cells,)``
    or, when ``time_bounds`` (``(T, 2)`` start and exclusive end) is given,
    ``(time, cells)``. The time axis follows CF: ``time`` is each step's start and
    its ``bounds`` attribute names ``time_bnds``. Every variable gets ``grid_mapping``.
    """
    import numpy as np
    import xarray as xr

    ids = np.asarray(cell_ids, dtype=np.uint64)
    order = np.argsort(ids, kind="stable")
    if np.any(np.diff(ids[order].astype(np.int64)) == 0):
        raise ValueError("cell_ids must be unique")
    data = {}
    for name, (v, va) in variables.items():
        v = np.asarray(v)
        if v.ndim == 2 and time_bounds is None:
            raise ValueError(f"{name} has a time axis but no time_bounds were given")
        data[name] = (("cells",) if v.ndim == 1 else ("time", "cells"), v[..., order],
                      {**va, "grid_mapping": "crs"})
    data["crs"] = ((), np.int8(0), cf_grid_mapping_attrs(depth))
    coords = {"cell_ids": ("cells", ids[order], {"standard_name": "healpix_index", "units": "1"})}
    if time_bounds is not None:
        tb = np.asarray(time_bounds, dtype="datetime64[ns]")
        coords["time"] = ("time", tb[:, 0], {"bounds": "time_bnds", "standard_name": "time"})
        data["time_bnds"] = (("time", "bnds"), tb)
    ds = xr.Dataset(data, coords=coords)
    ds.attrs.update(dggs_attrs(depth))
    ds.attrs.update(attrs)
    return ds


def concat_time(datasets):
    """Join datasets built by :func:`cell_dataset` (same cells) along time, in time order.

    Variables without a time axis (a per-cell flag, ``crs``) are taken from the first.
    """
    import numpy as np
    import xarray as xr

    first = datasets[0]["cell_ids"].values
    if any(not np.array_equal(d["cell_ids"].values, first) for d in datasets[1:]):
        raise ValueError("datasets cover different cells")
    out = xr.concat(datasets, dim="time", data_vars="minimal", coords="minimal",
                    compat="override", combine_attrs="override").sortby("time")
    if np.unique(out["time"].values).size != out.sizes["time"]:
        raise ValueError("the same time step appears twice")
    return out


def month_bounds(year: int, months) -> "np.ndarray":
    """``(len(months), 2)`` start and exclusive end of each month of ``year``."""
    import numpy as np

    starts = np.array([f"{year}-{m:02d}" for m in months], dtype="datetime64[M]")
    return np.stack([starts, starts + np.timedelta64(1, "M")], axis=1).astype("datetime64[s]")


def cf_grid_mapping_attrs(depth: int) -> dict:
    """Attributes of the scalar ``crs`` grid-mapping variable."""
    return {
        "grid_mapping_name": "healpix",
        "refinement_level": int(depth),
        "indexing_scheme": "nested",
        "reference_ellipsoid_name": "WGS84",
        "semi_major_axis": WGS84["semi_major_axis"],
        "inverse_flattening": WGS84["inverse_flattening"],
    }
