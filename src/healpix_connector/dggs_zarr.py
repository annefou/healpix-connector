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
