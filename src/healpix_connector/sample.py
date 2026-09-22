"""Attach environmental values to cells and records, with their support,
valid time and uncertainty components.

A value is never returned bare. Every sample states:

- **support**: the depth it was taken at, and whether it was aggregated from a
  finer dataset or inherited from a coarser one;
- **valid time**: whether the record's date falls inside the dataset's period
  (a climatology's reference period, or a time step), or is unknown;
- **uncertainty**, kept as separate components: the spread the source itself
  reports inside the cell, and the spread across the cells the record's
  positional uncertainty covers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from healpix_connector import __version__


@dataclass(frozen=True)
class Source:
    """A converted dataset: its grid, its period and what its numbers mean."""

    dataset: object  # xarray.Dataset
    depth: int
    variable: str
    title: str
    doi: str | None
    license: str | None
    valid_time: str | None
    period: tuple[str, str] | None
    resampling_method: str | None

    @property
    def cell_ids(self) -> np.ndarray:
        return np.asarray(self.dataset["cell_ids"].values, dtype=np.uint64)


def open_source(path, variable: str | None = None) -> Source:
    """Open a DGGS-Zarr dataset and read the metadata a sample must carry."""
    import xarray as xr

    ds = xr.open_zarr(path, consolidated=False)
    a = ds.attrs
    dggs = a.get("dggs") or {}
    if dggs.get("name") != "healpix" or dggs.get("indexing_scheme") != "nested":
        raise ValueError(f"{path} is not a NESTED HEALPix dataset (dggs={dggs})")
    if variable is None:
        candidates = [v for v in ds.data_vars
                      if not v.endswith(("_std", "_min", "_max")) and v not in ("crs", "pixel_count", "coverage")]
        if len(candidates) != 1:
            raise ValueError(f"specify variable; candidates: {candidates}")
        variable = candidates[0]
    period = a.get("climatology_bounds") or a.get("time_bounds")
    return Source(
        dataset=ds, depth=int(dggs["refinement_level"]), variable=variable,
        title=a.get("title", str(path)), doi=a.get("source_doi"), license=a.get("license"),
        valid_time=a.get("valid_time"), period=tuple(period) if period else None,
        resampling_method=a.get("resampling_method"),
    )


def _lookup(src: Source, cells: np.ndarray, name: str) -> np.ndarray:
    """Values of variable ``name`` for ``cells`` (NaN where the source has no cell)."""
    if name not in src.dataset:
        return np.full(cells.shape, np.nan)
    ids = src.cell_ids
    values = np.asarray(src.dataset[name].values, dtype=np.float64)
    pos = np.searchsorted(ids, cells)
    pos_clipped = np.clip(pos, 0, ids.size - 1)
    hit = ids[pos_clipped] == cells
    out = np.full(cells.shape, np.nan)
    out[hit] = values[pos_clipped[hit]]
    return out


def sample_cells(src: Source, cells, depth: int) -> dict:
    """Sample a source at ``cells`` given at ``depth``.

    Depths are reconciled explicitly: a coarser source is **inherited** (the
    parent's value, no new information), a finer source is **aggregated** (the
    mean of the children present). Both are reported in ``support``.
    """
    from healpix_geo import nested

    cells = np.asarray(cells, dtype=np.uint64)
    if depth == src.depth:
        value = _lookup(src, cells, src.variable)
        source_std = _lookup(src, cells, f"{src.variable}_std")
        how = "direct"
    elif depth > src.depth:
        parents = nested.zoom_to(cells, np.uint8(depth), np.uint8(src.depth))
        value = _lookup(src, parents, src.variable)
        source_std = _lookup(src, parents, f"{src.variable}_std")
        how = f"inherited from depth {src.depth}"
    else:
        children_per_cell = 4 ** (src.depth - depth)
        first = cells.astype(np.uint64) * np.uint64(children_per_cell)
        offsets = np.arange(children_per_cell, dtype=np.uint64)
        kids = (first[:, None] + offsets[None, :]).ravel()
        vals = _lookup(src, kids, src.variable).reshape(cells.size, children_per_cell)
        with np.errstate(invalid="ignore"):
            value = np.nanmean(vals, axis=1)
            source_std = np.nanstd(vals, axis=1)
        how = f"aggregated from depth {src.depth}"
    return {
        "value": value, "source_std": source_std,
        "support": {"depth": int(depth), "how": how, "source_depth": src.depth},
    }


def time_status(date: str | int | None, src: Source) -> str:
    """Where a record's date sits relative to the source's period."""
    if src.period is None:
        return "source period unknown"
    if date is None or (isinstance(date, float) and np.isnan(date)):
        return "record date unknown"
    year = int(str(date)[:4])
    start, end = int(str(src.period[0])[:4]), int(str(src.period[1])[:4])
    inside = start <= year <= end
    if src.valid_time == "climatology":
        return (f"inside the {start}-{end} climatology" if inside
                else f"outside the {start}-{end} climatology")
    return f"inside {start}-{end}" if inside else f"outside {start}-{end}"


def sample_records(records: list[dict], src: Source, depth: int) -> list[dict]:
    """Attach a source's values to records that already carry cells.

    Records must come from a connector's ``to_cells``: each needs ``cell`` and,
    for the positional component, ``footprint_cells`` and ``footprint_known``.
    """
    if not records:
        return []
    own = np.array([r["cell"] for r in records], dtype=np.uint64)
    at_own = sample_cells(src, own, depth)
    out = []
    for i, r in enumerate(records):
        fp = np.asarray(r.get("footprint_cells", [r["cell"]]), dtype=np.uint64)
        fp_values = sample_cells(src, fp, depth)["value"] if fp.size > 1 else np.array([at_own["value"][i]])
        finite = fp_values[np.isfinite(fp_values)]
        out.append({
            **r,
            "variable": src.variable,
            "value": float(at_own["value"][i]),
            "support": at_own["support"],
            "source_std": float(at_own["source_std"][i]),
            "positional_spread": float(finite.std()) if finite.size > 1 else 0.0,
            "positional_min": float(finite.min()) if finite.size else float("nan"),
            "positional_max": float(finite.max()) if finite.size else float("nan"),
            "positional_known": bool(r.get("footprint_known", False)),
            "time_status": time_status(r.get("eventDate") or r.get("year"), src),
        })
    return out


def provenance(src: Source, depth: int) -> dict:
    """What a result must carry to be citable and reproducible."""
    return {
        "dataset": src.title, "dataset_doi": src.doi, "license": src.license,
        "source_depth": src.depth, "sampled_at_depth": int(depth),
        "valid_time": src.valid_time, "period": list(src.period) if src.period else None,
        "resampling_method": src.resampling_method,
        "producer": f"healpix-connector {__version__}",
    }
