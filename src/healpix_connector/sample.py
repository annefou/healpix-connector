"""Attach environmental values to cells and records, with their support,
valid time and uncertainty components.

A value is never returned bare. Every sample states:

- **support**: the depth it was taken at, and whether it was aggregated from a
  finer dataset or inherited from a coarser one (and, for a source coarser than
  its cells, the width of the source box, ``support_m``);
- **valid time**: whether the record's date falls inside the dataset's period
  (a climatology's reference period) or, for a source with a time axis, which
  time step it falls in;
- **uncertainty**, kept as separate components: the spread the source itself
  reports inside the cell, and the spread across the cells the record's
  positional uncertainty covers.

Every source is a GRID4EARTH DGGS-Zarr dataset: one written to disk
(:func:`open_source`) or one built in memory by a live source's ``to_dataset``
(:func:`from_dataset`). A time axis follows CF: a ``time`` coordinate whose
``bounds`` attribute names a ``(time, 2)`` variable of start and (exclusive) end.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from healpix_connector import __version__

#: Per-variable companions written beside a variable; never the variable sampled.
COMPANION_SUFFIXES = ("_std", "_min", "_max", "_valid_share")
COMPANION_VARIABLES = ("crs", "pixel_count", "coverage", "filled_by_fallback")


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
    #: width of the source box when it is wider than a cell (a coarse grid), in metres
    support_m: float | None = None
    #: (T, 2) datetime64[s] start and exclusive end of each time step; None = no time axis
    time_bounds: np.ndarray | None = field(default=None, compare=False)

    @property
    def has_time(self) -> bool:
        return self.time_bounds is not None

    @cached_property
    def _index(self) -> tuple[np.ndarray, np.ndarray]:
        ids = np.asarray(self.dataset["cell_ids"].values, dtype=np.uint64)
        order = np.argsort(ids, kind="stable")
        return ids[order], order

    @property
    def cell_ids(self) -> np.ndarray:
        """The source's cells, sorted."""
        return self._index[0]

    def _array(self, name: str) -> np.ndarray | None:
        """``name`` as (cells,) or (time, cells), aligned on sorted cells; read once, then cached."""
        cache = self.__dict__.setdefault("_arrays", {})
        if name not in cache:
            if name not in self.dataset:
                cache[name] = None
            else:
                da = self.dataset[name]
                if "time" in da.dims:
                    da = da.transpose("time", ...)
                cache[name] = np.asarray(da.values, dtype=np.float64)[..., self._index[1]]
        return cache[name]

    def lookup(self, name: str, cells, steps=None) -> np.ndarray:
        """Values of ``name`` at ``cells`` (NaN where the source has no such cell).

        For a source with a time axis, ``steps`` gives the time-step index of each
        cell (one int, or one per cell); -1 means no step and yields NaN.
        """
        cells = np.asarray(cells, dtype=np.uint64)
        out = np.full(cells.shape, np.nan)
        arr = self._array(name)
        if arr is None or cells.size == 0:
            return out
        ids = self._index[0]
        pos = np.clip(np.searchsorted(ids, cells), 0, ids.size - 1)
        hit = ids[pos] == cells
        if arr.ndim == 1:
            out[hit] = arr[pos[hit]]
            return out
        if steps is None:
            raise ValueError(f"{self.title} has {arr.shape[0]} time steps; say which (steps=)")
        steps = np.broadcast_to(np.asarray(steps, dtype=np.int64), cells.shape)
        hit &= steps >= 0
        out[hit] = arr[steps[hit], pos[hit]]
        return out


def from_dataset(ds, variable: str | None = None, *, label: str | None = None) -> Source:
    """A :class:`Source` from a DGGS-Zarr dataset already open (or built in memory)."""
    a = ds.attrs
    dggs = a.get("dggs") or {}
    where = label or a.get("title", "dataset")
    if dggs.get("name") != "healpix" or dggs.get("indexing_scheme") != "nested":
        raise ValueError(f"{where} is not a NESTED HEALPix dataset (dggs={dggs})")
    cell_dim = dggs.get("spatial_dimension", "cells")
    if variable is None:
        candidates = [v for v in ds.data_vars
                      if cell_dim in ds[v].dims and not v.endswith(COMPANION_SUFFIXES)
                      and v not in COMPANION_VARIABLES]
        if len(candidates) != 1:
            raise ValueError(f"specify variable; candidates: {candidates}")
        variable = candidates[0]
    time_bounds = None
    period = a.get("climatology_bounds") or a.get("time_bounds")
    if "time" in ds[variable].dims:
        bname = ds["time"].attrs.get("bounds") or ds["time"].encoding.get("bounds")
        if not bname or bname not in ds:
            raise ValueError(f"{where}: {variable} has a time axis but no CF time bounds "
                             "(a 'bounds' attribute on 'time'); a step's extent cannot be guessed")
        time_bounds = np.asarray(ds[bname].transpose("time", ...).values).astype("datetime64[s]")
        if period is None:
            period = (_day(time_bounds[0, 0]), _day(time_bounds[-1, 1] - np.timedelta64(1, "s")))
    support_m = a.get("support_m")
    return Source(
        dataset=ds, depth=int(dggs["refinement_level"]), variable=variable,
        title=a.get("title", where), doi=a.get("source_doi"), license=a.get("license"),
        valid_time=a.get("valid_time"), period=tuple(period) if period else None,
        resampling_method=a.get("resampling_method"),
        support_m=None if support_m is None else float(support_m), time_bounds=time_bounds,
    )


def open_source(path, variable: str | None = None) -> Source:
    """Open a DGGS-Zarr dataset and read the metadata a sample must carry."""
    import xarray as xr

    return from_dataset(xr.open_zarr(path, consolidated=False), variable, label=str(path))


def _lookup(src: Source, cells: np.ndarray, name: str, steps=None) -> np.ndarray:
    """Values of variable ``name`` for ``cells`` (NaN where the source has no cell)."""
    return src.lookup(name, cells, steps)


def sample_cells(src: Source, cells, depth: int, steps=None) -> dict:
    """Sample a source at ``cells`` given at ``depth``.

    Depths are reconciled explicitly: a coarser source is **inherited** (the
    parent's value, no new information), a finer source is **aggregated** (the
    mean of the children present). Both are reported in ``support``.

    For a source with a time axis, ``steps`` is the time-step index of each cell
    (see :func:`time_step`).
    """
    from healpix_geo import nested

    cells = np.asarray(cells, dtype=np.uint64)
    if steps is not None:
        steps = np.broadcast_to(np.asarray(steps, dtype=np.int64), cells.shape)
    if depth == src.depth:
        value = src.lookup(src.variable, cells, steps)
        source_std = src.lookup(f"{src.variable}_std", cells, steps)
        how = "direct"
    elif depth > src.depth:
        parents = nested.zoom_to(cells, np.uint8(depth), np.uint8(src.depth))
        value = src.lookup(src.variable, parents, steps)
        source_std = src.lookup(f"{src.variable}_std", parents, steps)
        how = f"inherited from depth {src.depth}"
    else:
        children_per_cell = 4 ** (src.depth - depth)
        first = cells.astype(np.uint64) * np.uint64(children_per_cell)
        offsets = np.arange(children_per_cell, dtype=np.uint64)
        kids = (first[:, None] + offsets[None, :]).ravel()
        kid_steps = None if steps is None else np.repeat(steps, children_per_cell)
        vals = src.lookup(src.variable, kids, kid_steps).reshape(cells.size, children_per_cell)
        with np.errstate(invalid="ignore"):
            value = np.nanmean(vals, axis=1)
            source_std = np.nanstd(vals, axis=1)
        how = f"aggregated from depth {src.depth}"
    support = {"depth": int(depth), "how": how, "source_depth": src.depth}
    if src.support_m is not None:
        support["support_m"] = src.support_m
    return {"value": value, "source_std": source_std, "support": support}


# ------------------------------------------------------------------ time
def _day(t) -> str:
    return str(np.datetime64(t, "D"))


def _parse_date(s: str) -> tuple[np.datetime64, np.datetime64]:
    """Start and exclusive end of an ISO date at its own precision (year, month or day)."""
    s = s.strip().replace(" ", "T").split("T")[0]
    unit = {1: "Y", 2: "M", 3: "D"}.get(len(s.split("-")))
    if unit is None or (unit == "Y" and len(s) != 4):
        raise ValueError(s)
    start = np.datetime64(s, unit)
    return start.astype("datetime64[s]"), (start + np.timedelta64(1, unit)).astype("datetime64[s]")


def date_interval(date) -> tuple[np.datetime64, np.datetime64] | None:
    """A Darwin Core ``eventDate`` (or a year) as [start, end) seconds; None if unknown.

    Accepts ``2019``, ``2019-05``, ``2019-05-01``, a timestamp (read to the day),
    and ranges ``2019-05-01/2019-05-09`` or ``2019-05-01/09``. Raises ValueError
    on anything else.
    """
    if date is None or (isinstance(date, float) and np.isnan(date)):
        return None
    if isinstance(date, (int, float, np.integer, np.floating)):
        date = str(int(date))
    text = str(date).strip()
    if not text:
        return None
    if "/" in text:
        a, b = text.split("/", 1)
        a_parts, b_parts = a.split("T")[0].split("-"), b.split("T")[0].split("-")
        if len(b_parts) < len(a_parts):  # "2019-05-01/09": the end inherits the start's prefix
            b = "-".join(a_parts[:len(a_parts) - len(b_parts)] + b_parts)
        return _parse_date(a)[0], _parse_date(b)[1]
    return _parse_date(text)


def time_step(date, src: Source) -> tuple[int, str]:
    """For a source with a time axis: the step that contains the record's whole date.

    Returns ``(index, status)``; index -1 when there is none, because the date is
    unknown, outside the source's steps, or spans several of them (a year against
    monthly means). No step is chosen for a date that does not fit inside one.
    """
    if not src.has_time:
        raise ValueError(f"{src.title} has no time axis; use time_status")
    try:
        when = date_interval(date)
    except ValueError:
        return -1, f"record date not understood: {date!r}"
    if when is None:
        return -1, "record date unknown"
    start, end = when
    b = src.time_bounds
    inside = np.flatnonzero((b[:, 0] <= start) & (end <= b[:, 1]))
    if inside.size:
        i = int(inside[0])
        return i, f"inside time step {_day(b[i, 0])}..{_day(b[i, 1] - np.timedelta64(1, 's'))}"
    overlap = np.flatnonzero((b[:, 0] < end) & (start < b[:, 1]))
    if overlap.size:
        return -1, f"record date spans {overlap.size} time steps"
    return -1, f"outside {_day(b[0, 0])}..{_day(b[-1, 1] - np.timedelta64(1, 's'))}"


def time_status(date: str | int | None, src: Source) -> str:
    """Where a record's date sits relative to the source's period (or its time steps)."""
    if src.has_time:
        return time_step(date, src)[1]
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


# --------------------------------------------------------------- records
def sample_records(records: list[dict], src: Source, depth: int) -> list[dict]:
    """Attach a source's values to records that already carry cells.

    Records must come from a connector's ``to_cells``: each needs ``cell`` and,
    for the positional component, ``footprint_cells`` and ``footprint_known``.
    For a source with a time axis, each record is sampled at the time step that
    contains its ``eventDate`` (else its ``year``); the step is reported as
    ``time_step``, and a record whose date fits no single step gets NaN.
    """
    if not records:
        return []
    dates = [r.get("eventDate") or r.get("year") for r in records]
    if src.has_time:
        steps_status = [time_step(d, src) for d in dates]
        steps = np.array([s for s, _ in steps_status], dtype=np.int64)
        statuses = [s for _, s in steps_status]
    else:
        steps, statuses = None, [time_status(d, src) for d in dates]

    own = np.array([r["cell"] for r in records], dtype=np.uint64)
    at_own = sample_cells(src, own, depth, steps)

    # Every footprint in one lookup, then split back per record.
    fps = [np.asarray(r.get("footprint_cells", [r["cell"]]), dtype=np.uint64) for r in records]
    sizes = np.array([f.size for f in fps])
    multi = np.flatnonzero(sizes > 1)
    fp_values: dict[int, np.ndarray] = {}
    if multi.size:
        cells = np.concatenate([fps[i] for i in multi])
        fp_steps = None if steps is None else np.repeat(steps[multi], sizes[multi])
        flat = sample_cells(src, cells, depth, fp_steps)["value"]
        for i, chunk in zip(multi, np.split(flat, np.cumsum(sizes[multi])[:-1])):
            fp_values[int(i)] = chunk

    out = []
    for i, r in enumerate(records):
        vals = fp_values.get(i, np.array([at_own["value"][i]]))
        finite = vals[np.isfinite(vals)]
        row = {
            **r,
            "variable": src.variable,
            "value": float(at_own["value"][i]),
            "support": at_own["support"],
            "source_std": float(at_own["source_std"][i]),
            "positional_spread": float(finite.std()) if finite.size > 1 else 0.0,
            "positional_min": float(finite.min()) if finite.size else float("nan"),
            "positional_max": float(finite.max()) if finite.size else float("nan"),
            "positional_known": bool(r.get("footprint_known", False)),
            "time_status": statuses[i],
        }
        if src.has_time:
            k = int(steps[i])
            row["time_step"] = (None if k < 0 else
                                [_day(src.time_bounds[k, 0]),
                                 _day(src.time_bounds[k, 1] - np.timedelta64(1, "s"))])
        out.append(row)
    return out


def provenance(src: Source, depth: int) -> dict:
    """What a result must carry to be citable and reproducible."""
    p = {
        "dataset": src.title, "dataset_doi": src.doi, "license": src.license,
        "source_depth": src.depth, "sampled_at_depth": int(depth),
        "valid_time": src.valid_time, "period": list(src.period) if src.period else None,
        "resampling_method": src.resampling_method,
        "producer": f"healpix-connector {__version__}",
    }
    if src.has_time:
        p["time_steps"] = int(src.time_bounds.shape[0])
    if src.support_m is not None:
        p["support_m"] = src.support_m
    return p
