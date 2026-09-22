"""Climate DT (DestinE, IFS-NEMO) as a source of cell values.

The dataset definition, the MARS request, the Polytope endpoint and the
sphere -> WGS84 correction all come from GRID4EARTH's healpix-convert
(``ClimateDTConverter``); nothing about the source is redefined here. This
module adds only what a connector needs on top of it:

* a whole month per request instead of one timestep, so a decade is 120
  requests rather than ~90,000;
* aggregation to monthly means;
* the correction restricted to the cells the caller asked for.

Climate DT is delivered on HEALPix depth 7 (nside 128) NESTED, but on a sphere,
while our cell ids are WGS84. healpix-convert leaves that difference uncorrected
by default and declares ``nearest`` as its resampler for these fields, offering
``PSFResampler`` as an optional rigorous correction. Both are available here
through ``method``, restricted to the caller's cells (``out_cell_ids``), which
is what makes a decade affordable. ``nearest`` is the default because PSF is a
deconvolution: on precipitation it overshoots into negative values (114 of 480
Iberian cells for a January mean), which a non-negative quantity cannot take.
Resampling is linear, so correcting the monthly mean is the same as correcting
every timestep and then averaging.
"""

from __future__ import annotations

import calendar
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from healpix_convert.converters.climatedt import ClimateDTConverter
from healpix_convert.settings import climatedt as cfg

#: What healpix-convert says this source is. Kept as a reference, not a copy.
SOURCE = {
    "name": "Climate DT (DestinE)",
    "model": cfg.CDT_MODEL,
    "experiment": cfg.CDT_EXPERIMENT,
    "stream": cfg.CDT_STREAM,
    "type": cfg.CDT_TYPE,
    "resolution": cfg.CDT_RESOLUTION,
    "access": f"Polytope {cfg.POLYTOPE_ADDRESS}, collection {cfg.POLYTOPE_COLLECTION}",
    "grid": "HEALPix depth 7 (nside 128), NESTED, on a sphere",
    "converter": "healpix-convert ClimateDTConverter",
}

DEPTH = 7
#: Beyond healpix-convert's surface variables (its CDT_SFC_VARIABLE_META).
EXTRA_PARAMS = {"tprate": ("260048", "kg m-2 s-1", "Total precipitation rate")}
PARAM_ID = {name: meta[0] for name, meta in
            {**cfg.CDT_SFC_VARIABLE_META, **EXTRA_PARAMS}.items()}
#: PSF threshold healpix-convert uses for its optional ellipsoid correction.
PSF_THRESHOLD = 0.5
#: healpix-convert's declared resampler for these fields, and the alternative.
METHODS = ("nearest", "psf")


class MonthConverter(ClimateDTConverter):
    """ClimateDTConverter over a whole month, with a filename-safe cache path.

    Temporary (to move upstream): healpix-convert takes a single date and builds
    its cache name from it, so a MARS date range - which Polytope accepts and
    which the request itself already supports - cannot be cached. Only the file
    name is changed here; the request is the converter's own.
    """

    def _grib_path(self) -> Path:
        tag = f"cdt_{self.date}T{self.time}_{self.levtype}".replace("/", "-")
        return self.local_dir / f"{tag}.grib"

    def download(self) -> Path:
        """Fetch the GRIB if it is not already cached, and return its path."""
        path = self._grib_path()
        if not path.exists():
            self.local_dir.mkdir(parents=True, exist_ok=True)
            try:
                self._download(path)
            except SystemExit as exc:  # the polytope client exits on a refusal
                raise RuntimeError(f"Climate DT refused {self.date} {self.time}") from exc
        return path

    def dataset(self):
        """The GRIB as an xarray dataset, opened the way the converter opens it."""
        self._load_dataset(self.download())
        return self._ds


@dataclass(frozen=True)
class MonthlyMeans:
    """Monthly mean of each variable over ``cell_ids`` (depth-7 WGS84 cells)."""

    year: int
    month: int
    cell_ids: np.ndarray
    values: dict[str, np.ndarray]
    method: str
    support_m: float
    #: cells whose PSF kernel had too little support, filled by its own fallback
    filled_cells: np.ndarray


class MonthlySampler:
    """Monthly means of Climate DT variables over a fixed set of cells.

    The sphere -> WGS84 resampler is built once (the grid does not change from
    month to month) and reused, which is what makes a decade affordable.

    ``method="nearest"`` is healpix-convert's declared resampler for these
    fields. ``method="psf"`` uses its optional PSF correction instead; that one
    can leave cells with too little kernel support, which healpix-resample's own
    fallback fills, and it is unsuitable for non-negative fields such as
    precipitation (see the module docstring).
    """

    def __init__(self, cell_ids, variables: list[str], *, method: str = "nearest",
                 times: str = "0000/1200", local_dir: Path = Path("data/climatedt"),
                 keep_grib: bool = False):
        if method not in METHODS:
            raise ValueError(f"method must be one of {METHODS}, got {method!r}")
        self.cell_ids = np.asarray(cell_ids, dtype=np.int64)
        self.variables = list(variables)
        self.method = method
        self.times = times
        self.local_dir = Path(local_dir)
        self.keep_grib = keep_grib
        self._resampler = None
        self._filled: np.ndarray = np.array([], dtype=np.int64)

    def _converter(self, year: int, month: int) -> MonthConverter:
        last = calendar.monthrange(year, month)[1]
        return MonthConverter(
            date=f"{year}{month:02d}01/to/{year}{month:02d}{last:02d}",
            time=self.times,
            params="/".join(PARAM_ID[v] for v in self.variables),
            local_dir=self.local_dir,
        )

    def _build_resampler(self, ds) -> None:
        """The sphere -> WGS84 resampler, restricted to the cells asked for."""
        import healpix_resample

        lon = ds.longitude.values.ravel().astype(float)
        lat = ds.latitude.values.ravel().astype(float)
        common = dict(lon_deg=lon, lat_deg=lat, level=DEPTH, out_cell_ids=self.cell_ids,
                      ellipsoid="WGS84", verbose=False)
        if self.method == "nearest":
            self._resampler = healpix_resample.NearestResampler(**common)
            return
        strict = healpix_resample.PSFResampler(threshold=PSF_THRESHOLD, **common)
        # which cells the kernel cannot compute, before asking for the fallback
        probe = np.asarray(strict.resample(np.zeros(lon.size), lam=0.0).cell_data, dtype=float)
        self._filled = self.cell_ids[np.isnan(probe)]
        self._resampler = healpix_resample.PSFResampler(
            threshold=PSF_THRESHOLD, fill_missing_out_cells=bool(self._filled.size), **common)

    def month(self, year: int, month: int) -> MonthlyMeans:
        from healpix_connector.conventions import cell_size_m

        conv = self._converter(year, month)
        ds = conv.dataset()
        if self._resampler is None:
            self._build_resampler(ds)

        values = {}
        for name in self.variables:
            t_dim = next((d for d in ("valid_time", "time") if d in ds[name].dims), None)
            field = (ds[name].mean(t_dim) if t_dim else ds[name]).values.ravel().astype(float)
            res = self._resampler.resample(field)
            order = np.searchsorted(np.asarray(res.cell_ids), self.cell_ids)
            values[name] = np.asarray(res.cell_data, dtype=float)[order]

        if not self.keep_grib:
            conv._grib_path().unlink(missing_ok=True)
        return MonthlyMeans(year=year, month=month, cell_ids=self.cell_ids, values=values,
                            method=self.method, support_m=cell_size_m(DEPTH),
                            filled_cells=self._filled)


def month(year: int, month: int, variables: list[str], cell_ids, **kwargs) -> MonthlyMeans:
    """One month for one set of cells; see :class:`MonthlySampler` for a series."""
    return MonthlySampler(cell_ids, variables, **kwargs).month(year, month)
