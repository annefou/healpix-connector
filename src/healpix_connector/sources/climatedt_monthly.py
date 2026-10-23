"""Climate DT generation 2 monthly means (stream ``clmn``) as a source of cell values.

Generation 2 of the Climate DT (DestinE) publishes monthly means beside the hourly
``clte`` stream, at both resolutions: ``high`` = HEALPix depth 10 (nside 1024, ~6.4 km)
and ``standard`` = depth 7 (nside 128, ~51 km), NESTED, **on a sphere**; historical
1990-2014 (activity ``baseline``, experiment ``hist``) and SSP3-7.0 2015-2049
(activity ``projections``). Measured 2026-10-23 with IFS-FESOM (LUMI) and IFS-NEMO
(MareNostrum 5). healpix-convert's settings still describe generation 1 (``clte``,
depth 7, 2020-2039), which is why this source builds its own request.

What the monthly stream holds and what it does not (2026-10-23): ``avg_2t``
(228004), ``avg_tprate`` (235055), ``avg_2d`` (235168) and the radiation fluxes;
**no monthly maximum or minimum temperature**. Parameters are requested by numeric
id and every returned message is checked: a short name can resolve to another field
(``mx2t24`` came back as ``avg_snswrfcs``, clear-sky short-wave radiation).

Facts about the service (endpoints, years, parameter ids, transfer size) are read from
``climatedt_gen2.json``, which records when they were measured; the code holds none.

Sphere -> WGS84: never by matching ids (the same id is ~0.1 deg of latitude away at
mid-latitudes). Default ``bilinear_conservative`` (see ``healpix_connector.resampling``);
``nearest`` reproduces a spherical-cell lookup of the nearest centre. MareNostrum 5
drops long transfers, so a year is fetched as one request per parameter and
half-year.
"""

from __future__ import annotations

import calendar
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

METHODS = ("bilinear_conservative", "nearest")
#: Facts about the external service (endpoints, years, parameter ids, transfer size) live in
#: ``climatedt_gen2.json`` with the date they were measured, not in this code. Pass another
#: ``service`` dict (same shape) when the service changes or to try another experiment.
SERVICE_FILE = Path(__file__).with_name("climatedt_gen2.json")


def load_service(path: Path | None = None) -> dict:
    return json.loads(Path(path or SERVICE_FILE).read_text())


def _experiment(service: dict, year: int) -> dict:
    for e in service["experiments"]:
        if e["years"][0] <= year <= e["years"][1]:
            return e
    spans = ", ".join(f"{e['experiment']} {e['years'][0]}-{e['years'][1]}" for e in service["experiments"])
    raise ValueError(f"{year} is outside the stream ({spans}; measured {service['measured']})")


def request(year: int, param: str, months, *, model: str, resolution: str = "high", service: dict | None = None) -> dict:
    """The Polytope request for ``param`` (a key of the service's ``params``) over ``months`` of ``year``."""
    service = service or load_service()
    e = _experiment(service, year)
    if resolution not in service["resolutions"]:
        raise ValueError(f"resolution must be one of {tuple(service['resolutions'])}")
    return {**service["request"], "activity": e["activity"], "experiment": e["experiment"], "model": model,
            "resolution": resolution, "param": service["params"][param]["id"], "year": str(year),
            "month": "/".join(str(m) for m in months)}


def check_message(meta: dict, *, param: str, year: int, depth: int, service: dict | None = None) -> int:
    """Validate one GRIB message's metadata; return its month. Raises on any mismatch."""
    service = service or load_service()
    want = service["params"][param]["id"]
    if str(meta["paramId"]) != want:
        raise ValueError(f"asked for {param} ({want}), got {meta.get('shortName')} ({meta['paramId']})")
    if meta["gridType"] != "healpix" or int(meta["Nside"]) != 2 ** depth or meta["orderingConvention"] != "nested":
        raise ValueError(f"unexpected grid {meta['gridType']} Nside {meta['Nside']} {meta['orderingConvention']}")
    date = int(meta["dataDate"])
    if date // 10000 != year:
        raise ValueError(f"message for {date // 10000}, asked for {year}")
    return date // 100 % 100


def to_monthly_units(param: str, values: np.ndarray, year: int, month: int,
                     service: dict | None = None) -> tuple[np.ndarray, str]:
    """degC for temperatures; mm per month for the precipitation rate."""
    service = service or load_service()
    if service["params"][param]["units"] == "K":
        return values - 273.15, "degC"
    return values * 86400.0 * calendar.monthrange(year, month)[1], "mm month-1"


@dataclass(frozen=True)
class YearOfMonths:
    year: int
    cell_ids: np.ndarray
    values: dict          # param -> (12, n) array in the units of ``units``
    units: dict
    valid_share: dict     # param -> (12, n): share of each cell covered by source samples
    method: str
    depth: int
    requests: list


class MonthlyStreamSampler:
    """Monthly means of Climate DT generation 2 over fixed WGS84 cells, one year at a time.

    ``cell_ids`` are WGS84 NESTED cells at the stream's depth (10 for ``high``).
    The source samples used are the spherical cells covering those cells plus a
    one-cell ring, so the resampler sees every neighbour it needs.
    """

    def __init__(self, cell_ids, params: list[str], *, model: str = "IFS-FESOM", resolution: str = "high",
                 method: str = "bilinear_conservative", local_dir: Path = Path("data/climatedt_clmn"),
                 retries: int = 6, service: dict | None = None):
        import healpix_geo.nested as hn

        self.service = service or load_service()
        if method not in METHODS:
            raise ValueError(f"method must be one of {METHODS}, got {method!r}")
        if model not in self.service["endpoints"]:
            raise ValueError(f"model must be one of {tuple(self.service['endpoints'])}, got {model!r}")
        unknown = [p for p in params if p not in self.service["params"]]
        if unknown:
            raise ValueError(f"not in the monthly stream as measured {self.service['measured']}: {unknown}")
        self.depth = self.service["resolutions"][resolution]
        self.cell_ids = np.sort(np.asarray(cell_ids, dtype=np.uint64))
        self.params, self.model, self.resolution, self.method = list(params), model, resolution, method
        self.local_dir, self.retries = Path(local_dir), retries
        clon, clat = hn.healpix_to_lonlat(self.cell_ids, self.depth, ellipsoid="WGS84")
        seed = np.asarray(hn.lonlat_to_healpix(np.asarray(clon), np.asarray(clat), np.uint8(self.depth),
                                               ellipsoid="sphere"), np.uint64)
        ring = np.asarray(hn.kth_neighbourhood(seed, self.depth, 2), np.uint64).ravel()
        self.src = np.unique(ring[ring < np.uint64(12 * 4 ** self.depth)])
        slon, slat = hn.healpix_to_lonlat(self.src, self.depth, ellipsoid="sphere")
        if method == "bilinear_conservative":
            from healpix_connector.resampling import BilinearConservative
            self._op = BilinearConservative(slon, slat, self.cell_ids, self.depth)
        else:
            idx = np.asarray(hn.lonlat_to_healpix(np.asarray(clon), np.asarray(clat), np.uint8(self.depth),
                                                  ellipsoid="sphere"), np.uint64)
            self._nearest = np.searchsorted(self.src, idx)

    def _resample(self, field: np.ndarray):
        if self.method == "bilinear_conservative":
            return self._op(field)
        v = field[self._nearest]
        return v, np.isfinite(v).astype(float)

    def _retrieve(self, req: dict, target: Path) -> None:
        from polytope.api import Client

        client = Client(address=self.service["endpoints"][self.model], quiet=True)
        for attempt in range(self.retries):
            try:
                client.retrieve(self.service["collection"], req, str(target))
                return
            except (Exception, SystemExit) as exc:  # the polytope client exits on a refusal
                target.unlink(missing_ok=True)
                if attempt == self.retries - 1:
                    raise RuntimeError(f"Climate DT request failed: {json.dumps(req)}") from exc
                time.sleep(60 * (attempt + 1))

    def year(self, year: int) -> YearOfMonths:
        import eccodes

        self.local_dir.mkdir(parents=True, exist_ok=True)
        n = self.cell_ids.size
        values = {p: np.full((12, n), np.nan) for p in self.params}
        share = {p: np.full((12, n), np.nan) for p in self.params}
        units, reqs, seen = {}, [], set()
        step = int(self.service["transfer"]["months_per_request"])
        for p in self.params:
            for k, months in enumerate(range(m0, min(m0 + step, 13)) for m0 in range(1, 13, step)):
                req = request(year, p, months, model=self.model, resolution=self.resolution, service=self.service)
                path = self.local_dir / f"{self.model}_{self.resolution}_{p}_{year}_{k}.grib"
                self._retrieve(req, path)
                reqs.append(req)
                with open(path, "rb") as f:
                    while (h := eccodes.codes_grib_new_from_file(f)) is not None:
                        meta = {key: eccodes.codes_get(h, key) for key in
                                ("paramId", "shortName", "gridType", "Nside", "orderingConvention", "dataDate")}
                        month = check_message(meta, param=p, year=year, depth=self.depth, service=self.service)
                        field = eccodes.codes_get_values(h)[self.src.astype(np.int64)]
                        eccodes.codes_release(h)
                        field, units[p] = to_monthly_units(p, field, year, month, self.service)
                        values[p][month - 1], share[p][month - 1] = self._resample(field)
                        seen.add((p, month))
                path.unlink(missing_ok=True)
        missing = [(p, m) for p in self.params for m in range(1, 13) if (p, m) not in seen]
        if missing:
            raise RuntimeError(f"{year}: months missing from the stream: {missing}")
        return YearOfMonths(year, self.cell_ids, values, units, share, self.method, self.depth, reqs)
