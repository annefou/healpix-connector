"""Climate DT source: the parts that need no Polytope call.

The converter itself belongs to healpix-convert and is tested there; what is
ours is the month-long request, the cache name that survives a date range, and
the parameter table.
"""

import pytest

pytest.importorskip("healpix_convert")

from healpix_connector.sources import climatedt  # noqa: E402


def test_source_is_read_from_healpix_convert():
    from healpix_convert.settings import climatedt as cfg

    assert climatedt.SOURCE["experiment"] == cfg.CDT_EXPERIMENT
    assert climatedt.SOURCE["model"] == cfg.CDT_MODEL
    assert cfg.POLYTOPE_ADDRESS in climatedt.SOURCE["access"]


def test_param_ids_cover_the_converters_variables_and_precipitation():
    from healpix_convert.settings import climatedt as cfg

    for name, (param, *_ ) in cfg.CDT_SFC_VARIABLE_META.items():
        assert climatedt.PARAM_ID[name] == param
    assert climatedt.PARAM_ID["tprate"] == "260048"


def test_a_date_range_gives_a_usable_cache_name(tmp_path):
    conv = climatedt.MonthConverter(date="20300101/to/20300131", time="0000/1200",
                                    params="167", local_dir=tmp_path)
    path = conv._grib_path()
    assert "/" not in path.name and path.parent == tmp_path
    assert path.name == "cdt_20300101-to-20300131T0000-1200_sfc.grib"


def test_the_request_is_the_converters_own(tmp_path):
    from healpix_convert.settings import climatedt as cfg

    conv = climatedt.MonthConverter(date="20300101/to/20300131", time="0000/1200",
                                    params="167/260048", local_dir=tmp_path)
    assert (conv.experiment, conv.expver, conv.model, conv.resolution) == (
        cfg.CDT_EXPERIMENT, cfg.CDT_EXPVER, cfg.CDT_MODEL, cfg.CDT_RESOLUTION)
    assert conv.levtype == "sfc"


def test_sampler_declares_the_source_cell_as_support():
    from healpix_connector.conventions import cell_size_m

    sampler = climatedt.MonthlySampler([10786, 10792], ["t2m"])
    assert sampler._resampler is None  # built on the first month, not before
    assert round(cell_size_m(climatedt.DEPTH) / 1000) == 51


def test_a_month_becomes_a_one_step_dggs_source():
    import numpy as np

    from healpix_connector import sample as S

    cells = np.array([5, 3, 9], dtype=np.int64)  # in request order, not sorted
    m = climatedt.MonthlyMeans(year=2030, month=2, cell_ids=cells,
                               values={"t2m": np.array([280.0, 281.0, 282.0])}, method="psf",
                               support_m=0.0, filled_cells=np.array([9]))
    ds = m.to_dataset()
    assert ds["t2m"].attrs["units"] == "K" and ds.attrs["dggs"]["refinement_level"] == climatedt.DEPTH
    src = S.from_dataset(ds)
    assert src.variable == "t2m" and src.period == ("2030-02-01", "2030-02-28")
    step, _ = S.time_step("2030-02-14", src)
    got = S.sample_cells(src, np.array([3, 5, 9], dtype=np.uint64), climatedt.DEPTH, steps=step)
    assert got["value"].tolist() == [281.0, 280.0, 282.0]
    assert ds["filled_by_fallback"].sel(cells=ds.cell_ids == 9).item() == 1
