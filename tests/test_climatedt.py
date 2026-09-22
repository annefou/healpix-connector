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
