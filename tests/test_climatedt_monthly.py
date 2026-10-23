"""Climate DT generation 2 monthly stream: request, message checks and units (no Polytope call)."""

import numpy as np
import pytest

from healpix_connector.sources import climatedt_monthly as m


def test_request_switches_from_hist_to_ssp_after_2014():
    h = m.request(2014, "avg_2t", range(1, 7), model="IFS-FESOM")
    p = m.request(2015, "avg_2t", range(7, 13), model="IFS-FESOM")
    assert (h["activity"], h["experiment"]) == ("baseline", "hist")
    assert (p["activity"], p["experiment"]) == ("projections", "SSP3-7.0")
    assert h["stream"] == "clmn" and h["generation"] == "2" and h["resolution"] == "high"
    assert h["param"] == "228004" and p["month"] == "7/8/9/10/11/12"


def test_request_refuses_years_outside_generation_2():
    with pytest.raises(ValueError):
        m.request(1989, "avg_2t", [1], model="IFS-NEMO")
    with pytest.raises(ValueError):
        m.request(2050, "avg_2t", [1], model="IFS-NEMO")


def _meta(**kw):
    base = {"paramId": 228004, "shortName": "avg_2t", "gridType": "healpix", "Nside": 1024,
            "orderingConvention": "nested", "dataDate": 20300701}
    return {**base, **kw}


def test_a_message_for_another_parameter_is_refused():
    """The server answered mx2t24 with avg_snswrfcs (235051) on 2026-10-23."""
    with pytest.raises(ValueError, match="asked for avg_2t"):
        m.check_message(_meta(paramId=235051, shortName="avg_snswrfcs"), param="avg_2t", year=2030, depth=10)


def test_grid_and_year_are_checked():
    assert m.check_message(_meta(), param="avg_2t", year=2030, depth=10) == 7
    with pytest.raises(ValueError):
        m.check_message(_meta(Nside=128), param="avg_2t", year=2030, depth=10)
    with pytest.raises(ValueError):
        m.check_message(_meta(dataDate=20290701), param="avg_2t", year=2030, depth=10)


def test_service_facts_come_from_the_data_file_and_can_be_replaced():
    svc = m.load_service()
    assert svc["measured"] and svc["params"]["avg_2t"]["id"] == "228004"
    other = {**svc, "experiments": [{"activity": "projections", "experiment": "SSP2-4.5", "years": [2015, 2049]}]}
    assert m.request(2030, "avg_2t", [1], model="IFS-NEMO", service=other)["experiment"] == "SSP2-4.5"
    with pytest.raises(ValueError, match="outside the stream"):
        m.request(2000, "avg_2t", [1], model="IFS-NEMO", service=other)


def test_units():
    t, u = m.to_monthly_units("avg_2t", np.array([273.15, 300.15]), 2030, 1)
    np.testing.assert_allclose(t, [0.0, 27.0]); assert u == "degC"
    p, u = m.to_monthly_units("avg_tprate", np.array([1e-5]), 2030, 2)   # February 2030: 28 days
    np.testing.assert_allclose(p, [1e-5 * 86400 * 28]); assert u == "mm month-1"
