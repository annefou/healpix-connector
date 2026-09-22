"""Check CHELSA's published bioclimatic layers against its own monthly layers.

CHELSA's file specification and the GeoTIFFs' GDAL tags both describe bio4 as
degC with scale 0.1. Recomputing it from the monthly temperature layers of the
same climatology shows the published values are 100x the standard deviation -
WorldClim's convention ("standard deviation x100"), which CHELSA does not state.

This script recomputes bio4, bio12 and bio15 from CHELSA's own monthly tas and
pr layers over a small window, so the comparison needs no other dataset.

Usage: python examples/verify_chelsa_scaling.py
"""

import json
from pathlib import Path

import numpy as np

from healpix_connector.sources import chelsa
from healpix_connector.sources.geotiff import read_window

WINDOW = (-6.0, 39.0, -5.0, 40.0)  # 1 degree of inland Iberia
MONTHLY = ("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/"
           "1981-2010/{v}/CHELSA_{v}_{m:02d}_1981-2010_V.2.1.tif")


def monthly(var: str) -> np.ndarray:
    return np.stack([read_window(MONTHLY.format(v=var, m=m), WINDOW, pad_deg=0.0).values
                     for m in range(1, 13)])


def main():
    tas, pr = monthly("tas"), monthly("pr")
    published = {n: read_window(chelsa.url(n), WINDOW, pad_deg=0.0).values for n in (4, 12, 15)}

    sd_pop = float(np.nanmean(tas.std(axis=0, ddof=0)))
    sd_sample = float(np.nanmean(tas.std(axis=0, ddof=1)))
    cv = float(np.nanmean(100.0 * pr.std(axis=0, ddof=0) / pr.mean(axis=0)))
    annual = float(np.nanmean(pr.sum(axis=0)))
    pub = {n: float(np.nanmean(a)) for n, a in published.items()}

    report = {
        "window_lon_lat": WINDOW,
        "source": "CHELSA v2.1 1981-2010 monthly tas and pr, and the published bio layers",
        "bio4": {
            "published": round(pub[4], 2),
            "recomputed_sd_population_degC": round(sd_pop, 3),
            "recomputed_sd_sample_degC": round(sd_sample, 3),
            "ratio_to_population_sd": round(pub[4] / sd_pop, 2),
            "ratio_to_sample_sd": round(pub[4] / sd_sample, 2),
            "documented": "degC, scale 0.1, offset 0 (file specification and GDAL tags)",
            "actual": "100 x the population standard deviation, i.e. degC/100",
        },
        "bio12": {"published": round(pub[12], 1), "recomputed_annual_sum": round(annual, 1),
                  "ratio": round(pub[12] / annual, 3), "verdict": "matches; unit kg m-2 is right"},
        "bio15": {"published": round(pub[15], 2), "recomputed_cv_percent": round(cv, 2),
                  "ratio": round(pub[15] / cv, 3),
                  "documented": "kg m-2 (file specification)",
                  "actual": "a coefficient of variation in percent; the value is right, the unit is not"},
    }
    out = Path(__file__).parent / "results/chelsa_scaling_check.json"
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
