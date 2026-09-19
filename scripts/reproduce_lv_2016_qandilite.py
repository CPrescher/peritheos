#!/usr/bin/env python3
"""Independently reproduce Lv et al. (2016), pp. 303–305, doi:0794-1.

The source specifies least squares and EoS fit 5.2, but no input file or exact
weights. All-row equal-weight pressure fitting is the registered audit;
volume-weighted and propagated-error fits test the unspecified weighting.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from peritheos import get_eos_record

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/qandilite-lv-2016-table1-pv.csv"
OUTPUT = ROOT / "docs/data/lv-2016-qandilite.json"
PUBLISHED = {2: [603.21, 172.0], 3: [603.10, 175.0, 3.5]}


def load_data():
    """Read every printed coordinate and uncertainty without imputation."""
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        name: np.array([float(row[name]) if row[name] else np.nan for row in rows])
        for name in rows[0]
    }


def pressure(volume, parameters):
    """Implement the unnumbered p. 303 equation without Peritheos."""
    v0, k0 = parameters[:2]
    kp = parameters[2] if len(parameters) == 3 else 4.0
    eta = (v0 / np.asarray(volume)) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def slope(volume, parameters):
    """Analytic dP/dV for uncertainty propagation."""
    v0, k0 = parameters[:2]
    kp = parameters[2] if len(parameters) == 3 else 4.0
    eta = (v0 / volume) ** (1.0 / 3.0)
    c = 0.75 * (kp - 4)
    return (
        -0.5
        * k0
        / volume
        * (
            (7 * eta**7 - 5 * eta**5) * (1 + c * (eta**2 - 1))
            + (eta**7 - eta**5) * 2 * c * eta**2
        )
    )


def fit(data, order, mode="unweighted_pressure", start=0):
    """Fit simultaneously; no staged source fitting or hidden exclusions."""
    p, v, sv = (
        data[name][start:] for name in ("pressure_gpa", "volume_a3", "volume_sigma_a3")
    )
    sp = data["pressure_sigma_gpa"][start:].copy()
    # This assumption is confined to the effective-variance sensitivity. It is
    # never written into the source table as a measured zero uncertainty.
    sp[~np.isfinite(sp)] = 0.0

    def residual(parameters):
        if mode == "weighted_volume":
            predicted = np.array(
                [
                    brentq(
                        lambda vol: pressure(vol, parameters) - value,
                        parameters[0] * 0.7,
                        parameters[0] * 1.1,
                        xtol=1e-10,
                    )
                    for value in p
                ]
            )
            return (predicted - v) / sv
        delta = pressure(v, parameters) - p
        if mode == "effective_variance":
            return delta / np.sqrt(sp**2 + (slope(v, parameters) * sv) ** 2)
        return delta

    initial = np.array([600.0, 180.0] + ([4.0] if order == 3 else []))
    result = least_squares(
        residual,
        initial,
        bounds=(
            [580, 100] + ([0] if order == 3 else []),
            [630, 250] + ([10] if order == 3 else []),
        ),
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    if not result.success:
        raise RuntimeError(result.message)
    delta = pressure(v, result.x) - p
    return {
        "rows": len(p),
        "parameters": dict(zip(["V0", "K0", "K0_prime"][:order], result.x)),
        "pressure_rmse_gpa": float(np.sqrt(np.mean(delta**2))),
        "sum_squared_objective": float(np.sum(result.fun**2)),
    }


def reproduce():
    """Return source-observation, implementation-parity and refit diagnostics."""
    data = load_data()
    results = {}
    for order, parameters in PUBLISHED.items():
        record = get_eos_record(f"qandilite_lv_2016_bm{order}")
        independent = pressure(data["volume_a3"], parameters)
        residual = independent - data["pressure_gpa"]
        error = np.abs(independent - record.pressure(data["volume_a3"]))
        results[f"bm{order}"] = {
            "published_parameters": parameters,
            "published_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "published_max_abs_residual_gpa": float(np.max(np.abs(residual))),
            "high_pressure_checkpoint": {
                "volume_a3": 559.91,
                "observed_pressure_gpa": 14.9,
                "calculated_pressure_gpa": float(independent[-1]),
                "observed_volume_at_14_9_gpa_a3": 559.91,
                "calculated_volume_at_14_9_gpa_a3": float(record.volume(14.9)),
            },
            "implementation_max_error_gpa": float(np.max(error)),
            "all_rows_unweighted_pressure": fit(data, order),
            "all_rows_weighted_volume": fit(data, order, "weighted_volume"),
            "all_rows_effective_variance": fit(data, order, "effective_variance"),
            "dac_only_unweighted_pressure": fit(data, order, start=1),
        }
    return results


if __name__ == "__main__":
    artifact = {
        "doi": "10.1007/s00269-015-0794-1",
        "source": "Table 1, p. 303; all 18 printed rows at 300 K",
        "qualification": (
            "Source residual direction, weights, covariance and input file are not "
            "published. Equal-weight pressure fitting is the registered audit; "
            "weighted-volume and effective-variance fits are sensitivity checks. "
            "Only the latter assumes exact ambient pressure for propagation."
        ),
        "results": reproduce(),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2, allow_nan=False) + "\n")
    print(json.dumps(artifact, indent=2))
