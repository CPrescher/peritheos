"""Audit Crichton's printed coefficients against independently read Figure 2.

The nominal-temperature fit is a diagnostic, never an original-data refit.
No curve samples, recovered anchor, or cited Ming data enter that fit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from peritheos import get_eos_record, get_material_document

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/crichton-2016-vanadium-reproduction.json"
RECORD = "vanadium_bcc_crichton_2016_bm3_thermal"
POINTS = DATA / "vanadium-crichton-2016-figure2-points.csv"
CURVES = DATA / "vanadium-crichton-2016-figure2-curves.csv"
PUBLISHED = np.array([27.995, 150.4, 5.5, 4.8e-5, -2.4e-8, -0.0446])
NAMES = ("V0", "K0", "K0_prime", "alpha0", "alpha1", "dK_dT")
PLOT_DENOMINATOR = 27.995 / 1.005


def source_pressure(volume, temperature, parameters=PUBLISHED):
    """Independent BM3 + earlier-EosFit exponential integral, in cell A^3."""
    v0, k0, kp, a0, a1, dk = parameters
    temperature = np.asarray(temperature)
    v0t = v0 * np.exp(a0 * (temperature - 300.0) + a1 / 2 * (temperature**2 - 300.0**2))
    kt = k0 + dk * (temperature - 300.0)
    x = (v0t / np.asarray(volume)) ** (1 / 3)
    return 1.5 * kt * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))


def source_volume(pressure, temperature):
    return brentq(
        lambda v: float(source_pressure(v, temperature)) - pressure, 24.0, 30.0
    )


def _rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def reproduce():
    document = get_material_document("vanadium_bcc")
    for dataset in document["datasets"]:
        resource = dataset["resource"]
        payload = (ROOT / "peritheos/data" / resource["path"]).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == resource["sha256"]

    record = get_eos_record(RECORD)
    curves = _rows(CURVES)
    checkpoints = []
    for row in curves:
        p, t, observed = (
            float(row[k]) for k in ("pressure_gpa", "temperature_k", "relative_volume")
        )
        calculated = source_volume(p, t) / PLOT_DENOMINATOR
        checkpoints.append(
            dict(
                pressure_gpa=p,
                temperature_k=t,
                plotted_relative_volume=observed,
                calculated_relative_volume=calculated,
                difference=calculated - observed,
            )
        )
    max_curve_error = max(abs(row["difference"]) for row in checkpoints)
    assert max_curve_error < 0.001

    rows = [row for row in _rows(POINTS) if row["series"] == "filled"]
    volume = np.array([float(r["relative_volume"]) for r in rows]) * PLOT_DENOMINATOR
    pressure = np.array([float(r["pressure_gpa"]) for r in rows])
    temperature = np.array([float(r["temperature_nominal_k"]) for r in rows])
    # Dimensionless optimizer coordinates; all six source coefficients are free.
    scale = np.array([28.0, 150.0, 5.0, 1e-5, 1e-8, 0.01])

    def residual(x):
        return source_pressure(volume, temperature, x * scale) - pressure

    fit = least_squares(
        residual,
        PUBLISHED / scale,
        method="lm",
        max_nfev=3000,
        ftol=1e-11,
        xtol=1e-11,
        gtol=1e-11,
    )
    assert fit.success
    published_residual = residual(PUBLISHED / scale)
    native = record.pressure(volume, temperature=temperature, check_validity=False)
    parity = float(np.max(np.abs(native - source_pressure(volume, temperature))))
    assert parity < 1e-10
    refitted = fit.x * scale
    proxy = {
        "status": "diagnostic_only",
        "original_fit_reproduced": False,
        "observations": len(rows),
        "selection": "All 29 distinguishable filled Figure 2 symbols; exclude 10 cited open symbols and recovered star. No fitted-curve readings used.",
        "objective": "Unweighted pressure residuals at nominal color-bin temperatures; simultaneous six-parameter least squares. Source reports unweighted fitting, but exact row temperatures and 33 additional observations cannot be individually recovered.",
        "free_parameters": list(NAMES),
        "fixed_parameters": ["Tr"],
        "parameters": [
            dict(
                parameter=name,
                published=float(pub),
                refit=float(val),
                difference=float(val - pub),
            )
            for name, pub, val in zip(NAMES, PUBLISHED, refitted)
        ],
        "published_rmse_gpa": float(np.sqrt(np.mean(published_residual**2))),
        "refit_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        "jacobian_rank": int(np.linalg.matrix_rank(fit.jac)),
        "uncertainty": "No parameter errors inferred: temperature binning, digitization, incomplete sampling and unreported covariance preclude original-fit statistics.",
    }
    expansion = {
        str(p): 100 * (source_volume(p, 1000) / source_volume(p, 300) - 1)
        for p in (0, 10)
    }
    return {
        "record_identifier": RECORD,
        "curve_checkpoints": checkpoints,
        "curve_max_abs_relative_volume_difference": max_curve_error,
        "curve_tolerance": 0.001,
        "native_independent_pressure_max_difference_gpa": parity,
        "alpha_300k_per_k": 4.8e-5 - 2.4e-8 * 300,
        "alpha_k_300k_gpa_per_k": (4.8e-5 - 2.4e-8 * 300) * 150.4,
        "k0_1000k_gpa": 150.4 - 0.0446 * 700,
        "expansion_300_to_1000k_percent": expansion,
        "source_conclusion_expansion_percent": {"0": 2.7, "10": 1.2},
        "nominal_temperature_proxy_fit": proxy,
    }


def reports_close(actual, saved):
    """Compare reports while bounding cancellation-sensitive optimizer drift."""
    if isinstance(actual, dict):
        if actual.keys() != saved.keys():
            return False
        for key in actual:
            if actual.get("parameter") == "K0" and key == "difference":
                # Subtracting 150.4 GPa from the fitted modulus amplifies relative
                # drift: Python 3.9/SciPy 1.13 differs by 3.2e-5 GPa here.
                # Keep this absolute allowance confined to that derived value;
                # fitted coefficients, pressure RMSE and all other fields retain
                # their existing tolerances.
                if not np.isclose(actual[key], saved[key], rtol=2e-5, atol=5e-5):
                    return False
            elif not reports_close(actual[key], saved[key]):
                return False
        return True
    if isinstance(actual, list):
        return len(actual) == len(saved) and all(
            reports_close(x, y) for x, y in zip(actual, saved)
        )
    if isinstance(actual, (int, float)) and not isinstance(actual, bool):
        return bool(np.isclose(actual, saved, rtol=2e-5, atol=1e-10))
    return actual == saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    # Rounding avoids saving optimizer last-bit differences between BLAS builds.
    def rounded(value):
        if isinstance(value, float):
            return float(f"{value:.8g}")
        if isinstance(value, dict):
            return {k: rounded(v) for k, v in value.items()}
        if isinstance(value, list):
            return [rounded(v) for v in value]
        return value

    result = rounded(reproduce())
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        saved = json.loads(OUTPUT.read_text(encoding="utf-8"))

        if not reports_close(result, saved):
            raise SystemExit("Crichton reproduction report is stale")
    else:
        OUTPUT.write_text(text, encoding="utf-8")
        print(text)


if __name__ == "__main__":
    main()
