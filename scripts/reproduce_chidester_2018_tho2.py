"""Independently reproduce Chidester (2018) equations and refit official ThO2 rows."""

from __future__ import annotations

import argparse
import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos import get_eos_record

ROOT = Path(__file__).resolve().parents[1]
CELL_PER_MOLAR = 4e24 / 6.02214076e23
OUTPUT = ROOT / "docs/data/chidester-2018-tho2-reproduction.json"
PUBLISHED = {
    "thorianite": [26.379, 204.0, 0.0035],
    "tho2_cotunnite": [24.75, 190.0, 0.0037],
}
ERRORS = {"thorianite": [0.007, 2.0, 0.0003], "tho2_cotunnite": [0.06, 3.0, 0.0004]}


def source_pressure(volume, temperature, parameters):
    """Equations 1–3 in molar cm3/mol, GPa, K; independent of package evaluator."""
    v0, k0, alpha_kt = parameters
    f = ((v0 / np.asarray(volume)) ** (2 / 3) - 1) / 2
    return 3 * k0 * f * (1 + 2 * f) ** 2.5 + alpha_kt * (np.asarray(temperature) - 300)


@lru_cache(maxsize=1)
def reproduce():
    """Source-constrained fits, with explicit diagnostic weighting alternatives."""
    output = {}
    for material, published in PUBLISHED.items():
        path = ROOT / f"peritheos/data/datasets/chidester_2018_{material}_pvt.csv"
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        mask = np.array([r["included_in_fit"] == "1" for r in rows])

        def column(name):
            return np.array([float(r[name]) for r in rows])[mask]

        p, v, t = (
            column(n) for n in ("pressure_gpa", "molar_volume_cm3_mol", "temperature_k")
        )
        pe, ve, te = (
            column(n)
            for n in (
                "pressure_error_gpa",
                "molar_volume_error_cm3_mol",
                "temperature_error_k",
            )
        )
        # The ambient point fixes V0; its zero residual adds no constraint to K0/C.
        initial = np.array(published[1:] if material == "thorianite" else published)
        scale = np.array([200, 0.004] if material == "thorianite" else [25, 200, 0.004])

        def full(values):
            return np.r_[26.379, values] if material == "thorianite" else values

        pred = source_pressure(v, t, published)
        x = (published[0] / v) ** (1 / 3)
        dpdv = -0.5 * published[1] * (7 * x**7 - 5 * x**5) / v
        # Frozen at published coefficients: no parameter-dependent denominator bias.
        effective_error = np.sqrt(
            pe**2
            + (dpdv * ve) ** 2
            + (published[2] * te) ** 2
            + np.where(t == 300, (0.03 * p) ** 2, 0)
        )
        fits = {}
        for mode, denominator in (
            ("unweighted", np.ones(len(p))),
            ("pressure_errors", np.maximum(pe, 1e-8)),
            ("effective_errors", effective_error),
        ):
            fit = least_squares(
                lambda values: (
                    (source_pressure(v, t, full(values * scale)) - p) / denominator
                ),
                initial / scale,
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
            )
            fitted = full(fit.x * scale)
            fits[mode] = {
                "parameters_molar": fitted.tolist(),
                "rmse_gpa": float(
                    np.sqrt(np.mean((source_pressure(v, t, fitted) - p) ** 2))
                ),
                "solver_success": bool(fit.success),
            }
        record = get_eos_record(f"{material}_chidester_2018_bm3_linear_thermal")
        native = record.pressure(v * CELL_PER_MOLAR, temperature=t)
        # Independent measured, off-reference states: original supplement IDs.
        ids = (
            ["B5_209", "B5_253"] if material == "thorianite" else ["B25_446", "B70_562"]
        )
        benchmarks = []
        for sample in ids:
            row = next(r for r in rows if r["sample"] == sample)
            bv, bt, bp = (
                float(row[n])
                for n in ("molar_volume_cm3_mol", "temperature_k", "pressure_gpa")
            )
            calc = float(source_pressure(bv, bt, published))
            tolerance = 2 * float(row["pressure_error_gpa"])
            benchmarks.append(
                {
                    "sample": sample,
                    "molar_volume_cm3_mol": bv,
                    "temperature_k": bt,
                    "measured_pressure_gpa": bp,
                    "calculated_pressure_gpa": calc,
                    "absolute_tolerance_gpa": tolerance,
                    "within_tolerance": abs(calc - bp) <= tolerance,
                }
            )
        output[material] = {
            "total_rows": len(rows),
            "selected_rows": int(mask.sum()),
            "excluded_samples": [
                r["sample"] for r, selected in zip(rows, mask) if not selected
            ],
            "published_rmse_gpa": float(np.sqrt(np.mean((pred - p) ** 2))),
            "native_max_difference_gpa": float(np.max(np.abs(native - pred))),
            "refits": fits,
            "benchmarks": benchmarks,
            "selection": "All high-T rows plus 300 K rows below 15 GPa; measured ambient V0 fixed."
            if material == "thorianite"
            else "All high-T rows; no 300 K observations exist.",
            "qualification": "Published objective, weights, covariance, confidence levels and fit staging are unspecified. Each diagnostic jointly varies the source-free parameters once. Effective-error weights propagate P,V,T errors at published coefficients and add the stated Ar 3% accuracy in quadrature as a sensitivity assumption; cross-row scale correlation is not known. No strict statistical parity claimed.",
        }
    return output


def ledger_outcome(record):
    """Register numerical similarity without claiming recovered source weights."""
    material = record["identifier"].split("_chidester_")[0]
    result = reproduce()[material]
    fitted = result["refits"]["effective_errors"]["parameters_molar"]
    comparisons = []
    for index, name in enumerate(("V0", "K0", "alpha_KT")):
        if material == "thorianite" and index == 0:
            continue
        factor = CELL_PER_MOLAR if index == 0 else 1
        pub, value, err = (
            a[index] * factor for a in (PUBLISHED[material], fitted, ERRORS[material])
        )
        comparisons.append(
            dict(
                parameter=name,
                published=pub,
                refit=value,
                published_error=err,
                refit_error=None,
                relative_difference=abs(value - pub) / abs(pub),
                within_combined_2sigma=None,
                similar=abs(value - pub) <= err,
            )
        )
    return {
        "status": "similar"
        if all(c["similar"] for c in comparisons)
        else "parity_not_achieved",
        "dataset_identifiers": record["fit_datasets"],
        "observations": result["selected_rows"],
        "selection": result["selection"],
        "fit_kind": "source_constrained_bm3_linear_thermal_diagnostic",
        "objective": "joint pressure residuals divided by fixed propagated coordinate errors including Ar 3% accuracy",
        "free_parameters": [c["parameter"] for c in comparisons],
        "parameters": comparisons,
        "rmse_gpa": result["refits"]["effective_errors"]["rmse_gpa"],
        "published_rmse_gpa": result["published_rmse_gpa"],
        "solver_success": result["refits"]["effective_errors"]["solver_success"],
        "reason": result["qualification"],
        "reproduction": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = reproduce()
    if any(
        r["native_max_difference_gpa"] > 1e-10
        or not all(b["within_tolerance"] for b in r["benchmarks"])
        for r in output.values()
    ):
        raise SystemExit("source reproduction failed")
    payload = json.dumps(output, indent=2, sort_keys=True) + "\n"
    if args.check:
        # Re-evaluate all evidence, allowing insignificant solver variation.
        saved = json.loads(OUTPUT.read_text())

        def compare(actual, expected):
            if isinstance(actual, dict):
                assert actual.keys() == expected.keys()
                for key, value in actual.items():
                    compare(value, expected[key])
            elif isinstance(actual, list):
                assert len(actual) == len(expected)
                for value, other in zip(actual, expected):
                    compare(value, other)
            elif isinstance(actual, float):
                np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=1e-10)
            else:
                assert actual == expected

        compare(output, saved)
    else:
        OUTPUT.write_text(payload)
    print(payload)


if __name__ == "__main__":
    main()
