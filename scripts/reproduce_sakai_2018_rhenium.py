"""Refit Sakai (2018) Figure 10 and check its maximum-pressure evaluation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos import get_eos_record
from peritheos.eos.rt import Vinet

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/sakai-2018-rhenium-reproduction.json"
RECORD = "rhenium_sakai_2018_yokoo_pt_vinet"
DATASET = "rhenium_sakai_2018_figure10"
DATA = ROOT / "peritheos/data/datasets/rhenium-sakai-2018-figure10.csv"


def load_points() -> list[dict]:
    provenance = json.loads(
        DATA.with_name("rhenium-sakai-2018-figure10-source.json").read_text()
    )
    if hashlib.sha256(DATA.read_bytes()).hexdigest() != provenance["csv_sha256"]:
        raise ValueError("Sakai Figure 10 CSV checksum mismatch")
    with DATA.open() as stream:
        return list(csv.DictReader(stream))


def fit_points(pressure, volume, *, objective="pressure", initial=(300.0, 4.0)) -> dict:
    """Fit K0 and K0' with fixed V0; source errors do not enter the objective."""
    pressure, volume = np.asarray(pressure), np.asarray(volume)

    def residual(parameters):
        k0, kp = parameters
        if objective == "volume":
            return Vinet(V0=29.47, K0=k0, K0_prime=kp).volume(pressure) - volume
        x = (volume / 29.47) ** (1 / 3)
        return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x)) - pressure

    fit = least_squares(
        residual,
        initial,
        bounds=([1, 0], [1000, 15]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    covariance = (
        np.linalg.inv(fit.jac.T @ fit.jac) * np.sum(fit.fun**2) / (len(volume) - 2)
    )
    errors = np.sqrt(np.diag(covariance))
    predicted = Vinet(V0=29.47, K0=fit.x[0], K0_prime=fit.x[1]).pressure(volume)
    published = get_eos_record(RECORD).pressure(volume)
    return {
        "observations": len(volume),
        "objective": objective + "_residuals",
        "parameters": dict(zip(("K0", "K0_prime"), fit.x.tolist())),
        "conditional_standard_errors": dict(zip(("K0", "K0_prime"), errors.tolist())),
        "conditional_parameter_correlation": float(covariance[0, 1] / np.prod(errors)),
        "pressure_rmse_gpa": float(np.sqrt(np.mean((predicted - pressure) ** 2))),
        "published_pressure_rmse_gpa": float(
            np.sqrt(np.mean((published - pressure) ** 2))
        ),
        "observed_pressure_range_gpa": [float(min(pressure)), float(max(pressure))],
        "solver_success": bool(fit.success),
        "within_published_quoted_errors": bool(
            abs(fit.x[0] - 358) <= 10 and abs(fit.x[1] - 4.8) <= 0.2
        ),
    }


def refit_figure10() -> dict:
    rows = load_points()
    fits = {}
    for name, series in {
        "rp01_yokoo_pt": {"rp01_yokoo_pt"},
        "all_yokoo_pt": {"rp01_yokoo_pt", "micro17_yokoo_pt"},
        "rp01_dewaele_pt_comparison": {"rp01_dewaele_pt"},
    }.items():
        selected = [r for r in rows if r["series"] in series]
        p = np.array([float(r["pressure_gpa"]) for r in selected])
        v = np.array([float(r["volume_a3"]) for r in selected])
        fits[name] = fit_points(p, v)
        if name == "rp01_yokoo_pt":
            fits["rp01_yokoo_pt_volume_objective"] = fit_points(
                p, v, objective="volume"
            )
            source = json.loads(
                DATA.with_name("rhenium-sakai-2018-figure10-source.json").read_text()
            )
            # Deliberately larger than vector/tick rounding: systematic +/-0.1 pt
            # center offsets. This is a sensitivity exercise, not a measured sigma.
            offsets = [
                fit_points(
                    p + dx * source["pressure_affine_map"][0],
                    v + dy * source["volume_affine_map"][0],
                )["parameters"]
                for dx, dy in product((-0.1, 0.1), repeat=2)
            ]
            fits["coordinate_offset_sensitivity"] = {
                "assumed_systematic_offset_pt": 0.1,
                "K0_range_gpa": [
                    min(r["K0"] for r in offsets),
                    max(r["K0"] for r in offsets),
                ],
                "K0_prime_range": [
                    min(r["K0_prime"] for r in offsets),
                    max(r["K0_prime"] for r in offsets),
                ],
            }
    return fits


def ledger_outcome(record: dict) -> dict:
    fits = refit_figure10()
    fit = fits["rp01_yokoo_pt"]
    parameters = []
    for name, source, error in (("K0", 358.0, 10.0), ("K0_prime", 4.8, 0.2)):
        value = fit["parameters"][name]
        difference = value - source
        similar = abs(difference) / source <= (0.15 if name == "K0" else 0.20)
        if name == "K0_prime":
            similar |= abs(difference) <= 1.0
        parameters.append(
            {
                "parameter": name,
                "published": source,
                "published_error": error,
                "refit": value,
                "refit_error": fit["conditional_standard_errors"][name],
                "difference": difference,
                "relative_difference": abs(difference) / source,
                "within_combined_2sigma": None,
                "similar": bool(similar),
                "within_reported_error": abs(difference) <= error,
            }
        )
    return {
        "status": "parity"
        if all(p["within_reported_error"] and p["similar"] for p in parameters)
        else "similar"
        if all(p["similar"] for p in parameters)
        else "parity_not_achieved",
        "parity_basis": "within_reported_parameter_errors",
        "dataset_identifiers": [DATASET],
        "observations": fit["observations"],
        "fit_kind": "diagnostic_vector_digitized_vinet",
        "objective": "pressure_residuals",
        "selection": "26 RP01 Yokoo-Pt gray diamonds only; source membership inferred, not explicitly enumerated",
        "free_parameters": ["K0", "K0_prime"],
        "fixed_parameters": ["V0"],
        "parameters": parameters,
        "rmse_gpa": fit["pressure_rmse_gpa"],
        "published_rmse_gpa": fit["published_pressure_rmse_gpa"],
        "observed_pressure_range_gpa": fit["observed_pressure_range_gpa"],
        "sensitivity": fits,
        "qualification": (
            "Numerical parity: 26 RP01 Yokoo-Pt Figure 10 points reproduce both "
            "coefficients within the quoted errors, with comparable conditional "
            "refit standard errors (10.23 GPa and 0.180). Parity is based on the "
            "reported parameter error widths, not a formal combined-two-sigma claim. "
            "Exact row selection, weights and source error confidence are unspecified. "
            "Including six Micro17 points shifts K0/K0' to about 327.1/5.373; this "
            "selection sensitivity is retained in the dedicated Sakai 2018 audit."
        ),
    }


def reference_pressure(volume: float, v0: float, k0: float, kp: float) -> float:
    """Standard Vinet expression evaluated independently of Peritheos."""
    x = (volume / v0) ** (1.0 / 3.0)
    return 3.0 * k0 * (1.0 - x) / x**2 * math.exp(1.5 * (kp - 1) * (1 - x))


def reproduce() -> dict:
    # Article pp. 10-11 (final pp. 116-117): coefficients and a separate
    # reported pressure evaluation. These are source constants, not fit data.
    volume, v0, k0, kp = 18.65, 29.47, 358.0, 4.8
    reported = 464.0
    independent = reference_pressure(volume, v0, k0, kp)
    actual = float(get_eos_record(RECORD).pressure(volume))
    # Half a last printed digit in each rounded quantity, not a confidence
    # interval and not independent sampling of the published parameter errors.
    rounded = [
        reference_pressure(v, a, b, c)
        for v, a, b, c in product(
            (volume - 0.005, volume + 0.005),
            (v0 - 0.005, v0 + 0.005),
            (k0 - 0.5, k0 + 0.5),
            (kp - 0.05, kp + 0.05),
        )
    ]
    return {
        "record_identifier": RECORD,
        "doi": "10.1080/08957959.2018.1448082",
        "source_locations": ["Article p. 10, below Figure 10", "Article p. 11"],
        "classification": "qualified_digitized_refit",
        "figure10_refits": refit_figure10(),
        "volume_angstrom3": volume,
        "reported_pressure_gpa": reported,
        "independent_vinet_pressure_gpa": independent,
        "catalog_pressure_gpa": actual,
        "catalog_minus_reported_gpa": actual - reported,
        "printed_rounding_pressure_bounds_gpa": [min(rounded), max(rounded)],
        "reported_value_within_rounding_bounds": min(rounded)
        <= reported
        <= max(rounded),
        "catalog_matches_independent_expression": math.isclose(
            actual, independent, rel_tol=1e-12
        ),
        "qualification": (
            "The reported maximum pressure lies outside the co-compression fit "
            "range. Rounding bounds are a compatibility check, not an uncertainty "
            "interval or a recovery of the authors' unrounded coefficients. "
            "Original fit observations and weights are unavailable."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    assert result["reported_value_within_rounding_bounds"]
    assert result["catalog_matches_independent_expression"]
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        return int(not OUTPUT.exists() or OUTPUT.read_text() != text)
    OUTPUT.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
