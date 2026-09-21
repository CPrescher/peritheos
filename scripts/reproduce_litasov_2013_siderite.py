"""Audit Litasov's siderite Table 2 without replacing published coefficients.

The sample calculations and regressions implement the equations independently
of Peritheos. The separately labelled gold comparison uses the catalog only to
test whether its corrected calibration is interchangeable with the 2013 scale.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATASET = "siderite_fe095mn005_litasov_2013_table2_pvt"
RECORD = "siderite_fe095mn005_litasov_2013_bm3"
THERMAL_RECORD = "siderite_fe095mn005_litasov_2013_joint_thermal_refit"
OUTPUT = ROOT / "docs/data/litasov-2013-siderite-reproduction.json"
NAMES = ["K0", "K0_prime", "alpha0", "alpha1", "dK_dT"]
SCALE = np.array([1.0, 1.0, 1e-5, 1e-8, 1.0])
TABLE = np.array([120.0, 3.57, 3.77, 0.06, -0.015])
ABSTRACT = np.array([120.0, 3.57, 3.57, 0.06, -0.015])


def observations():
    return np.genfromtxt(
        ROOT
        / "peritheos/data/datasets/siderite-fe095mn005-litasov-2013-table2-pvt.csv",
        delimiter=",",
        names=True,
    )


def bm3(volume, v0, k0, kp):
    """BM3 in conventional hexagonal cell angstrom cubed and GPa."""
    x = (v0 / np.asarray(volume)) ** (1 / 3)
    return 1.5 * k0 * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))


def thermal_pressure(parameters, volume, temperature):
    """Litasov (2007) Eqs. (1)-(4), cited by the 2013 Methods.

    Solver coordinates scale alpha0 by 1e-5 and alpha1 by 1e-8; alpha0 is
    the absolute-temperature intercept, not the expansivity at 300 K.
    """
    k0, kp, a0, a1, dk = parameters * SCALE
    dt = temperature - 300
    v0 = 293.4 * np.exp(a0 * dt + a1 * (temperature**2 - 300**2) / 2)
    return bm3(volume, v0, k0 + dk * dt, kp)


def metrics(residual):
    return {
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual_gpa": float(max(abs(residual))),
        "residuals_gpa": residual.tolist(),
    }


def fixed_volume_fit(volume, pressure, v0):
    """Exact unweighted pressure least squares in K0 and K0*(Kprime-4)."""
    x = (v0 / volume) ** (1 / 3)
    first = 1.5 * (x**7 - x**5)
    design = np.column_stack([first, first * 0.75 * (x**2 - 1)])
    linear = np.linalg.lstsq(design, pressure, rcond=None)[0]
    k0, kp = linear[0], 4 + linear[1] / linear[0]
    residual = bm3(volume, v0, k0, kp) - pressure
    covariance = np.linalg.inv(design.T @ design) * sum(residual**2) / (len(volume) - 2)
    transform = np.array([[1, 0], [-linear[1] / k0**2, 1 / k0]])
    covariance = transform @ covariance @ transform.T
    return {
        "parameters": {"K0": float(k0), "K0_prime": float(kp)},
        "estimated_standard_errors": dict(
            zip(NAMES[:2], np.sqrt(np.diag(covariance)).tolist())
        ),
        "estimated_covariance": covariance.tolist(),
        "covariance_scope": "OLS residual-variance estimate with n-2 degrees of freedom; not source covariance or absolute measurement errors.",
        **metrics(residual),
    }


def rt_reproduction():
    data = observations()
    data = data[data["temperature_k"] == 300]
    volume, pressure = data["volume_a3"], data["pressure_gpa"]
    result = fixed_volume_fit(volume, pressure, 293.4)
    return {
        "observations": len(data),
        "selected_rows": data["row"].astype(int).tolist(),
        "objective": "unweighted_pressure_residuals",
        "weights": "equal; the 0.1 GPa pressure bound is not a row standard deviation",
        "fixed_parameters": {"V0": 293.4, "Tr": 300},
        "published": metrics(bm3(volume, 293.4, 120, 3.57) - pressure),
        **result,
    }


def reproduce():
    data = observations()
    volume, pressure, temperature = (
        data[n] for n in ["volume_a3", "pressure_gpa", "temperature_k"]
    )
    rt = rt_reproduction()

    def fit(mask, cold=None, weights=None):
        v, p, t = volume[mask], pressure[mask], temperature[mask]
        weights = np.ones(len(v)) if weights is None else weights

        def unpack(x):
            return x if cold is None else np.r_[cold, x]

        def objective(x):
            return (thermal_pressure(unpack(x), v, t) - p) * weights

        start = TABLE if cold is None else TABLE[2:]
        result = least_squares(
            objective, start, xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=5000
        )
        covariance = (
            np.linalg.pinv(result.jac.T @ result.jac)
            * sum(result.fun**2)
            / (len(v) - len(start))
        )
        scale = SCALE if cold is None else SCALE[2:]
        covariance *= np.outer(scale, scale)
        return {
            "parameters": dict(zip(NAMES, (unpack(result.x) * SCALE).tolist())),
            "free_parameters": NAMES if cold is None else NAMES[2:],
            "estimated_covariance_free_parameters": covariance.tolist(),
            "estimated_standard_errors": np.sqrt(np.diag(covariance)).tolist(),
            "covariance_scope": "Residual-scaled local least-squares estimate; no absolute-sigma or source-covariance claim.",
            "solver_success": bool(result.success),
            "observations": len(v),
            **metrics(thermal_pressure(unpack(result.x), v, t) - p),
        }

    all_rows = temperature > 0
    hot = temperature > 300
    # Sensitivity only: ignores unknown row pressure errors; not the source weights.
    dv = 1e-3
    derivative = (
        thermal_pressure(TABLE, volume + dv, temperature)
        - thermal_pressure(TABLE, volume - dv, temperature)
    ) / (2 * dv)
    sigma_from_volume = abs(derivative) * data["volume_sigma_a3"]
    axes = {}
    cold = data[temperature == 300]
    for axis, v0, k0, kp in [("a", 103.4, 166, 14), ("c", 3633.8, 59, 2.7)]:
        fictive = cold[f"{axis}_a"] ** 3
        axes[axis] = {
            "fictive_reference_volume_a3": v0,
            "published_parameters": {"K0": k0, "K0_prime": kp},
            "published": metrics(bm3(fictive, v0, k0, kp) - cold["pressure_gpa"]),
            "refit": fixed_volume_fit(fictive, cold["pressure_gpa"], v0),
            "scope": "Axial fictive volume; not a physical unit-cell volume EOS. All 27 room-temperature rows, no exclusions.",
        }

    from peritheos import get_eos_record

    gold = get_eos_record("gold_sokolova_2013_holzapfel_4")
    gold_residual = gold.pressure(data["gold_volume_a3"], temperature) - pressure
    highest = data[np.argmax(pressure)]
    benchmark = float(bm3(highest["volume_a3"], 293.4, 120, 3.57))
    return {
        "format": "peritheos.literature-reproduction",
        "doi": "10.1016/j.pepi.2013.07.011",
        "observations": len(data),
        "run_counts": {str(i): int(sum(data["run"] == i)) for i in range(1, 5)},
        "max_abs_hexagonal_volume_rounding_difference_a3": float(
            max(abs(np.sqrt(3) / 2 * data["a_a"] ** 2 * data["c_a"] - volume))
        ),
        "rt_refit": rt,
        "independent_source_benchmark": {
            "selection": "highest tabulated pressure, specified before comparison",
            "row": int(highest["row"]),
            "temperature_k": float(highest["temperature_k"]),
            "volume_a3": float(highest["volume_a3"]),
            "reported_pressure_gpa": float(highest["pressure_gpa"]),
            "predicted_pressure_gpa": benchmark,
            "residual_gpa": benchmark - highest["pressure_gpa"],
            "tolerance_gpa": 0.1,
            "tolerance_basis": "Table 2 pressure-uncertainty upper bound; not a fitted tolerance or a row sigma.",
        },
        "thermal_published_table": metrics(
            thermal_pressure(TABLE, volume, temperature) - pressure
        ),
        "thermal_published_abstract": metrics(
            thermal_pressure(ABSTRACT, volume, temperature) - pressure
        ),
        "thermal_joint_refit": fit(all_rows),
        "thermal_volume_uncertainty_weight_sensitivity": fit(
            all_rows, weights=1 / sigma_from_volume
        ),
        "thermal_staged_published_cold_sensitivity": fit(hot, cold=TABLE[:2]),
        "thermal_staged_refitted_cold_sensitivity": fit(
            hot, cold=np.array(list(rt["parameters"].values()))
        ),
        "fitting_qualification": "Section 3 explicitly varies all five coefficients jointly. The abstract suggests a staged calculation, included as sensitivity only. Source objective, weights, covariance and any unprinted exclusions are unknown; these diagnostics do not establish original-fit parity.",
        "axial_diagnostics": axes,
        "gold_calibration_diagnostic": {
            "related_record": gold.identifier,
            "scope": "Catalog uses the corrected 2016 realization of Sokolova 2013; this is a compatibility diagnostic, not an independent reconstruction of the exact 2013 pressure scale.",
            "recalibration_applied": False,
            **metrics(gold_residual),
        },
    }


def ledger_outcome(record):
    if record["identifier"] == THERMAL_RECORD:
        result = reproduce()["thermal_joint_refit"]
        stored = {
            **record["eos"]["parameters"],
            **record["thermal"]["parameters"],
        }
        parameters = []
        for index, name in enumerate(NAMES):
            expected, fitted = stored[name], result["parameters"][name]
            error = result["estimated_standard_errors"][index]
            parameters.append(
                {
                    "parameter": f"rt_eos.{name}" if index < 2 else name,
                    "published": expected,
                    "published_error": (
                        record["parameter_errors"][name]
                        if index < 2
                        else record["thermal"]["parameter_errors"][name]
                    ),
                    "refit": fitted,
                    "refit_error": error,
                    "difference": fitted - expected,
                    "relative_difference": abs(fitted - expected) / abs(expected),
                    "within_combined_2sigma": bool(abs(fitted - expected) <= 2 * error),
                    "similar": bool(
                        np.isclose(
                            fitted,
                            expected,
                            # The small alpha1 coefficient is especially
                            # sensitive to finite-difference optimizer drift.
                            rtol=2e-4 if name == "alpha1" else 5e-6,
                            atol=SCALE[index] * 1e-9,
                        )
                    ),
                }
            )
        data = observations()
        stored_scaled = np.array([stored[name] for name in NAMES]) / SCALE
        stored_pressure = thermal_pressure(
            stored_scaled, data["volume_a3"], data["temperature_k"]
        )
        residual = stored_pressure - data["pressure_gpa"]
        fitted_scaled = np.array([result["parameters"][name] for name in NAMES]) / SCALE
        curve_matches = np.allclose(
            stored_pressure,
            thermal_pressure(fitted_scaled, data["volume_a3"], data["temperature_k"]),
            rtol=0,
            atol=5e-6,
        )
        return {
            "status": "parity"
            if all(p["similar"] for p in parameters) and curve_matches
            else "parity_not_achieved",
            "parity_basis": "reproduction_of_stored_peritheos_refit",
            "qualification": "Reproduces the separately labelled Peritheos refit, not either published thermal coefficient set. All 111 printed states enter an unweighted joint pressure-residual fit; V0 and Tr are fixed. The ledger's published column contains the stored refit values for this comparison. Formal covariance excludes fixed V0 uncertainty and calibration systematics. Both source thermal variants remain non-executable evidence.",
            "dataset_identifiers": [DATASET],
            "observations": 111,
            "selected_rows": list(range(1, 112)),
            "fit_kind": "joint_thermal_reference_state_refit",
            "objective": "unweighted_pressure_residuals",
            "absolute_sigma": False,
            "free_parameters": [p["parameter"] for p in parameters],
            "parameters": parameters,
            "rmse_gpa": result["rmse_gpa"],
            "published_rmse_gpa": metrics(residual)["rmse_gpa"],
            "observed_pressure_range_gpa": [0, 33.01],
            "observed_temperature_range_k": [300, 1673],
            "solver_success": result["solver_success"],
        }
    result = rt_reproduction()
    parameters = []
    for name, fitted in result["parameters"].items():
        published = record["eos"]["parameters"][name]
        error = record["parameter_errors"][name]
        refit_error = result["estimated_standard_errors"][name]
        parameters.append(
            {
                "parameter": name,
                "published": published,
                "published_error": error,
                "refit": fitted,
                "refit_error": refit_error,
                "difference": fitted - published,
                "relative_difference": abs(fitted - published) / abs(published),
                "within_combined_2sigma": None,
                "similar": abs(fitted - published) <= error,
            }
        )
    return {
        "status": "similar"
        if all(p["similar"] for p in parameters)
        else "parity_not_achieved",
        "qualification": "Point estimates fall within the published errors, but source weighting and covariance are unspecified. All 27 printed 300 K rows are retained; the source RMS range is not recovered. The conflicting thermal coefficient is withheld. See literature-reproductions/litasov-2013-siderite.md.",
        "dataset_identifiers": [DATASET],
        "observations": 27,
        "selected_rows": result["selected_rows"],
        "fit_kind": "fixed_V0_linearized_bm3",
        "objective": "unweighted_pressure_residuals",
        "absolute_sigma": False,
        "free_parameters": NAMES[:2],
        "parameters": parameters,
        "rmse_gpa": result["rmse_gpa"],
        "published_rmse_gpa": result["published"]["rmse_gpa"],
        "observed_pressure_range_gpa": [0, 33.01],
        "observed_temperature_range_k": [300, 300],
        "solver_success": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = json.dumps(reproduce(), indent=2, allow_nan=False) + "\n"
    if args.check:
        if OUTPUT.read_text() != text:
            raise SystemExit("Siderite reproduction is stale")
    else:
        OUTPUT.write_text(text)
