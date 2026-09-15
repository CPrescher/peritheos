"""Reproduce Zha's isochore fit and audit Xian's printed thermal surrogate.

Reference calculations use the printed expressions directly, independently of
Peritheos. Derived isotherms are never treated as primary observations.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/zha-2004-xian-2022-rhenium-reproduction.json"
RECORD = "rhenium_zha_2004_bm3_log_thermal"
DATASET = "rhenium_zha_2004_table2_paired_pvt"


def rows(table: str) -> list[dict[str, str]]:
    with (DATA / f"rhenium-zha-2004-{table}.csv").open(newline="") as stream:
        return list(csv.DictReader(stream))


def bm3(volume, v0=29.4087, k0=360.0, kp=4.5):
    """Zha Equation (3), conventional-cell volume and GPa."""
    x = (v0 / np.asarray(volume)) ** (1 / 3)
    return 1.5 * k0 * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))


def pressure(volume, temperature):
    """Zha Equation (6) with Table III coefficients."""
    return bm3(volume) + (0.00776 - 0.00815 * np.log(29.4087 / volume)) * (
        temperature - 300
    )


def gold_calibration_diagnostic(observations: list[dict[str, str]]) -> dict:
    """Compare Zha Table IIa with Anderson's final and exploratory choices.

    The alternative pairs are Section III thermal-consistency trials, not new
    calibrated EOS records. Rounding envelopes are distinct from measurement
    uncertainties and are never used as regression weights.
    """
    volume = np.array([float(r["au_volume_a3"]) for r in observations])
    temperature = np.array([float(r["temperature_k"]) for r in observations])
    observed = np.array([float(r["pressure_gpa"]) for r in observations])
    mean_a = np.array([float(r["au_a_mean_angstrom"]) for r in observations])
    averaged_a = np.array(
        [
            np.mean([float(r[f"au_a_{hkl}_angstrom"]) for hkl in [111, 200, 220, 311]])
            for r in observations
        ]
    )

    def au_pressure(
        v=volume, t=temperature, va=67.847, k=166.65, kp=5.4823, a=0.00714, b=-0.0115
    ):
        return bm3(v, va, k, kp) + (a + b * np.log(va / v)) * (t - 300)

    canonical = au_pressure()
    residual = canonical - observed
    alternatives = []
    for kp, b in [(6.39, -0.0052), (6.12, -0.0071), (5.5, -0.0115), (5.21, -0.0135)]:
        delta = au_pressure(kp=kp, b=b) - observed
        alternatives.append(
            {
                "K0_prime": kp,
                "dK_dT_V_gpa_per_k": b,
                "pressure_residuals_gpa": delta.tolist(),
                "rmse_gpa": float(np.sqrt(np.mean(delta**2))),
                "max_abs_residual_gpa": float(max(abs(delta))),
            }
        )
    # Corners of half-last-digit intervals; use the lattice-parameter precision
    # (larger volume interval than the printed volume precision).
    corners = []
    for sv, st, sva, sk, skp, sa, sb in product([-1, 1], repeat=7):
        corners.append(
            au_pressure(
                v=(mean_a + sv * 0.00005) ** 3,
                t=temperature + st * 0.05,
                va=67.847 + sva * 0.0005,
                k=166.65 + sk * 0.005,
                kp=5.4823 + skp * 0.00005,
                a=0.00714 + sa * 0.000005,
                b=-0.0115 + sb * 0.00005,
            )
        )
    lower = np.min(corners, axis=0) - (observed + 0.005)
    upper = np.max(corners, axis=0) - (observed - 0.005)
    return {
        "related_record": "gold_anderson_1989_bm3_1",
        "source_reference_volume_a3": 67.847,
        "temperature_k": temperature.tolist(),
        "printed_pressure_gpa": observed.tolist(),
        "calculated_pressure_gpa": canonical.tolist(),
        "pressure_residuals_gpa": residual.tolist(),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual_gpa": float(max(abs(residual))),
        "catalog_default_V0_pressure_residuals_gpa": (
            au_pressure(va=67.79) - observed
        ).tolist(),
        "max_abs_volume_minus_printed_mean_a_cubed_a3": float(
            max(abs(volume - mean_a**3))
        ),
        "max_abs_pressure_change_using_mean_of_four_a_gpa": float(
            max(abs(au_pressure(v=averaged_a**3) - canonical))
        ),
        "rounding_corner_residual_interval_gpa": np.column_stack(
            [lower, upper]
        ).tolist(),
        "rounding_scope": "Half last printed digit for a, T, Va, K0, Kprime, alpha_KT, dK_dT_V and tabulated P; corner evaluations only, not experimental uncertainty.",
        "anderson_section3_exploratory_pairs": alternatives,
        "finding": "The final Anderson Equation (29) with Zha's Va does not reproduce all printed pressures. Source rounding and lattice averaging do not explain the difference. The Section III pair Kprime=6.39, dK_dT_V=-0.0052 is closer but still not an exact reproduction and is not identified by Zha as its chosen realization. Calibration remains partially resolved.",
    }


def reproduce() -> dict:
    observations = rows("table2-paired-pvt")
    volume = np.array([float(r["volume_a3"]) for r in observations])
    temperature = np.array([float(r["temperature_k"]) for r in observations])
    observed = np.array([float(r["pressure_gpa"]) for r in observations])
    dt = temperature - 300
    design = np.column_stack([np.ones(8), np.log(29.4087 / volume)])
    slopes = (observed - bm3(volume)) / dt
    coefficients = np.linalg.lstsq(design, slopes, rcond=None)[0]
    direct = np.linalg.lstsq(design * dt[:, None], observed - bm3(volume), rcond=None)[
        0
    ]
    residual = pressure(volume, temperature) - observed
    fitted_residual = bm3(volume) + (design @ coefficients) * dt - observed
    source = np.array([0.00776, -0.00815])
    names = ["alpha_KT_ref", "dK_dT_V"]
    comparison_volume = 29.4087 / np.array([1.0, 1.1, 1.2])
    free_checks = []
    for r in rows("table4-isotherms"):
        if r["branch"] != "free":
            continue
        t = float(r["temperature_k"])
        difference = pressure(comparison_volume, t) - bm3(
            comparison_volume,
            float(r["v0_a3"]),
            float(r["k0_gpa"]),
            float(r["k0_prime"]),
        )
        free_checks.append(
            {"temperature_k": t, "max_abs_difference_gpa": float(max(abs(difference)))}
        )
    fixed = {
        float(r["temperature_k"]): r
        for r in rows("table4-isotherms")
        if r["branch"] == "fixed_V0"
    }
    grid_errors = []
    continuous_errors = []
    for r in rows("table5-pressure-grid"):
        t, p = float(r["temperature_k"]), float(r["pressure_gpa"])
        v = 29.4087 * (1 - float(r["compression"]))
        fit = fixed.get(t)
        predicted = (
            bm3(v)
            if fit is None
            else bm3(
                v, float(fit["v0_a3"]), float(fit["k0_gpa"]), float(fit["k0_prime"])
            )
        )
        grid_errors.append(float(predicted - p))
        continuous_errors.append(float(pressure(v, t) - p))
    source_xian = json.loads(
        (ROOT / "docs/data/xian-2022-rhenium-source-audit.json").read_text()
    )
    polynomial = source_xian["printed_surrogate"]
    xian_checks = []
    for t in [300.0, 1000.0, 2000.0, 3200.0]:
        v0 = float(
            np.polynomial.polynomial.polyval(
                t, polynomial["specific_volume_cm3_g_coefficients"]
            )
        )
        k0 = float(
            np.polynomial.polynomial.polyval(t, polynomial["K0_gpa_coefficients"])
        )
        kp = float(
            np.polynomial.polynomial.polyval(t, polynomial["K0_prime_coefficients"])
        )
        _, v1, v2 = polynomial["specific_volume_cm3_g_coefficients"]
        x = (0.047893 / v0) ** (1 / 3)
        p = 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))
        xian_checks.append(
            {
                "temperature_k": t,
                "specific_volume_cm3_g": v0,
                "zero_pressure_expansivity_per_k": (v1 + 2 * v2 * t) / v0,
                "pressure_at_printed_cold_V0_gpa": float(p),
            }
        )
    return {
        "zha_isochore_refit": {
            "observations": 8,
            "objective": "unweighted residuals of (P-P300)/(T-300), the staged isochore construction in Section III",
            "parameters": dict(zip(names, coefficients.tolist())),
            "relative_coefficient_differences": dict(
                zip(names, (abs(coefficients - source) / abs(source)).tolist())
            ),
            "published_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "published_pressure_max_abs_residual_gpa": float(max(abs(residual))),
            "refit_pressure_rmse_gpa": float(np.sqrt(np.mean(fitted_residual**2))),
            "pressure_residuals_gpa": residual.tolist(),
            "direct_pressure_objective_sensitivity": dict(zip(names, direct.tolist())),
        },
        "zha_table4_free_branch_checks": free_checks,
        "zha_gold_calibration_diagnostic": gold_calibration_diagnostic(observations),
        "zha_table5": {
            "states": 147,
            "fixed_V0_branch_max_abs_residual_gpa": max(abs(v) for v in grid_errors),
            "continuous_eq6_max_abs_difference_gpa": max(
                abs(v) for v in continuous_errors
            ),
        },
        "xian_literal_equations_37_40": xian_checks,
    }


def ledger_outcome(record: dict) -> dict:
    result = reproduce()["zha_isochore_refit"]
    parameters = []
    for name, value in result["parameters"].items():
        published = record["thermal"]["parameters"][name]
        parameters.append(
            {
                "parameter": name,
                "published": published,
                "published_error": None,
                "refit": value,
                "refit_error": None,
                "difference": value - published,
                "relative_difference": abs(value - published) / abs(published),
                "within_combined_2sigma": None,
                "similar": abs(value - published) / abs(published) < 0.0005,
            }
        )
    return {
        "status": "similar"
        if all(p["similar"] for p in parameters)
        else "parity_not_achieved",
        "reason": "Both thermal point estimates recover within 0.05%; source errors are unavailable. The source staged isochore objective differs from direct pressure least squares. See literature-reproductions/zha-2004-rhenium.md.",
        "dataset_identifiers": [DATASET],
        "observations": 8,
        "fit_kind": "isochore_slope_linear_regression",
        "objective": result["objective"],
        "free_parameters": list(result["parameters"]),
        "parameters": parameters,
        "rmse_gpa": result["refit_pressure_rmse_gpa"],
        "published_rmse_gpa": result["published_pressure_rmse_gpa"],
        "observed_pressure_range_gpa": [6.41, 8.47],
        "observed_temperature_range_k": [1380.3, 1914.5],
        "solver_success": True,
    }


if __name__ == "__main__":
    OUTPUT.write_text(json.dumps(reproduce(), indent=2, allow_nan=False) + "\n")
