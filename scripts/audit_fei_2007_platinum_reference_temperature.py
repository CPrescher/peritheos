"""Test adding absolute Debye pressure to the unchanged Fei RT cold curve."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import least_squares

from scripts.reproduce_fei_2007_platinum import (
    PARAMETERS,
    observations,
    source_pressure,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data"
OUTPUT = DATA / "fei-2007-platinum-reference-temperature-sensitivity.json"


def reference_pressure(volume, coefficients):
    """Unreferenced 300 K vibrational pressure, excluding zero-point energy."""
    volume = np.asarray(volume, dtype=float)
    v0, _, _, gamma0, q, theta0 = coefficients
    gamma = gamma0 * (volume / v0) ** q
    theta = theta0 * (volume / v0) ** (-gamma)

    def energy(th):
        y = th / 300.0
        integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-11, epsrel=1e-11)[
            0
        ]
        return 9 * R * 300.0 * integral / y**3

    energy300 = np.array([energy(th) for th in theta.flat]).reshape(volume.shape)
    return gamma * energy300 / (volume * N_A * 1e-30 / 4) / 1e9


def fit(volume, temperature, target, cold, *, omit_reference, q_min=None):
    def unpack(x):
        c = cold.copy()
        c[3:5] = x
        return c

    def evaluate(c):
        p = source_pressure(volume, temperature, coefficients=c)
        return p + reference_pressure(volume, c) if omit_reference else p

    solutions = []
    for q_start in (0.05, 0.5, 2.0):
        result = least_squares(
            lambda x: evaluate(unpack(x)) - target,
            [cold[3], q_start],
            bounds=([-np.inf, q_min], [np.inf, np.inf])
            if q_min is not None
            else (-np.inf, np.inf),
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
        )
        assert result.success
        solutions.append(result)
    rmses = [float(np.sqrt(np.mean(r.fun**2))) for r in solutions]
    assert max(rmses) - min(rmses) < 1e-8
    result = solutions[1]
    c = unpack(result.x)
    assert np.array_equal(c[[0, 1, 2, 5]], cold[[0, 1, 2, 5]])
    return {
        "gamma0": float(c[3]),
        "q": float(c[4]),
        "q_lower_bound": q_min,
        "hot_rmse_gpa": float(np.sqrt(np.mean(result.fun**2))),
        "q_at_lower_bound": bool(q_min is not None and c[4] - q_min < 1e-7),
        "coefficients": c.tolist(),
        "solver_success": bool(result.success),
        "q_start_sensitivity_rmse_gpa": rmses,
        "added_pressure_at_V0_300K_gpa": float(reference_pressure(c[0], c))
        if omit_reference
        else 0.0,
    }


def main():
    staged = json.loads(
        (DATA / "fei-2007-platinum-fixed-rt-thermal-subset.json").read_text()
    )
    cold = np.array(PARAMETERS["platinum"])
    cold[2] = staged["variants"]["all_43_rt_rows"][
        "cold_parameters_fixed_in_thermal_stage"
    ]["K0_prime"]
    volume, temperature, target, paired, reduced = observations()
    hot = temperature > 300
    gold_volume = np.array([float(row["gold_volume_a3"]) for row in paired])
    gold_offset = reference_pressure(gold_volume, np.array(PARAMETERS["gold"]))
    wrong_au_target = target.copy()
    wrong_au_target[36:] = reduced + gold_offset
    cases = {}
    for label, omit, pressure_target in (
        ("correct_300K_reference", False, target),
        ("Pt_omits_300K_subtraction", True, target),
        ("Pt_and_Au_omit_300K_subtraction_in_hot_calibration", True, wrong_au_target),
    ):
        cases[label] = {
            "unconstrained_q": fit(
                volume[hot],
                temperature[hot],
                pressure_target[hot],
                cold,
                omit_reference=omit,
            ),
            "q_nonnegative": fit(
                volume[hot],
                temperature[hot],
                pressure_target[hot],
                cold,
                omit_reference=omit,
                q_min=0.0,
            ),
        }
    for name, result in cases["correct_300K_reference"].items():
        previous = staged["variants"]["all_43_rt_rows"]["thermal_fits"][name]
        assert abs(result["hot_rmse_gpa"] - previous["hot_rmse_gpa"]) < 1e-9

    source = json.loads(
        (DATA / "fei-2007-platinum-source-curve-check.json").read_text()
    )
    figure_checks = {}
    published = np.array(PARAMETERS["platinum"])
    for t in (300, 1473, 1873):
        rows = [r for r in source["points"] if r["temperature_k"] == t]
        v = np.array([r["volume_a3"] for r in rows])
        figure_p = np.array([r["figure_pressure_gpa"] for r in rows])
        wrong = source_pressure(v, t) + reference_pressure(v, published)
        figure_checks[str(t)] = {
            "rmse_gpa": float(np.sqrt(np.mean((wrong - figure_p) ** 2))),
            "max_absolute_difference_gpa": float(np.max(abs(wrong - figure_p))),
        }
    # A consistent conversion to a 0 K baseline leaves the original pressure
    # predictions unchanged: subtract E(300) from the cold curve first.
    cold_pressure = source_pressure(volume[hot], 300, coefficients=cold)
    absolute_thermal = (
        source_pressure(volume[hot], temperature[hot], coefficients=cold)
        - cold_pressure
        + reference_pressure(volume[hot], cold)
    )
    converted_cold = cold_pressure - reference_pressure(volume[hot], cold)
    identity_error = float(
        np.max(
            abs(
                converted_cold
                + absolute_thermal
                - source_pressure(volume[hot], temperature[hot], coefficients=cold)
            )
        )
    )
    assert identity_error < 1e-10
    report = {
        "hypothesis": "Treat the unchanged RT Vinet curve as a 0 K baseline: P=P_RT+gamma*E_D(T)/V_m, omitting subtraction of E_D(300). Debye energy excludes zero-point energy.",
        "correct_expression": "P=P_RT+gamma*(E_D(T)-E_D(300))/V_m",
        "fixed_cold_coefficients": cold.tolist(),
        "observations": int(sum(hot)),
        "objective": "Same 35 hot Fei (2004) rows and equal pressure weights as the fixed-RT diagnostic. No RT residuals enter these fits.",
        "cases": cases,
        "published_coefficients_added_pressure_V0_300K_gpa": float(
            reference_pressure(published[0], published)
        ),
        "published_coefficients_omitted_subtraction_vs_figure2": figure_checks,
        "consistent_0K_baseline_conversion_max_identity_error_gpa": identity_error,
        "qualification": "Counterfactual diagnostics, not evidence about the author's actual code. The combined Au/Pt case changes hot calibration targets while retaining the previously frozen RT cold coefficients. Omitting the subtraction changes RT predictions despite unchanged cold coefficients. Canonical EOS coefficients, source measurements and normal calibration targets remain unchanged.",
        "input_sha256": staged["input_sha256"],
    }
    OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "cases": cases,
                "figure_checks": figure_checks,
                "published_offset_at_V0": report[
                    "published_coefficients_added_pressure_V0_300K_gpa"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
