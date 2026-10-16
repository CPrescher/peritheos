#!/usr/bin/env python3
"""Resolve printed-data discrepancies without modifying observations or EOS fits.

Benchmark inversions below reconstruct calculated tables only. Excluding a
checkpoint diagnoses consistency; it never changes an observational regression.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, R
from scipy.optimize import brentq, minimize

from scripts.investigate_rubidium_halides_farla_2025 import (
    AUTHOR_THETA,
    calibration_replay,
    cscl_pressure,
)
from scripts.reproduce_rubidium_halides_farla_2025 import (
    PREFIX,
    SOURCE,
    read_rows,
    reference_pressure,
)

ROOT = Path(__file__).resolve().parents[1]
BURNMAN_REPLAY = ROOT / "docs/data" / f"{PREFIX}-burnman-replay.json"


def cscl_input_intervals():
    """Invert pressure/volume rounding boxes, keeping uncertainty separate."""
    rows = [r for r in read_rows("rbcl", "observations") if r["source_table"] != "S1"]
    outliers = calibration_replay()["gui_defaults_legacy_na"][
        "outliers_beyond_rounding"
    ]
    intervals = []
    temperature_offsets = []
    for row in rows:
        v, t, p = [
            float(row[k])
            for k in ("cscl_volume_a3", "temperature_k", "cscl_pressure_gpa")
        ]
        # P increases with a0 and T, and decreases with V over this domain.
        low = brentq(
            lambda a: cscl_pressure(v - 0.0005, t, reference_a0=a) - (p - 0.005),
            4.10,
            4.14,
        )
        high = brentq(
            lambda a: cscl_pressure(v + 0.0005, t, reference_a0=a) - (p + 0.005),
            4.10,
            4.14,
        )
        intervals.append((low, high))
        low_t = brentq(lambda x: cscl_pressure(v - 0.0005, x) - (p - 0.005), 100, 2200)
        high_t = brentq(lambda x: cscl_pressure(v + 0.0005, x) - (p + 0.005), 100, 2200)
        temperature_offsets.append((low_t - t, high_t - t))
    details = []
    for outlier in outliers:
        row = rows[outlier["observation_index"] - 1]
        v, t, p = [
            float(row[k])
            for k in ("cscl_volume_a3", "temperature_k", "cscl_pressure_gpa")
        ]
        required_v = brentq(lambda x: cscl_pressure(x, t) - p, v - 1, v + 1)
        t_low = brentq(lambda x: cscl_pressure(v - 0.0005, x) - (p - 0.005), 100, 2200)
        t_high = brentq(lambda x: cscl_pressure(v + 0.0005, x) - (p + 0.005), 100, 2200)
        minimum_shift = max(0, t_low - t, t - t_high)
        a0_range = intervals[outlier["observation_index"] - 1]
        details.append(
            {
                **outlier,
                "inverse_volume_a3": required_v,
                "inverse_volume_shift_a3": required_v - v,
                "inverse_volume_shift_over_reported_esd": (required_v - v)
                / float(row["cscl_volume_uncertainty_a3"]),
                "inverse_lattice_shift_angstrom": required_v ** (1 / 3) - v ** (1 / 3),
                "temperature_shift_interval_after_PV_rounding_k": [
                    t_low - t,
                    t_high - t,
                ],
                "minimum_absolute_temperature_shift_k": minimum_shift,
                "compatible_with_target_control_plus_minus_1k": minimum_shift <= 1,
                "reference_a0_interval_at_printed_T_angstrom": list(a0_range),
            }
        )
    lower = max(low for low, _ in intervals)
    upper = min(high for _, high in intervals)
    return {
        "distinct_rows": len(rows),
        "source_control_tolerance_k": 1,
        "source_control_location": "Official supplement Figure S1 caption, rendered page 4, lines 71-72",
        "nominal_temperature_half_integer_rounding_k": 0.5,
        "outliers": details,
        "outliers_remaining_with_half_integer_T_rounding": sum(
            d["minimum_absolute_temperature_shift_k"] > 0.5 for d in details
        ),
        "outliers_remaining_with_target_control_tolerance": sum(
            not d["compatible_with_target_control_plus_minus_1k"] for d in details
        ),
        "common_a0_intersection_angstrom": [lower, upper],
        "common_a0_possible_at_printed_T": lower <= upper,
        "common_temperature_offset_intersection_k": [
            max(a for a, _ in temperature_offsets),
            min(b for _, b in temperature_offsets),
        ],
        "common_temperature_offset_possible": max(a for a, _ in temperature_offsets)
        <= min(b for _, b in temperature_offsets),
        "a0_lower_bound_row": int(np.argmax([a for a, _ in intervals]) + 1),
        "a0_upper_bound_row": int(np.argmin([b for _, b in intervals]) + 1),
        "qualification": (
            "Inverse inputs are diagnostic requirements, not recovered measurements. "
            "Reported 30 K temperature uncertainty is distinct from rounding and the stated "
            "1 K target-control tolerance. No source temperature, volume or pressure is corrected. "
            "The common-a0 test holds the recovered model coefficients and printed temperatures fixed."
        ),
    }


def rbcl_minimax(
    *,
    exclude=(),
    theta=177.04,
    reference_rounding=False,
    normalize_to_fitted_v0=True,
    eosfit_r=False,
    replacement=None,
):
    """Minimize maximum ratio discrepancy using exact pressure interval constraints.

    Four coefficient bounds are half the published last digit. Optionally
    allow the printed V0/theta half digits too. Successful convergence is
    numerical evidence, not a global infeasibility certificate.
    """
    molar, _, coefs = SOURCE["rbcl"]
    rows = read_rows("rbcl", "benchmarks")
    ratio, temperature, pressure = [
        np.array([float(r[k]) for r in rows])
        for k in ("volume_ratio", "temperature_k", "pressure_gpa")
    ]
    if replacement:
        for index, value in replacement.items():
            ratio[index] = value
    center = np.r_[coefs, molar, theta]
    widths = np.r_[
        [0.005] * 4,
        0.005 if reference_rounding else 0,
        0.5 if reference_rounding else 0,
    ]
    variable = widths > 0
    keep = np.array([i not in exclude for i in range(len(ratio))])

    def parameters(z):
        x = center.copy()
        x[variable] += z[:-1] * widths[variable]
        return x

    def evaluate(z, r, t=temperature):
        x = parameters(z)
        denominator = x[4] if normalize_to_fitted_v0 else molar
        volume = r * denominator * 1e24 / Avogadro
        value = reference_pressure(volume, t, x[4], x[5], x[:4])
        if eosfit_r:
            cold = reference_pressure(volume, np.full_like(t, 300), x[4], x[5], x[:4])
            value = cold + (value - cold) * 8.314 / R
        return value

    start = np.zeros(sum(variable) + 1)
    start[-1] = 2
    scale = (evaluate(start, ratio - 5e-5) - evaluate(start, ratio + 5e-5)) / 2

    def constraints(z):
        half_width = z[-1] * 5e-5
        return np.r_[
            ((evaluate(z, ratio - half_width) - pressure) / scale)[keep],
            ((pressure - evaluate(z, ratio + half_width)) / scale)[keep],
        ]

    fit = minimize(
        lambda z: float(z[-1]),
        start,
        method="SLSQP",
        bounds=[(-1, 1)] * sum(variable) + [(0, 20)],
        constraints=[{"type": "ineq", "fun": constraints}],
        options={"ftol": 1e-11, "maxiter": 500},
    )
    candidate = parameters(fit.x)
    calculated = np.array(
        [
            brentq(lambda r: evaluate(fit.x, r, t) - p, 0.2, 1.2, xtol=2e-14)
            for p, t in zip(pressure, temperature)
        ]
    )
    difference = calculated - ratio
    return {
        "retained_rows": int(sum(keep)),
        "excluded_zero_based_indices": list(exclude),
        "reference_rounding_allowed": reference_rounding,
        "normalization": "same_fitted_v0" if normalize_to_fitted_v0 else "printed_v0",
        "gas_constant_j_mol_k": 8.314 if eosfit_r else R,
        "candidate_parameter_order": [
            "K0_gpa",
            "K0_prime",
            "gamma0",
            "q",
            "molar_v0_cm3_mol",
            "theta0_k",
        ],
        "candidate_parameters": candidate.tolist(),
        "optimizer_success": bool(fit.success),
        "optimizer_message": fit.message,
        "minimum_scaled_constraint": float(min(constraints(fit.x))),
        "minimax_half_last_digits": float(fit.x[-1]),
        "maximum_retained_absolute_ratio_difference": float(max(abs(difference[keep]))),
        "all_retained_within_ratio_rounding": bool(
            max(abs(difference[keep])) <= 5e-5 + 1e-9
        ),
        "states": [
            {
                "index": i,
                "pressure_gpa": float(pressure[i]),
                "temperature_k": float(temperature[i]),
                "source_ratio": float(rows[i]["volume_ratio"]),
                "diagnostic_target_ratio": float(ratio[i]),
                "calculated_ratio": float(calculated[i]),
                "calculated_ratio_rounded": f"{calculated[i]:.4f}",
                "difference_from_source": float(
                    calculated[i] - float(rows[i]["volume_ratio"])
                ),
                "retained": bool(keep[i]),
            }
            for i in range(len(rows))
        ],
        "qualification": "Calculated-table consistency diagnostic only; no original coefficient recovery, observational refit, catalog replacement or global infeasibility proof.",
    }


def rbcl_omitted_prediction_range(localized):
    """Numerically bound the omitted ratio while the other 53 obey rounding."""
    index = localized["excluded_zero_based_indices"][0]
    rows = read_rows("rbcl", "benchmarks")
    ratio, temperature, pressure = [
        np.array([float(r[k]) for r in rows])
        for k in ("volume_ratio", "temperature_k", "pressure_gpa")
    ]
    molar, _, center = SOURCE["rbcl"]
    center = np.array(center)
    volume0 = molar * 1e24 / Avogadro
    keep = np.arange(len(rows)) != index

    def pressure_at(z, r, t=temperature):
        return reference_pressure(
            r * volume0, t, molar, AUTHOR_THETA["rbcl"], center + z * 0.005
        )

    scale = (
        pressure_at(np.zeros(4), ratio - 5e-5) - pressure_at(np.zeros(4), ratio + 5e-5)
    ) / 2

    def constraints(z):
        return np.r_[
            ((pressure_at(z, ratio - 5e-5) - pressure) / scale)[keep],
            ((pressure - pressure_at(z, ratio + 5e-5)) / scale)[keep],
        ]

    def predicted(z):
        return brentq(
            lambda r: pressure_at(z, r, temperature[index]) - pressure[index],
            0.2,
            1.2,
            xtol=2e-14,
        )

    start = (np.array(localized["candidate_parameters"][:4]) - center) / 0.005
    extrema = []
    for sign in (1, -1):
        fit = minimize(
            lambda z: sign * (predicted(z) - ratio[index]) / 5e-5,
            start,
            method="SLSQP",
            bounds=[(-1, 1)] * 4,
            constraints=[{"type": "ineq", "fun": constraints}],
            options={"ftol": 1e-11, "maxiter": 500},
        )
        extrema.append(
            {
                "optimizer_success": bool(fit.success),
                "minimum_scaled_constraint": float(min(constraints(fit.x))),
                "ratio": float(predicted(fit.x)),
                "candidate_coefficients": (center + fit.x * 0.005).tolist(),
            }
        )
    return {
        "minimization": extrema[0],
        "maximization": extrema[1],
        "printed_ratio_interval": [
            float(ratio[index] - 5e-5),
            float(ratio[index] + 5e-5),
        ],
        "qualification": "Numerical extrema under the four coefficient half-digit bounds, fixed V0 and precise author theta; not certified global interval bounds.",
    }


def rbcl_localization():
    cases = {
        "fixed_author_theta": rbcl_minimax(),
        "fixed_printed_theta": rbcl_minimax(theta=177),
        "eosfit_gas_constant": rbcl_minimax(eosfit_r=True),
        "reference_half_digits_same_normalization": rbcl_minimax(
            theta=177, reference_rounding=True
        ),
        "reference_half_digits_printed_normalization": rbcl_minimax(
            theta=177, reference_rounding=True, normalize_to_fitted_v0=False
        ),
    }
    omitted = [rbcl_minimax(exclude=(i,)) for i in range(54)]
    compatible = [
        i for i, fit in enumerate(omitted) if fit["all_retained_within_ratio_rounding"]
    ]
    candidates = [
        {
            "omitted_index": i,
            "half_last_digits": f["minimax_half_last_digits"],
            "optimizer_success": f["optimizer_success"],
            "retained_rounding_match": f["all_retained_within_ratio_rounding"],
        }
        for i, f in enumerate(omitted)
    ]
    result = {
        "full_table_cases": cases,
        "leave_one_out": candidates,
        "single_omissions_with_compatible_remaining_table": compatible,
    }
    if compatible:
        index = compatible[0]
        result["localized_candidate"] = omitted[index]
        result["omitted_prediction_range"] = rbcl_omitted_prediction_range(
            omitted[index]
        )
        changed = float(omitted[index]["states"][index]["calculated_ratio_rounded"])
        result["one_cell_replacement_diagnostic"] = rbcl_minimax(
            replacement={index: changed}
        )
    result["literal_published_table_parity"] = False
    result["source_table_modified"] = False
    result["qualification"] = (
        "A compatible 53-row reconstruction localizes a discrepancy. A one-cell alternative is a hypothesis, not an author correction or original-fit reproduction."
    )
    return result


def reproduce():
    replay = json.loads(BURNMAN_REPLAY.read_text())
    rows = [r for r in read_rows("rbcl", "observations") if r["source_table"] != "S1"]
    reference = cscl_pressure(
        np.array([float(r["cscl_volume_a3"]) for r in rows]),
        np.array([float(r["temperature_k"]) for r in rows]),
    )
    replay["maximum_absolute_reference_difference_gpa"] = float(
        max(abs(reference - replay["pressures_gpa"]))
    )
    return {
        "audit_date": "2026-10-16",
        "source_doi": "10.1063/5.0248905",
        "author_contact": "none",
        "author_fit_parity": False,
        "original_input_files_recovered": False,
        "historical_software_execution": replay,
        "cscl_input_diagnostics": cscl_input_intervals(),
        "rbcl_benchmark_diagnostics": rbcl_localization(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(reproduce(), indent=2) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result)


if __name__ == "__main__":
    main()
