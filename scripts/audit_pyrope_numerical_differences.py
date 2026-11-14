"""Explain pyrope replay differences using frozen native numerical evidence.

Repository-only: python -m scripts.audit_pyrope_numerical_differences.
The compatibility controls are diagnostics, not changes to library defaults.
Native kernel values were captured from the version-locked external console;
see numerics-2026-11-14-evidence.json for source, binary and capture provenance.
"""

from __future__ import annotations

import csv
import json
from itertools import product
from unittest.mock import patch

import numpy as np

import scripts.refit_pyrope_angel_2022_subset as model
from scripts.audit_pyrope_author_inputs import DATA, NAMES, PREFIX

STEM = "numerics-2026-11-14"
SCALE = np.array([170, 4.5, 770, 1.2, 1])
PRECISE = model.pressure_at_observations
LOCAL_SIGMA = model.effective_sigma


def inputs():
    audit = json.loads((DATA / (PREFIX + "-audit.json")).read_text())
    # Prepared input coordinates, not calculated checkpoints from the old fit.
    keys = ["dataset", "row", "kind", "P", "T", "Y", "sP", "sT", "sY", "scale_index"]
    rows = [{k: row[k] for k in keys} for row in audit["independent_fits"][0]["rows"]]
    with (DATA / (STEM + "-native-kernel.csv")).open() as f:
        captured = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]
    evidence = json.loads((DATA / (STEM + "-evidence.json")).read_text())
    return rows, model.arrays(rows), captured, evidence, audit


def coarse_pressure(c, d, *, y=None, t=None, seed="observed", tight=False):
    """Double-precision structural replay of CFML Get_Pressure_K stopping.

    K(P,T) here uses a tightly solved volume, whereas CFML uses its native
    spline/precision path. This is deliberately not called a bitwise clone.
    """
    y = d["Y"] if y is None else y
    t = d["T"] if t is None else t
    pth, cv = model.thermal(t, c)
    scales = np.array([c[i] if i is not None else 1 for i in d["scale_index"]])
    x = y / scales
    pressure = model.cold(x, c)[0] + pth
    k = d["KS"]
    if not np.any(k):
        return pressure
    target = x[k]
    slope = c[3] ** 2 * cv[k] * t[k] / (1000 * model.V0**2)
    pt = pth[k]

    def modulus(pp):
        v = model.V0 * (1 + c[1] * (pp - pt) / c[0]) ** (-1 / c[1])
        for _ in range(20):
            pc, kt = model.cold(v, c)
            delta = (pc + pt - pp) * v / kt
            v += delta
            if np.max(np.abs(delta)) < 1e-12:
                break
        return model.cold(v, c)[1] + slope * v

    if isinstance(seed, str):
        pp = d["P"][k].copy() if seed == "observed" else (target - c[0]) / c[1]
    else:
        pp = np.asarray(seed)[k].copy()
    step = np.full(len(target), c[0] / 100)
    previous = np.zeros(len(target))
    previous_delta = np.zeros(len(target))
    previous_step = step.copy()
    active = np.ones(len(target), dtype=bool)
    ktol, prel, pmin = (1e-12, 1e-12, 1e-12) if tight else (0.001, 0.0005, 0.0001)
    for cycle in range(100):
        delta = target - modulus(pp)
        if cycle:
            stop = (
                (abs(delta) < ktol * target)
                | (abs(pp - previous) < prel * pp)
                | (abs(step) < pmin)
            )
            active &= ~stop
            if not np.any(active):
                break
            test = np.divide(
                delta * (previous - pp),
                delta - previous_delta,
                out=step.copy(),
                where=abs(delta - previous_delta) > 10 * np.finfo(np.float32).tiny,
            )
            step = np.where(test < step, test, np.copysign(step, test))
        previous = pp.copy()
        previous_delta = delta
        if cycle > 1:
            step = np.where(
                abs(step) > abs(previous_step),
                np.copysign(0.9 * previous_step, step),
                step,
            )
        pp[active] += step[active]
        previous_step = step.copy()
    else:
        raise ValueError("Diagnostic modulus inversion did not converge")
    pressure[k] = pp
    return pressure


def forward_sigma(c, d, function=PRECISE):
    """EosFit's positive-one-error pressure increments; fallback only if zero."""
    p = function(c, d)
    terms = []
    for key, sigma, argument in [("Y", "sY", "y"), ("T", "sT", "t")]:
        increment = function(c, d, **{argument: d[key] + d[sigma]}) - p
        zero = abs(increment) < np.finfo(np.float32).tiny
        if np.any(zero):
            wider = (function(c, d, **{argument: d[key] + 3 * d[sigma]}) - p) / 3
            increment = np.where(zero, wider, increment)
        terms.append(increment)
    return np.sqrt(d["sP"] ** 2 + terms[0] ** 2 + terms[1] ** 2)


def jacobian(c, d):
    columns = []
    for i, scale in enumerate(SCALE):
        h = scale * 1e-5
        plus, minus = c.copy(), c.copy()
        plus[i] += h
        minus[i] -= h
        columns.append((PRECISE(plus, d) - PRECISE(minus, d)) / (2 * h))
    return np.array(columns).T


def covariance(j, w, chi):
    a = j * np.sqrt(w[:, None]) * SCALE
    _, r = np.linalg.qr(a, mode="reduced")
    inverse = np.linalg.solve(r, np.eye(5))
    return (inverse @ inverse.T) * np.outer(SCALE, SCALE) * chi


def step(c, j, w, observed, calculated):
    a = j * np.sqrt(w[:, None]) * SCALE
    shift = (
        np.linalg.lstsq(a, (observed - calculated) * np.sqrt(w), rcond=None)[0] * SCALE
    )
    return c + shift


def control_fit(rows, initial, function):
    with (
        patch.object(model, "pressure_at_observations", function),
        patch.object(
            model, "effective_sigma", lambda c, d: forward_sigma(c, d, function)
        ),
    ):
        result = model.solve(
            rows,
            initial,
            list(range(5)),
            names=NAMES,
            scaling=SCALE,
            lower=np.array([50, 2, 50, 0.05, 0.9]),
            upper=np.array([350, 10, 3000, 4, 1.1]),
        )
    # Retain the scientific result without duplicating another 194 input rows.
    return {k: v for k, v in result.items() if k != "rows"}


def run():
    rows, d, captured, evidence, old = inputs()
    c = np.array(list(evidence["starting_parameters"].values()))
    pn = np.array([r["native_pressure_before_step_gpa"] for r in captured])
    wn = np.array([r["native_weight"] for r in captured])
    jn = np.array([[r["dP_d_" + name] for name in NAMES] for r in captured])
    p, wp, jp = PRECISE(c, d), 1 / LOCAL_SIGMA(c, d) ** 2, jacobian(c, d)
    observed = d["P"].astype(np.float32).astype(float)
    chi = evidence["reduced_chi_square"]
    ix = evidence["covariance_parameter_indices_zero_based"]
    native_parameters = np.array(evidence["parameters59_after_step"])[ix]
    native_errors = np.array(evidence["errors59_after_step"])[ix]
    width_cases = []
    for use_j, use_w in product([False, True], repeat=2):
        checked = covariance(jn if use_j else jp, wn if use_w else wp, chi)
        width = np.sqrt(checked.diagonal())
        width_cases.append(
            dict(
                native_jacobian=use_j,
                native_weights=use_w,
                widths=dict(zip(NAMES, width.tolist())),
                relative_to_native_errors=dict(
                    zip(NAMES, (width / native_errors - 1).tolist())
                ),
            )
        )
    objective_cases = []
    for use_p, use_w in product([False, True], repeat=2):
        residual = (pn if use_p else p) - observed
        value = float(np.sum(residual**2 * (wn if use_w else wp)) / 189)
        objective_cases.append(
            dict(native_pressure=use_p, native_weights=use_w, reduced_sum_squares=value)
        )
    transitions = []
    for use_p, use_w, use_j in product([False, True], repeat=3):
        fitted = step(
            c, jn if use_j else jp, wn if use_w else wp, observed, pn if use_p else p
        )
        transitions.append(
            dict(
                native_pressure=use_p,
                native_weights=use_w,
                native_jacobian=use_j,
                after_one_step=dict(zip(NAMES, fitted.tolist())),
            )
        )
    reconstructed = step(c, jn, wn, observed, pn)
    groups = []
    comparison_rows = []
    coarse = coarse_pressure(c, d)
    wf = 1 / forward_sigma(c, d) ** 2
    wc = 1 / forward_sigma(c, d, coarse_pressure) ** 2
    for dataset in dict.fromkeys(r["dataset"] for r in rows):
        mask = np.array([r["dataset"] == dataset for r in rows])
        groups.append(
            dict(
                dataset=dataset,
                rows=int(mask.sum()),
                python_weighted_sum=float(
                    np.sum((p[mask] - observed[mask]) ** 2 * wp[mask])
                ),
                native_weighted_sum=float(
                    np.sum((pn[mask] - observed[mask]) ** 2 * wn[mask])
                ),
                maximum_pressure_difference_gpa=float(np.max(abs(p[mask] - pn[mask]))),
                maximum_relative_weight_difference=float(
                    np.max(abs(wp[mask] / wn[mask] - 1))
                ),
            )
        )
    for i, r in enumerate(rows):
        comparison_rows.append(
            dict(
                fit_row=i + 1,
                dataset=r["dataset"],
                kind=r["kind"],
                python_pressure_gpa=float(p[i]),
                native_pressure_gpa=float(pn[i]),
                double_precision_coarse_pressure_gpa=float(coarse[i]),
                python_derivative_weight=float(wp[i]),
                python_forward_weight=float(wf[i]),
                coarse_forward_weight=float(wc[i]),
                native_weight=float(wn[i]),
            )
        )
    with (DATA / (STEM + "-row-comparison.csv")).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=comparison_rows[0])
        writer.writeheader()
        writer.writerows(comparison_rows)
    fits = {
        "precise_forward": control_fit(rows, c, PRECISE),
        "coarse_forward": control_fit(rows, c, coarse_pressure),
    }
    report = dict(
        checked_on="2026-11-14",
        scope="Attribution of small numerical replay differences, not a replacement EOS or new experimental errors",
        starting_parameters=dict(zip(NAMES, c.tolist())),
        comparison_degrees_of_freedom=189,
        native_exact_reduced_chi_square=chi,
        original_python_reduced_chi_square=old["independent_fits"][0][
            "conditional_reduced_chi_square"
        ],
        objective_component_swaps=objective_cases,
        dataset_contributions=groups,
        covariance_component_swaps=width_cases,
        one_step_component_swaps=transitions,
        native_step_reconstruction=dict(
            parameters=dict(zip(NAMES, reconstructed.tolist())),
            actual_parameters=dict(zip(NAMES, native_parameters.tolist())),
            maximum_difference_in_float32_ulps=float(
                np.max(
                    abs(reconstructed - native_parameters)
                    / np.spacing(native_parameters.astype(np.float32))
                )
            ),
        ),
        covariance_reconstruction_maximum_relative_width_error=max(
            abs(v) for v in width_cases[-1]["relative_to_native_errors"].values()
        ),
        native_report_minus_reconstructed_cycle_chi_square=chi
        - objective_cases[-1]["reduced_sum_squares"],
        propagation_and_inversion_controls=fits,
        native_pressure_rounding_and_inversion_qualification="Primary attribution uses exact captured pre-step native pressures, weights and Jacobian, not rounded printed columns. The structural coarse-inversion control uses double precision and tight internal K(P,T), not a bitwise emulation of CFML's single-precision spline evaluation.",
        comparison_state_qualification="Common native cycle starting coefficients and identical supplied decimal observations; native observations enter its arithmetic as float32. Python's original high-precision V0 and pressure law are retained, so representation differences remain in the numerical-kernel term.",
        error_scaling_qualification="Covariance swaps hold the native residual-scale factor fixed to separate weights/Jacobian. The original Python factor differs by about0.026percent in variance; all errors remain local conditional widths.",
        catalog_changed=False,
        library_default_changed=False,
        author_contacted=False,
    )
    (DATA / (STEM + "-analysis.json")).write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in [
                    "native_exact_reduced_chi_square",
                    "original_python_reduced_chi_square",
                    "covariance_reconstruction_maximum_relative_width_error",
                    "native_step_reconstruction",
                ]
            },
            indent=2,
        )
    )
    return report


if __name__ == "__main__":
    run()
