"""Analytic-derivative pyrope replay; no change to the public fitting API.

Repository-only prototype for the supplied 194-row, pressure-residual IRLS
protocol. The public Peritheos fitter has a separate errors-in-variables
objective; this routine must not be presented as its implementation.
"""

from __future__ import annotations

import json
from time import perf_counter

import numpy as np
from scipy.optimize import least_squares

from scripts.audit_pyrope_author_inputs import DATA, NAMES, PREFIX
from scripts.refit_pyrope_angel_2022_subset import (
    GAS_CONSTANT,
    NODES,
    TR,
    V0,
    WEIGHTS,
    N,
    arrays,
    cold,
    ks_volume,
    solve,
    thermal,
)

SCALING = np.array([170, 4.5, 770, 1.2, 1])
LOWER = np.array([50, 2, 50, 0.05, 0.9])
UPPER = np.array([350, 10, 3000, 4, 1.1])


def energy_theta_derivative(theta, temperature):
    y = theta / np.atleast_1d(temperature)
    x = 0.5 * y[:, None] * (NODES + 1)
    integral = 0.5 * y * np.sum(x**3 / np.expm1(x) * WEIGHTS, axis=1) / y**3
    return 9 * N * GAS_CONSTANT * (1 / np.expm1(y) - 3 * integral / y)


def evaluate(c, d):
    """Pressure and exact local parameter/coordinate derivatives.

    Modulus rows use implicit differentiation of KS(V,T,c)=KS_observed/scale.
    Debye energy and Cv retain the independently checked 64-node quadrature.
    """
    t = d["T"]
    pth, cv = thermal(t, c)
    scales = np.array([c[i] if i is not None else 1 for i in d["scale_index"]])
    v = d["Y"] / scales
    k = d["KS"]
    if np.any(k):
        v[k] = ks_volume(v[k], t[k], c, cv[k])
    pc, kt = cold(v, c)
    h = (V0 / v) ** (2 / 3)
    f = 0.5 * (h - 1)
    b = 1.5 * (c[1] - 4)
    shape = 1 + 7 * f + b * f * (2 + 9 * f)
    kt_v = -c[0] * h**2.5 * (5 * shape + h * (7 + b * (2 + 18 * f))) / (3 * v)
    pc_kp = 4.5 * c[0] * f**2 * h**2.5
    kt_kp = 1.5 * c[0] * h**2.5 * f * (2 + 9 * f)
    p_v = -kt / v
    p_t_fixed_v = c[3] * cv / (1000 * V0)
    p_theta = (
        c[3]
        / (1000 * V0)
        * (energy_theta_derivative(c[2], t) - energy_theta_derivative(c[2], TR)[0])
    )
    j = np.column_stack([pc / c[0], pc_kp, p_theta, pth / c[3], np.zeros_like(pc)])
    p_y = p_v / scales
    p_t = p_t_fixed_v.copy()
    scaled = np.array([i is not None for i in d["scale_index"]])
    scaled_volumes = scaled & ~k
    j[scaled_volumes, 4] = (
        -p_v[scaled_volumes] * d["Y"][scaled_volumes] / scales[scaled_volumes] ** 2
    )
    if np.any(k):
        y = c[2] / t
        cv_theta = (
            9
            * N
            * GAS_CONSTANT
            / t
            * (y * np.exp(y) / np.expm1(y) ** 2 - 3 * (cv / (9 * N * GAS_CONSTANT)) / y)
        )
        cv_t = -c[2] / t * cv_theta
        slope = c[3] ** 2 * cv * t / (1000 * V0**2)
        slope_theta = c[3] ** 2 * cv_theta * t / (1000 * V0**2)
        slope_t = c[3] ** 2 * (cv + t * cv_t) / (1000 * V0**2)
        ks_v = kt_v + slope
        ks_c = np.column_stack(
            [kt / c[0], kt_kp, slope_theta * v, 2 * slope * v / c[3], np.zeros_like(pc)]
        )
        ks_c[scaled, 4] = d["Y"][scaled] / scales[scaled] ** 2
        j[k] -= p_v[k, None] * ks_c[k] / ks_v[k, None]
        p_y[k] = p_v[k] / (scales[k] * ks_v[k])
        p_t[k] -= p_v[k] * slope_t[k] * v[k] / ks_v[k]
    return pc + pth, j, p_y, p_t


def fit(rows, initial):
    d = arrays(rows)
    c = np.array(initial, dtype=float)
    trace = []
    calls = 0
    for iteration in range(30):
        _, _, py, pt = evaluate(c, d)
        sigma = np.sqrt(d["sP"] ** 2 + (py * d["sY"]) ** 2 + (pt * d["sT"]) ** 2)
        cache_c, cache_result = None, None

        def cached(values):
            nonlocal cache_c, cache_result, calls
            if cache_c is None or not np.array_equal(values, cache_c):
                cache_c = values.copy()
                cache_result = evaluate(values, d)
                calls += 1
            return cache_result

        result = least_squares(
            lambda values: (cached(values)[0] - d["P"]) / sigma,
            c,
            jac=lambda values: cached(values)[1] / sigma[:, None],
            bounds=(LOWER, UPPER),
            x_scale=SCALING,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-9,
            max_nfev=400,
        )
        shift = float(np.max(abs((result.x - c) / SCALING)))
        c = result.x
        trace.append(
            dict(
                iteration=iteration + 1,
                max_scaled_shift=shift,
                optimizer_success=bool(result.success),
                nfev=result.nfev,
                njev=result.njev,
            )
        )
        if shift < 1e-8 and result.success:
            break
    pressure, jac, py, pt = evaluate(c, d)
    sigma = np.sqrt(d["sP"] ** 2 + (py * d["sY"]) ** 2 + (pt * d["sT"]) ** 2)
    z = (pressure - d["P"]) / sigma
    j = jac / sigma[:, None]
    score = np.abs(j.T @ z) / (np.linalg.norm(j, axis=0) * np.linalg.norm(z))
    _, r = np.linalg.qr(j * SCALING, mode="reduced")
    inverse = np.linalg.solve(r, np.eye(5))
    covariance = (inverse @ inverse.T) * np.outer(SCALING, SCALING)
    chi = float(z @ z / 189)
    return dict(
        parameters=dict(zip(NAMES, c.tolist())),
        numerically_converged=bool(
            trace[-1]["max_scaled_shift"] < 1e-8
            and result.success
            and np.max(score) < 1e-6
        ),
        conditional_reduced_chi_square=chi,
        conditional_unscaled_covariance=covariance.tolist(),
        conditional_residual_scaled_widths=dict(
            zip(NAMES, (np.sqrt(covariance.diagonal() * chi)).tolist())
        ),
        maximum_normalized_stationarity_score=float(np.max(score)),
        scaled_jacobian_condition_number=float(np.linalg.cond(j * SCALING)),
        at_parameter_bound=[
            NAMES[i]
            for i in range(5)
            if min(c[i] - LOWER[i], UPPER[i] - c[i]) < 1e-5 * SCALING[i]
        ],
        model_and_derivative_evaluations=calls,
        irls_trace=trace,
    )


def run():
    original = json.loads((DATA / (PREFIX + "-audit.json")).read_text())
    keys = ["dataset", "row", "kind", "P", "T", "Y", "sP", "sT", "sY", "scale_index"]
    rows = [{k: r[k] for k in keys} for r in original["independent_fits"][0]["rows"]]
    starts = [[169.339, 4.55314, 771.096, 1.1847, 0.97828], [170, 4.5, 800, 1.2, 0.98]]
    fits = [fit(rows, initial) for initial in starts]
    timings = {"finite_difference_seconds": [], "analytic_seconds": []}
    for i in range(5):
        # Alternate the order to reduce simple warm-state/order bias.
        names = ["finite_difference_seconds", "analytic_seconds"]
        if i % 2:
            names.reverse()
        for name in names:
            start = perf_counter()
            if name.startswith("finite"):
                solve(
                    rows,
                    starts[0],
                    list(range(5)),
                    names=NAMES,
                    scaling=SCALING,
                    lower=LOWER,
                    upper=UPPER,
                )
            else:
                fit(rows, starts[0])
            timings[name].append(perf_counter() - start)
    report = dict(
        scope="Analytic-derivative improvement to the repository pyrope replay only; public errors-in-variables fitter and weighting policy unchanged",
        fits=fits,
        timing_repeats=timings,
        benchmark_qualification="Five local paired alternating-order repeats of this complete194rowfit, not a general library performance claim",
        library_default_changed=False,
        source_observations_changed=False,
    )
    (DATA / "numerics-2026-11-14-analytic-fit.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(fits[0], indent=2))
    print(
        "median speed ratio",
        np.median(timings["finite_difference_seconds"])
        / np.median(timings["analytic_seconds"]),
    )
    return report


if __name__ == "__main__":
    run()
