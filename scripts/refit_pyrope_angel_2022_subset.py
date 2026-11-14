"""Conditional source-model refits of recovered pyrope data, not catalog EOSs.

Repository-only: ``python -m scripts.refit_pyrope_angel_2022_subset``.
The predeclared cases/assumptions are in pyrope/refit-plan.json. This does not
recover missing author inputs or optimize unknown error/weight multipliers.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
from scipy.constants import Avogadro
from scipy.optimize import least_squares

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_pyrope_angel_2022_observations import DATA, observations

V0, TR, N = 113.13, 298.0, 20
GAS_CONSTANT = float(np.float32(8.314))
NODES, WEIGHTS = np.polynomial.legendre.leggauss(64)
NAMES = ["K0_gpa", "K0_prime", "theta0_k", "gamma0", "scale_milani", "scale_hartwig"]
SCALING = np.array([170, 4.5, 770, 1.2, 1, 1])
LOWER = np.array([50, 2, 50, 0.05, 0.98, 0.98])
UPPER = np.array([350, 10, 3000, 4, 1.02, 1.02])


def thermal(t, c):
    """Gauss-Legendre Debye energy and Cv, checked against independent quad."""
    t = np.atleast_1d(t)

    def integral(temp, heat_capacity=False):
        y = c[2] / temp
        x = 0.5 * y[:, None] * (NODES + 1)
        f = x**3 / np.expm1(x)
        if heat_capacity:
            f = x**4 * np.exp(x) / np.expm1(x) ** 2
        return 0.5 * y * np.sum(f * WEIGHTS, axis=1) / y**3

    energy = 9 * N * GAS_CONSTANT * t * integral(t)
    reference = 9 * N * GAS_CONSTANT * TR * integral(np.array([TR]))[0]
    cv = 9 * N * GAS_CONSTANT * integral(t, heat_capacity=True)
    return c[3] * (energy - reference) / (1000 * V0), cv


def cold(v, c):
    eta = (V0 / np.asarray(v)) ** (1 / 3)
    a = 0.75 * (c[1] - 4)
    pressure = 1.5 * c[0] * (eta**7 - eta**5) * (1 + a * (eta**2 - 1))
    kt = (
        0.5
        * c[0]
        * eta
        * (
            (7 * eta**6 - 5 * eta**4) * (1 + a * (eta**2 - 1))
            + (eta**7 - eta**5) * 2 * a * eta
        )
    )
    return pressure, kt


def ks_volume(ks, t, c, cv):
    """Invert the source's KS(V,T) convention; no observed XRD V used."""
    v = V0 * (c[0] / ks) ** (1 / c[1])
    slope = c[3] ** 2 * cv * t / (1000 * V0**2)
    for _ in range(30):
        _, kt = cold(v, c)
        h = 0.0001
        derivative = (cold(v + h, c)[1] - cold(v - h, c)[1]) / (2 * h) + slope
        step = (kt + slope * v - ks) / derivative
        v = np.clip(v - step, V0 * 0.60, V0 * 1.12)
        if np.max(np.abs(step)) < 1e-10:
            break
    if np.max(np.abs(cold(v, c)[1] + slope * v - ks)) > 1e-7:
        raise ValueError("KS inversion failed within declared limits")
    return v


def prepared_rows():
    rows = []
    for row in observations("milani-2015-table-a1.csv"):
        rows.append(
            {
                "dataset": "milani_vt",
                "row": int(row["source_row"]),
                "kind": "V",
                "T": float(row["temperature_k"]),
                "P": 0.0001,
                "Y": float(row["unit_cell_volume_a3"]) * Avogadro / 8e24,
                "sY": float(row["unit_cell_volume_sd_a3"]) * Avogadro / 8e24,
                "sT": None,
                "sP": None,
                "scale_index": 4,
            }
        )
    for row in observations("hartwig-galkin-2021-table-3.csv"):
        t = float(row["temperature_k"])
        rows.append(
            {
                "dataset": "hartwig_vt",
                "row": int(row["source_row"]),
                "kind": "V",
                "T": t,
                "P": 0.0001,
                "Y": float(row["measured_molar_volume_cm3_mol"]),
                "sY": 0.03,
                "sT": 0.01 * t,
                "sP": None,
                "scale_index": 5,
            }
        )
    for row in observations("zou-2012b-table-i.csv"):
        t, p = float(row["temperature_k"]), float(row["pressure_gpa"])
        base = {
            "row": int(row["source_row"]),
            "T": t,
            "P": p,
            "sT": 0.02 * t,
            "sP": float(row["pressure_gpa_error"]),
            "scale_index": None,
        }
        if t == 300:
            rows.append(
                {
                    **base,
                    "dataset": "zou_pv",
                    "kind": "V",
                    "Y": V0 * float(row["relative_volume"]),
                    "sY": V0 * float(row["relative_volume_error"]),
                }
            )
        rows.append(
            {
                **base,
                "dataset": "zou_ks",
                "kind": "KS",
                "Y": float(row["ks_gpa"]),
                "sY": float(row["ks_gpa_error"]),
            }
        )
    return rows


def arrays(rows):
    return {
        key: np.array([row[key] if row[key] is not None else 0 for row in rows])
        for key in ("T", "P", "Y", "sT", "sP", "sY")
    } | {
        "KS": np.array([row["kind"] == "KS" for row in rows]),
        "scale_index": [row["scale_index"] for row in rows],
    }


def pressure_at_observations(c, d, *, y=None, t=None):
    y = d["Y"] if y is None else y
    t = d["T"] if t is None else t
    pth, cv = thermal(t, c)
    scales = np.array([c[i] if i is not None else 1 for i in d["scale_index"]])
    v = y / scales
    k = d["KS"]
    if np.any(k):
        v[k] = ks_volume(v[k], t[k], c, cv[k])
    return cold(v, c)[0] + pth


def effective_sigma(c, d):
    h, ht = 0.0001, 0.01
    py = (
        pressure_at_observations(c, d, y=d["Y"] + h)
        - pressure_at_observations(c, d, y=d["Y"] - h)
    ) / (2 * h)
    pt = (
        pressure_at_observations(c, d, t=d["T"] + ht)
        - pressure_at_observations(c, d, t=d["T"] - ht)
    ) / (2 * ht)
    # None in raw rows means that term is omitted, not known zero error.
    return np.sqrt(d["sP"] ** 2 + (py * d["sY"]) ** 2 + (pt * d["sT"]) ** 2)


def solve(
    rows, initial, free, *, names=NAMES, scaling=SCALING, lower=LOWER, upper=UPPER
):
    d, c, trace = arrays(rows), np.array(initial, dtype=float), []
    converged = False
    for iteration in range(30):
        sigma = effective_sigma(c, d)
        fixed = c.copy()

        def residual(x):
            cc = fixed.copy()
            cc[free] = x
            return (pressure_at_observations(cc, d) - d["P"]) / sigma

        result = least_squares(
            residual,
            c[free],
            bounds=(lower[free], upper[free]),
            jac="3-point",
            x_scale=scaling[free],
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-9,
            max_nfev=400,
        )
        new = c.copy()
        new[free] = result.x
        shift = float(np.max(np.abs((new - c) / scaling)))
        trace.append(
            {
                "iteration": iteration + 1,
                "max_scaled_shift": shift,
                "optimizer_success": bool(result.success),
                "nfev": result.nfev,
            }
        )
        c = new
        if shift < 1e-8 and result.success:
            converged = True
            break
    sigma = effective_sigma(c, d)
    residual = pressure_at_observations(c, d) - d["P"]
    normalized = residual / sigma
    jacobian = []
    for i in free:
        h = scaling[i] * 1e-5
        plus, minus = c.copy(), c.copy()
        plus[i] += h
        minus[i] -= h
        jacobian.append(
            (pressure_at_observations(plus, d) - pressure_at_observations(minus, d))
            / (2 * h * sigma)
        )
    j = np.array(jacobian).T
    scales = scaling[free]
    _, singular, vectors = np.linalg.svd(j * scales, full_matrices=False)
    rank = int(np.sum(singular > singular[0] * 1e-12))
    if rank != len(free):
        raise ValueError("Conditional covariance is rank deficient")
    inverse = vectors.T / singular
    covariance = np.einsum("ik,jk->ij", inverse, inverse) * np.outer(scales, scales)
    dof = len(rows) - len(free)
    chisq = float(normalized @ normalized)
    widths = np.sqrt(np.diag(covariance))
    return {
        "parameters": dict(zip(names, c.tolist())),
        "free_parameters": [names[i] for i in free],
        "numerically_converged": converged,
        "irls_trace": trace,
        "at_parameter_bound": [
            names[i]
            for i in free
            if min(c[i] - lower[i], upper[i] - c[i]) < 1e-5 * scaling[i]
        ],
        "objective_sum_squares": chisq,
        "degrees_of_freedom": dof,
        "conditional_reduced_chi_square": chisq / dof,
        "covariance_order": [names[i] for i in free],
        "scaled_jacobian_rank": rank,
        "scaled_jacobian_singular_values": singular.tolist(),
        "conditional_unscaled_covariance": covariance.tolist(),
        "conditional_unscaled_widths": dict(
            zip([names[i] for i in free], widths.tolist())
        ),
        "conditional_residual_scaled_widths": dict(
            zip(
                [names[i] for i in free],
                (widths * np.sqrt(max(1, chisq / dof))).tolist(),
            )
        ),
        "uncertainty_qualification": "Local inverse weighted-Jacobian covariance conditional on this subset/protocol and independent initial coordinate errors. Residual-scaled widths use max(1, chi_square/dof), matching documented EosFit convention; neither is the original author covariance/confidence or total systematic uncertainty.",
        "rows": [
            {
                **row,
                "calculated_pressure_gpa": float(row["P"] + r),
                "calculated_minus_observed_pressure_gpa": float(r),
                "effective_pressure_error_gpa": float(s),
                "standardized_pressure_residual": float(z),
            }
            for row, r, s, z in zip(rows, residual, sigma, normalized)
        ],
    }


def run():
    p, source_cov = source("pyrope")
    initial = np.array([p[2], p[3], p[11], p[18], 1.001, 1.0003])
    rows = prepared_rows()
    cases = {
        "milani_vt": ([r for r in rows if r["dataset"] == "milani_vt"], [2, 3, 4]),
        "hartwig_vt": ([r for r in rows if r["dataset"] == "hartwig_vt"], [2, 3, 5]),
        "combined_vt": ([r for r in rows if r["dataset"].endswith("vt")], [2, 3, 4, 5]),
        "joint_recovered": (rows, [0, 1, 2, 3, 4, 5]),
    }
    report = {
        "scope": "Conditional source-model fits of available subsets; not original 2022 fit reproduction or replacement models",
        "author_workflow_reproduced": False,
        "source_parameters": dict(zip(NAMES[:4], initial[:4].tolist())),
        "source_covariance_marginal_widths": dict(
            zip(NAMES[:4], np.sqrt(np.diag(source_cov)[[1, 2, 10, 17]]).tolist())
        ),
        "gas_constant": GAS_CONSTANT,
        "fixed": {"V0_cm3_mol": V0, "Tref_k": TR, "atoms_per_formula": N},
        "plan_sha256": hashlib.sha256(
            (DATA / "refit-plan.json").read_bytes()
        ).hexdigest(),
        "source_data_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in (
                "milani-2015-table-a1.csv",
                "hartwig-galkin-2021-table-3.csv",
                "zou-2012b-table-i.csv",
            )
        },
        "cases": {},
    }
    for name, (selected, free) in cases.items():
        baseline = solve(selected, initial, [i for i in free if i >= 4])
        seed = np.array([baseline["parameters"][key] for key in NAMES])
        fitted = solve(selected, seed, free)
        fitted["source_coefficient_differences"] = {
            NAMES[i]: fitted["parameters"][NAMES[i]] - initial[i] for i in free if i < 4
        }
        report["cases"][name] = {"source_eos_scales_only": baseline, "refit": fitted}
        print(
            name,
            fitted["numerically_converged"],
            fitted["parameters"],
            fitted["conditional_reduced_chi_square"],
            flush=True,
        )
    (DATA / "subset-refits.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    run()
