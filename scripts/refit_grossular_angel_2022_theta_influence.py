"""Theta influence with the retained independent2%Conrad diagnostic errors.

Repository-only: python -m scripts.refit_grossular_angel_2022_theta_influence.
Every case is predeclared; this never writes a catalog or source record.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter

import numpy as np
from scipy.optimize import least_squares

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, GAS_CONSTANT
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_conrad_processing import prepare as conrad_rows
from scripts.refit_grossular_angel_2022_dataset_influence import derivatives, group_of
from scripts.refit_grossular_angel_2022_full_sensitivity import write_csv
from scripts.refit_grossular_angel_2022_subset import (
    arrays,
    effective_sigma,
    pressure,
    solve,
)

PLAN = DATA / "theta-influence-plan.json"


def plan():
    return json.loads(PLAN.read_text())


def cases():
    return [
        dict(id=f"{n}__{a}", normalization=n, action=a)
        for n in plan()["normalizations"]
        for a in ["baseline", "published_theta"]
        + ["omit_" + g["id"] for g in plan()["omissions"]]
    ]


def prepare(case):
    rows, starts, config = conrad_rows(
        dict(normalization=case["normalization"], variant="extra_020")
    )
    if case["action"].startswith("omit_"):
        omitted = next(
            g["groups"]
            for g in plan()["omissions"]
            if case["action"] == "omit_" + g["id"]
        )
        rows = [r for r in rows if group_of(r) not in omitted]
    assert sum(r["dataset"] == "conrad_ks" for r in rows) == 3
    used = set(range(4)) | {r["scale"] for r in rows if r["scale"] is not None}
    indices = [i for i in range(len(config["names"])) if i in used]
    remap = {old: new for new, old in enumerate(indices)}
    for r in rows:
        if r["scale"] is not None:
            r["scale"] = remap[r["scale"]]
    config = {
        k: [v[i] for i in indices] if k == "names" else v[indices]
        for k, v in config.items()
    }
    return rows, [[s[i] for i in indices] for s in starts], config


def solve_fixed_theta(rows, initial, config):
    """Same IRLS objective and tolerances, with theta fixed at the source value."""
    data, c = arrays(rows), np.array(initial, dtype=float)
    c[2] = plan()["published_theta_k"]
    active = [i for i in range(len(c)) if i != 2]
    scales = config["scaling"][active]
    trace, converged = [], False
    for iteration in range(30):
        sigma = effective_sigma(c, data)

        def expand(values):
            full = c.copy()
            full[active] = values
            return full

        fit = least_squares(
            lambda values: (pressure(expand(values), data) - data["P"]) / sigma,
            c[active],
            bounds=(config["lower"][active], config["upper"][active]),
            jac="3-point",
            x_scale=scales,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-9,
            max_nfev=400,
        )
        shift = float(np.max(np.abs((fit.x - c[active]) / scales)))
        c[active] = fit.x
        trace.append(
            dict(
                iteration=iteration + 1,
                max_scaled_shift=shift,
                success=bool(fit.success),
                nfev=fit.nfev,
            )
        )
        if shift < 1e-8 and fit.success:
            converged = True
            break
    sigma = effective_sigma(c, data)
    residual = pressure(c, data) - data["P"]
    jac = derivatives(c, data, sigma, config["scaling"])[:, active]
    _, singular, vectors = np.linalg.svd(jac * scales, full_matrices=False)
    if np.sum(singular > singular[0] * 1e-12) != len(active):
        raise ValueError("Fixed-theta covariance is rank deficient")
    inverse = vectors.T / singular
    covariance = np.einsum("ik,jk->ij", inverse, inverse) * np.outer(scales, scales)
    names = [config["names"][i] for i in active]
    objective = float(np.sum((residual / sigma) ** 2))
    return dict(
        parameters=dict(zip(config["names"], c.tolist())),
        numerically_converged=converged,
        at_parameter_bound=[
            config["names"][i]
            for i in active
            if min(c[i] - config["lower"][i], config["upper"][i] - c[i])
            < 1e-5 * config["scaling"][i]
        ],
        irls_trace=trace,
        objective_sum_squares=objective,
        degrees_of_freedom=len(rows) - len(active),
        conditional_reduced_chi_square=objective / (len(rows) - len(active)),
        covariance_order=names,
        conditional_unscaled_covariance=covariance.tolist(),
        conditional_unscaled_widths=dict(
            zip(names, np.sqrt(np.diag(covariance)).tolist())
        ),
        scaled_jacobian_singular_values=singular.tolist(),
        rows=[
            r
            | dict(
                calculated_pressure_gpa=float(r["P"] + v),
                pressure_residual_gpa=float(v),
                effective_error_gpa=float(s),
            )
            for r, v, s in zip(rows, residual, sigma)
        ],
    )


def temperature_bin(temperature):
    return (
        "below298"
        if temperature < 298
        else "298to699"
        if temperature < 700
        else "700andabove"
    )


def profile_diagnostics(normalization, by_id):
    baseline = by_id[f"{normalization}__baseline"]
    fixed = by_id[f"{normalization}__published_theta"]
    rows, _, config = prepare(fixed)
    data = arrays(rows)
    c = np.array(list(fixed["fits"][0]["parameters"].values()))
    free = np.array(list(baseline["fits"][0]["parameters"].values()))
    p, _ = source("grossular")
    published = np.array([p[2], p[3], p[11], p[18], p[53], p[54]])
    sigma = effective_sigma(c, data)
    residuals = {
        label: (pressure(coefficients, data) - data["P"]) / sigma
        for label, coefficients in [
            ("published_all", published),
            ("fixed_theta", c),
            ("free", free),
        ]
    }
    jac = derivatives(c, data, sigma, config["scaling"])
    active = [0, 1, 3, 4, 5]
    nuisance = jac[:, active] * config["scaling"][active]
    coefficients = np.linalg.lstsq(nuisance, jac[:, 2], rcond=None)[0]
    efficient = jac[:, 2] - np.einsum("ij,j->i", nuisance, coefficients)
    information = float(efficient @ efficient)
    row_receipt = [
        dict(
            normalization=normalization,
            group=group_of(r),
            dataset=r["dataset"],
            source_row=r["row"],
            P_gpa=r["P"],
            T_k=r["T"],
            temperature_bin=temperature_bin(r["T"]),
            common_profile_error_gpa=float(sigma[i]),
            **{
                f"{label}_standardized_residual": float(values[i])
                for label, values in residuals.items()
            },
            efficient_theta_derivative_per_k=float(efficient[i]),
            information_fraction=float(efficient[i] ** 2 / information),
            local_upward_theta_contribution_k=float(
                -efficient[i] * residuals["fixed_theta"][i] / information
            ),
        )
        for i, r in enumerate(rows)
    ]

    def aggregate(selected, **labels):
        costs = {
            label: sum(r[f"{label}_standardized_residual"] ** 2 for r in selected)
            for label in residuals
        }
        return dict(
            normalization=normalization,
            **labels,
            rows=len(selected),
            **{f"{label}_sum_squares": value for label, value in costs.items()},
            fixed_minus_free_sum_squares=costs["fixed_theta"] - costs["free"],
            information_fraction=sum(r["information_fraction"] for r in selected),
            local_upward_theta_contribution_k=sum(
                r["local_upward_theta_contribution_k"] for r in selected
            ),
        )

    groups = [
        aggregate([r for r in row_receipt if r["group"] == g], group=g)
        for g in Counter(group_of(r) for r in rows)
    ]
    bins = [
        aggregate(
            [r for r in row_receipt if (r["group"], r["temperature_bin"]) == (g, b)],
            group=g,
            temperature_bin=b,
        )
        for g, b in dict.fromkeys(
            (r["group"], r["temperature_bin"]) for r in row_receipt
        )
    ]
    covariance = np.array(baseline["fits"][0]["conditional_unscaled_covariance"])
    diagnostic = dict(
        common_errors="Effective errors frozen at fixed-theta fit for all curves and local derivatives",
        efficient_information_per_k2=information,
        nuisance_orthogonality_max=float(
            np.max(np.abs(np.einsum("ij,i->j", nuisance, efficient)))
        ),
        local_predicted_theta_shift_k=float(
            -efficient @ residuals["fixed_theta"] / information
        ),
        actual_theta_shift_k=float(free[2] - c[2]),
        baseline_theta_gamma_correlation=float(
            covariance[2, 3] / np.sqrt(covariance[2, 2] * covariance[3, 3])
        ),
        group_sums=groups,
    )
    return diagnostic, row_receipt, groups, bins


def native_profile_rows(normalization, by_id):
    rows, _, _ = prepare(by_id[f"{normalization}__baseline"])
    p, _ = source("grossular")
    curves = {"published_all": np.array([p[2], p[3], p[11], p[18], p[53], p[54]])}
    curves.update(
        {
            label: np.array(
                list(
                    by_id[f"{normalization}__{action}"]["fits"][0][
                        "parameters"
                    ].values()
                )
            )
            for label, action in [
                ("baseline", "baseline"),
                ("fixed_theta", "published_theta"),
            ]
        }
    )
    output = []
    for r in rows:
        row = dict(
            normalization=normalization,
            group=group_of(r),
            dataset=r["dataset"],
            source_row=r["row"],
            P_gpa=r["P"],
            T_k=r["T"],
            coordinate="KS" if r["KS"] else "V",
            unit="GPa" if r["KS"] else "cm3/mol",
            observed=r["Y"],
        )
        for label, coefficients in curves.items():
            model = p | dict(zip([2, 3, 11, 18], coefficients[:4]))
            scale = coefficients[r["scale"]] if r["scale"] is not None else 1
            predicted = (
                scale
                * properties(model, r["P"], r["T"], GAS_CONSTANT)[
                    "ks_gpa" if r["KS"] else "volume_cm3_mol"
                ]
            )
            row[label] = predicted
            row[f"observed_minus_{label}"] = r["Y"] - predicted
            row[f"{label}_multiplier"] = scale
        output.append(row)
    return output


def normalization_effects(by_id):
    """Paired normalization contrasts; no additional fitting or causal claim."""
    output = []
    baseline_gap = (
        by_id["linear_298__baseline"]["fits"][0]["parameters"]["theta0_k"]
        - by_id["plotted__baseline"]["fits"][0]["parameters"]["theta0_k"]
    )
    for case in cases():
        if case["normalization"] != "plotted":
            continue
        action = case["action"]
        old = by_id[f"plotted__{action}"]["fits"][0]["parameters"]
        new = by_id[f"linear_298__{action}"]["fits"][0]["parameters"]
        gap = new["theta0_k"] - old["theta0_k"]
        output.append(
            dict(
                action=action,
                linear_minus_plotted_theta_k=gap,
                normalization_interaction_vs_baseline_k=gap - baseline_gap,
                linear_minus_plotted_kprime=new["K0_prime"] - old["K0_prime"],
            )
        )
    return output


def run():
    results, residuals, summaries = [], [], []
    for case in cases():
        rows, starts, config = prepare(case)
        fixed = case["action"] == "published_theta"
        fits = [
            solve_fixed_theta(rows, start, config)
            if fixed
            else solve(rows, start, **config)
            for start in starts
        ]
        residuals.extend(dict(case=case["id"]) | r for r in fits[0]["rows"])
        for fit in fits:
            del fit["rows"]
        difference = float(
            np.max(
                np.abs(
                    (
                        np.array(list(fits[0]["parameters"].values()))
                        - np.array(list(fits[1]["parameters"].values()))
                    )
                    / config["scaling"]
                )
            )
        )
        results.append(
            case
            | dict(
                rows=len(rows),
                dataset_counts=dict(Counter(r["dataset"] for r in rows)),
                free_parameters=fits[0]["covariance_order"],
                fixed_parameters={"theta0_k": plan()["published_theta_k"]}
                if fixed
                else {},
                fits=fits,
                maximum_scaled_start_difference=difference,
            )
        )
        print(
            f"{case['id']}: theta={fits[0]['parameters']['theta0_k']:.6f}, Kprime={fits[0]['parameters']['K0_prime']:.6f}, converged={all(f['numerically_converged'] for f in fits)}",
            flush=True,
        )
    by_id = {r["id"]: r for r in results}
    for r in results:
        values = r["fits"][0]["parameters"]
        baseline = by_id[f"{r['normalization']}__baseline"]["fits"][0]["parameters"]
        r["effects"] = {
            k: dict(minus_baseline=v - baseline[k]) for k, v in values.items()
        }
        summaries.append(
            {k: r[k] for k in ["id", "normalization", "action", "rows"]}
            | values
            | dict(
                theta_minus_baseline_k=values["theta0_k"] - baseline["theta0_k"],
                theta_minus_source_k=values["theta0_k"] - plan()["published_theta_k"],
                kprime_minus_source=values["K0_prime"] - plan()["published_kprime"],
                objective=r["fits"][0]["objective_sum_squares"],
                free_parameters=len(r["free_parameters"]),
            )
        )
    profiles, profile_rows, groups, bins, native = {}, [], [], [], []
    for n in plan()["normalizations"]:
        profiles[n], row, group, temperature = profile_diagnostics(n, by_id)
        profile_rows.extend(row)
        groups.extend(group)
        bins.extend(temperature)
        native.extend(native_profile_rows(n, by_id))
    report = dict(
        audit_date=plan()["declared_on"],
        scope="theta_dataset_influence_with_independent2percentConrad_diagnostic",
        author_workflow_reproduced=False,
        physical_validation=False,
        catalog_changes=False,
        author_contacted=False,
        conrad_error_model="Formal upper width plus additional independent2%KS in quadrature; diagnostic, not measured errors",
        all_three_conrad_points_retained=True,
        cases=results,
        profiles=profiles,
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        all_starts_converged=all(
            f["numerically_converged"] for r in results for f in r["fits"]
        ),
        any_start_at_bound=any(
            bool(f["at_parameter_bound"]) for r in results for f in r["fits"]
        ),
        maximum_scaled_start_difference=max(
            r["maximum_scaled_start_difference"] for r in results
        ),
        input_sha256=plan()["input_sha256"],
        normalization_comparisons=normalization_effects(by_id),
    )
    keys = list(dict.fromkeys(k for r in summaries for k in r))
    write_csv(
        "theta-influence-summary.csv", [{k: r.get(k) for k in keys} for r in summaries]
    )
    write_csv("theta-influence-residuals.csv", residuals)
    write_csv("theta-influence-profile-rows.csv", profile_rows)
    write_csv("theta-influence-groups.csv", groups)
    write_csv("theta-influence-temperature-bins.csv", bins)
    write_csv("theta-influence-native-comparisons.csv", native)
    write_csv(
        "theta-influence-normalization-interactions.csv", normalization_effects(by_id)
    )
    (DATA / "theta-influence-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    run()
