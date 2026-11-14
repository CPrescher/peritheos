"""Declared dataset influence diagnostics, never a replacement EOS.

Repository-only: python -m scripts.refit_grossular_angel_2022_dataset_influence.
See dataset-influence-plan.json for the fixed design and comparison limits.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter

import numpy as np
from scipy.optimize import least_squares

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import (
    DATA,
    GAS_CONSTANT,
    observations,
)
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_full_sensitivity import write_csv
from scripts.refit_grossular_angel_2022_normalization_specimen import prepared_case
from scripts.refit_grossular_angel_2022_subset import (
    arrays,
    effective_sigma,
    pressure,
    solve,
)

PLAN = DATA / "dataset-influence-plan.json"


def plan():
    return json.loads(PLAN.read_text())


def group_of(row):
    return "isaak_ks" if row["dataset"] == "isaak_ks_estimated" else row["dataset"]


def cases():
    return [
        dict(id=f"{norm}__{action}", normalization=norm, action=action)
        for norm in plan()["normalizations"]
        for action in ["baseline", "published_kprime"]
        + ["omit_" + g["id"] for g in plan()["omissions"]]
    ]


def prepare(case):
    rows, starts, config = prepared_case(
        dict(normalization=case["normalization"], isaak="as_plotted")
    )
    if case["action"].startswith("omit_"):
        omission = next(
            g for g in plan()["omissions"] if "omit_" + g["id"] == case["action"]
        )
        rows = [r for r in rows if group_of(r) not in omission["groups"]]
    used = set(range(4)) | {r["scale"] for r in rows if r["scale"] is not None}
    indices = [i for i in range(len(config["names"])) if i in used]
    remap = {old: new for new, old in enumerate(indices)}
    for row in rows:
        if row["scale"] is not None:
            row["scale"] = remap[row["scale"]]
    config = {
        key: [value[i] for i in indices] if key == "names" else value[indices]
        for key, value in config.items()
    }
    starts = [[s[i] for i in indices] for s in starts]
    return rows, starts, config


def derivatives(c, data, sigma, scaling):
    columns = []
    for i, scale in enumerate(scaling):
        step = scale * 1e-5
        plus, minus = c.copy(), c.copy()
        plus[i] += step
        minus[i] -= step
        columns.append(
            (pressure(plus, data) - pressure(minus, data)) / (2 * step * sigma)
        )
    return np.array(columns).T


def solve_fixed_kprime(rows, initial, config):
    """Same IRLS protocol, optimizing five coefficients with Kprime fixed."""
    data = arrays(rows)
    c = np.array(initial, dtype=float)
    c[1] = plan()["published_kprime"]
    active = [i for i in range(len(c)) if i != 1]
    names = [config["names"][i] for i in active]
    scales = config["scaling"][active]
    trace, converged = [], False
    for iteration in range(30):
        sigma = effective_sigma(c, data)

        def expand(x):
            full = c.copy()
            full[active] = x
            return full

        fit = least_squares(
            lambda x: (pressure(expand(x), data) - data["P"]) / sigma,
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
        raise ValueError("Fixed-Kprime covariance is rank deficient")
    inverse = vectors.T / singular
    covariance = (inverse @ inverse.T) * np.outer(scales, scales)
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


def profile_diagnostics(normalization, by_id):
    """Fixed-error local score and information, not a likelihood attribution."""
    baseline = by_id[f"{normalization}__baseline"]
    constrained = by_id[f"{normalization}__published_kprime"]
    rows, _, config = prepare(constrained)
    data = arrays(rows)
    c = np.array(list(constrained["fits"][0]["parameters"].values()))
    free_c = np.array(list(baseline["fits"][0]["parameters"].values()))
    p, _ = source("grossular")
    source_c = np.array([p[2], p[3], p[11], p[18], p[53], p[54]])
    sigma = effective_sigma(c, data)
    residuals = {
        "published_all": (pressure(source_c, data) - data["P"]) / sigma,
        "fixed_kprime": (pressure(c, data) - data["P"]) / sigma,
        "free": (pressure(free_c, data) - data["P"]) / sigma,
    }
    jac = derivatives(c, data, sigma, config["scaling"])
    nuisance = jac[:, [0, 2, 3, 4, 5]]
    coefficients = np.linalg.lstsq(nuisance, jac[:, 1], rcond=None)[0]
    efficient = jac[:, 1] - np.einsum("ij,j->i", nuisance, coefficients)
    information = float(efficient @ efficient)
    output = []
    for i, row in enumerate(rows):
        output.append(
            dict(
                normalization=normalization,
                group=group_of(row),
                dataset=row["dataset"],
                source_row=row["row"],
                P_gpa=row["P"],
                T_k=row["T"],
                common_profile_error_gpa=float(sigma[i]),
                **{
                    f"{name}_standardized_residual": float(value[i])
                    for name, value in residuals.items()
                },
                efficient_kprime_derivative=float(efficient[i]),
                information_fraction=float(efficient[i] ** 2 / information),
                local_upward_kprime_contribution=float(
                    -efficient[i] * residuals["fixed_kprime"][i] / information
                ),
            )
        )
    grouped = []
    for group in Counter(group_of(r) for r in rows):
        selected = [r for r in output if r["group"] == group]
        costs = {
            name: sum(r[f"{name}_standardized_residual"] ** 2 for r in selected)
            for name in residuals
        }
        grouped.append(
            dict(
                normalization=normalization,
                group=group,
                rows=len(selected),
                **{f"{name}_sum_squares": value for name, value in costs.items()},
                fixed_minus_free_sum_squares=costs["fixed_kprime"] - costs["free"],
                information_fraction=sum(r["information_fraction"] for r in selected),
                local_upward_kprime_contribution=sum(
                    r["local_upward_kprime_contribution"] for r in selected
                ),
            )
        )
    diagnostics = dict(
        common_errors="Effective errors at constrained fit; frozen for all three curves and local derivatives",
        efficient_information=information,
        nuisance_orthogonality_max=float(
            np.max(np.abs(np.einsum("ij,i->j", nuisance, efficient)))
        ),
        local_predicted_kprime_shift=float(
            -efficient @ residuals["fixed_kprime"] / information
        ),
        actual_kprime_shift=float(free_c[1] - c[1]),
        group_sums=grouped,
    )
    return diagnostics, output, grouped


def conrad_native_audit(by_id):
    """Forward checks on all five literal rows; excluded rows stay excluded."""
    p, _ = source("grossular")
    output = []
    for normalization in plan()["normalizations"]:
        baseline = by_id[f"{normalization}__baseline"]["fits"][0]["parameters"]
        omitted = by_id[f"{normalization}__omit_conrad_ks"]["fits"][0]["parameters"]
        for r in observations("conrad-1999-table-6.csv"):
            k = (float(r["c11_gpa"]) + 2 * float(r["c12_gpa"])) / 3
            row = dict(
                normalization=normalization,
                source_row=int(r["source_row"]),
                selection=r["selection_angel_2022"],
                P_gpa=float(r["pressure_gpa"]),
                adopted_T_k=298,
                c11_token=r["c11_gpa_token"],
                c12_token=r["c12_gpa_token"],
                derived_ks_gpa=k,
            )
            for name, parameters in [
                ("published", None),
                ("baseline", baseline),
                ("omit_conrad", omitted),
            ]:
                model = (
                    p
                    if parameters is None
                    else p | dict(zip([2, 3, 11, 18], list(parameters.values())[:4]))
                )
                predicted = properties(model, row["P_gpa"], 298, GAS_CONSTANT)["ks_gpa"]
                row[f"{name}_ks_gpa"] = predicted
                row[f"observed_minus_{name}_ks_gpa"] = k - predicted
            output.append(row)
    return output


def run():
    results, residuals, summary = [], [], []
    for case in cases():
        rows, starts, config = prepare(case)
        fixed = case["action"] == "published_kprime"
        fits = [
            solve_fixed_kprime(rows, start, config)
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
        result = case | dict(
            rows=len(rows),
            dataset_counts=dict(Counter(r["dataset"] for r in rows)),
            free_parameters=fits[0]["covariance_order"],
            fixed_parameters={"K0_prime": plan()["published_kprime"]} if fixed else {},
            fits=fits,
            maximum_scaled_start_difference=difference,
        )
        results.append(result)
        print(
            f"{case['id']}: Kprime={fits[0]['parameters']['K0_prime']:.6f}, converged={all(f['numerically_converged'] for f in fits)}",
            flush=True,
        )
    by_id = {r["id"]: r for r in results}
    for result in results:
        fitted = result["fits"][0]
        baseline = by_id[f"{result['normalization']}__baseline"]["fits"][0][
            "parameters"
        ]
        result["effects"] = {
            name: dict(minus_baseline=value - baseline[name])
            for name, value in fitted["parameters"].items()
        }
        summary.append(
            {k: result[k] for k in ["id", "normalization", "action", "rows"]}
            | fitted["parameters"]
            | dict(
                kprime_minus_baseline=fitted["parameters"]["K0_prime"]
                - baseline["K0_prime"],
                objective_sum_squares=fitted["objective_sum_squares"],
                free_parameters=len(result["free_parameters"]),
                both_starts_converged=all(
                    f["numerically_converged"] for f in result["fits"]
                ),
            )
        )
    profile, profile_rows, groups = {}, [], []
    for normalization in plan()["normalizations"]:
        profile[normalization], rows, grouped = profile_diagnostics(
            normalization, by_id
        )
        profile_rows.extend(rows)
        groups.extend(grouped)
    report = dict(
        audit_date=plan()["declared_on"],
        scope="conditional_dataset_influence_diagnostic",
        author_workflow_reproduced=False,
        physical_validation=False,
        source_coefficients_changed=False,
        catalog_changes=False,
        author_contacted=False,
        plan=PLAN.name,
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        cases=results,
        profile=profile,
        all_starts_converged=all(
            f["numerically_converged"] for r in results for f in r["fits"]
        ),
        any_start_at_bound=any(
            bool(f["at_parameter_bound"]) for r in results for f in r["fits"]
        ),
        maximum_scaled_start_difference=max(
            r["maximum_scaled_start_difference"] for r in results
        ),
        input_sha256={
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [
                "normalization-specimen-plan.json",
                "normalization-specimen-refits.json",
                "full-sensitivity-plan.json",
                "full-sensitivity-refits.json",
            ]
        },
    )
    # CSV columns must be shared even when an omitted scale is absent.
    keys = list(dict.fromkeys(k for row in summary for k in row))
    write_csv(
        "dataset-influence-summary.csv",
        [{k: row.get(k) for k in keys} for row in summary],
    )
    write_csv("dataset-influence-residuals.csv", residuals)
    write_csv("dataset-influence-profile-rows.csv", profile_rows)
    write_csv("dataset-influence-groups.csv", groups)
    write_csv("dataset-influence-conrad-native.csv", conrad_native_audit(by_id))
    (DATA / "dataset-influence-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    run()
