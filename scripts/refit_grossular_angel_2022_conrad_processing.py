"""Fixed Conrad input/error controls, with every selected point retained.

Repository-only: python -m scripts.refit_grossular_angel_2022_conrad_processing.
All scenarios are predeclared in conrad-processing-plan.json.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import (
    DATA,
    GAS_CONSTANT,
    observations,
)
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_dataset_influence import derivatives
from scripts.refit_grossular_angel_2022_full_sensitivity import write_csv
from scripts.refit_grossular_angel_2022_normalization_specimen import prepared_case
from scripts.refit_grossular_angel_2022_subset import (
    arrays,
    effective_sigma,
    pressure,
    solve,
)

PLAN = DATA / "conrad-processing-plan.json"


def plan():
    return json.loads(PLAN.read_text())


def cases():
    return [
        dict(id=f"{n}__{v['id']}", normalization=n, variant=v["id"])
        for n in plan()["normalizations"]
        for v in plan()["variants"]
    ]


def density_controls():
    p, _ = source("grossular")
    original = {
        int(r["source_row"]): r for r in observations("conrad-1999-table-6.csv")
    }
    result = []
    for point in plan()["density_inputs"]:
        r = original[point["source_row"]]
        new_rho = (
            3.60
            * 125.35
            / properties(p, point["P_gpa"], 298, GAS_CONSTANT)["volume_cm3_mol"]
        )
        for channel, modulus in [("vp", "c11_gpa"), ("vs", "c44_gpa")]:
            velocity = point[f"{channel}_calc_km_s"]
            cij = float(r[modulus])
            rho = cij / velocity**2
            result.append(
                dict(
                    source_row=point["source_row"],
                    pressure_gpa=point["P_gpa"],
                    channel=channel,
                    calculated_velocity_km_s=velocity,
                    inferred_density_g_cm3=rho,
                    rounding_density_lower=cij / (velocity + 0.005) ** 2,
                    rounding_density_upper=cij / (velocity - 0.005) ** 2,
                    published_eos_calculated_density=new_rho,
                    calculated_modulus_multiplier=new_rho / rho,
                )
            )
    return result


def prepare(case):
    rows, starts, config = prepared_case(
        dict(normalization=case["normalization"], isaak="as_plotted")
    )
    literal = {int(r["source_row"]): r for r in observations("conrad-1999-table-6.csv")}
    densities = {(r["source_row"], r["channel"]): r for r in density_controls()}
    for r in rows:
        if r["dataset"] != "conrad_ks":
            continue
        variant = case["variant"]
        raw = literal[r["row"]]
        a, b = float(raw["c11_sd_gpa"]), 2 * float(raw["c12_sd_gpa"])
        if variant == "cij_negative_correlation":
            r["sY"] = abs(a - b) / 3
        elif variant == "cij_independent":
            r["sY"] = np.hypot(a, b) / 3
        elif variant.startswith("extra_"):
            fraction = int(variant.split("_")[1]) / 1000
            r["sY"] = float(np.hypot(r["sY"], fraction * r["Y"]))
        elif variant in ["pressure_02", "pressure_05"]:
            r["sP"] = 0.2 if variant == "pressure_02" else 0.5
        elif variant == "round_01":
            r["Y"] = round(r["Y"], 1)
        elif variant.startswith("density_"):
            factor = densities[(r["row"], variant.split("_")[1])][
                "calculated_modulus_multiplier"
            ]
            r["Y"] *= factor
            r["sY"] *= factor
        elif variant == "free_conrad_scale":
            r["scale"] = 6
        elif variant not in ["baseline", "shared_020"]:
            raise ValueError(f"Unknown variant {variant}")
    if case["variant"] == "free_conrad_scale":
        config["names"].append("scale_conrad_ks")
        for key, value in [("scaling", 1), ("lower", 0.95), ("upper", 1.05)]:
            config[key] = np.append(config[key], value)
        for start in starts:
            start.append(1)
    return rows, starts, config


def weight_state(c, data, *, shared=False):
    sigma = effective_sigma(c, data)
    indices = np.array([i for i, r in enumerate(data["dataset"]) if r == "conrad_ks"])
    block = np.diag(sigma[indices] ** 2)
    if shared:
        py = (
            pressure(c, data, y=data["Y"] + 0.0001)
            - pressure(c, data, y=data["Y"] - 0.0001)
        ) / 0.0002
        u = py[indices] * 0.02 * data["Y"][indices]
        block += np.outer(u, u)
    return sigma, indices, block, np.linalg.cholesky(block)


def whiten(values, state):
    sigma, indices, _, factor = state
    result = values / sigma if values.ndim == 1 else values / sigma[:, None]
    result[indices] = solve_triangular(factor, values[indices], lower=True)
    return result


def solve_shared(rows, initial, config):
    """IRLS with a three-row shared calibration covariance block."""
    data = arrays(rows) | {"dataset": [r["dataset"] for r in rows]}
    c = np.array(initial, dtype=float)
    trace, converged = [], False
    for iteration in range(30):
        weights = weight_state(c, data, shared=True)
        fit = least_squares(
            lambda cc: whiten(pressure(cc, data) - data["P"], weights),
            c,
            bounds=(config["lower"], config["upper"]),
            jac="3-point",
            x_scale=config["scaling"],
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-9,
            max_nfev=400,
        )
        shift = float(np.max(np.abs((fit.x - c) / config["scaling"])))
        c = fit.x
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
    weights = weight_state(c, data, shared=True)
    # Derivatives() divides by its argument; ones give unweighted dP/dc.
    jac = whiten(derivatives(c, data, np.ones(len(rows)), config["scaling"]), weights)
    _, singular, vectors = np.linalg.svd(jac * config["scaling"], full_matrices=False)
    if np.sum(singular > singular[0] * 1e-12) != len(c):
        raise ValueError("Shared-error covariance is rank deficient")
    inverse = vectors.T / singular
    covariance = np.einsum("ik,jk->ij", inverse, inverse) * np.outer(
        config["scaling"], config["scaling"]
    )
    objective = float(np.sum(whiten(pressure(c, data) - data["P"], weights) ** 2))
    return dict(
        parameters=dict(zip(config["names"], c.tolist())),
        numerically_converged=converged,
        at_parameter_bound=[
            config["names"][i]
            for i in range(len(c))
            if min(c[i] - config["lower"][i], config["upper"][i] - c[i])
            < 1e-5 * config["scaling"][i]
        ],
        irls_trace=trace,
        objective_sum_squares=objective,
        degrees_of_freedom=len(rows) - len(c),
        conditional_reduced_chi_square=objective / (len(rows) - len(c)),
        covariance_order=config["names"],
        conditional_unscaled_covariance=covariance.tolist(),
        conditional_unscaled_widths=dict(
            zip(config["names"], np.sqrt(np.diag(covariance)).tolist())
        ),
        scaled_jacobian_singular_values=singular.tolist(),
    )


def residual_receipt(rows, fit, shared):
    data = arrays(rows) | {"dataset": [r["dataset"] for r in rows]}
    c = np.array(list(fit["parameters"].values()))
    state = weight_state(c, data, shared=shared)
    residual = pressure(c, data) - data["P"]
    weighted = whiten(residual, state)
    widths = state[0].copy()
    widths[state[1]] = np.sqrt(np.diag(state[2]))
    exported = [
        r
        | dict(
            pressure_residual_gpa=float(v),
            marginal_effective_error_gpa=float(s),
            whitened_pressure_residual=float(w),
        )
        for r, v, s, w in zip(rows, residual, widths, weighted)
    ]
    return exported, dict(
        source_rows=[rows[i]["row"] for i in state[1]],
        covariance_pressure_gpa2=state[2].tolist(),
        lower_cholesky=state[3].tolist(),
    )


def run():
    results, residuals, inputs, summaries = [], [], [], []
    for case in cases():
        rows, starts, config = prepare(case)
        shared = case["variant"] == "shared_020"
        fits = [
            solve_shared(rows, start, config)
            if shared
            else solve(rows, start, **config)
            for start in starts
        ]
        for fit in fits:
            fit.pop("rows", None)
        exported, block = residual_receipt(rows, fits[0], shared)
        residuals.extend(dict(case=case["id"]) | r for r in exported)
        inputs.extend(
            dict(
                case=case["id"],
                input_kind="calculated_density_reduction"
                if case["variant"].startswith("density_")
                else "diagnostic_rounding"
                if case["variant"] == "round_01"
                else "literal_derived_modulus",
            )
            | r
            for r in rows
            if r["dataset"] == "conrad_ks"
        )
        a, b = (np.array(list(f["parameters"].values())) for f in fits)
        difference = float(np.max(np.abs((a - b) / config["scaling"])))
        results.append(
            case
            | dict(
                rows=len(rows),
                free_parameters=config["names"],
                fits=fits,
                conrad_error_block=block,
                maximum_scaled_start_difference=difference,
            )
        )
        print(
            f"{case['id']}: Kprime={a[1]:.7f}, converged={all(f['numerically_converged'] for f in fits)}",
            flush=True,
        )
    p, _ = source("grossular")
    source_values = dict(
        zip(
            [
                "K0_gpa",
                "K0_prime",
                "theta0_k",
                "gamma0",
                "scale_kono_ks",
                "scale_gwanmesia_ks",
            ],
            [p[2], p[3], p[11], p[18], p[53], p[54]],
        )
    )
    baseline = {
        n: next(
            r["fits"][0]["parameters"]
            for r in results
            if r["normalization"] == n and r["variant"] == "baseline"
        )
        for n in plan()["normalizations"]
    }
    for r in results:
        values = r["fits"][0]["parameters"]
        r["effects"] = {
            name: dict(
                minus_baseline=values[name] - baseline[r["normalization"]][name],
                minus_source=values[name] - value,
            )
            for name, value in source_values.items()
        }
        summaries.append(
            {k: r[k] for k in ["id", "normalization", "variant", "rows"]}
            | values
            | dict(
                kprime_minus_source=values["K0_prime"] - p[3],
                kprime_minus_baseline=values["K0_prime"]
                - baseline[r["normalization"]]["K0_prime"],
                objective=r["fits"][0]["objective_sum_squares"],
                free_parameters=len(r["free_parameters"]),
            )
        )
    report = dict(
        audit_date=plan()["declared_on"],
        scope="conrad_preprocessing_and_error_diagnostic",
        author_workflow_reproduced=False,
        physical_validation=False,
        catalog_changes=False,
        author_contacted=False,
        all_three_conrad_rows_retained=True,
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        cases=results,
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
            f: hashlib.sha256((DATA / f).read_bytes()).hexdigest()
            for f in [
                "conrad-1999-table-6.csv",
                "normalization-specimen-plan.json",
                "normalization-specimen-refits.json",
                "dataset-influence-refits.json",
            ]
        },
    )
    keys = list(dict.fromkeys(k for r in summaries for k in r))
    write_csv(
        "conrad-processing-summary.csv",
        [{k: r.get(k) for k in keys} for r in summaries],
    )
    write_csv("conrad-processing-residuals.csv", residuals)
    write_csv("conrad-processing-inputs.csv", inputs)
    write_csv("conrad-processing-density-controls.csv", density_controls())
    (DATA / "conrad-processing-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    run()
