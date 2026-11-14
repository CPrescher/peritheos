"""Declared Du/Bosenick processing controls; never writes source records.

Repository-only: python -m scripts.refit_grossular_du_bosenick.
See du-bosenick-plan.json for assumptions, covariance limits and compute cap.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, observations
from scripts.refit_grossular_angel_2022_conrad_processing import (
    prepare as reference_rows,
)
from scripts.refit_grossular_angel_2022_dataset_influence import derivatives
from scripts.refit_grossular_angel_2022_full_sensitivity import write_csv
from scripts.refit_grossular_angel_2022_normalization_specimen import (
    denominators,
    line_at,
)
from scripts.refit_grossular_angel_2022_subset import arrays, effective_sigma, pressure

PLAN = DATA / "du-bosenick-plan.json"


def plan():
    return json.loads(PLAN.read_text())


def cases():
    result = []
    for background in plan()["backgrounds"]:
        for anchor in plan()["bosenick_anchors"]:
            combinations = [(c, "table_nacl") for c in plan()["bosenick_covariances"]]
            combinations += [
                ("full_independent_T", t)
                for t in plan()["du_treatments"]
                if t != "table_nacl"
            ]
            for covariance, du in combinations:
                result.append(
                    dict(
                        id=f"{background}__{anchor}__{covariance}__{du}",
                        background=background,
                        anchor=anchor,
                        covariance=covariance,
                        du=du,
                    )
                )
    return result


def bosenick_anchor(anchor):
    raw = observations("bosenick-geiger-1997-table-2.csv")
    x = np.array([float(r["molar_volume_cm3_mol"]) for r in raw])
    t = np.array([float(r["temperature_k"]) for r in raw])
    sx = np.array([float(r["molar_volume_sd_cm3_mol"]) for r in raw])
    st = np.where(t >= 100, 3.0, 5.0)
    w, dt = np.zeros(len(raw)), np.zeros(len(raw))
    if anchor == "measured_295":
        w[0] = 1
    elif anchor == "linear_298":
        w[1] = -(298 - t[0]) / (t[0] - t[1])
        w[0] = 1 - w[1]
    else:
        raise ValueError(anchor)
    denominator = float(w @ x)
    if anchor == "linear_298":
        dt[0] = (x[1] - denominator) / (t[0] - t[1])
        dt[1] = (denominator - x[0]) / (t[0] - t[1])
    y = 125.35 * x[1:] / denominator
    jy = np.zeros((19, 20))
    jy[np.arange(19), np.arange(1, 20)] = 125.35 / denominator
    jy -= np.outer(y / denominator, w)
    yt = -np.outer(y / denominator, dt)
    return dict(
        x=x,
        t=t,
        sx=sx,
        st=st,
        weights=w,
        denominator=denominator,
        denominator_temperature_derivative=dt,
        y=y,
        volume_jacobian=jy,
        temperature_jacobian=yt,
    )


def du_temperature(temperature, treatment):
    celsius = temperature - 273.15
    if treatment == "tc_reconstructed":
        # Stable small root of 3.6e-5*T1^2-.982*T1+T2+.565=0.
        discriminant = 0.982**2 - 4 * 3.6e-5 * (celsius + 0.565)
        value = 2 * (celsius + 0.565) / (0.982 + np.sqrt(discriminant))
        return float(value + 273.15), float(1 / (0.982 - 7.2e-5 * value))
    if treatment == "double_corrected":
        return (
            -3.6e-5 * celsius**2 + 0.982 * celsius - 0.565 + 273.15,
            0.982 - 7.2e-5 * celsius,
        )
    return temperature, 1.0


def prepare(case):
    rows, starts, config = reference_rows(
        dict(normalization=case["background"], variant="extra_020")
    )
    bose = bosenick_anchor(case["anchor"])
    du_raw = observations("du-2015-preprint-table-2.csv")
    du = {int(r["source_row"]): r for r in du_raw}
    du_t = [du_temperature(float(r["temperature_k"]), case["du"])[0] for r in du_raw]
    du_x = [
        float(r["lattice_a_angstrom"]) ** 3
        if case["du"] == "lattice"
        else float(r["cell_volume_a3"])
        for r in du_raw
    ]
    d = denominators(case["background"])["du_vt"]
    if case["background"] == "linear_298":
        d = line_at(298, du_t[0], du_x[0], du_t[1], du_x[1])
    elif case["du"] == "lattice":
        d *= du_x[0] / float(du_raw[0]["cell_volume_a3"])
    for r in rows:
        if r["dataset"] == "bosenick_vt":
            i = r["row"] - 1
            r["Y"] = float(bose["y"][i - 1])
            r["sY"] = float(125.35 * bose["sx"][i] / bose["denominator"])
        elif r["dataset"] == "du_vt":
            raw = du[r["row"]]
            index = r["row"] - 1
            r["T"], multiplier = du_temperature(float(raw["temperature_k"]), case["du"])
            r["sT"] = 0.2 * multiplier
            r["Y"] = 125.35 * du_x[index] / d
            sy = (
                (
                    3
                    * float(raw["lattice_a_angstrom"]) ** 2
                    * float(raw["lattice_a_sd_angstrom"])
                )
                if case["du"] == "lattice"
                else float(raw["cell_volume_sd_a3"])
            )
            r["sY"] = 125.35 * sy / d
    return (
        rows,
        starts,
        config,
        dict(bosenick_denominator=bose["denominator"], du_denominator=d),
    )


def weight_state(c, data, case):
    sigma = effective_sigma(c, data)
    py = (
        pressure(c, data, y=data["Y"] + 0.0001)
        - pressure(c, data, y=data["Y"] - 0.0001)
    ) / 0.0002
    pt = (
        pressure(c, data, t=data["T"] + 0.01) - pressure(c, data, t=data["T"] - 0.01)
    ) / 0.02
    blocks = []
    for dataset in ["bosenick_vt", "du_vt"]:
        idx = np.array([i for i, name in enumerate(data["dataset"]) if name == dataset])
        covariance = np.diag(sigma[idx] ** 2)
        if dataset == "bosenick_vt" and case["covariance"] != "diagonal":
            a = bosenick_anchor(case["anchor"])
            jv = py[idx, None] * a["volume_jacobian"]
            weighted_jv = jv * a["sx"]
            covariance = np.einsum("ik,jk->ij", weighted_jv, weighted_jv)
            if case["covariance"] == "volume_anchor":
                covariance += np.diag(
                    data["sP"][idx] ** 2 + (pt[idx] * data["sT"][idx]) ** 2
                )
            else:
                jt = py[idx, None] * a["temperature_jacobian"]
                jt[np.arange(19), np.arange(1, 20)] += pt[idx]
                weighted_jt = jt * a["st"]
                if case["covariance"] == "full_shared_T":
                    vector = np.sum(weighted_jt, axis=1)
                    covariance += np.outer(vector, vector)
                else:
                    covariance += np.einsum("ik,jk->ij", weighted_jt, weighted_jt)
                covariance += np.diag(data["sP"][idx] ** 2)
        elif dataset == "du_vt" and case["du"] == "shared_T":
            covariance = np.diag(
                data["sP"][idx] ** 2 + (py[idx] * data["sY"][idx]) ** 2
            )
            covariance += np.outer(pt[idx] * data["sT"][idx], pt[idx] * data["sT"][idx])
        factor = np.linalg.cholesky(covariance)  # Never add jitter or clip eigenvalues.
        blocks.append(
            dict(dataset=dataset, indices=idx, covariance=covariance, factor=factor)
        )
    return sigma, blocks


def whiten(values, state):
    sigma, blocks = state
    result = values / sigma if values.ndim == 1 else values / sigma[:, None]
    for block in blocks:
        idx = block["indices"]
        result[idx] = solve_triangular(block["factor"], values[idx], lower=True)
    return result


def solve_case(rows, initial, config, case):
    data = arrays(rows) | {"dataset": [r["dataset"] for r in rows]}
    c = np.array(initial, dtype=float)
    trace, converged = [], False
    for iteration in range(30):
        state = weight_state(c, data, case)
        fit = least_squares(
            lambda cc: whiten(pressure(cc, data) - data["P"], state),
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
                success=bool(fit.success),
                max_scaled_shift=shift,
                nfev=fit.nfev,
            )
        )
        if shift < 1e-8 and fit.success:
            converged = True
            break
    state = weight_state(c, data, case)
    jac = whiten(derivatives(c, data, np.ones(len(rows)), config["scaling"]), state)
    _, singular, vectors = np.linalg.svd(jac * config["scaling"], full_matrices=False)
    if np.sum(singular > singular[0] * 1e-12) != len(c):
        raise ValueError("Covariance is rank deficient")
    inverse = vectors.T / singular
    covariance = np.einsum("ik,jk->ij", inverse, inverse) * np.outer(
        config["scaling"], config["scaling"]
    )
    residual = pressure(c, data) - data["P"]
    z = whiten(residual, state)
    objective = float(z @ z)
    widths = state[0].copy()
    blocks = []
    for b in state[1]:
        widths[b["indices"]] = np.sqrt(np.diag(b["covariance"]))
        blocks.append(
            dict(
                dataset=b["dataset"],
                source_rows=[rows[i]["row"] for i in b["indices"]],
                pressure_covariance_gpa2=b["covariance"].tolist(),
                cholesky=b["factor"].tolist(),
            )
        )
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
        pressure_covariance_blocks=blocks,
        rows=[
            r
            | dict(
                pressure_residual_gpa=float(v),
                marginal_effective_error_gpa=float(s),
                whitened_pressure_residual=float(w),
            )
            for r, v, s, w in zip(rows, residual, widths, z)
        ],
    )


def run():
    results, residuals, inputs, summary = [], [], [], []
    for case in cases():
        rows, starts, config, reduction = prepare(case)
        fits = [solve_case(rows, start, config, case) for start in starts]
        residuals.extend(dict(case=case["id"]) | r for r in fits[0]["rows"])
        inputs.extend(
            dict(case=case["id"]) | r
            for r in rows
            if r["dataset"] in ["du_vt", "bosenick_vt", "conrad_ks"]
        )
        for f in fits:
            del f["rows"]
        a, b = (np.array(list(f["parameters"].values())) for f in fits)
        difference = float(np.max(np.abs((a - b) / config["scaling"])))
        results.append(
            case
            | dict(
                reduction=reduction,
                rows=len(rows),
                fits=fits,
                maximum_scaled_start_difference=difference,
            )
        )
        print(
            f"{case['id']}: theta={a[2]:.5f}, Kprime={a[1]:.6f}, converged={all(f['numerically_converged'] for f in fits)}",
            flush=True,
        )
    prior = json.loads((DATA / "theta-influence-refits.json").read_text())
    p, _ = source("grossular")
    for r in results:
        values = r["fits"][0]["parameters"]
        baseline = next(
            c["fits"][0]["parameters"]
            for c in prior["cases"]
            if c["id"] == r["background"] + "__baseline"
        )
        r["effects"] = {
            k: dict(minus_baseline=v - baseline[k]) for k, v in values.items()
        }
        summary.append(
            {k: r[k] for k in ["id", "background", "anchor", "covariance", "du"]}
            | values
            | dict(
                theta_minus_source_k=values["theta0_k"] - p[11],
                kprime_minus_source=values["K0_prime"] - p[3],
                theta_minus_baseline_k=values["theta0_k"] - baseline["theta0_k"],
                objective=r["fits"][0]["objective_sum_squares"],
            )
        )
    report = dict(
        audit_date=plan()["declared_on"],
        scope="conditional_Du_preprocessing_Bosenick_anchor_covariance",
        author_workflow_reproduced=False,
        physical_validation=False,
        catalog_changes=False,
        author_contacted=False,
        all219rows_retained=True,
        conrad_diagnostic="Formal upper widths plus independent2%KS, unchanged",
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        input_sha256=plan()["input_sha256"],
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
    )
    write_csv("du-bosenick-summary.csv", summary)
    write_csv("du-bosenick-residuals.csv", residuals)
    write_csv("du-bosenick-inputs.csv", inputs)
    anchor_receipt = {
        a: {
            k: v.tolist() if isinstance(v, np.ndarray) else v
            for k, v in bosenick_anchor(a).items()
        }
        for a in plan()["bosenick_anchors"]
    }
    (DATA / "du-bosenick-anchor-jacobians.json").write_text(
        json.dumps(anchor_receipt, indent=2, allow_nan=False) + "\n"
    )
    (DATA / "du-bosenick-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    run()
