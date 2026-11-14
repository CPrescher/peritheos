"""Check anchor/error algebra and retained-data covariance diagnostics."""

import csv
import hashlib
import json

import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA
from scripts.refit_grossular_angel_2022_conrad_processing import (
    prepare as reference_rows,
)
from scripts.refit_grossular_angel_2022_dataset_influence import derivatives
from scripts.refit_grossular_angel_2022_subset import arrays, pressure
from scripts.refit_grossular_du_bosenick import (
    PLAN,
    bosenick_anchor,
    cases,
    du_temperature,
    prepare,
    solve_case,
    weight_state,
    whiten,
)


def report():
    return json.loads((DATA / "du-bosenick-refits.json").read_text())


def test_all219rows_conrad_widths_and_background_are_preserved():
    assert len(cases()) == 32
    for case in cases():
        rows, starts, config, reduction = prepare(case)
        old, _, _ = reference_rows(
            dict(normalization=case["background"], variant="extra_020")
        )
        assert len(rows) == 219 and len(starts) == 2 and len(config["names"]) == 6
        assert [r for r in rows if r["dataset"] not in ["du_vt", "bosenick_vt"]] == [
            r for r in old if r["dataset"] not in ["du_vt", "bosenick_vt"]
        ]
        assert sum(r["dataset"] == "conrad_ks" for r in rows) == 3
        assert sum(r["dataset"] == "bosenick_vt" for r in rows) == 19
        assert sum(r["dataset"] == "du_vt" for r in rows) == 13
        anchor = bosenick_anchor(case["anchor"])
        np.testing.assert_allclose(
            reduction["bosenick_denominator"], anchor["denominator"]
        )
        for row in rows:
            if row["dataset"] == "bosenick_vt":
                i = row["row"] - 1
                np.testing.assert_allclose(
                    row["Y"], 125.35 * anchor["x"][i] / anchor["denominator"]
                )
                np.testing.assert_allclose(
                    row["sY"], 125.35 * anchor["sx"][i] / anchor["denominator"]
                )


def test_celsius_calibration_round_trip_derivatives_and_anchor_jacobians():
    for temperature in [296.1, 393.7, 447.1, 859.6]:
        tc, derivative = du_temperature(temperature, "tc_reconstructed")
        corrected, forward = du_temperature(tc, "double_corrected")
        np.testing.assert_allclose(corrected, temperature, atol=1e-10)
        np.testing.assert_allclose(derivative * forward, 1)
        h = 0.01
        finite = (
            du_temperature(temperature + h, "tc_reconstructed")[0]
            - du_temperature(temperature - h, "tc_reconstructed")[0]
        ) / (2 * h)
        np.testing.assert_allclose(finite, derivative, rtol=1e-7)
    for name in ["measured_295", "linear_298"]:
        a = bosenick_anchor(name)
        np.testing.assert_allclose(
            np.einsum("ij,j->i", a["volume_jacobian"], a["x"]), 0, atol=1e-10
        )

        # Independently differentiate the complete ratio and its nominal-T
        # extrapolation, including the retained270K numerator/denominator.
        def values(x, t):
            d = (
                x[0]
                if name == "measured_295"
                else x[0] + (298 - t[0]) * (x[0] - x[1]) / (t[0] - t[1])
            )
            return 125.35 * x[1:] / d

        for variable, jacobian, step in [
            ("x", "volume_jacobian", 1e-4),
            ("t", "temperature_jacobian", 1e-3),
        ]:
            finite = []
            for j in range(20):
                plus, minus = a[variable].copy(), a[variable].copy()
                plus[j] += step
                minus[j] -= step
                if variable == "x":
                    finite.append(
                        (values(plus, a["t"]) - values(minus, a["t"])) / (2 * step)
                    )
                else:
                    finite.append(
                        (values(a["x"], plus) - values(a["x"], minus)) / (2 * step)
                    )
            np.testing.assert_allclose(
                np.array(finite).T, a[jacobian], rtol=2e-5, atol=1e-8
            )
    linear = bosenick_anchor("linear_298")
    np.testing.assert_allclose(linear["denominator"], 125.30104)
    assert linear["weights"][1] < 0  # This row must not be treated as independent of D.


def test_dense_gls_complete_pressure_jacobian_and_inverse_information():
    case = next(
        c
        for c in report()["cases"]
        if c["id"] == "linear_298__linear_298__full_independent_T__table_nacl"
    )
    rows, _, config, _ = prepare(case)
    data = arrays(rows) | {"dataset": [r["dataset"] for r in rows]}
    c = np.array(list(case["fits"][0]["parameters"].values()))
    state = weight_state(c, data, case)
    residual = pressure(c, data) - data["P"]
    z = whiten(residual, state)
    for b in state[1]:
        idx, cov, factor = b["indices"], b["covariance"], b["factor"]
        np.testing.assert_allclose(np.einsum("ik,jk->ij", factor, factor), cov)
        dense = np.einsum("i,i->", residual[idx], np.linalg.solve(cov, residual[idx]))
        np.testing.assert_allclose(np.sum(z[idx] ** 2), dense)
    b = state[1][0]
    a = bosenick_anchor("linear_298")
    idx = b["indices"]
    # Numerical pressure Jacobians include direct T effects and anchor-T
    # effects; compare their propagated covariance with the production block.
    columns = []
    widths = []
    for variable, step in [("x", 1e-4), ("t", 1e-3)]:
        for j in range(20):
            predictions = []
            for sign in [1, -1]:
                x, t = a["x"].copy(), a["t"].copy()
                (x if variable == "x" else t)[j] += sign * step
                d = x[0] + (298 - t[0]) * (x[0] - x[1]) / (t[0] - t[1])
                ynew, tnew = data["Y"].copy(), data["T"].copy()
                ynew[idx] = 125.35 * x[1:] / d
                tnew[idx] = t[1:]
                predictions.append(pressure(c, data, y=ynew, t=tnew)[idx])
            columns.append((predictions[0] - predictions[1]) / (2 * step))
            widths.append((a["sx"] if variable == "x" else a["st"])[j])
    weighted = np.array(columns).T * np.array(widths)
    np.testing.assert_allclose(
        np.einsum("ik,jk->ij", weighted, weighted),
        b["covariance"],
        rtol=2e-5,
        atol=1e-9,
    )
    jac = whiten(derivatives(c, data, np.ones(len(rows)), config["scaling"]), state)
    score = np.abs(np.einsum("ij,i->j", jac, z)) / (
        np.linalg.norm(jac, axis=0) * np.linalg.norm(z)
    )
    assert np.max(score) < 2e-5
    scaled = jac * config["scaling"]
    info = np.einsum("ij,ik->jk", scaled, scaled)
    covariance = np.array(
        case["fits"][0]["conditional_unscaled_covariance"]
    ) / np.outer(config["scaling"], config["scaling"])
    np.testing.assert_allclose(
        np.einsum("ij,jk->ik", info, covariance), np.eye(6), atol=1e-6
    )


def test_baseline_parity_replays_covariance_and_all_costs():
    r = report()
    assert r["plan_sha256"] == hashlib.sha256(PLAN.read_bytes()).hexdigest()
    assert r["all_starts_converged"] and not r["any_start_at_bound"]
    prior = json.loads((DATA / "theta-influence-refits.json").read_text())
    with (DATA / "du-bosenick-residuals.csv").open() as f:
        residuals = list(csv.DictReader(f))
    assert len(residuals) == 32 * 219
    for case in r["cases"]:
        a, b = case["fits"]
        np.testing.assert_allclose(
            list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-6
        )
        for fit in case["fits"]:
            assert np.linalg.eigvalsh(fit["conditional_unscaled_covariance"]).min() > 0
            assert fit["degrees_of_freedom"] == 213
            for block in fit["pressure_covariance_blocks"]:
                assert np.linalg.eigvalsh(block["pressure_covariance_gpa2"]).min() > 0
        selected = [x for x in residuals if x["case"] == case["id"]]
        np.testing.assert_allclose(
            sum(float(x["whitened_pressure_residual"]) ** 2 for x in selected),
            a["objective_sum_squares"],
        )
        if (
            case["covariance"] == "diagonal"
            and case["du"] == "table_nacl"
            and case["anchor"]
            == ("measured_295" if case["background"] == "plotted" else "linear_298")
        ):
            old = next(
                x
                for x in prior["cases"]
                if x["id"] == case["background"] + "__baseline"
            )
            np.testing.assert_allclose(
                list(a["parameters"].values()),
                list(old["fits"][0]["parameters"].values()),
                rtol=1e-6,
            )
    for ident in [
        "plotted__measured_295__volume_anchor__table_nacl",
        "linear_298__linear_298__full_independent_T__table_nacl",
        "linear_298__linear_298__full_shared_T__table_nacl",
        "linear_298__linear_298__full_independent_T__shared_T",
    ]:
        case = next(c for c in r["cases"] if c["id"] == ident)
        rows, starts, config, _ = prepare(case)
        replay = solve_case(rows, starts[1], config, case)
        np.testing.assert_allclose(
            list(replay["parameters"].values()),
            list(case["fits"][0]["parameters"].values()),
            rtol=1e-6,
        )


def test_previous_artifacts_untouched_and_inputs_exported():
    with (DATA / "du-bosenick-inputs.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 32 * (19 + 13 + 3)
    assert sum(r["dataset"] == "conrad_ks" for r in rows) == 32 * 3
    for name in [
        "artifact-manifest.json",
        "followup-artifact-manifest.json",
        "expanded-artifact-manifest.json",
        "normalized-artifact-manifest.json",
        "full-sensitivity-artifact-manifest.json",
        "normalization-specimen-artifact-manifest.json",
        "dataset-influence-artifact-manifest.json",
        "conrad-processing-artifact-manifest.json",
        "theta-influence-artifact-manifest.json",
    ]:
        for item in json.loads((DATA / name).read_text())["files"]:
            assert (
                hashlib.sha256((DATA / item["file"]).read_bytes()).hexdigest()
                == item["sha256"]
            )
