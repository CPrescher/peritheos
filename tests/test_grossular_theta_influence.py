"""Retained Conrad diagnostic and theta influence with controlled deletions."""

import csv
import hashlib
import json

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, GAS_CONSTANT
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_conrad_processing import prepare as conrad_rows
from scripts.refit_grossular_angel_2022_dataset_influence import derivatives, group_of
from scripts.refit_grossular_angel_2022_subset import (
    arrays,
    effective_sigma,
    pressure,
    solve,
)
from scripts.refit_grossular_angel_2022_theta_influence import (
    PLAN,
    cases,
    normalization_effects,
    plan,
    prepare,
    profile_diagnostics,
    solve_fixed_theta,
    temperature_bin,
)


def report():
    return json.loads((DATA / "theta-influence-refits.json").read_text())


def test_controlled_selection_retains_conrad_errors_and_scale_identities():
    assert len(cases()) == 30
    assert len(plan()["omissions"]) == 13
    for case in cases():
        rows, starts, config = prepare(case)
        old, _, old_config = conrad_rows(
            dict(normalization=case["normalization"], variant="extra_020")
        )
        lookup = {(r["dataset"], r["row"]): r for r in old}
        deleted = (
            []
            if not case["action"].startswith("omit_")
            else next(
                g["groups"]
                for g in plan()["omissions"]
                if case["action"] == "omit_" + g["id"]
            )
        )
        assert "conrad_ks" not in deleted
        assert len(rows) == sum(group_of(r) not in deleted for r in old)
        assert len(starts) == 2
        for r in rows:
            prior = lookup[(r["dataset"], r["row"])]
            assert {k: v for k, v in r.items() if k != "scale"} == {
                k: v for k, v in prior.items() if k != "scale"
            }
            if r["scale"] is not None:
                assert (
                    config["names"][r["scale"]] == old_config["names"][prior["scale"]]
                )
        assert [r for r in rows if r["dataset"] == "conrad_ks"] == [
            r for r in old if r["dataset"] == "conrad_ks"
        ]
        assert len([r for r in rows if r["dataset"] == "conrad_ks"]) == 3
    assert [temperature_bin(t) for t in [20, 297.9, 298, 699.9, 700, 1350]] == [
        "below298",
        "below298",
        "298to699",
        "298to699",
        "700andabove",
        "700andabove",
    ]


def test_baseline_parity_all_starts_and_conditional_covariance_dimensions():
    r = report()
    assert r["plan_sha256"] == hashlib.sha256(PLAN.read_bytes()).hexdigest()
    assert r["all_starts_converged"] and not r["any_start_at_bound"]
    assert r["all_three_conrad_points_retained"] and not r["author_workflow_reproduced"]
    prior = json.loads((DATA / "conrad-processing-refits.json").read_text())
    for case in r["cases"]:
        a, b = case["fits"]
        np.testing.assert_allclose(
            list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-6
        )
        for fit in case["fits"]:
            covariance = np.array(fit["conditional_unscaled_covariance"])
            assert covariance.shape == (len(case["free_parameters"]),) * 2
            assert np.linalg.eigvalsh(covariance).min() > 0
            assert fit["degrees_of_freedom"] == case["rows"] - len(
                case["free_parameters"]
            )
        if case["action"] == "baseline":
            old = next(
                x
                for x in prior["cases"]
                if x["id"] == case["normalization"] + "__extra_020"
            )
            assert a["parameters"] == old["fits"][0]["parameters"]
        if case["action"] == "published_theta":
            assert a["parameters"]["theta0_k"] == 749.825
            assert "theta0_k" not in case["free_parameters"]
            assert case["fixed_parameters"] == {"theta0_k": 749.825}


def test_constrained_theta_and_deletions_replay_with_stationary_free_scores():
    for action in ["published_theta", "omit_du_vt", "omit_bosenick_vt", "omit_kono_ks"]:
        case = next(c for c in report()["cases"] if c["id"] == "plotted__" + action)
        rows, starts, config = prepare(case)
        fit = (
            solve_fixed_theta(rows, starts[1], config)
            if action == "published_theta"
            else solve(rows, starts[1], **config)
        )
        np.testing.assert_allclose(
            list(fit["parameters"].values()),
            list(case["fits"][0]["parameters"].values()),
            rtol=1e-6,
        )
        if action == "published_theta":
            c = np.array(list(fit["parameters"].values()))
            data = arrays(rows)
            sigma = effective_sigma(c, data)
            residual = (pressure(c, data) - data["P"]) / sigma
            for i in [0, 1, 3, 4, 5]:
                h = config["scaling"][i] * 2e-6
                plus, minus = c.copy(), c.copy()
                plus[i] += h
                minus[i] -= h
                column = (pressure(plus, data) - pressure(minus, data)) / (
                    2 * h * sigma
                )
                score = abs(np.einsum("i,i->", column, residual)) / (
                    np.linalg.norm(column) * np.linalg.norm(residual)
                )
                assert score < 2e-5
            jac = derivatives(c, data, sigma, config["scaling"])[:, [0, 1, 3, 4, 5]]
            scaling = config["scaling"][[0, 1, 3, 4, 5]]
            information = np.einsum("ij,ik->jk", jac * scaling, jac * scaling)
            covariance = np.array(fit["conditional_unscaled_covariance"]) / np.outer(
                scaling, scaling
            )
            np.testing.assert_allclose(
                np.einsum("ij,jk->ik", information, covariance), np.eye(5), atol=1e-6
            )


def test_profile_costs_information_bins_and_normalization_interactions_reconcile():
    r = report()
    by_id = {c["id"]: c for c in r["cases"]}
    for normalization in plan()["normalizations"]:
        profile, rows, groups, bins = profile_diagnostics(normalization, by_id)
        assert len(rows) == 219 and len(groups) == 12
        np.testing.assert_allclose(sum(x["information_fraction"] for x in groups), 1)
        np.testing.assert_allclose(sum(x["information_fraction"] for x in bins), 1)
        np.testing.assert_allclose(
            sum(x["local_upward_theta_contribution_k"] for x in groups),
            profile["local_predicted_theta_shift_k"],
        )
        assert profile["nuisance_orthogonality_max"] < 1e-8
        np.testing.assert_allclose(
            profile["efficient_information_per_k2"],
            r["profiles"][normalization]["efficient_information_per_k2"],
        )
        for g in groups:
            selected = [x for x in rows if x["group"] == g["group"]]
            selected_bins = [x for x in bins if x["group"] == g["group"]]
            assert sum(x["rows"] for x in selected_bins) == g["rows"]
            for curve in ["published_all", "fixed_theta", "free"]:
                cost = sum(x[f"{curve}_standardized_residual"] ** 2 for x in selected)
                np.testing.assert_allclose(cost, g[f"{curve}_sum_squares"])
                np.testing.assert_allclose(
                    sum(x[f"{curve}_sum_squares"] for x in selected_bins), cost
                )
    expected = normalization_effects(by_id)
    assert expected == r["normalization_comparisons"]
    baseline = expected[0]["linear_minus_plotted_theta_k"]
    for x in expected:
        np.testing.assert_allclose(
            x["normalization_interaction_vs_baseline_k"],
            x["linear_minus_plotted_theta_k"] - baseline,
        )
    with (DATA / "theta-influence-residuals.csv").open() as f:
        residuals = list(csv.DictReader(f))
    assert len(residuals) == sum(c["rows"] for c in r["cases"])
    for case in r["cases"]:
        selected = [x for x in residuals if x["case"] == case["id"]]
        cost = sum(
            (float(x["pressure_residual_gpa"]) / float(x["effective_error_gpa"])) ** 2
            for x in selected
        )
        np.testing.assert_allclose(cost, case["fits"][0]["objective_sum_squares"])


def test_native_du_checks_and_all_prior_artifacts_are_preserved():
    with (DATA / "theta-influence-native-comparisons.csv").open() as f:
        native = list(csv.DictReader(f))
    assert len(native) == 438
    r = report()
    by_id = {c["id"]: c for c in r["cases"]}
    p, _ = source("grossular")
    for normalization in plan()["normalizations"]:
        # Check both native coordinates against fast pressure inversion for
        # every Du state and every modulus state of the constrained fit.
        case = by_id[f"{normalization}__published_theta"]
        rows, _, _ = prepare(case)
        c = np.array(list(case["fits"][0]["parameters"].values()))
        selected = [x for x in native if x["normalization"] == normalization]
        synthetic = []
        for row, exported in zip(rows, selected):
            assert row["dataset"] == exported["dataset"] and row["row"] == int(
                exported["source_row"]
            )
            for curve in ["published_all", "baseline", "fixed_theta"]:
                np.testing.assert_allclose(
                    float(exported[f"observed_minus_{curve}"]),
                    row["Y"] - float(exported[curve]),
                )
            if row["dataset"] == "du_vt" or row["KS"]:
                synthetic.append(row | dict(Y=float(exported["fixed_theta"])))
        np.testing.assert_allclose(
            pressure(c, arrays(synthetic)), [x["P"] for x in synthetic], atol=1e-7
        )
        first_du = next(
            x for x in selected if x["dataset"] == "du_vt" and int(x["source_row"]) == 2
        )
        model = p | dict(zip([2, 3, 11, 18], c[:4]))
        independent = properties(model, 0.0001, 393.7, GAS_CONSTANT)["volume_cm3_mol"]
        np.testing.assert_allclose(float(first_du["fixed_theta"]), independent)
    for filename in [
        "artifact-manifest.json",
        "followup-artifact-manifest.json",
        "expanded-artifact-manifest.json",
        "normalized-artifact-manifest.json",
        "full-sensitivity-artifact-manifest.json",
        "normalization-specimen-artifact-manifest.json",
        "dataset-influence-artifact-manifest.json",
        "conrad-processing-artifact-manifest.json",
    ]:
        for entry in json.loads((DATA / filename).read_text())["files"]:
            assert (
                hashlib.sha256((DATA / entry["file"]).read_bytes()).hexdigest()
                == entry["sha256"]
            )
