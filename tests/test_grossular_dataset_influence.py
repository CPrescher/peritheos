"""Check the scientific controls, constrained fit and influence accounting."""

import csv
import hashlib
import json

import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA
from scripts.refit_grossular_angel_2022_dataset_influence import (
    PLAN,
    cases,
    conrad_native_audit,
    group_of,
    plan,
    prepare,
    profile_diagnostics,
    solve_fixed_kprime,
)
from scripts.refit_grossular_angel_2022_normalization_specimen import prepared_case
from scripts.refit_grossular_angel_2022_subset import (
    arrays,
    effective_sigma,
    pressure,
    solve,
)


def report():
    return json.loads((DATA / "dataset-influence-refits.json").read_text())


def test_deletions_preserve_all_retained_values_and_remove_unused_scales():
    assert len(cases()) == 32
    omissions = plan()["omissions"]
    single = omissions[:12]
    groups = [g["groups"][0] for g in single]
    assert len(set(groups)) == 12 and all(len(g["groups"]) == 1 for g in single)
    for case in cases():
        rows, starts, config = prepare(case)
        original, _, old_config = prepared_case(
            dict(normalization=case["normalization"], isaak="as_plotted")
        )
        lookup = {(r["dataset"], r["row"]): r for r in original}
        assert len(starts) == 2
        removed = (
            []
            if not case["action"].startswith("omit_")
            else next(
                g["groups"] for g in omissions if case["action"] == "omit_" + g["id"]
            )
        )
        assert len(rows) == sum(group_of(r) not in removed for r in original)
        for r in rows:
            old = lookup[(r["dataset"], r["row"])]
            assert {k: v for k, v in r.items() if k != "scale"} == {
                k: v for k, v in old.items() if k != "scale"
            }
            if old["scale"] is not None:
                assert config["names"][r["scale"]] == old_config["names"][old["scale"]]
            else:
                assert r["scale"] is None
        for dataset in ["kono", "gwanmesia"]:
            assert (f"scale_{dataset}_ks" in config["names"]) == any(
                r["dataset"] == f"{dataset}_ks" for r in rows
            )


def test_all_fits_and_baselines_and_covariance_have_declared_scope():
    r = report()
    assert r["plan_sha256"] == hashlib.sha256(PLAN.read_bytes()).hexdigest()
    assert r["all_starts_converged"] and not r["any_start_at_bound"]
    assert not r["author_workflow_reproduced"] and not r["catalog_changes"]
    prior = json.loads((DATA / "normalization-specimen-refits.json").read_text())
    for case in r["cases"]:
        a, b = case["fits"]
        assert case["maximum_scaled_start_difference"] < 1e-6
        np.testing.assert_allclose(
            list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-6
        )
        if case["action"] == "baseline":
            old = next(
                c
                for c in prior["cases"]
                if c["id"] == case["normalization"] + "__as_plotted"
            )
            assert a["parameters"] == old["fits"][0]["parameters"]
        for fit in case["fits"]:
            covariance = np.array(fit["conditional_unscaled_covariance"])
            assert covariance.shape == (len(case["free_parameters"]),) * 2
            assert np.linalg.eigvalsh(covariance).min() > 0
            assert fit["degrees_of_freedom"] == case["rows"] - len(
                case["free_parameters"]
            )
        if case["action"] == "published_kprime":
            assert a["parameters"]["K0_prime"] == 5.06534
            assert "K0_prime" not in case["free_parameters"]
            assert len(case["free_parameters"]) == 5


def test_constrained_fit_and_remapped_scale_deletion_replay():
    for action in ["published_kprime", "omit_kono_ks", "omit_conrad_ks"]:
        case = next(c for c in report()["cases"] if c["id"] == "plotted__" + action)
        rows, starts, config = prepare(case)
        fit = (
            solve_fixed_kprime(rows, starts[1], config)
            if action == "published_kprime"
            else solve(rows, starts[1], **config)
        )
        np.testing.assert_allclose(
            list(fit["parameters"].values()),
            list(case["fits"][0]["parameters"].values()),
            rtol=1e-6,
        )
        if action == "published_kprime":
            # Independent finite difference of the frozen-error objective:
            # every free coefficient must be stationary at the constrained fit.
            c = np.array(list(fit["parameters"].values()))
            d = arrays(rows)
            sigma = effective_sigma(c, d)
            residual = (pressure(c, d) - d["P"]) / sigma
            for i in [0, 2, 3, 4, 5]:
                h = config["scaling"][i] * 2e-6
                plus, minus = c.copy(), c.copy()
                plus[i] += h
                minus[i] -= h
                derivative = (pressure(plus, d) - pressure(minus, d)) / (2 * h * sigma)
                relative_score = abs(derivative @ residual) / (
                    np.linalg.norm(derivative) * np.linalg.norm(residual)
                )
                assert relative_score < 2e-5


def test_profile_partitions_match_exports_and_local_finite_difference():
    r = report()
    by_id = {c["id"]: c for c in r["cases"]}
    with (DATA / "dataset-influence-profile-rows.csv").open() as f:
        exported = list(csv.DictReader(f))
    for normalization in plan()["normalizations"]:
        diagnostic, rows, groups = profile_diagnostics(normalization, by_id)
        assert len(rows) == 219 and len(groups) == 12
        np.testing.assert_allclose(sum(g["information_fraction"] for g in groups), 1)
        np.testing.assert_allclose(
            sum(g["local_upward_kprime_contribution"] for g in groups),
            diagnostic["local_predicted_kprime_shift"],
        )
        np.testing.assert_allclose(
            diagnostic["local_predicted_kprime_shift"],
            diagnostic["actual_kprime_shift"],
            atol=1e-4,
        )
        stored = r["profile"][normalization]
        np.testing.assert_allclose(
            diagnostic["efficient_information"], stored["efficient_information"]
        )
        selected = [x for x in exported if x["normalization"] == normalization]
        np.testing.assert_allclose(
            [float(x["information_fraction"]) for x in selected],
            [x["information_fraction"] for x in rows],
        )
        # Row exports partition the fixed-weight curve costs exactly.
        for g in groups:
            subset = [x for x in selected if x["group"] == g["group"]]
            for curve in ["published_all", "fixed_kprime", "free"]:
                cost = sum(
                    float(x[f"{curve}_standardized_residual"]) ** 2 for x in subset
                )
                np.testing.assert_allclose(cost, g[f"{curve}_sum_squares"])
    with (DATA / "dataset-influence-residuals.csv").open() as f:
        residuals = list(csv.DictReader(f))
    assert len(residuals) == sum(c["rows"] for c in r["cases"])
    for case in r["cases"]:
        selected = [x for x in residuals if x["case"] == case["id"]]
        cost = sum(
            (float(x["pressure_residual_gpa"]) / float(x["effective_error_gpa"])) ** 2
            for x in selected
        )
        np.testing.assert_allclose(cost, case["fits"][0]["objective_sum_squares"])


def test_conrad_derived_moduli_and_forward_residuals_and_preservation():
    r = report()
    rows = conrad_native_audit({c["id"]: c for c in r["cases"]})
    assert len(rows) == 10 and sum(x["selection"] == "selected" for x in rows) == 6
    np.testing.assert_allclose(
        [x["derived_ks_gpa"] for x in rows[:5]],
        [174.5366666667, 187.0933333333, 205.3933333333, 214.27, 225.23],
    )
    for x in rows:
        for curve in ["published", "baseline", "omit_conrad"]:
            np.testing.assert_allclose(
                x[f"observed_minus_{curve}_ks_gpa"],
                x["derived_ks_gpa"] - x[f"{curve}_ks_gpa"],
            )
    # An independent forward model also inverts each selected observed KS at
    # the pressure predicted by the fast fit equation.
    from scripts.audit_garnet_angel_2022 import source
    from scripts.audit_grossular_angel_2022_observations import GAS_CONSTANT
    from scripts.audit_pyrope_angel_2022_observations import properties

    baseline = r["cases"][0]
    retained, _, _ = prepare(baseline)
    selected = [x for x in retained if x["dataset"] == "conrad_ks"]
    c = np.array(list(baseline["fits"][0]["parameters"].values()))
    p, _ = source("grossular")
    p.update(dict(zip([2, 3, 11, 18], c[:4])))
    predicted_pressure = pressure(c, arrays(selected))
    for row, pp in zip(selected, predicted_pressure):
        np.testing.assert_allclose(
            properties(p, pp, 298, GAS_CONSTANT)["ks_gpa"], row["Y"], atol=1e-6
        )
    for filename in [
        "artifact-manifest.json",
        "followup-artifact-manifest.json",
        "expanded-artifact-manifest.json",
        "normalized-artifact-manifest.json",
        "full-sensitivity-artifact-manifest.json",
        "normalization-specimen-artifact-manifest.json",
    ]:
        for entry in json.loads((DATA / filename).read_text())["files"]:
            assert (
                hashlib.sha256((DATA / entry["file"]).read_bytes()).hexdigest()
                == entry["sha256"]
            )
