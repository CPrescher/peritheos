"""Conrad uncertainty controls retain data and distinguish shared errors."""

import csv
import hashlib
import json

import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA
from scripts.refit_grossular_angel_2022_conrad_processing import (
    PLAN,
    cases,
    density_controls,
    prepare,
    residual_receipt,
    solve_shared,
    weight_state,
    whiten,
)
from scripts.refit_grossular_angel_2022_dataset_influence import derivatives
from scripts.refit_grossular_angel_2022_normalization_specimen import prepared_case
from scripts.refit_grossular_angel_2022_subset import arrays, pressure, solve


def report():
    return json.loads((DATA / "conrad-processing-refits.json").read_text())


def test_all_cases_retain_three_conrad_rows_and_unchanged_other_data():
    assert len(cases()) == 26
    for case in cases():
        rows, starts, config = prepare(case)
        original, _, _ = prepared_case(
            dict(normalization=case["normalization"], isaak="as_plotted")
        )
        assert len(rows) == len(original) == 219
        assert len(starts) == 2
        assert [r for r in rows if r["dataset"] != "conrad_ks"] == [
            r for r in original if r["dataset"] != "conrad_ks"
        ]
        conrad = [r for r in rows if r["dataset"] == "conrad_ks"]
        assert [r["row"] for r in conrad] == [1, 2, 3]
        assert [r["P"] for r in conrad] == [1.6, 3.5, 6.95]
        assert all(r["T"] == 298 and r["KS"] and r["sT"] is None for r in conrad)
        for r in conrad:
            old = next(
                x
                for x in original
                if x["dataset"] == "conrad_ks" and x["row"] == r["row"]
            )
            if (
                not case["variant"].startswith("density_")
                and case["variant"] != "round_01"
            ):
                assert r["Y"] == old["Y"]
            if case["variant"].startswith("extra_"):
                fraction = int(case["variant"].split("_")[1]) / 1000
                np.testing.assert_allclose(
                    r["sY"] ** 2, old["sY"] ** 2 + (fraction * r["Y"]) ** 2
                )
            if case["variant"] == "free_conrad_scale":
                assert r["scale"] == 6 and config["names"][6] == "scale_conrad_ks"
        assert len(config["names"]) == (
            7 if case["variant"] == "free_conrad_scale" else 6
        )


def test_density_reduction_and_formal_covariance_bounds_are_explicit():
    density = density_controls()
    assert len(density) == 6
    for r in density:
        # The computed density change is inside printed-velocity rounding
        # bounds; it does not resolve a mismatch at the printed precision.
        assert (
            r["rounding_density_lower"]
            <= r["published_eos_calculated_density"]
            <= r["rounding_density_upper"]
        )
        np.testing.assert_allclose(
            r["calculated_modulus_multiplier"],
            r["published_eos_calculated_density"] / r["inferred_density_g_cm3"],
        )
    cases_by_variant = {
        c["variant"]: c for c in cases() if c["normalization"] == "plotted"
    }
    widths = {}
    for variant in ["baseline", "cij_negative_correlation", "cij_independent"]:
        rows, _, _ = prepare(cases_by_variant[variant])
        widths[variant] = np.array(
            [r["sY"] for r in rows if r["dataset"] == "conrad_ks"]
        )
    np.testing.assert_allclose(
        widths["cij_negative_correlation"], [0.01, 0.17 / 3, 0.02]
    )
    assert np.all(widths["cij_negative_correlation"] <= widths["cij_independent"])
    assert np.all(widths["cij_independent"] <= widths["baseline"])
    for variant in ["density_vp", "density_vs"]:
        rows, _, _ = prepare(cases_by_variant[variant])
        baseline, _, _ = prepare(cases_by_variant["baseline"])
        for row, old in zip(rows, baseline):
            if row["dataset"] == "conrad_ks":
                np.testing.assert_allclose(row["Y"] / old["Y"], row["sY"] / old["sY"])


def test_shared_whitening_matches_dense_gls_and_rank_one_identity():
    case = next(c for c in report()["cases"] if c["id"] == "plotted__shared_020")
    rows, _, config = prepare(case)
    data = arrays(rows) | {"dataset": [r["dataset"] for r in rows]}
    c = np.array(list(case["fits"][0]["parameters"].values()))
    state = weight_state(c, data, shared=True)
    sigma, idx, block, factor = state
    residual = pressure(c, data) - data["P"]
    z = whiten(residual, state)
    np.testing.assert_allclose(np.einsum("ik,jk->ij", factor, factor), block)
    dense = np.einsum("i,i->", residual[idx], np.linalg.solve(block, residual[idx]))
    np.testing.assert_allclose(np.sum(z[idx] ** 2), dense)
    u = np.sqrt(np.diag(block) - sigma[idx] ** 2)
    # u has positive entries because dP/dKS is positive on these states.
    rank_one = np.sum((residual[idx] / sigma[idx]) ** 2) - np.sum(
        u * residual[idx] / sigma[idx] ** 2
    ) ** 2 / (1 + np.sum((u / sigma[idx]) ** 2))
    np.testing.assert_allclose(dense, rank_one)
    jac = derivatives(c, data, np.ones(len(rows)), config["scaling"])
    weighted_jac = whiten(jac, state)
    np.testing.assert_allclose(weighted_jac[idx], np.linalg.solve(factor, jac[idx]))
    scores = np.einsum("ij,i->j", weighted_jac, z)
    relative = np.abs(scores) / (
        np.linalg.norm(weighted_jac, axis=0) * np.linalg.norm(z)
    )
    assert np.max(relative) < 2e-5
    scaled_jac = weighted_jac * config["scaling"]
    information = np.einsum("ij,ik->jk", scaled_jac, scaled_jac)
    covariance = np.array(case["fits"][0]["conditional_unscaled_covariance"])
    scaled_covariance = covariance / np.outer(config["scaling"], config["scaling"])
    np.testing.assert_allclose(
        np.einsum("ij,jk->ik", information, scaled_covariance), np.eye(6), atol=1e-6
    )


def test_fits_replay_and_all_conditional_covariances_match_free_dimensions():
    r = report()
    assert r["plan_sha256"] == hashlib.sha256(PLAN.read_bytes()).hexdigest()
    assert r["all_starts_converged"] and not r["any_start_at_bound"]
    assert r["all_three_conrad_rows_retained"] and not r["author_workflow_reproduced"]
    old = json.loads((DATA / "normalization-specimen-refits.json").read_text())
    for case in r["cases"]:
        a, b = case["fits"]
        np.testing.assert_allclose(
            list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-6
        )
        for fit in case["fits"]:
            covariance = np.array(fit["conditional_unscaled_covariance"])
            count = len(case["free_parameters"])
            assert covariance.shape == (count, count)
            assert np.linalg.eigvalsh(covariance).min() > 0
            assert fit["degrees_of_freedom"] == 219 - count
        if case["variant"] == "baseline":
            prior = next(
                c
                for c in old["cases"]
                if c["id"] == case["normalization"] + "__as_plotted"
            )
            assert a["parameters"] == prior["fits"][0]["parameters"]
    for variant in ["extra_020", "shared_020", "free_conrad_scale"]:
        case = next(c for c in r["cases"] if c["id"] == "plotted__" + variant)
        rows, starts, config = prepare(case)
        fit = (
            solve_shared(rows, starts[1], config)
            if variant == "shared_020"
            else solve(rows, starts[1], **config)
        )
        np.testing.assert_allclose(
            list(fit["parameters"].values()),
            list(case["fits"][0]["parameters"].values()),
            rtol=1e-6,
        )


def test_exported_whitened_costs_inputs_and_all_previous_artifacts():
    r = report()
    with (DATA / "conrad-processing-residuals.csv").open() as f:
        exported = list(csv.DictReader(f))
    assert len(exported) == 26 * 219
    with (DATA / "conrad-processing-inputs.csv").open() as f:
        inputs = list(csv.DictReader(f))
    assert len(inputs) == 26 * 3
    for case in r["cases"]:
        selected = [x for x in exported if x["case"] == case["id"]]
        objective = sum(float(x["whitened_pressure_residual"]) ** 2 for x in selected)
        np.testing.assert_allclose(objective, case["fits"][0]["objective_sum_squares"])
        rows, _, _ = prepare(case)
        computed, block = residual_receipt(
            rows, case["fits"][0], case["variant"] == "shared_020"
        )
        np.testing.assert_allclose(
            [float(x["whitened_pressure_residual"]) for x in selected],
            [x["whitened_pressure_residual"] for x in computed],
        )
        np.testing.assert_allclose(
            block["covariance_pressure_gpa2"],
            case["conrad_error_block"]["covariance_pressure_gpa2"],
        )
        if case["variant"] != "shared_020":
            np.testing.assert_allclose(
                [float(x["whitened_pressure_residual"]) for x in selected],
                [
                    float(x["pressure_residual_gpa"])
                    / float(x["marginal_effective_error_gpa"])
                    for x in selected
                ],
            )
    for name in [
        "artifact-manifest.json",
        "followup-artifact-manifest.json",
        "expanded-artifact-manifest.json",
        "normalized-artifact-manifest.json",
        "full-sensitivity-artifact-manifest.json",
        "normalization-specimen-artifact-manifest.json",
        "dataset-influence-artifact-manifest.json",
    ]:
        for entry in json.loads((DATA / name).read_text())["files"]:
            assert (
                hashlib.sha256((DATA / entry["file"]).read_bytes()).hexdigest()
                == entry["sha256"]
            )
