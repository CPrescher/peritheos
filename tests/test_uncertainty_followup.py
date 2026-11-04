"""Independent uncertainty cases and scientific-preservation checks for the audit."""

import copy
import hashlib
import json
from collections import Counter

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import convert_eosmat, load_eosmat, validate_eosmat_document
from scripts.reconstruct_similar_eos_uncertainties import (
    INPUT,
    OUTPUT,
    ROOT,
    augment_outcome,
    covariance_from_residual,
    reconstruct,
    saved_followup,
    verify_preserved_sources,
)


def test_covariance_agrees_with_closed_weighted_linear_regression_and_units():
    t = np.linspace(-2, 3, 21)
    design = np.column_stack([np.ones(len(t)), 1e6 * t])
    y = design @ [2, 3e-6] + 0.1 * np.sin(3 * t)
    sigma = 1 + 0.1 * t
    xw, yw = design / sigma[:, None], y / sigma
    coefficients = np.linalg.solve(xw.T @ xw, xw.T @ yw)

    def residual(x):
        return (design @ x - y) / sigma

    expected = (
        np.linalg.inv(xw.T @ xw)
        * (residual(coefficients) @ residual(coefficients))
        / 19
    )
    result = covariance_from_residual(residual, coefficients, ["intercept", "slope"])
    np.testing.assert_allclose(result["matrix"], expected, rtol=1e-7, atol=1e-25)
    assert result["degrees_of_freedom"] == 19
    # Rescaling a coordinate's units scales both covariance axes.
    transform = np.diag([1, 1e6])
    rescaled = covariance_from_residual(
        lambda z: residual(z / [1, 1e6]),
        coefficients * [1, 1e6],
        ["intercept", "slope"],
    )
    np.testing.assert_allclose(
        rescaled["matrix"], transform @ expected @ transform, rtol=1e-7
    )
    with pytest.raises(ValueError, match="not identifiable"):
        covariance_from_residual(
            lambda x: np.arange(10) - (x[0] + x[1]), [1.0, 2.0], ["a", "b"]
        )


def test_all_cases_have_errors_but_holds_and_source_values_survive():
    inputs = json.loads(INPUT.read_text())
    result = json.loads(OUTPUT.read_text())
    assert len(inputs["records"]) == len(result["records"]) == 63
    assert Counter(r["status"] for r in result["records"]) == {
        "parity": 63,
    }
    original = {r["record_identifier"]: r for r in inputs["records"]}
    assert (
        sum(
            not r.get("parameters")
            or any(p.get("refit_error") is None for p in r["parameters"])
            for r in inputs["records"]
        )
        == 51
    )
    for row in result["records"]:
        old = original[row["record_identifier"]]
        for p in row["parameters"]:
            assert p["refit_error"] > 0 and np.isfinite(p["refit_error"])
        for before, after in zip(old.get("parameters", []), row["parameters"]):
            for key in ("parameter", "published", "published_error", "refit"):
                assert before[key] == after[key]
            if before.get("refit_error") is not None:
                assert before["refit_error"] == after["refit_error"]
            assert before.get("within_combined_2sigma") == after.get(
                "within_combined_2sigma"
            )
        if row["status"] == "parity" and row["record_identifier"] not in {
            "ca_perovskite_tetragonal_chen_2018_vinet",
            "superhydrous_phase_b_lt_litasov_2007_htbm_t03",
        }:
            assert all(
                abs(p["refit"] - p["published"]) <= p["published_error"]
                for p in row["parameters"]
            )
        uncertainty = row["covariance"]
        if "matrix" in uncertainty:
            c = np.array(uncertainty["matrix"])
            errors = np.array([p["refit_error"] for p in row["parameters"]])
            np.testing.assert_allclose(np.diag(c), errors**2, rtol=1e-12)
            correlation = c / np.outer(errors, errors)
            assert np.linalg.eigvalsh(correlation).min() > 0
            np.testing.assert_allclose(c, c.T, rtol=1e-12, atol=1e-30)
        assert uncertainty["previous_assessment"]["status"] == "similar"
    verify_preserved_sources()
    for path, digest in result["source_sha256"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest


@pytest.mark.parametrize(
    "rid,source_width_divisor",
    [
        ("ca_perovskite_tetragonal_chen_2018_vinet", 2),
        ("superhydrous_phase_b_lt_litasov_2007_htbm_t03", 1),
    ],
)
def test_combined_error_parity_retains_source_confidence_and_scope(
    rid, source_width_divisor
):
    row = saved_followup()[rid]
    assert row["status"] == "parity"
    assert any(abs(p["difference"]) > p["published_error"] for p in row["parameters"])
    for p in row["parameters"]:
        combined = np.sqrt(
            (p["published_error"] / source_width_divisor) ** 2 + p["refit_error"] ** 2
        )
        assert abs(p["difference"]) <= 2 * combined
        assert p["similar"] is True
    assert "author workflow remain unrecovered" in row["decision"]
    assert "shared-data cross-fit covariance is unknown" in row["decision"]


def test_farla_parameter_parity_preserves_separate_calibration_and_benchmark_findings():
    rid = "rbcl_b2_farla_2025_bm3_mgd_cscl"
    row = saved_followup()[rid]
    assert row["status"] == "parity"
    assert {p["parameter"] for p in row["parameters"]} == {
        "K0",
        "K0_prime",
        "gamma0",
        "q",
    }
    assert (
        max(
            abs(p["refit"] - p["published"]) / p["published_error"]
            for p in row["parameters"]
        )
        < 0.40
    )
    native = load_eosmat(ROOT / "peritheos/data/materials/rbcl.eosmat")
    record = next(r for r in native["eos_records"] if r["identifier"] == rid)
    assessment = record["evidence"]["reproduction"]["assessment"]
    assert assessment["status"] == "parity"
    assert assessment["scope"] == "parameter_agreement"
    assert assessment["basis_code"] == "within_reported_parameter_errors"
    assert (
        "CsCl calibration replay and source benchmark discrepancies remain unresolved as separate findings"
        in row["decision"]
    )
    validation = record["evidence"]["scientific_validation"]
    check = validation["primary_data_check"]["reconstruction_assessment"]
    assert check["author_fit_parity"] is False
    assert check["author_workflow_parity"] is False
    assert "53/54" in check["published_curve_consistency"]
    verify_preserved_sources()


@pytest.mark.parametrize(
    "rid,material,max_source_width,required_qualification",
    [
        (
            "diamond_datchi_2007_vinet_1",
            "diamond",
            0.45,
            "source-unconfirmed run-2/equal-volume hypothesis",
        ),
        (
            "platinum_fei_2004_bm3_300k",
            "platinum",
            0.67,
            "source's fixed/free mask is unrecovered",
        ),
        (
            "siderite_fe095mn005_litasov_2013_bm3",
            "siderite_fe095mn005",
            0.23,
            "source-data agreement finding is unresolved separately",
        ),
    ],
)
def test_remaining_parameter_matches_retain_their_distinct_limits(
    rid, material, max_source_width, required_qualification
):
    row = saved_followup()[rid]
    assert row["status"] == "parity"
    assert (
        max(
            abs(p["refit"] - p["published"]) / p["published_error"]
            for p in row["parameters"]
        )
        < max_source_width
    )
    assert required_qualification in row["decision"]
    native = load_eosmat(ROOT / "peritheos/data/materials" / (material + ".eosmat"))
    record = next(r for r in native["eos_records"] if r["identifier"] == rid)
    assessment = record["evidence"]["reproduction"]["assessment"]
    assert assessment["status"] == "parity"
    assert assessment["scope"] == "parameter_agreement"
    assert assessment["basis_code"] == "within_reported_parameter_errors"
    assert required_qualification in " ".join(assessment["qualifications"])
    if material == "diamond":
        assert "14/24" in assessment["fit"]["selection"]
        assert (
            record["evidence"]["scientific_validation"]["primary_data_check"]["status"]
            == "plot_only"
        )
    elif material == "platinum":
        assert assessment["fit"]["validation_fixed_parameters"] == ["V0"]
        assert assessment["fit"]["fixed_values"] == {"V0": 60.38}
        assert assessment["fit"]["free_parameters"] == ["K0", "K0_prime"]
        assert "not the preferred Pt pressure scale" in row["decision"]
    else:
        assert assessment["fit"]["observations"] == 27
        assert "0.22958 GPa" in row["decision"]
        calibration = record["evidence"]["scientific_validation"][
            "primary_source_check"
        ]["finding"]
        assert "V0 fixed at 293.4" in calibration
    verify_preserved_sources()


def test_komabayashi_propagation_matches_independent_small_noise_refits():
    from scripts.reproduce_wang_zhou_komabayashi_mgsio3 import observations

    rid = "mgsio3_post_perovskite_komabayashi_2008_thermal_bm3"
    case = saved_followup()[rid]
    rows = np.array(observations()[rid])
    v, t = rows[:, 1], rows[:, 3]
    center = np.array([p["refit"] for p in case["parameters"]])

    # An independent explicit BM2+thermal expression, not the native evaluator.
    def pressure(k0, beta):
        alpha1, dkdt = beta * [1e-9, 0.01]
        alpha0 = 1.7e-5 - 300 * alpha1
        vt = 163.813 * np.exp(alpha0 * (t - 300) + 0.5 * alpha1 * (t * t - 300**2))
        ratio = vt / v
        return 1.5 * (k0 + dkdt * (t - 300)) * (ratio ** (7 / 3) - ratio ** (5 / 3))

    cold = t == 300
    shape = 1.5 * ((163.813 / v[cold]) ** (7 / 3) - (163.813 / v[cold]) ** (5 / 3))
    beta0 = center[1:] / [1e-9, 0.01]
    truth = pressure(center[0], beta0)
    rng = np.random.default_rng(251025)
    sigma = 0.01
    draws = []
    for _ in range(600):
        perturbed = truth + rng.normal(0, sigma, len(t))
        fitted_k = shape @ perturbed[cold] / (shape @ shape)
        fit = least_squares(
            lambda b: pressure(fitted_k, b) - perturbed,
            beta0,
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
        )
        assert fit.success
        draws.append([fitted_k, *(fit.x * [1e-9, 0.01])])
    empirical = np.cov(draws, rowvar=False) / sigma**2
    propagated = (
        np.array(case["covariance"]["matrix"]) / case["covariance"]["residual_variance"]
    )
    np.testing.assert_allclose(empirical, propagated, rtol=0.15, atol=1e-23)
    assert 18 < case["covariance"]["degrees_of_freedom"] < 19
    assert case["covariance"]["stage_uncertainty_propagated"] is True
    assert case["covariance"]["correlation_matrix"][0][1] > 0.7


def test_ma_existing_covariance_converts_both_axes():
    source = json.loads((ROOT / "docs/data/ma-2024-kcl-reproduction.json").read_text())[
        "fits"
    ]["reported_errors"]
    case = saved_followup()["kcl_b2_ma_2024_bm3_mgd"]
    d = np.diag([1 / 0.602214076, 1, 1])
    expected = d @ source["diagnostic_covariance"] @ d
    np.testing.assert_allclose(case["covariance"]["matrix"], expected, rtol=1e-14)
    np.testing.assert_allclose(
        [p["refit_error"] for p in case["parameters"]], np.sqrt(np.diag(expected))
    )


def test_every_assessment_preserves_native_source_and_compatibility_exports():
    verify_preserved_sources()
    cases = saved_followup()
    ledger = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    ledger_rows = {r["record_identifier"]: r for r in ledger["records"]}
    for material in {r["material_identifier"] for r in cases.values()}:
        native = load_eosmat(ROOT / "peritheos/data/materials" / (material + ".eosmat"))
        restored = convert_eosmat(convert_eosmat(native, version=3), version=4)
        validate_eosmat_document(native)
        validate_eosmat_document(restored)
        for record, roundtrip in zip(native["eos_records"], restored["eos_records"]):
            if record["identifier"] not in cases:
                continue
            assessment = record["evidence"]["reproduction"]["assessment"]
            assert assessment == roundtrip["evidence"]["reproduction"]["assessment"]
            assert assessment["status"] == cases[record["identifier"]]["status"]
            assert assessment["scope"] == "parameter_agreement"
            assert (
                assessment["fit"]["covariance"]
                == cases[record["identifier"]]["covariance"]
            )
            row = ledger_rows[record["identifier"]]
            assert row["parameters"] == assessment["parameters"]
            assert "Earlier audit statement" in row["qualification"]
            invalid = copy.deepcopy(row)
            invalid["parameters"][0]["published"] *= 1.01
            with pytest.raises(ValueError, match="Changed source"):
                augment_outcome(record, invalid)
            invalid = copy.deepcopy(row)
            invalid["status"] = "refit_failed"
            with pytest.raises(ValueError, match="Original reconstruction failed"):
                augment_outcome(record, invalid)


def test_complete_numerical_reconstruction_is_repeatable():
    inputs = json.loads(INPUT.read_text())["records"]
    replay = [reconstruct(row) for row in inputs]
    stored = saved_followup()
    for case in replay:
        original = stored[case["record_identifier"]]
        assert case["status"] == original["status"]
        np.testing.assert_allclose(
            [p["refit_error"] for p in case["parameters"]],
            [p["refit_error"] for p in original["parameters"]],
            rtol=3e-4,
        )
