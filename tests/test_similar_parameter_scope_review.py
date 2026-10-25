"""Numerical scope, missing-error semantics and complete source preservation."""

import copy
import csv
import json
from collections import Counter

import numpy as np
import pytest

from peritheos import convert_eosmat, load_eosmat, validate_eosmat_document
from scripts.review_similar_parameter_scope import (
    INPUT,
    OUTPUT,
    ROOT,
    TANGE_ERRORS,
    augment_outcome,
    calculate,
    check_report,
    review,
    saved_review,
    verify_sources,
)


def test_all_96_decisions_are_reproducible_and_sources_preserved():
    expected = json.loads(OUTPUT.read_text())
    check_report(expected, calculate())
    assert len(expected["records"]) == 96
    assert Counter(r["status"] for r in expected["records"]) == {
        "parity": 24,
        "similar": 72,
    }
    assert Counter(
        r["scope"] for r in expected["records"] if r["status"] == "parity"
    ) == {"parameter_agreement": 21, "independent_refit": 2, "source_equation": 1}
    verify_sources()


def test_all_promotions_meet_their_named_criterion_without_treating_missing_as_zero():
    for case in saved_review().values():
        if case["status"] != "parity":
            continue
        free = [
            p
            for p in case["parameters"]
            if p["parameter"] in case["free_comparison_parameters"]
        ]
        assert free
        for p in free:
            delta = abs(p["refit"] - p["published"])
            if case["parity_basis"] == "within_reported_parameter_errors":
                assert p["published_error"] > 0
                assert delta <= p["published_error"]
            elif (
                case["parity_basis"] == "conditional_combined_error_parameter_agreement"
            ):
                assert p["published_error"] > 0 and p["refit_error"] > 0
                assert delta <= 2 * np.hypot(p["published_error"], p["refit_error"])
            elif (
                case["parity_basis"]
                == "conditional_reported_widths_and_print_precision"
            ):
                width = p["published_error"]
                if width is None:
                    width = case["printed_half_units"][p["parameter"]]
                assert width > 0 and delta <= width
            else:
                assert case["scope"] == "independent_refit"
                assert delta <= max(1e-8, 1e-6 * abs(p["published"]))


def test_native_assessments_roundtrip_and_unpromoted_ledger_findings_are_unchanged():
    original = {
        r["record_identifier"]: r for r in json.loads(INPUT.read_text())["records"]
    }
    ledger = {
        r["record_identifier"]: r
        for r in json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())[
            "records"
        ]
    }
    for material in {r["material_identifier"] for r in saved_review().values()}:
        native = load_eosmat(ROOT / "peritheos/data/materials" / (material + ".eosmat"))
        restored = convert_eosmat(convert_eosmat(native, version=3), version=4)
        validate_eosmat_document(native)
        validate_eosmat_document(restored)
        for record, roundtrip in zip(native["eos_records"], restored["eos_records"]):
            rid = record["identifier"]
            if rid not in saved_review():
                continue
            case = saved_review()[rid]
            if case["status"] == "similar":
                assert ledger[rid] == original[rid]
                continue
            assessment = record["evidence"]["reproduction"]["assessment"]
            assert assessment == roundtrip["evidence"]["reproduction"]["assessment"]
            assert assessment["status"] == "parity"
            assert assessment["scope"] == case["scope"]
            assert assessment["basis_code"] == case["parity_basis"]
            assert "Earlier audit wording" in ledger[rid]["qualification"]
            before = {p["parameter"]: p for p in original[rid].get("parameters", [])}
            for p in assessment["parameters"]:
                if p["parameter"] in before:
                    assert p["published"] == before[p["parameter"]]["published"]
                    assert p["refit"] == before[p["parameter"]]["refit"]


@pytest.mark.parametrize("compound,n", [("rbbr", 123), ("rbi", 121)])
def test_farla_covariance_uses_each_compounds_actual_inputs_and_fixed_weights(
    compound, n
):
    case = saved_review()[compound + "_b2_farla_2025_bm3_mgd_cscl"]
    covariance = case["covariance"]
    assert covariance["residual_count"] == n
    assert covariance["degrees_of_freedom"] == n - 4
    assert covariance["max_abs_residual_jacobian_cosine"] < 0.01
    assert covariance["derivative_step_halving_relative_change"] < 2e-4
    matrix = np.array(covariance["matrix"])
    np.testing.assert_allclose(
        np.diag(matrix), [p["refit_error"] ** 2 for p in case["parameters"]], rtol=1e-12
    )
    assert np.linalg.eigvalsh(covariance["correlation_matrix"]).min() > 0
    assert "fixed V0/theta0" in " ".join(case["qualifications"])


def test_tange_recovers_only_primary_table_2_errors_with_unchanged_estimates():
    original = {
        r["record_identifier"]: r for r in json.loads(INPUT.read_text())["records"]
    }
    for rid, errors in TANGE_ERRORS.items():
        case = saved_review()[rid]
        before = {p["parameter"]: p for p in original[rid]["parameters"]}
        for p in case["parameters"]:
            if p["parameter"] in errors:
                assert before[p["parameter"]]["published_error"] is None
                assert p["published_error"] == errors[p["parameter"]]
            assert p["refit_error"] == before[p["parameter"]]["refit_error"]
        assert "standard deviations" in " ".join(case["qualifications"])


def test_marquardt_printed_precision_does_not_become_a_source_standard_error():
    with (
        ROOT
        / "peritheos/data/datasets/mg090fe010o-marquardt-2009-table-s1-ls-models.csv"
    ).open() as stream:
        tokens = {
            int(r["source_row"]): r["k0_gpa"]
            for r in csv.DictReader(stream)
            if r["source_row"].isdigit()
        }
    matches = []
    for index in range(1, 13):
        case = saved_review()[f"mg090fe010o_marquardt_2009b_ls_bm3_s1_{index:02d}"]
        p = case["parameters"][0]
        assert tokens[index].isdigit()
        assert p["published_error"] is None
        assert case["printed_half_units"] == {"K0": 0.5}
        matched = abs(p["refit"] - int(tokens[index])) <= 0.5
        assert (case["status"] == "parity") == matched
        if matched:
            matches.append(index)
    assert matches == [1, 2, 3, 6, 8, 9, 10, 11, 12]


def test_solomatova_priors_and_fixed_constraints_are_not_independent_recovery():
    cases = [
        r
        for r in saved_review().values()
        if "_solomatova_2016_" in r["record_identifier"] and r["status"] == "parity"
    ]
    assert len(cases) == 5
    for case in cases:
        if "_ls_" in case["record_identifier"]:
            assert "K0_prime" in case["fixed_comparison_parameters"]
            assert "K0_prime" not in case["free_comparison_parameters"]
        assert all(p["refit_error"] is None for p in case["parameters"])
        assert "does not independently validate" in " ".join(case["qualifications"])
        assert case["parity_basis"] == "within_reported_parameter_errors"


def test_missing_error_and_covariance_convention_holds_are_not_waived():
    cases = saved_review()
    for rid in ["palladium_frost_2023_vinet_1", "palladium_frost_2023_bm3_2"]:
        assert cases[rid]["status"] == "similar"
        assert "no favorable rescaling" in " ".join(cases[rid]["qualifications"])
    assert cases["iceviii_besson_1994_bm3_1"]["status"] == "similar"
    original = next(
        r
        for r in json.loads(INPUT.read_text())["records"]
        if r["record_identifier"] == "naalsi2o6_zhao_1997_bm3_1"
    )
    altered = copy.deepcopy(original)
    altered["parameters"][0]["refit"] *= 1.5
    assert review(altered)["status"] == "similar"
    # Row 10 is near the half-last-digit boundary: numerical closeness to a
    # saved fit must not hide a newly failed printed-precision comparison.
    row = next(
        r
        for r in json.loads(INPUT.read_text())["records"]
        if r["record_identifier"].endswith("marquardt_2009b_ls_bm3_s1_10")
    )
    changed = copy.deepcopy(row)
    changed["parameters"][0]["refit"] = 176.4999
    with pytest.raises(ValueError, match="no longer meets its scoped criterion"):
        augment_outcome({"identifier": row["record_identifier"]}, changed)


def test_independent_refit_scope_replays_actual_solvers():
    from scripts.refit_chen_2010_argon import fit, observations
    from scripts.reproduce_aip_handbook_cu_vinet_refit import refit

    argon = fit(observations())["parameters"]
    case = saved_review()["argon_fcc_chen_2010_bm3_digitized_refit"]
    for p in case["parameters"]:
        assert argon[p["parameter"]] == pytest.approx(p["published"], rel=1e-6)
    cu = refit()["refit"]
    case = saved_review()["cu_sun_2010_low_vn_refit"]
    assert cu["k0_gpa"] == pytest.approx(case["parameters"][0]["published"], rel=1e-6)
    assert cu["k0_prime"] == pytest.approx(case["parameters"][1]["published"], rel=1e-6)
