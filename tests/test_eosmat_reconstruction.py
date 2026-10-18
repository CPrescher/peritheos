"""Named reconstruction outcomes retain their actual values and scientific scope."""

import copy
import hashlib
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    convert_eosmat,
    create_eosmat,
    eosmat_schema,
    get_material_document,
    load_eosmat,
    validate_eosmat_document,
)
from peritheos._eosmat_format import format_eosmat
from peritheos.errors import EosmatError
from scripts.generate_material_summaries import reconstruction
from scripts.run_eosmat_audit import historical_reconstruction_view
from scripts.standardize_reconstruction_metadata import (
    REPORT,
    STATUS,
    assessment_for_record,
    standardize_document,
)

ROOT = Path(__file__).resolve().parents[1]


def sample():
    document = create_eosmat(
        "example",
        "X",
        "BM3",
        {"V0": 10, "K0": 160, "K0_prime": 4},
        temperature_ref=300,
    )
    record = document["eos_records"][0]
    row = {
        "record_identifier": record["identifier"],
        "status": "parity_not_achieved",
        "reason": "Synthetic test: the free K0 comparison fails.",
        "parameters": [
            {
                "parameter": "K0",
                "published": 160,
                "published_error": None,
                "refit": 180.12345678901235,
                "refit_error": 2.5,
                "within_combined_2sigma": None,
                "similar": False,
            }
        ],
        "fixed_parameters": ["V0", "K0_prime"],
        "objective": "equal_weight_pressure_residuals",
        "qualification": "Synthetic conditional diagnostic, not author workflow.",
    }
    return document, row


def test_migration_is_additive_idempotent_and_does_not_promote_source_validation():
    document, row = sample()
    record = document["eos_records"][0]
    record["evidence"] = {
        "scientific_validation": {"status": "primary_source_validated"},
        "reproduction": {"historical_finding": [None, -0.0, "source finding"]},
    }
    before = copy.deepcopy(document)
    result = standardize_document(document, {row["record_identifier"]: row})
    assert document == before
    assessment = result["eos_records"][0]["evidence"]["reproduction"].pop("assessment")
    assert result == before
    assert assessment["status"] == "not_reproduced"
    assert assessment["scope"] == "parameter_agreement"
    assert assessment["parameters"] == row["parameters"]
    assert assessment["fit"]["fixed_parameters"] == ["V0", "K0_prime"]
    result["eos_records"][0]["evidence"]["reproduction"]["assessment"] = assessment
    assert standardize_document(result, {row["record_identifier"]: row}) == result
    with pytest.raises(ValueError, match="No existing reconstruction ledger"):
        standardize_document(before, {})
    result["eos_records"][0]["evidence"]["reproduction"]["assessment"]["basis"] = (
        "User finding"
    )
    with pytest.raises(ValueError, match="Existing assessment differs"):
        standardize_document(result, {row["record_identifier"]: row})


def test_bundled_outcomes_and_actual_reconstructed_values_match_named_ledger():
    rows = {
        r["record_identifier"]: r
        for r in json.loads((ROOT / REPORT).read_text())["records"]
    }
    checked = failed = 0
    for path in sorted((ROOT / "peritheos/data/materials").glob("*.eosmat")):
        document = load_eosmat(path)
        for record in document["eos_records"]:
            row = rows[record["identifier"]]
            assessment = record["evidence"]["reproduction"]["assessment"]
            assert assessment == assessment_for_record(record, row)
            assert assessment["status"] == STATUS[row["status"]]
            assert assessment["source_status"] == row["status"]
            assert assessment.get("parameters", []) == row.get("parameters", [])
            assert assessment["scope"] != "author_workflow"
            if row["status"] == "parity_not_achieved":
                assert assessment["parameters"]
                failed += 1
            checked += 1
    assert checked == len(rows) and failed > 0


def test_scope_distinguishes_independent_refit_and_source_equation():
    document, row = sample()
    record = document["eos_records"][0]
    record["record_kind"] = "refit"
    row.update(
        status="parity", original_publication_reproduction_status="not_reproduced"
    )
    assessment = assessment_for_record(record, row)
    assert assessment["scope"] == "independent_refit"
    assert any("not_reproduced" in v for v in assessment["qualifications"])
    assert "not recovery of the original publication" in assessment["qualifications"][0]
    row["status"] = "source_reconstruction"
    assert assessment_for_record(record, row)["scope"] == "source_equation"
    del row["parameters"]
    assert "parameters" not in assessment_for_record(record, row)


def test_comparison_survives_exports_without_replacing_calculation_parameters():
    document, row = sample()
    result = standardize_document(document, {row["record_identifier"]: row})
    legacy = convert_eosmat(result, version=3)
    restored = convert_eosmat(legacy, version=4)
    assert (
        restored["eos_records"][0]["evidence"]["reproduction"]
        == result["eos_records"][0]["evidence"]["reproduction"]
    )
    assert legacy["eos_records"][0]["eos"]["parameters"]["K0"] == 160
    assert (
        legacy["eos_records"][0]["reproduction"]["assessment"]["parameters"][0]["refit"]
        == 180.12345678901235
    )
    before = Material.from_eosmat(
        document, require_primary_validation=False
    ).eos_records[0]
    after = Material.from_eosmat(result, require_primary_validation=False).eos_records[
        0
    ]
    assert before.pressure(9) == after.pressure(9)
    rendered = json.loads(format_eosmat(result))
    assessment = rendered["eos_records"][0]["evidence"]["reproduction"]["assessment"]
    assert list(assessment)[:4] == ["status", "scope", "source_status", "basis"]
    assert list(assessment["parameters"][0])[:3] == [
        "parameter",
        "published",
        "published_error",
    ]
    sheet = "\n".join(reconstruction(legacy["eos_records"][0]))
    assert "180.12345678901235" in sheet and "160" in sheet
    assert "Fixed/adopted coefficients" in sheet and "do not replace" in sheet
    assessment = legacy["eos_records"][0]["reproduction"]["assessment"]
    assessment["parameters"].append(
        {"parameter": "K0_prime", "published": 4, "refit": 4}
    )
    sheet = "\n".join(reconstruction(legacy["eos_records"][0]))
    assert "| 4 | Not specified | 4 | Not specified | Fixed/adopted |" in sheet
    assessment["fit"] = {
        "fixed_parameters": "K0_prime was constrained",
        "free_parameters": None,
    }
    sheet = "\n".join(reconstruction(legacy["eos_records"][0]))
    assert "| 4 | Not specified | 4 | Not specified | See protocol |" in sheet


@pytest.mark.parametrize(
    "invalid",
    [
        {"parameter": "K0", "published": 160},
        {"parameter": "", "published": 160, "refit": 180},
        {"parameter": "K0", "published": True, "refit": 180},
        {"parameter": "K0", "published": 160, "refit": "180"},
        {"parameter": "K0", "published": 160, "refit": 180, "refit_error": -1},
    ],
)
def test_invalid_comparison_values_fail_both_schema_and_semantic_validation(invalid):
    document, row = sample()
    result = standardize_document(document, {row["record_identifier"]: row})
    result["eos_records"][0]["evidence"]["reproduction"]["assessment"]["parameters"] = [
        invalid
    ]
    assert not Draft202012Validator(eosmat_schema(4)).is_valid(result)
    with pytest.raises(EosmatError):
        validate_eosmat_document(result)


def test_duplicate_names_nonfinite_values_and_invalid_protocol_are_rejected():
    for invalid in ("duplicates", "nan", "fit"):
        document, row = sample()
        result = standardize_document(document, {row["record_identifier"]: row})
        assessment = result["eos_records"][0]["evidence"]["reproduction"]["assessment"]
        if invalid == "duplicates":
            assessment["parameters"] *= 2
        elif invalid == "nan":
            assessment["parameters"][0]["refit"] = float("nan")
        else:
            assessment["fit"] = []
        with pytest.raises(EosmatError):
            validate_eosmat_document(result)


def test_historical_calibration_fingerprint_ignores_only_added_assessment():
    document = get_material_document("mgo")
    projected = historical_reconstruction_view(document)
    record = next(
        r
        for r in projected["eos_records"]
        if r["identifier"] == "mgo_b1_tange_2009_vinet"
    )
    report = json.loads((ROOT / "docs/data/miozzi-2020-tange-refit.json").read_text())
    assert (
        hashlib.sha256(
            json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        == report["calibration_record_sha256"]
    )
    assert "assessment" in document["eos_records"][0]["reproduction"]
    record["eos"]["parameters"]["K0"] += 1
    assert (
        hashlib.sha256(
            json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        != report["calibration_record_sha256"]
    )
    authored, row = sample()
    authored = convert_eosmat(
        standardize_document(authored, {row["record_identifier"]: row}), version=3
    )
    authored["eos_records"][0]["reproduction"]["assessment"].pop("source_status")
    assert historical_reconstruction_view(authored) == authored
