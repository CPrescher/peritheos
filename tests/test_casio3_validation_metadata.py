"""Source verification must not imply original-fit or measurement validation."""

import hashlib
import json
from pathlib import Path

import jsonschema
import pytest

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]
SUN = "ca_perovskite_sun_2016_bm3_3"
FU = "ca_perovskite_fu_2023_bm3_mgd_refit"
REFIT = "ca_perovskite_fu_2023_candidate_data_unweighted_bm3_mgd_refit"
CHEN = "ca_perovskite_tetragonal_chen_2018_vinet"


def records():
    return {
        row["identifier"]: row
        for material in ("ca_perovskite", "ca_perovskite_tetragonal")
        for row in get_material_document(material)["eos_records"]
    }


@pytest.mark.parametrize(
    "material, scientific_payload_sha256",
    [
        (
            "ca_perovskite",
            "a2a224de2f134b3dc7b80aaaa1bea913bc6818ae24ec72f5313cf8559b92230d",
        ),
        (
            "ca_perovskite_tetragonal",
            "f482eba6ddf5d35f0c637104debaebdbc07dad529f8cf2dd954bcc1ba9242b34",
        ),
    ],
)
def test_qualification_preserves_all_scientific_payloads_and_schema(
    material, scientific_payload_sha256
):
    document = get_material_document(material)
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    jsonschema.validate(document, schema)
    # Reviewed scientific payload checkpoint, including the 2026-10-08 Wang
    # source recovery, excluded-row corrections and optional 64-point thermal refit. Freezes
    # every other record field, dataset, checksum, provenance and classification.
    for row in document["eos_records"]:
        row.pop("scientific_validation")
    payload = json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(payload).hexdigest() == scientific_payload_sha256


def test_fu_source_verification_and_independent_refit_are_distinct():
    by_id = records()
    published, refit = by_id[FU], by_id[REFIT]
    check = published["scientific_validation"]
    assert published["record_kind"] == "published"
    assert check["status"] == "primary_source_validated"
    assert check["reproduction_status"] == "not_reproduced"
    assert "NOT REPRODUCED" in check["note"]
    assert "every prediction is invalid" in check["note"]
    assert "fit_datasets" not in published
    assert check["primary_data_check"]["observation_count"] == 174
    assert refit["record_kind"] == "refit"
    assert refit["derived_from_record"] == FU
    assert refit["scientific_validation"]["reproduction_status"] == "refit_reproduced"
    assert (
        refit["scientific_validation"]["original_publication_reproduction_status"]
        == "not_reproduced"
    )


def test_sun_recovery_closure_is_a_report_and_keeps_workbook_uncertainties():
    check = records()[SUN]["scientific_validation"]
    assert (
        check["reproduction_status"] == "parameters_reproduced_at_published_precision"
    )
    recovery = check["original_data_recovery"]
    assert recovery["recovery_status"] == "closed"
    assert recovery["evidence_type"] == "user_reported_author_communication"
    assert recovery["author"] == "Ningyu Sun"
    assert "independently archived" in recovery["finding"]
    workbook = recovery["supplied_workbook"]
    assert workbook["filename"] == "Sun et al. 2016 data.xlsx"
    assert len(workbook["unresolved"]) == 4
    assert "no workbook values or blanket corrections" in workbook["disposition"]
    assert check["primary_data_check"]["status"] == "bundled"


def test_chen_rounded_table_check_retains_missing_weighting_and_scope():
    record = records()[CHEN]
    check = record["scientific_validation"]
    assert record["record_kind"] == "published"
    assert check["reproduction_status"] == "qualified_numerical_check"
    assert "objective weights" in check["note"]
    assert "not additional independent" in check["usage_recommendation"]
    assert "2-sigma" in check["usage_recommendation"]


def test_aggregate_ledger_carries_qualifications_without_rewriting_sources():
    ledger = json.loads(
        (ROOT / "peritheos/data/primary-source-audit.json").read_text(encoding="utf-8")
    )
    by_id = {row["record"]: row for row in ledger["records"]}
    for identifier in (SUN, FU, REFIT, CHEN):
        check = records()[identifier]["scientific_validation"]
        for field in ("status", "note", "reproduction_status", "usage_recommendation"):
            assert by_id[identifier][field] == check[field]
