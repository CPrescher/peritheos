"""Human-facing examples are complete, honest and executable."""

import copy
import json
import re
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, convert_eosmat, eosmat_schema, load_eosmat
from peritheos._eosmat_audit import check_block, derived_report_links
from peritheos._eosmat_metadata import SCOPES, STATUSES
from peritheos.errors import EosmatError

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "name",
    ["simple-material", "thermal-material", "curated-material", "shared-material"],
)
def test_complete_teaching_files_are_native_unreviewed_and_numerically_compatible(name):
    document = load_eosmat(ROOT / "docs/examples" / f"{name}.eosmat")
    Draft202012Validator(eosmat_schema(5)).validate(document)
    assert document["units"]["pressure"] == "GPa"
    assert document["units"]["volume"] == "angstrom^3/conventional_unit_cell"
    native = Material.from_eosmat(document, require_primary_validation=False)
    old = Material.from_eosmat(
        convert_eosmat(document, version=3), require_primary_validation=False
    )
    for a, b in zip(native.eos_records, old.eos_records):
        assert a.scientific_validation_status == "pending_primary_source_check"
        assert float(a.pressure(0.99 * a.reference_volume)) == float(
            b.pressure(0.99 * b.reference_volume)
        )
    for dataset in document.get("datasets", []):
        assert dataset["kind"] == "calculated_checkpoint"
        assert not dataset["used_by_eos_records"]


def test_documented_record_file_and_audit_sidecar_load_as_described(tmp_path):
    """The split example in docs/eosmat-v5.md is valid, consistent and complete."""
    reference = (ROOT / "docs/eosmat-v5.md").read_text(encoding="utf-8")
    blocks = {
        name: json.loads(
            re.search(
                rf"<!-- example: {name} -->\n```json\n(.*?)\n```", reference, re.DOTALL
            )[1]
        )
        for name in ("record file", "audit sidecar")
    }
    record_file, sidecar = blocks["record file"], blocks["audit sidecar"]
    (tmp_path / "materials").mkdir()
    (tmp_path / "audits").mkdir()
    path = tmp_path / "materials/example_oxide.eosmat"
    path.write_text(json.dumps(record_file), encoding="utf-8")
    (tmp_path / "audits/example_oxide.audit.json").write_text(
        json.dumps(sidecar), encoding="utf-8"
    )
    for document in (record_file, load_eosmat(path)):
        Draft202012Validator(eosmat_schema(5)).validate(document)
    assert load_eosmat(path, include_audit=False) == record_file
    merged = load_eosmat(path)
    record = merged["eos_records"][0]
    expected = copy.deepcopy(sidecar["records"]["example_oxide_bm3"]["evidence"])
    # As documented, absent report links are filled in from the identifiers.
    expected["reproduction"]["assessment"].update(
        derived_report_links("example_oxide", "example_oxide_bm3")
    )
    assert record["evidence"] == expected
    stored = {k: v for k, v in record["check"].items() if k != "audit"}
    assert stored == check_block(record, merged)
    loaded = Material.from_eosmat(merged, require_primary_validation=False)
    assert (loaded.eos_records[0].reproduction, loaded.eos_records[0].warning) == (
        "reproduced",
        None,
    )
    assert loaded.eos_records[0].scientific_validation_status == (
        "pending_primary_source_check"
    )
    unmerged = Material.from_eosmat(record_file, require_primary_validation=False)
    assert unmerged.eos_records[0].reproduction == "reproduced"
    with pytest.raises(EosmatError) as error:
        unmerged.eos_records[0].audit
    assert error.value.code == "eosmat.audit_not_loaded"
    legacy = convert_eosmat(merged, version=3)["eos_records"][0]
    assert "check" not in legacy and legacy["reproduction"]["assessment"]
    (tmp_path / "audits/example_oxide.audit.json").unlink()
    with pytest.raises(EosmatError) as error:
        load_eosmat(path)
    assert error.value.code == "eosmat.audit_missing"


def test_label_guide_covers_supported_assessment_terms_and_actual_generic_thresholds():
    guide = (ROOT / "docs/eos-labels.md").read_text()
    for value in [
        *STATUSES,
        *SCOPES,
        "parity_not_achieved",
        "primary_source_validated",
        "deferred",
        "pending_primary_source_check",
    ]:
        assert f"`{value}`" in guide
    assert "Fixed/adopted coefficients" in guide
    assert "5%" in guide and "15%" in guide and "25%" in guide and "30%" in guide


def test_dataset_citations_cannot_hide_a_shared_metadata_dependency():
    document = load_eosmat(ROOT / "docs/examples/curated-material.eosmat")
    document["datasets"][0]["reference"] = {"shared_ref": "source"}
    with pytest.raises(ValueError, match="Dataset references must remain inline"):
        convert_eosmat(document)
