"""Human-facing examples are complete, honest and executable."""

from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, convert_eosmat, eosmat_schema, load_eosmat
from peritheos._eosmat_metadata import SCOPES, STATUSES

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "name",
    ["simple-material", "thermal-material", "curated-material", "shared-material"],
)
def test_complete_teaching_files_are_native_unreviewed_and_numerically_compatible(name):
    document = load_eosmat(ROOT / "docs/examples" / f"{name}.eosmat")
    Draft202012Validator(eosmat_schema(4)).validate(document)
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
