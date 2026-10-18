"""Compact native records preserve source semantics and reject hidden overrides."""

import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    convert_eosmat,
    create_eosmat,
    eosmat_schema,
    load_eosmat,
    save_eosmat,
    validate_eosmat_document,
)
from peritheos.__main__ import main
from peritheos._eosmat_format import format_eosmat
from peritheos.observations import observation_material_fingerprint

ROOT = Path(__file__).resolve().parents[1]


def equivalent(value):
    if isinstance(value, dict):
        return {
            k: equivalent(v)
            for k, v in value.items()
            if not (k in {"parameter_errors", "fixed_parameters"} and not v)
        }
    if isinstance(value, list):
        return [equivalent(v) for v in value]
    return value


def sample():
    document = create_eosmat(
        "example", "X", "BM3", {"V0": 10, "K0": 160, "K0_prime": 4}, temperature_ref=300
    )
    citation = {
        "authors": ["Example"],
        "year": 2026,
        "source": "Fictional teaching source",
        "title": "Fictional citation for parser tests, not scientific evidence. " * 4,
    }
    document["eos_records"][0]["reference"] = citation
    sibling = copy.deepcopy(document["eos_records"][0])
    sibling["identifier"] = "example_vinet"
    sibling["model"] = "Vinet"
    document["eos_records"].append(sibling)
    document["extension"] = {
        "shared_ref": "literal_application_extension",
        "matrix": [[3, None], [1, 2]],
    }
    return document


def test_entire_catalog_compacts_without_changing_any_original_scientific_field():
    schema = Draft202012Validator(eosmat_schema(4))
    for path in sorted((ROOT / "peritheos/data/materials").glob("*.eosmat")):
        original = json.loads(path.read_text())
        before = copy.deepcopy(original)
        compact = convert_eosmat(original, compact=True)
        schema.validate(compact)
        restored = convert_eosmat(compact, version=3)
        # A generated sharing pool is additional, duplicated metadata. All
        # original fields, source arrays and unknown extensions remain exact.
        restored.pop("shared_metadata", None)
        assert equivalent(restored) == equivalent(
            convert_eosmat(original, version=3)
        ), path
        assert compact.get("datasets") == original.get("datasets"), path
        assert observation_material_fingerprint(
            compact
        ) == observation_material_fingerprint(original)
        assert convert_eosmat(compact, compact=True) == compact
        assert original == before


def test_shared_file_is_native_lossless_and_numerically_equivalent(tmp_path):
    original = sample()
    compact = convert_eosmat(original, compact=True)
    assert len(format_eosmat(compact)) < len(format_eosmat(original))
    assert compact["eos_records"][0]["reference"] == {"shared_ref": "example_2026"}
    assert compact["extension"] == original["extension"]
    before = Material.from_eosmat(original, require_primary_validation=False)
    after = Material.from_eosmat(compact, require_primary_validation=False)
    for a, b in zip(before.eos_records, after.eos_records):
        assert float(a.pressure(9)) == float(b.pressure(9))
        assert a.reference == b.reference
    path = tmp_path / "compact.eosmat"
    save_eosmat(path, compact)
    assert load_eosmat(path) == compact
    assert after.to_eosmat()["shared_metadata"] == compact["shared_metadata"]
    assert (
        after.to_eosmat()["eos_records"][0]["reference"]
        == original["eos_records"][0]["reference"]
    )
    inline = convert_eosmat(compact, inline=True)
    assert inline["eos_records"] == original["eos_records"]
    assert inline["extension"] == original["extension"]


def test_unused_definitions_and_nested_references_survive_compatibility_exports():
    document = sample()
    document["shared_metadata"] = {
        "citation": document["eos_records"][0]["reference"],
        "nested": {"shared_ref": "citation"},
        "unused": {
            "source_note": "Retain supporting evidence even if no record cites it."
        },
    }
    document["eos_records"][0]["reference"] = {"shared_ref": "nested"}
    for version in [3, 4]:
        inline = convert_eosmat(document, version=version, inline=True)
        assert (
            inline["shared_metadata"]["unused"] == document["shared_metadata"]["unused"]
        )
        assert (
            inline["eos_records"][0]["reference"]
            == document["shared_metadata"]["citation"]
        )
        assert "shared_ref" not in json.dumps(inline["shared_metadata"])


def test_inline_citation_extensions_are_not_mistaken_for_reference_overrides():
    document = sample()
    document["eos_records"][0]["reference"]["shared_ref"] = "application_annotation"
    validate_eosmat_document(document)
    assert (
        convert_eosmat(document, version=3)["eos_records"][0]["reference"]
        == document["eos_records"][0]["reference"]
    )


@pytest.mark.parametrize(
    "pool,reference",
    [
        ({}, {"shared_ref": "missing"}),
        ({"a": {"shared_ref": "a"}}, {"shared_ref": "a"}),
        ({"a": {"shared_ref": "b"}, "b": {"shared_ref": "a"}}, {"shared_ref": "a"}),
        ({"a": {}}, {"shared_ref": "a", "title": "Hidden override"}),
        ({"a": {}}, {"shared_ref": ["a"]}),
        ({"bad-name": {}}, "Plain reference"),
        ({"a": "not an object"}, "Plain reference"),
    ],
)
def test_invalid_references_are_rejected(pool, reference):
    document = sample()
    document["shared_metadata"] = pool
    document["eos_records"][0]["reference"] = reference
    with pytest.raises(ValueError):
        validate_eosmat_document(document)


def test_expansion_budget_rejects_exponential_metadata():
    document = sample()
    pool = {"a": {"note": "source"}}
    name = "a"
    for i in range(18):
        following = f"layer_{i}"
        pool[following] = {"left": {"shared_ref": name}, "right": {"shared_ref": name}}
        name = following
    document["shared_metadata"] = pool
    with pytest.raises(ValueError, match="expansion"):
        validate_eosmat_document(document)


@pytest.mark.parametrize(
    "key,value",
    [("eos_records", None), ("eos_records", 1), ("datasets", None), ("datasets", 1)],
)
def test_invalid_collection_types_keep_structural_validation_errors(key, value):
    document = sample()
    document[key] = value
    with pytest.raises(ValueError):
        validate_eosmat_document(document)


def test_explicit_reconstruction_assessment_is_scoped_and_never_inferred():
    document = sample()
    record = document["eos_records"][0]
    assert "evidence" not in record
    assessment = {
        "status": "parity",
        "scope": "parameter_agreement",
        "basis": "Within reported parameter error widths under stated fixed constraints",
        "qualifications": [
            "Author weights, covariance and original workflow are unrecovered."
        ],
    }
    record["evidence"] = {
        "reproduction": {
            "assessment": assessment,
            "legacy_note": "Keep historical findings",
        }
    }
    validate_eosmat_document(document)
    Draft202012Validator(eosmat_schema(4)).validate(document)
    assert (
        convert_eosmat(document, version=3)["eos_records"][0]["reproduction"][
            "assessment"
        ]
        == assessment
    )
    for key, value in [
        ("status", "primary_source_validated"),
        ("scope", "unknown"),
        ("basis", ""),
        ("qualifications", [None]),
    ]:
        invalid = copy.deepcopy(document)
        invalid["eos_records"][0]["evidence"]["reproduction"]["assessment"][key] = value
        with pytest.raises(ValueError):
            validate_eosmat_document(invalid)
        assert list(Draft202012Validator(eosmat_schema(4)).iter_errors(invalid))


def test_installed_commands_work_outside_checkout_without_repository_scripts(tmp_path):
    path = tmp_path / "source.eosmat"
    path.write_text(json.dumps(sample()))
    assert main(["validate", str(path)]) == 0
    assert main(["format", "--check", str(path)]) == 1
    assert main(["format", str(path)]) == 0
    assert main(["format", "--check", str(path)]) == 0
    assert load_eosmat(path) == sample()
    target = tmp_path / "compact.eosmat"
    assert main(["convert", str(path), str(target), "--compact"]) == 0
    inline = tmp_path / "inline.eosmat"
    assert main(["convert", str(target), str(inline), "--inline"]) == 0
    result = subprocess.run(
        [sys.executable, "-m", "peritheos", "validate", str(inline)],
        cwd=tmp_path,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    original = path.read_bytes()
    bad = tmp_path / "bad.eosmat"
    bad.write_text("{}")
    with pytest.raises(SystemExit):
        main(["format", str(path), str(bad)])
    assert path.read_bytes() == original


def test_native_formatter_puts_all_eos_before_long_structure_and_does_not_change_values():
    document = sample()
    document["atom_sites"] = [{"element": "X", "x": 0, "y": 0, "z": 0}]
    rendered = format_eosmat(document)
    assert rendered.index('"eos_records"') < rendered.index('"atom_sites"')
    assert json.loads(rendered) == document
