"""Bundled v4 files stay native, scientifically complete and API-compatible."""

import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

from peritheos import (
    get_eos_record_document,
    get_material_document,
    list_material_documents,
    load_eosmat,
)
from peritheos.observations import observation_material_fingerprint
from scripts.eosmat_io import read_material, render_material, write_material
from scripts.migrate_eosmat_v4 import migrate_file, scientific_content

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("identifier", list_material_documents())
def test_native_catalog_and_compatibility_views_preserve_complete_content(identifier):
    path = ROOT / "peritheos/data/materials" / (identifier + ".eosmat")
    native = load_eosmat(path)
    assert native["format_version"] == 4
    assert "format" not in native
    assert "shared_metadata" not in native
    assert native["units"]["pressure"] == "GPa"
    assert native["units"]["temperature"] == "K"
    assert native["units"]["volume"] == "angstrom^3/conventional_unit_cell"
    assert native == get_material_document(identifier, format_version=None)
    assert native == get_material_document(identifier, format_version=4)
    legacy = get_material_document(identifier)
    assert legacy["format_version"] == 3
    assert scientific_content(native) == scientific_content(legacy)
    assert observation_material_fingerprint(native) == observation_material_fingerprint(
        legacy
    )
    assert native.get("datasets") == legacy.get("datasets")
    for record, old in zip(native["eos_records"], legacy["eos_records"]):
        assert record["model"] == old["eos"]["type"]
        assert (
            record["evidence"]["scientific_validation"] == old["scientific_validation"]
        )
        assert "eos" not in record
    assert render_material(native) == path.read_text(encoding="utf-8")
    assert not migrate_file(path, check=True)


def test_record_lookup_supports_native_fields_and_returns_defensive_copies():
    record = get_material_document("gold")["eos_records"][0]
    identifier = record["identifier"]
    assert get_eos_record_document(identifier) == record
    native = get_eos_record_document(identifier, format_version=4)
    assert native == get_eos_record_document(identifier, format_version=None)
    native["parameters"].clear()
    assert get_eos_record_document(identifier, format_version=4)["parameters"]


def test_migration_check_and_repeat_are_read_only_and_maintenance_keeps_v4(tmp_path):
    original = get_material_document("gold")
    original["extension"] = {"ordered": [3, 1, 2], "coefficient": 1.2345678901234567}
    before = copy.deepcopy(original)
    path = tmp_path / "gold.eosmat"
    payload = json.dumps(original).encode()
    path.write_bytes(payload)
    assert migrate_file(path, check=True)
    assert path.read_bytes() == payload
    assert migrate_file(path)
    native_bytes = path.read_bytes()
    assert scientific_content(load_eosmat(path)) == scientific_content(before)
    assert not migrate_file(path)
    assert path.read_bytes() == native_bytes
    maintenance = read_material(path)
    assert maintenance["format_version"] == 3
    write_material(path, maintenance)
    assert path.read_bytes() == native_bytes
    assert original == before
    result = subprocess.run(
        [sys.executable, "-m", "scripts.migrate_eosmat_v4", "--check", str(path)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_invalid_migration_preserves_existing_bytes(tmp_path):
    path = tmp_path / "invalid.eosmat"
    payload = b'{"format_version": 4, "identifier": "incomplete"}'
    path.write_bytes(payload)
    with pytest.raises(ValueError):
        migrate_file(path)
    assert path.read_bytes() == payload


def test_preservation_check_detects_number_spelling_and_extension_changes():
    original = get_material_document("gold")
    original["extension"] = {"value": 0.0, "parameter_errors": {}}
    changed = copy.deepcopy(original)
    changed["extension"]["value"] = -0.0
    assert scientific_content(original) != scientific_content(changed)
    changed["extension"]["value"] = 0
    assert scientific_content(original) != scientific_content(changed)
    changed["extension"] = {"value": 0.0}
    assert scientific_content(original) != scientific_content(changed)
