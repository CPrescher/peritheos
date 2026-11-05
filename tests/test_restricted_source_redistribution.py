"""Admission changes selection, preserving source evidence and qualifications."""

import csv
import hashlib
import json
from collections import Counter

import pytest

from peritheos_build import ROOT, digest, stage_datasets

PROOF_PATH = ROOT / "docs/data/restricted-source-numerical-review.json"
PROOF = json.loads(PROOF_PATH.read_bytes())
HELD_KEY = "gold/gold_sakai_2025_table_s8_lineage"


@pytest.mark.parametrize("key", sorted(PROOF["cases"]))
def test_original_declarations_resources_and_scientific_kinds_remain_exact(
    key,
    independent_generation_baseline,
):
    material, identifier = key.split("/")
    relative = f"peritheos/data/materials/{material}.eosmat"
    payload = independent_generation_baseline["original_material_bytes"].get(
        material, (ROOT / relative).read_bytes()
    )
    assert (
        hashlib.sha256(payload).hexdigest() == PROOF["material_file_bindings"][relative]
    )
    original = json.loads(payload)
    declaration = next(d for d in original["datasets"] if d["identifier"] == identifier)
    assert digest(declaration) == PROOF["declaration_bindings"][key]
    relative = "peritheos/data/" + declaration["resource"]["path"]
    assert (
        hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        == PROOF["resource_bindings"][relative]
    )
    with (ROOT / relative).open() as stream:
        reader = csv.DictReader(stream)
        profile = PROOF["cases"][key]["profile"]
        assert reader.fieldnames == profile["columns"]
        assert len(list(reader)) == profile["rows"]
    if "sakai_2025_table_s9" in identifier:
        assert declaration["kind"] == "derived_pressure_volume_grid"
    elif "figure4" in identifier:
        assert declaration["kind"] == "pressure_volume_temperature_digitized"
    elif "ismailova" in identifier:
        assert declaration["license"] == "CC BY-NC 4.0"
        assert declaration["kind"] == "selected_pressure_lattice_volume_states"


def test_optional_package_keeps_restricted_source_flags_and_excludes_lineage(
    tmp_path,
    independent_generation_baseline,
):
    assert (
        independent_generation_baseline["proof"]["before_catalog_scientific_sha256"]
        == PROOF["catalog_scientific_sha256"]
    )
    package = stage_datasets(tmp_path / "optional") / "peritheos_datasets"
    manifest = json.loads((package / "manifest.json").read_bytes())
    review = json.loads((package / "redistribution-review.json").read_bytes())
    assert review["restricted_source_review"] == PROOF
    assert manifest["dataset_count"] == 451
    entries = {e["review_key"]: e for e in review["entries"]}
    selected = {
        e["review_key"]: e
        for e in manifest["policy"]["datasets"]
        if e.get("review_key")
    }
    assert len(PROOF["cases"]) == 22
    assert Counter(c["decision"] for c in PROOF["cases"].values()) == {
        "admit": 21,
        "hold": 1,
    }
    for key, case in PROOF["cases"].items():
        entry = entries[key]
        assert case["source_status"] == entry["source_status"] == "restricted_source"
        assert entry["source_status"] == entry["prior_decision"]["source_status"]
        assert entry["prior_decision"]["commit"] == "294a651f"
        assert entry["prior_decision"]["decision"] == "hold"
        assert (
            entry["followup_evidence"]["sha256"]
            == hashlib.sha256(PROOF_PATH.read_bytes()).hexdigest()
        )
        relative = entry["resource"]["path"]
        if key == HELD_KEY:
            assert entry["decision"] == "hold"
            assert key not in selected
            assert not (package / "data" / relative).exists()
        else:
            assert entry["decision"] == "admit"
            assert selected[key]["source_status"] == "restricted_source"
            assert (
                selected[key]["redistribution_basis"]
                == "restricted_source_numerical_assessment"
            )
            assert (package / "data" / relative).read_bytes() == (
                ROOT / "peritheos/data" / relative
            ).read_bytes()
    assert len([e for e in review["entries"] if e["decision"] == "hold"]) == 5
    assert (
        sum(e["source_status"] == "restricted_source" for e in review["entries"]) == 22
    )
    assert (
        package / "notices/restricted-source-numerical.LICENSE.md"
    ).read_bytes() == (
        ROOT / "peritheos/data/datasets/restricted-source-numerical.LICENSE.md"
    ).read_bytes()
    assert not any(
        p.suffix.lower() in {".pdf", ".xlsx", ".xls", ".png", ".jpg", ".zip"}
        for p in package.rglob("*")
    )
