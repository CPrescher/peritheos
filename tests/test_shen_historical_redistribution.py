"""Admission preserves numerical evidence and separates facts from source grants."""

import csv
import hashlib
import io
import json
from collections import Counter

import numpy as np
import pytest

from peritheos import Dataset, get_material_document
from peritheos_build import ROOT, digest, stage_datasets

PROOF_PATH = ROOT / "docs/data/shen-historical-redistribution-review.json"
PROOF = json.loads(PROOF_PATH.read_bytes())
DATA = ROOT / "peritheos/data"


@pytest.mark.parametrize("key", sorted(PROOF["declaration_bindings"]))
def test_admitted_original_declarations_and_source_bytes_remain_exact(key):
    material, identifier = key.split("/")
    original = json.loads((DATA / "materials" / f"{material}.eosmat").read_bytes())
    declaration = next(d for d in original["datasets"] if d["identifier"] == identifier)
    assert digest(declaration) == PROOF["declaration_bindings"][key]
    relative = "peritheos/data/" + declaration["resource"]["path"]
    assert (
        hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        == PROOF["resource_bindings"][relative]
    )
    dataset = Dataset.from_mapping(declaration)
    expected = (
        3511
        if identifier.startswith("shen_smith")
        else 20
        if material == "gold"
        else 25
    )
    assert len(dataset) == expected
    if material == "cu_handbook_1972_legacy":
        assert declaration["license"].startswith("All rights reserved source")
    else:
        assert "license" not in declaration


def test_shen_pair_identities_errors_and_missing_readout_are_preserved():
    document = get_material_document("gold")
    declaration = next(
        d
        for d in document["datasets"]
        if d["identifier"] == "shen_smith_2026_table_s1_simultaneous_volumes"
    )
    dataset = Dataset.from_mapping(declaration)
    profile = PROOF["source_evidence"]["shen_smith"]["profile"]
    assert Counter(dataset["phase"].tolist()) == profile["phase_counts"]
    sequences = []
    for value in dataset["measurement_sequence"].tolist():
        if isinstance(value, str):
            sequences.append(value)
        else:
            assert np.isnan(value)  # Dataset represents the original blank as missing.
            sequences.append("")
    assert Counter(sequences) == profile["sequence_counts"]
    assert np.all(dataset["unit_cell_volume_standard_error_a3"] > 0)
    assert np.all(dataset["unit_cell_volume_a3"] > 0)
    rows = list(
        csv.DictReader(
            io.StringIO((DATA / declaration["resource"]["path"]).read_text())
        )
    )
    keys = {
        (r["experiment"], r["run_number"], r["phase"], r["measurement_sequence"])
        for r in rows
    }
    assert len(keys) == 3511
    assert ("DAC-2", "194", "Pt", "first") not in keys
    assert ("DAC-2", "194", "Pt", "last") in keys
    for experiment, count in [("DAC-1", 228), ("DAC-2", 194)]:
        assert {
            int(r["run_number"]) for r in rows if r["experiment"] == experiment
        } == set(range(1, count + 1))
    assert set(profile["columns"]) == set(rows[0])
    assert "pressure_gpa" not in rows[0]


def test_historical_source_interpretations_blanks_and_branches_are_not_replaced():
    copper = next(
        d
        for d in get_material_document("cu_handbook_1972_legacy")["datasets"]
        if d["identifier"] == "aip_handbook_table4d12_cu"
    )
    rows = list(
        csv.DictReader(io.StringIO((DATA / copper["resource"]["path"]).read_text()))
    )
    assert rows[0] == {
        "pressure_kbar": "15",
        "pressure_gpa": "1.5",
        "v_over_v0": "0.9900",
        "note": "printed as 09900; corrected printing error",
    }
    assert len(rows) == 25
    assert not {"0", "5", "10"}.intersection(r["pressure_kbar"] for r in rows)
    gold = next(
        d
        for d in get_material_document("gold")["datasets"]
        if d["identifier"] == "gold_anderson_1989_table1_thermodynamic_inputs"
    )
    rows = list(
        csv.DictReader(io.StringIO((DATA / gold["resource"]["path"]).read_text()))
    )
    branches = [r for r in rows if r["temperature_k"] == "300"]
    assert len(branches) == 2
    assert branches[0]["cp_printed_0p1_j_g_k"] == "1.1288"
    assert branches[0]["cp_reconstructed_0p1_j_g_k"] == "1.2888"
    assert branches[1]["cp_printed_0p1_j_g_k"] == ""
    assert branches[0]["c_prime_gpa"] == "14.60"
    assert branches[1]["c_prime_gpa"] == "14.85"
    assert gold["kind"] == "heterogeneous_thermodynamic_compilation"


def test_installed_review_and_notices_bind_admissions_without_bundling_originals(
    tmp_path,
    independent_generation_baseline,
):
    assert (
        independent_generation_baseline["proof"]["before_catalog_scientific_sha256"]
        == PROOF["catalog_scientific_sha256"]
    )
    prior = json.loads(
        (ROOT / "docs/data/morard-campbell-redistribution-followup.json").read_bytes()
    )
    assert (
        PROOF["catalog_scientific_sha256"] == prior["after_catalog_scientific_sha256"]
    )
    package = stage_datasets(tmp_path / "optional") / "peritheos_datasets"
    manifest = json.loads((package / "manifest.json").read_bytes())
    review = json.loads((package / "redistribution-review.json").read_bytes())
    assert manifest["dataset_count"] == 451
    assert review["shen_historical_review"] == PROOF
    selected = {
        e["review_key"]: e
        for e in manifest["policy"]["datasets"]
        if e.get("review_key")
    }
    decisions = {e["review_key"]: e for e in review["entries"]}
    for key in PROOF["declaration_bindings"]:
        assert decisions[key]["decision"] == "admit"
        assert decisions[key]["prior_decision"]["decision"] == "hold"
        assert (
            decisions[key]["source_status"] == decisions[key]["previous_source_status"]
        )
        assert (
            decisions[key]["followup_evidence"]["sha256"]
            == hashlib.sha256(PROOF_PATH.read_bytes()).hexdigest()
        )
        assert selected[key]["declaration_sha256"] == PROOF["declaration_bindings"][key]
    assert (
        package / "notices/shen-historical-collections.LICENSE.md"
    ).read_bytes() == (
        DATA / "datasets/shen-historical-collections.LICENSE.md"
    ).read_bytes()
    assert not any(
        p.suffix.lower() in {".pdf", ".xlsx", ".xls"} for p in package.rglob("*")
    )
    held = [e for e in review["entries"] if e["decision"] == "hold"]
    assert len(held) == 5
    assert {e["basis"] for e in held} == {"source_nc_nd", "large_collection"}
