import json
from collections import Counter
from pathlib import Path

import pytest

from peritheos import Material, get_material_document, validate_eosmat_document
from scripts.reproduce_cohen_lin_2014_fesio3 import reproduce

ROOT = Path(__file__).parents[1]
DOI = "10.1103/physrevb.90.140102"
MATERIALS = (
    "fesio3_bridgmanite",
    "fesio3_post_perovskite",
    "fesio3_post_perovskite_ii",
)


def test_cohen_lin_materials_are_phase_specific_and_normatively_valid():
    documents = [get_material_document(identifier) for identifier in MATERIALS]
    assert {document["formula"] for document in documents} == {"FeSiO3"}
    assert {document["formula_units_per_cell"] for document in documents} == {4}
    assert len({document["phase"] for document in documents}) == 3
    for document in documents:
        validate_eosmat_document(document)
        records = [
            record
            for record in document["eos_records"]
            if record["reference"]["doi"].lower() == DOI
        ]
        assert len(records) == 1
        record = records[0]
        assert record["reference"]["doi"].lower() == DOI
        assert record["temperature_ref"] == 0.0
        assert record["pressure_range_status"] == "theoretical"
        assert record["scientific_validation"]["status"] == ("primary_source_validated")


def test_cohen_lin_converts_table_iii_formula_volumes_to_z4_cells():
    expected = {
        "fesio3_bridgmanite": (44.31, 225.0, 4.42),
        "fesio3_post_perovskite": (44.90, 189.0, 4.73),
        "fesio3_post_perovskite_ii": (45.45, 195.0, 4.67),
    }
    for identifier, (v0_formula, k0, kp) in expected.items():
        record = next(
            record
            for record in get_material_document(identifier)["eos_records"]
            if record["reference"]["doi"].lower() == DOI
        )
        parameters = record["eos"]["parameters"]
        assert parameters == {"V0": 4.0 * v0_formula, "K0": k0, "K0_prime": kp}
        executable = Material.from_eosmat(get_material_document(identifier))
        eos = executable.get_eos_record(record["identifier"])
        assert eos.pressure(parameters["V0"]) == pytest.approx(0.0, abs=1.0e-12)


def test_cohen_lin_table_i_structures_preserve_composition_and_cell_volume():
    ppv = get_material_document("fesio3_post_perovskite")
    ppv_ii = get_material_document("fesio3_post_perovskite_ii")
    assert ppv["space_group"] == "Cmcm"
    assert ppv_ii["space_group"] == "Cmmm"
    for document in (ppv, ppv_ii):
        contents = Counter()
        for site in document["atom_sites"]:
            contents[site["element"]] += site["site_multiplicity"] * site["occupancy"]
        assert contents == {"Fe": 4.0, "Si": 4.0, "O": 12.0}
    ppv_volume = ppv["lattice"]["a"] * ppv["lattice"]["b"] * ppv["lattice"]["c"]
    ppv_ii_volume = (
        ppv_ii["lattice"]["a"] * ppv_ii["lattice"]["b"] * ppv_ii["lattice"]["c"]
    )
    assert ppv_volume / 4.0 == pytest.approx(33.98, abs=0.05)
    assert ppv_ii_volume / 4.0 == pytest.approx(34.49, abs=0.05)


def test_cohen_lin_reproduction_matches_independent_100_gpa_checkpoints():
    result = reproduce()
    assert result["pressure_gpa"] == 100.0
    for phase in result["phases"].values():
        differences = phase["absolute_differences"]
        assert differences["V100_a3_per_formula"] < 0.04
        assert differences["K100_gpa"] < 1.2
        assert differences["K100_prime"] < 0.004


def test_cohen_lin_audit_disposes_all_three_litcurate_candidates_once():
    audit = (
        ROOT / "docs" / "literature-reproductions" / "cohen-lin-2014-fesio3.md"
    ).read_text(encoding="utf-8")
    candidates = json.loads(
        (ROOT / "docs" / "data" / "litcurate-eos-candidates.json").read_text(
            encoding="utf-8"
        )
    )["records"]
    same_doi = [row for row in candidates if row["publication"]["doi"].lower() == DOI]
    assert len(same_doi) == 3
    for row in same_doi:
        assert audit.count(row["identifier"]) == 1
