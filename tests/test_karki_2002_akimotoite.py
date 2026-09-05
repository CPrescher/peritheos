import json
from pathlib import Path

import pytest

from peritheos import Material, get_material_document, validate_eosmat_document
from scripts.reproduce_karki_2002_akimotoite import reproduce

ROOT = Path(__file__).parents[1]
DOI = "10.1029/2001jb000702"


def _records():
    document = get_material_document("akimotoite")
    return [
        record
        for record in document["eos_records"]
        if record["reference"]["doi"].lower() == DOI
    ]


def test_karki_four_temperature_specific_bm4_records_are_stored():
    document = get_material_document("akimotoite")
    validate_eosmat_document(document)
    records = _records()
    assert len(records) == 4
    assert [record["temperature_ref"] for record in records] == [
        0.0,
        300.0,
        1000.0,
        2000.0,
    ]
    assert all(record["eos"]["type"] == "BM4" for record in records)
    assert all("default_for" not in record for record in records)


def test_karki_table_three_coefficients_and_cell_basis_are_exact():
    expected = [
        (261.66, 210.0, 4.57, -0.041),
        (265.20, 201.0, 4.64, -0.042),
        (270.18, 182.0, 4.86, -0.051),
        (291.54, 153.0, 5.20, -0.067),
    ]
    for record, values in zip(_records(), expected):
        parameters = record["eos"]["parameters"]
        assert tuple(parameters.values()) == pytest.approx(values)
        assert record["volume_basis"]["formula_units"] == 6.0


def test_karki_bm4_curves_execute_and_recover_reference_values():
    document = get_material_document("akimotoite")
    material = Material.from_eosmat(document)
    for stored in _records():
        parameters = stored["eos"]["parameters"]
        eos = material.get_eos_record(stored["identifier"])
        assert eos.pressure(parameters["V0"]) == pytest.approx(0.0, abs=1.0e-12)
        assert eos.eos.bulk_modulus(parameters["V0"]) == pytest.approx(
            parameters["K0"], abs=1.0e-10
        )


def test_karki_reproduction_checks_conversion_and_temperature_trends():
    result = reproduce()
    assert result["formula_units_per_cell"] == 6.0
    assert result["records"]["static"]["V0_a3_conventional_z6"] == pytest.approx(261.66)
    assert result["records"]["2000_K"]["V0_a3_conventional_z6"] == pytest.approx(291.54)
    assert result["temperature_trends"] == {
        "volume_increases_monotonically": True,
        "bulk_modulus_decreases_monotonically": True,
    }
    for row in result["records"].values():
        assert row["pressure_at_v0_gpa"] == pytest.approx(0.0, abs=1.0e-12)
        assert row["bulk_modulus_at_v0_gpa"] == pytest.approx(row["K0_gpa"])
        assert row["pressure_at_0.90_v0_gpa"] > 0.0


def test_karki_audit_disposes_all_five_litcurate_candidates_once():
    audit = (
        ROOT / "docs" / "literature-reproductions" / "karki-2002-akimotoite.md"
    ).read_text(encoding="utf-8")
    candidates = json.loads(
        (ROOT / "docs" / "data" / "litcurate-eos-candidates.json").read_text(
            encoding="utf-8"
        )
    )["records"]
    same_doi = [row for row in candidates if row["publication"]["doi"].lower() == DOI]
    assert len(same_doi) == 5
    for row in same_doi:
        assert audit.count(row["identifier"]) == 1
