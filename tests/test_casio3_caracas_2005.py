import math

import pytest

from peritheos import Material, get_material_document
from scripts.apply_primary_source_audit import audit_record

RECORD_ID = "ca_perovskite_caracas_2005_bm3_3"


def test_caracas_2005_adds_exactly_one_cubic_literature_record():
    document = get_material_document("ca_perovskite")
    records = [
        record
        for record in document["eos_records"]
        if record["reference"].get("doi", "").lower() == "10.1029/2004gl022144"
    ]

    assert len(records) == 1
    record = records[0]
    assert record["identifier"] == RECORD_ID
    assert record["eos"] == {
        "type": "BM3",
        "parameters": {"V0": 44.579, "K0": 250.0, "K0_prime": 4.098},
        "model": "birch_murnaghan_3",
    }
    assert record["fixed_parameters"] == []
    assert record["parameter_errors"] == {
        "V0": None,
        "K0": None,
        "K0_prime": None,
    }
    assert record["temperature_ref"] == 0.0
    assert record["pressure_range_status"] == "theoretical"
    assert record["experimental_pressure_range_gpa"] == [0.0, 160.0]
    assert record["pressure_calibration"]["status"] == "not_applicable"


def test_caracas_2005_audit_distinguishes_all_bm3_and_bm4_rows():
    document = get_material_document("ca_perovskite")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    rows = record["scientific_validation"]["reported_parameterizations"]

    assert len(rows) == 9
    assert len({row["structure"] for row in rows}) == 9
    assert all(set(row) == {"structure", "glazer", "bm3", "bm4"} for row in rows)
    assert all("K0_double_prime_gpa_inverse" in row["bm4"] for row in rows)
    assert next(row for row in rows if row["structure"] == "P42/nmc")["bm4"] == {
        "V0_a3_per_formula_unit": 44.566,
        "K0_gpa": 251.0,
        "K0_prime": 3.977,
        "K0_double_prime_gpa_inverse": -0.001,
    }

    for row in rows:
        bm3 = row["bm3"]
        implied = (
            -((bm3["K0_prime"] - 4.0) * (bm3["K0_prime"] - 3.0) + 35.0 / 9.0)
            / bm3["K0_gpa"]
        )
        assert bm3["implied_K0_double_prime_gpa_inverse"] == pytest.approx(
            implied, abs=5.0e-10
        )


def test_caracas_2005_audit_regeneration_preserves_parameterization_resolution():
    document = get_material_document("ca_perovskite")
    source = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )

    regenerated = audit_record(source, "ca_perovskite.eosmat")

    assert (
        regenerated["scientific_validation"]["reported_parameterizations"]
        == (source["scientific_validation"]["reported_parameterizations"])
    )
    assert (
        regenerated["scientific_validation"]["parameterization_resolution"]
        == (source["scientific_validation"]["parameterization_resolution"])
    )
    assert regenerated["scientific_validation"]["audit_date"] == "2026-09-04"


def test_caracas_2005_bm3_reproduces_published_130_gpa_density():
    material = Material.from_eosmat(
        get_material_document("ca_perovskite"), record_identifiers=[RECORD_ID]
    )
    record = material.get_eos_record(RECORD_ID)
    volume = record.volume(130.0)

    formula_mass_g_mol = 40.078 + 28.0855 + 3.0 * 15.999
    avogadro = 6.02214076e23
    density_g_cm3 = formula_mass_g_mol / (volume * 1.0e-24 * avogadro)

    assert volume == pytest.approx(33.4363170991, rel=1.0e-10)
    assert density_g_cm3 == pytest.approx(5.77, abs=0.01)
    assert record.pressure(volume) == pytest.approx(130.0, rel=1.0e-11)
    assert math.isfinite(record.eos.bulk_modulus(volume))
