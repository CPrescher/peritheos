import json
from pathlib import Path

import pytest

from peritheos import Material, get_material_document, validate_eosmat_document
from scripts.reproduce_metsue_2012_bridgmanite import reproduce

ROOT = Path(__file__).parents[1]
DOI = "10.1111/j.1365-246x.2012.05511.x"
FE_MATERIAL = "mg09375fe00625sio3_bridgmanite"


def test_metsue_fe_material_preserves_six_spin_configuration_fits():
    document = get_material_document(FE_MATERIAL)
    validate_eosmat_document(document)
    assert document["formula"] == "Mg0.9375Fe0.0625SiO3"
    assert document["formula_units_per_cell"] == 4
    records = document["eos_records"]
    assert len(records) == 6
    assert sum(record.get("default_for") == "equilibrium" for record in records) == 1
    assert sum("high-spin" in record["label"] for record in records) == 3
    assert sum("low-spin" in record["label"] for record in records) == 3
    for record in records:
        assert record["fixed_parameters"] == ["K0_prime"]
        assert record["eos"]["parameters"]["K0_prime"] == 3.94
        assert record["reference"]["doi"].lower() == DOI
        assert record["scientific_validation"]["status"] == ("primary_source_validated")


def test_metsue_pure_control_is_a_separate_nondefault_source_record():
    document = get_material_document("bridgmanite")
    record = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == "bridgmanite_metsue_2012_static_bm3_1"
    )
    assert "default_for" not in record
    assert record["eos"]["parameters"] == {
        "V0": pytest.approx(161.216424307, abs=1.0e-9),
        "K0": 258.7857,
        "K0_prime": 3.94,
    }
    assert record["fixed_parameters"] == ["K0_prime"]
    assert record["reference"]["doi"].lower() == DOI


def test_metsue_all_seven_records_execute_and_anchor_zero_pressure():
    documents = [
        get_material_document(FE_MATERIAL),
        get_material_document("bridgmanite"),
    ]
    identifiers = {
        record["identifier"]
        for document in documents
        for record in document["eos_records"]
        if record["reference"]["doi"].lower() == DOI
    }
    assert len(identifiers) == 7
    for document in documents:
        material = Material.from_eosmat(document)
        for stored in document["eos_records"]:
            if stored["identifier"] not in identifiers:
                continue
            parameters = stored["eos"]["parameters"]
            eos = material.get_eos_record(stored["identifier"])
            assert eos.pressure(parameters["V0"]) == pytest.approx(0.0, abs=1.0e-12)


def test_metsue_reproduction_confirms_conversion_and_small_model_spread():
    result = reproduce()
    assert result["molar_to_z4_a3"] == pytest.approx(6.642156268695387)
    assert result["records"]["HS_model_1"]["V0_a3_z4"] == pytest.approx(
        161.620267408, abs=1.0e-9
    )
    assert result["records"]["pure_MgSiO3"]["V0_a3_z4"] == pytest.approx(
        161.216424307, abs=1.0e-9
    )
    for group in result["configuration_spread"].values():
        assert group["maximum_relative_volume_spread_percent"] < 0.01
        assert group["maximum_relative_modulus_spread_percent"] < 0.07
    assert result["model1_composition_trends"]["HS"]["dlnK0_dXFe"] == pytest.approx(
        0.035, abs=0.002
    )
    assert result["model1_composition_trends"]["LS"]["dlnK0_dXFe"] == pytest.approx(
        0.123, abs=0.003
    )


def test_metsue_audit_disposes_all_eleven_litcurate_rows_once():
    audit = (
        ROOT / "docs" / "literature-reproductions" / "metsue-2012-bridgmanite.md"
    ).read_text(encoding="utf-8")
    candidates = json.loads(
        (ROOT / "docs" / "data" / "litcurate-eos-candidates.json").read_text(
            encoding="utf-8"
        )
    )["records"]
    same_doi = [row for row in candidates if row["publication"]["doi"].lower() == DOI]
    assert len(same_doi) == 11
    for row in same_doi:
        assert audit.count(row["identifier"]) == 1
