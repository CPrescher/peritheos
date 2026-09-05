import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document, validate_eosmat_document
from scripts.reproduce_oganov_wang_chizmeshya_eos import reproduce

ROOT = Path(__file__).parents[1]
OGANOV_DOI = "10.1103/physrevb.67.224110"
WANG_DOI = "10.1029/95jb03254"
CHIZMESHYA_DOI = "10.1029/96gl02624"
KARKI_DOI = "10.1029/1999jb900069"


def _doi_records(document: dict, doi: str) -> list[dict]:
    return [
        record
        for record in document["eos_records"]
        if record["reference"].get("doi", "").lower() == doi
    ]


def test_oganov_preserves_four_static_and_seven_corrected_bm3_fits():
    document = get_material_document("mgo")
    validate_eosmat_document(document)
    records = _doi_records(document, OGANOV_DOI)
    assert len(records) == 11
    assert sum("pressure_corrected" in record["identifier"] for record in records) == 7
    observed = {
        record["identifier"]: tuple(record["eos"]["parameters"].values())
        for record in records
    }
    assert observed == {
        "mgo_oganov_2003_ecp_large_core_static_bm3": (77.629, 151.707, 4.212),
        "mgo_oganov_2003_ecp_small_core_static_bm3": (76.595, 150.839, 4.052),
        "mgo_oganov_2003_paw_large_core_static_bm3": (76.049, 154.183, 4.141),
        "mgo_oganov_2003_paw_small_core_static_bm3": (76.947, 150.597, 4.103),
        "mgo_oganov_2003_pressure_corrected_static_bm3": (73.425, 181.24, 3.997),
        "mgo_oganov_2003_pressure_corrected_0k_bm3": (74.439, 173.48, 4.014),
        "mgo_oganov_2003_pressure_corrected_298k_bm3": (74.67, 170.53, 4.036),
        "mgo_oganov_2003_pressure_corrected_1000k_bm3": (76.549, 152.595, 4.13),
        "mgo_oganov_2003_pressure_corrected_2000k_bm3": (79.915, 127.719, 4.244),
        "mgo_oganov_2003_pressure_corrected_3000k_bm3": (83.772, 106.11, 4.331),
        "mgo_oganov_2003_pressure_corrected_4000k_bm3": (88.006, 88.473, 4.385),
    }


def test_wang_and_chizmeshya_preserve_only_source_owned_complete_fits():
    document = get_material_document("ca_perovskite")
    validate_eosmat_document(document)
    wang = _doi_records(document, WANG_DOI)
    chizmeshya = _doi_records(document, CHIZMESHYA_DOI)
    assert len(wang) == 4
    assert len(chizmeshya) == 5
    assert [tuple(record["eos"]["parameters"].values()) for record in wang] == [
        (45.58, 232.0, 4.8),
        (45.71, 244.0, 4.8),
        (45.3, 282.0, 4.0),
        (45.47, 268.0, 4.3),
    ]
    assert [tuple(record["eos"]["parameters"].values()) for record in chizmeshya] == [
        (45.04, 241.8, 4.15),
        (45.02, 241.0, 4.16),
        (45.06, 238.2, 4.18),
        (45.55, 237.9, 4.0),
        (45.62, 227.0, 4.29),
    ]
    assert wang[0]["fixed_parameters"] == ["K0_prime"]
    assert chizmeshya[3]["fixed_parameters"] == ["K0_prime"]


def test_all_twenty_accepted_curves_are_executable_and_anchor_zero_pressure():
    for material_name, dois in [
        ("mgo", {OGANOV_DOI}),
        ("ca_perovskite", {WANG_DOI, CHIZMESHYA_DOI}),
    ]:
        document = get_material_document(material_name)
        executable = Material.from_eosmat(document)
        for record in document["eos_records"]:
            if record["reference"]["doi"].lower() not in dois:
                continue
            parameters = record["eos"]["parameters"]
            model = executable.get_eos_record(record["identifier"])
            assert model.pressure(parameters["V0"]) == pytest.approx(0.0, abs=1e-12)
            assert np.isfinite(model.pressure(0.9 * parameters["V0"]))


def test_primary_table_resources_are_complete_and_unchanged():
    resources = {
        "mgo-oganov-2003-table3-pv.csv": (
            8,
            "bd8327b505436092cf2b36d66aff26da0950918659d5c7847781f62ffbb2a213",
        ),
        "ca-perovskite-wang-1996-table1-room-temperature.csv": (
            14,
            "71f8f01e5149e140e66251dbc3ff2c0e2a7d02aede82deae8d1f79cd2823d446",
        ),
    }
    for filename, (expected_rows, expected_sha256) in resources.items():
        path = ROOT / "peritheos" / "data" / "datasets" / filename
        with path.open(newline="", encoding="utf-8") as stream:
            assert len(list(csv.DictReader(stream))) == expected_rows
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_sha256


def test_numerical_reproduction_matches_available_primary_data():
    result = reproduce()
    assert result["oganov"]["record_count"] == 11
    assert result["oganov"]["table3_observations"] == 8
    assert result["oganov"]["paw_small_core_table3_pressure_rmse_gpa"] < 0.15
    assert result["oganov"]["pressure_corrected_static_table3_rmse_gpa"] < 0.15
    assert result["wang"]["room_temperature"]["observations"] == 12
    assert result["wang"]["room_temperature"]["V0"] == pytest.approx(45.58, abs=0.01)
    assert result["wang"]["room_temperature"]["fixed_V0_K0"] == pytest.approx(
        232.0, abs=1.0
    )
    mao = result["wang"]["mao_1989_reanalyses"]
    assert mao["equal_weight_excluding_below_1_gpa"]["V0"] == pytest.approx(
        45.71, abs=0.01
    )
    assert mao["equal_weight_excluding_below_1_gpa"]["K0"] == pytest.approx(
        244.0, abs=0.5
    )
    assert result["chizmeshya"]["record_count"] == 5


def test_karki_elastic_moduli_are_not_stored_as_volume_eos_records():
    found = []
    for path in (ROOT / "peritheos" / "data" / "materials").glob("*.eosmat"):
        document = json.loads(path.read_text(encoding="utf-8"))
        found.extend(_doi_records(document, KARKI_DOI))
    assert found == []


def test_audit_disposes_every_candidate_under_the_four_dois_once():
    audit = (
        ROOT
        / "docs"
        / "literature-reproductions"
        / "oganov-wang-chizmeshya-karki-eos-audit.md"
    ).read_text(encoding="utf-8")
    candidates = json.loads(
        (ROOT / "docs" / "data" / "litcurate-eos-candidates.json").read_text(
            encoding="utf-8"
        )
    )["records"]
    dois = {OGANOV_DOI, WANG_DOI, CHIZMESHYA_DOI, KARKI_DOI}
    same_doi = [row for row in candidates if row["publication"]["doi"].lower() in dois]
    assert len(same_doi) == 44
    assert all(audit.count(row["identifier"]) == 1 for row in same_doi)
