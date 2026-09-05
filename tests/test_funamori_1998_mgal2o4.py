import csv
import hashlib
import math
from pathlib import Path

import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_funamori_1998_mgal2o4 import reproduce

ROOT = Path(__file__).parents[1]
DOI = "10.1029/98JB01575"
CASES = {
    "mgal2o4_cafe2o4": {
        "record": "mgal2o4_cafe2o4_funamori_1998_bm2_1",
        "dataset": "mgal2o4_cafe2o4_funamori_1998_text_pv",
        "resource": "mgal2o4-cafe2o4-funamori-1998-text-pv.csv",
        "sha256": "72c4d48219e5be683fdf507178155dc9bb2fba5313ef40c71ce2bc021a69826d",
        "space_group": "Pbnm",
        "parameters": {"V0": 240.3, "K0": 211.0},
        "errors": {"V0": 0.2, "K0": 6.0},
        "high_pressure": 33.4,
        "high_volume": 212.406277266,
        "diagnostic": "cafe2o4_type",
        "refit_k0": 211.436871484,
    },
    "mgal2o4_cati2o4": {
        "record": "mgal2o4_cati2o4_funamori_1998_bm2_1",
        "dataset": "mgal2o4_cati2o4_funamori_1998_text_pv",
        "resource": "mgal2o4-cati2o4-funamori-1998-text-pv.csv",
        "sha256": "a0a299230ab34aae5acea5ed50d917d363637d6dacd1fc4b525a5a3b32de6257",
        "space_group": "Cmcm",
        "parameters": {"V0": 240.3, "K0": 206.0},
        "errors": {"V0": 0.4, "K0": 3.0},
        "high_pressure": 64.3,
        "high_volume": 195.540421098,
        "diagnostic": "cati2o4_type",
        "refit_k0": 206.401375420,
    },
}


def _rows(resource):
    path = ROOT / "peritheos" / "data" / "datasets" / resource
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_funamori_phase_identity_and_published_bm2(material_identifier, case):
    document = get_material_document(material_identifier)
    assert document["formula"] == "MgAl2O4"
    assert document["formula_units_per_cell"] == 4
    assert document["space_group"] == case["space_group"]
    assert document["atom_sites"] == []
    assert len(document["eos_records"]) == 1

    stored = document["eos_records"][0]
    assert stored["identifier"] == case["record"]
    assert stored["reference"]["doi"].lower() == DOI.lower()
    assert stored["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": case["parameters"],
    }
    assert stored["parameter_errors"] == case["errors"]
    assert stored["scientific_validation"]["status"] == "primary_source_validated"

    executable = Material.from_eosmat(document).get_eos_record(case["record"])
    assert executable.pressure(case["parameters"]["V0"]) == pytest.approx(0.0)


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_funamori_text_data_transcription_and_checksum(material_identifier, case):
    dataset = get_material_document(material_identifier)["datasets"][0]
    assert dataset["identifier"] == case["dataset"]
    assert dataset["used_by_eos_records"] == [case["record"]]
    path = ROOT / "peritheos" / "data" / "datasets" / case["resource"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == case["sha256"]
    assert dataset["resource"]["sha256"] == case["sha256"]

    rows = _rows(case["resource"])
    assert len(rows) == 2
    assert rows[0]["observation_kind"] == "recovered_ambient"
    assert rows[1]["observation_kind"] == "in_situ_after_heating"
    assert float(rows[1]["pressure_gpa"]) == case["high_pressure"]
    lattice_volume = math.prod(
        float(rows[1][axis]) for axis in ("a_angstrom", "b_angstrom", "c_angstrom")
    )
    assert lattice_volume == pytest.approx(case["high_volume"], abs=5.0e-10)


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_funamori_published_curve_and_fixed_v0_refit_have_parity(
    material_identifier, case
):
    diagnostic = reproduce()[case["diagnostic"]]
    assert diagnostic["fixed_v0_refit_k0_gpa"] == pytest.approx(
        case["refit_k0"], abs=1.0e-9
    )
    assert (
        abs(diagnostic["fixed_v0_refit_k0_gpa"] - case["parameters"]["K0"])
        < case["errors"]["K0"]
    )

    document = get_material_document(material_identifier)
    stored = document["eos_records"][0]
    executable = Material.from_eosmat(document).get_eos_record(case["record"])
    assert executable.pressure(case["high_volume"]) == pytest.approx(
        diagnostic["published_curve_pressure_gpa"], abs=1.0e-10
    )
    stored_refit = stored["scientific_validation"]["independent_refit"]["parameters"]
    assert stored_refit["K0"] == pytest.approx(case["refit_k0"], abs=1.0e-8)


def test_funamori_audit_disposes_all_same_doi_candidates_once():
    audit = (
        ROOT
        / "docs"
        / "literature-reproductions"
        / "funamori-1998-mgal2o4-transformations.md"
    ).read_text(encoding="utf-8")
    for candidate in (
        "litcurate_ca9a93f28c82d5bc",
        "litcurate_45e8795108756c17",
        "litcurate_cc1efc12c6e1c6e0",
    ):
        assert audit.count(candidate) == 1
