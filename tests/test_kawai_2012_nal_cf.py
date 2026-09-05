import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, load_eosmat
from scripts.reproduce_kawai_2012_nal_cf import MOLAR_A3, reproduce_phase

ROOT = Path(__file__).resolve().parents[1]
MATERIAL_ROOT = ROOT / "peritheos" / "data" / "materials"
DATA_ROOT = ROOT / "peritheos" / "data"
DOI = "10.2138/am.2012.3915"
CASES = {
    "nal": {
        "material": "namg2al5sio12_nal",
        "record": "namg2al5sio12_nal_kawai_2012_bm3_1",
        "dataset": "namg2al5sio12_nal_kawai_2012_figure2c_vector_digitized",
        "csv": "datasets/namg2al5sio12-nal-kawai-2012-figure2c-vector-digitized.csv",
        "sha256": "25507f00d6de30eab57a3b3de1af1763feca09bf09f2287d676d41e489243bbd",
        "space_group": "P63/m",
        "number": 176,
        "z": 1.0,
        "rows": 6,
        "pressure_range": [0.0, 50.0],
        "parameters": {"V0": 35.8 * 3.0 / MOLAR_A3, "K0": 217.7, "K0_prime": 4.08},
        "proxy_formula": {"Ca": 1.0, "Mg": 2.0, "Al": 6.0, "O": 12.0},
        "rmse": 0.100424430486,
        "maximum": 0.112884772937,
        "refit": [178.233853876982, 217.562095743266, 4.088521901863],
        "refit_rmse": 0.003452426412,
    },
    "cf": {
        "material": "namg2al5sio12_cf",
        "record": "namg2al5sio12_cf_kawai_2012_bm3_1",
        "dataset": "namg2al5sio12_cf_kawai_2012_figure2d_vector_digitized",
        "csv": "datasets/namg2al5sio12-cf-kawai-2012-figure2d-vector-digitized.csv",
        "sha256": "76e05b388ce7ad11c73f7a216e5bf8ce009d4d47e3f94409d75489b80a1ffa41",
        "space_group": "Pbnm",
        "number": 62,
        "z": 4.0 / 3.0,
        "rows": 8,
        "pressure_range": [0.0, 150.0],
        "parameters": {"V0": 35.2 * 4.0 / MOLAR_A3, "K0": 213.2, "K0_prime": 4.12},
        "proxy_formula": {"Na": 4.0, "Al": 4.0, "Si": 4.0, "O": 16.0},
        "rmse": 0.151787368860,
        "maximum": 0.198344317864,
        "refit": [233.944765268382, 213.205584668798, 4.132190010213],
        "refit_rmse": 0.030775628572,
    },
}


def _document(case):
    return load_eosmat(MATERIAL_ROOT / f"{case['material']}.eosmat")


@pytest.mark.parametrize("phase", CASES)
def test_kawai_records_are_composition_phase_specific_and_executable(phase):
    case = CASES[phase]
    document = _document(case)
    source = document["eos_records"][0]
    record = Material.from_eosmat(document).eos_records[0]

    assert document["formula"] == "NaMg2Al5SiO12"
    assert document["space_group"] == case["space_group"]
    assert document["space_group_number"] == case["number"]
    assert document["formula_units_per_cell"] == pytest.approx(case["z"])
    assert source["identifier"] == case["record"]
    assert source["reference"]["doi"] == DOI
    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": case["parameters"],
    }
    assert source["temperature_ref"] == 0.0
    assert source["experimental_pressure_range_gpa"] == case["pressure_range"]
    assert source["pressure_range_status"] == "theoretical"
    assert source["parameter_errors"] == {"V0": None, "K0": None, "K0_prime": None}
    assert source["fixed_parameters"] == []
    assert source["pressure_calibration"]["status"] == "not_applicable"
    assert source["pressure_calibration"]["recalculation"]["status"] == "not_applicable"

    v0 = case["parameters"]["V0"]
    assert record.pressure(v0, 0.0) == pytest.approx(0.0, abs=1.0e-13)
    pressure = np.linspace(
        case["pressure_range"][0] + 0.1,
        case["pressure_range"][1] - 0.1,
        9,
    )
    volume = record.volume(pressure, 0.0, check_validity=True)
    assert record.pressure(volume, 0.0, check_validity=True) == pytest.approx(
        pressure, rel=1.0e-10
    )


@pytest.mark.parametrize("phase", CASES)
def test_kawai_structure_models_are_explicit_unmodified_proxies(phase):
    case = CASES[phase]
    document = _document(case)
    contents = {}
    for site in document["atom_sites"]:
        assert "topology proxy" in site["source_label"]
        contents[site["element"]] = contents.get(site["element"], 0.0) + (
            site["site_multiplicity"] * site["occupancy"]
        )
    assert contents == case["proxy_formula"]
    assert (
        "No target occupancies" in document["notes"]
        or "No Mg substitution" in document["notes"]
    )
    assert "not NaMg2Al5SiO12 coordinates" in document["notes"]


@pytest.mark.parametrize("phase", CASES)
def test_kawai_vector_digitizations_are_checksummed_and_volume_consistent(phase):
    case = CASES[phase]
    document = _document(case)
    dataset = document["datasets"][0]
    path = DATA_ROOT / case["csv"]
    assert dataset["identifier"] == case["dataset"]
    assert dataset["used_by_eos_records"] == [case["record"]]
    assert dataset["provenance"]["type"] == "digitized_from_figure"
    assert dataset["provenance"]["source_pdf_sha256"] == (
        "8273dc07c0ba3dca62724d490e035e945ee6c4123d7c6319f140afc44c7ec089"
    )
    assert hashlib.sha256(path.read_bytes()).hexdigest() == case["sha256"]
    assert dataset["resource"]["sha256"] == case["sha256"]

    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == case["rows"]
    assert {row["used_in_published_fit"] for row in rows} == {"1"}
    assert [float(row["pressure_gpa"]) for row in rows] == sorted(
        float(row["pressure_gpa"]) for row in rows
    )
    for row in rows:
        full_formula = float(row["molar_volume_cm3_per_mol_namg2al5sio12"])
        o4 = float(row["molar_volume_cm3_per_mol_o4"])
        cell = float(row["volume_a3_conventional_cell"])
        assert full_formula == pytest.approx(3.0 * o4, abs=2.0e-9)
        assert cell == pytest.approx(
            o4 * (3.0 if phase == "nal" else 4.0) / MOLAR_A3,
            abs=5.0e-9,
        )


@pytest.mark.parametrize("phase", CASES)
def test_kawai_digitized_curves_reproduce_published_bm3(phase):
    case = CASES[phase]
    result = reproduce_phase(phase)
    assert result["rows"] == case["rows"]
    assert result["published_curve_volume_rmse_a3"] == pytest.approx(
        case["rmse"], abs=5.0e-12
    )
    assert result["published_curve_max_abs_volume_residual_a3"] == pytest.approx(
        case["maximum"], abs=5.0e-12
    )
    assert result["volume_residual_refit"] == pytest.approx(case["refit"], abs=5.0e-6)
    assert result["volume_residual_refit_rmse_a3"] == pytest.approx(
        case["refit_rmse"], abs=5.0e-12
    )


def test_kawai_audit_disposes_every_same_doi_litcurate_candidate():
    path = ROOT / "docs" / "literature-reproductions" / "kawai-2012-nal-cf.md"
    content = path.read_text(encoding="utf-8")
    identifiers = {
        "litcurate_c42ee5edc21ddfe5",
        "litcurate_3b528f9e781f555b",
        "litcurate_d289c4c72b25a892",
        "litcurate_b8fa7310ff9994db",
        "litcurate_f6c0b448f52c70d9",
        "litcurate_527df255bad7681e",
        "litcurate_dba4bb2c50bd7a50",
        "litcurate_74a5d6778933b7d6",
        "litcurate_6217113a1e9daf50",
    }
    assert all(identifier in content for identifier in identifiers)
    assert content.count("**ACCEPT — production record**") == 2
    assert content.count("**REJECT as a Kawai") == 7
    assert "Final production count" in content and "**2 records**" in content
