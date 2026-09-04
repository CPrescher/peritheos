import csv
import hashlib
import json
import math
from pathlib import Path

import pytest

from peritheos import Material, load_eosmat
from scripts.reproduce_qin_2023_naalsio4 import reproduce_sample

ROOT = Path(__file__).resolve().parents[1]
MATERIAL_ROOT = ROOT / "peritheos" / "data" / "materials"
DATA_ROOT = ROOT / "peritheos" / "data" / "datasets"

CASES = (
    {
        "material": "na093al102si100o4_calcium_ferrite",
        "record": "na093al102si100o4_calcium_ferrite_qin_2023_bm3_1",
        "dataset": "na093al102si100o4_qin_2023_table_s3_pv",
        "csv": "na093al102si100o4-qin-2023-table-s3-pv.csv",
        "sha256": "50acd8887e5b9215ad0e47528f81c8fba4db6a8bc1e1a5276dc1169c36ad4a0a",
        "formula": "Na0.93Al1.02Si1.00O4",
        "parameters": {"V0": 241.6, "K0": 220.0, "K0_prime": 2.6},
        "errors": {"V0": 0.1, "K0": 4.0, "K0_prime": 0.3},
        "rows": 22,
        "pressure_range": [0.0001, 41.0],
        "last_volume": 207.3,
        "last_calculated_pressure": 40.60961206892379,
        "rmse": 0.32627826994900294,
        "max_residual": 0.6438583996006795,
        "unweighted_parameters": [
            242.01380086577035,
            209.58508880982552,
            3.103588024098269,
        ],
        "weighted_parameters": [
            242.0273016839734,
            208.87612577656526,
            3.136837395053564,
        ],
        "weighted_reduced_chi_square": 1.1189465362793038,
    },
    {
        "material": "na088al099fe013si094o4_calcium_ferrite",
        "record": "na088al099fe013si094o4_calcium_ferrite_qin_2023_bm3_1",
        "dataset": "na088al099fe013si094o4_qin_2023_table_s4_pv",
        "csv": "na088al099fe013si094o4-qin-2023-table-s4-pv.csv",
        "sha256": "44ed46a6c00a7ea5715af87fcf1e55c5ae1c80acaecf3230c5205d8754e6ab0c",
        "formula": "Na0.88Al0.99Fe0.13Si0.94O4",
        "parameters": {"V0": 244.2, "K0": 211.0, "K0_prime": 2.6},
        "errors": {"V0": 0.2, "K0": 6.0, "K0_prime": 0.3},
        "rows": 10,
        "pressure_range": [0.0001, 44.0],
        "last_volume": 206.6,
        "last_calculated_pressure": 43.2093788098509,
        "rmse": 0.3887326932802591,
        "max_residual": 0.7906211901491034,
        "unweighted_parameters": [
            244.32647528735754,
            205.7755681835045,
            2.9166477472839603,
        ],
        "weighted_parameters": [
            244.05732254870946,
            201.29958178891874,
            3.107153595241988,
        ],
        "weighted_reduced_chi_square": 0.536114208923388,
    },
)


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["material"])
def test_qin_material_records_load_execute_and_invert(case):
    document = load_eosmat(MATERIAL_ROOT / f"{case['material']}.eosmat")
    material = Material.from_eosmat(document)
    source = document["eos_records"][0]
    record = material.eos_records[0]

    assert document["formula"] == case["formula"]
    assert document["space_group"] == "Pbnm"
    assert document["space_group_number"] == 62
    assert document["formula_units_per_cell"] == 4
    assert source["identifier"] == case["record"]
    assert source["reference"]["doi"] == "10.2138/am-2022-8432"
    assert source["eos"]["model"] == "birch_murnaghan_3"
    assert source["eos"]["parameters"] == case["parameters"]
    assert source["parameter_errors"] == case["errors"]
    assert source["fixed_parameters"] == []
    assert source["parameter_covariance"] is None
    assert source["parameter_error_confidence"] is None
    assert source["temperature_ref"] == 293.0
    assert source["experimental_pressure_range_gpa"] == case["pressure_range"]
    assert source["fit_datasets"] == [case["dataset"]]
    assert source["scientific_validation"]["status"] == "primary_source_validated"

    assert record.pressure(case["parameters"]["V0"], 293.0) == pytest.approx(
        0.0, abs=1.0e-13
    )
    pressure = record.pressure(case["last_volume"], 293.0, check_validity=True)
    assert pressure == pytest.approx(case["last_calculated_pressure"], abs=1.0e-10)
    assert record.volume(pressure, 293.0, check_validity=True) == pytest.approx(
        case["last_volume"], abs=1.0e-9
    )


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["dataset"])
def test_qin_primary_table_transcription_and_checksum(case):
    document = json.loads(
        (MATERIAL_ROOT / f"{case['material']}.eosmat").read_text(encoding="utf-8")
    )
    dataset = document["datasets"][0]
    path = DATA_ROOT / case["csv"]

    assert dataset["identifier"] == case["dataset"]
    assert dataset["used_by_eos_records"] == [case["record"]]
    assert dataset["resource"]["path"] == f"datasets/{case['csv']}"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == case["sha256"]
    assert dataset["resource"]["sha256"] == case["sha256"]

    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == case["rows"]
    assert float(rows[0]["pressure_gpa"]) == pytest.approx(0.0001)
    assert float(rows[-1]["pressure_gpa"]) == pytest.approx(case["pressure_range"][1])
    assert float(rows[-1]["volume_a3_conventional_cell"]) == pytest.approx(
        case["last_volume"]
    )
    assert {row["used_in_published_fit"] for row in rows} == {"1"}
    assert sum(int(row["source_note_a"]) for row in rows) == 1


def test_qin_fe_free_source_inconsistencies_are_preserved():
    document = json.loads(
        (MATERIAL_ROOT / "na093al102si100o4_calcium_ferrite.eosmat").read_text(
            encoding="utf-8"
        )
    )
    with (DATA_ROOT / "na093al102si100o4-qin-2023-table-s3-pv.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        rows = list(csv.DictReader(stream))
    anomalous = next(row for row in rows if float(row["pressure_gpa"]) == 38.7)

    lattice_volume = math.prod(
        float(anomalous[name]) for name in ("a_angstrom", "b_angstrom", "c_angstrom")
    )
    assert lattice_volume == pytest.approx(210.1498988016)
    assert float(anomalous["volume_a3_conventional_cell"]) == pytest.approx(208.8)
    inconsistencies = document["eos_records"][0]["scientific_validation"][
        "reported_inconsistencies"
    ]
    assert len(inconsistencies) == 3
    assert any("38.7 GPa" in item["field"] for item in inconsistencies)
    assert any("ambient volume" in item["field"] for item in inconsistencies)
    assert any("superscript a" in item["field"] for item in inconsistencies)


def test_qin_fe_bearing_structure_proxy_is_explicit_and_does_not_infer_fe_sites():
    document = json.loads(
        (MATERIAL_ROOT / "na088al099fe013si094o4_calcium_ferrite.eosmat").read_text(
            encoding="utf-8"
        )
    )

    assert {site["element"] for site in document["atom_sites"]} == {
        "Na",
        "Al",
        "Si",
        "O",
    }
    assert "isostructural topology proxy" in document["notes"]
    assert "not an exact refinement" in document["notes"]
    assert (
        "site-specific Fe occupancies are not given and are not inferred"
        in document["notes"]
    )


@pytest.mark.parametrize(
    ("sample", "case"),
    (("fe_free", CASES[0]), ("fe_bearing", CASES[1])),
)
def test_qin_independent_bm3_regression(sample, case):
    result = reproduce_sample(sample)

    assert result["rows"] == case["rows"]
    assert result["published_curve_pressure_rmse_gpa"] == pytest.approx(
        case["rmse"], abs=1.0e-12
    )
    assert result["published_curve_max_abs_pressure_residual_gpa"] == pytest.approx(
        case["max_residual"], abs=1.0e-12
    )
    high_pressure = result["high_pressure_state"]
    assert high_pressure["volume_a3"] == pytest.approx(case["last_volume"])
    assert high_pressure["calculated_pressure_gpa"] == pytest.approx(
        case["last_calculated_pressure"], abs=1.0e-12
    )

    unweighted = result["unweighted_pressure_refit"]
    weighted = result["volume_uncertainty_weighted_refit"]
    assert unweighted["parameters"] == pytest.approx(
        case["unweighted_parameters"], rel=1.0e-5, abs=1.0e-8
    )
    assert weighted["parameters"] == pytest.approx(
        case["weighted_parameters"], rel=1.0e-5, abs=1.0e-8
    )
    assert unweighted["rmse_gpa"] < result["published_curve_pressure_rmse_gpa"]
    assert weighted["reduced_chi_square"] == pytest.approx(
        case["weighted_reduced_chi_square"], rel=1.0e-5
    )
