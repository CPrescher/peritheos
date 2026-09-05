import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import pytest

from peritheos import Material, load_eosmat
from scripts.reproduce_fu_2024_bridgmanite import reproduce_model

ROOT = Path(__file__).resolve().parents[1]
MATERIAL_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "materials"
    / "mg088fe010al014si090o3_bridgmanite.eosmat"
)
DATA_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "mg088fe010al014si090o3-bridgmanite-fu-2024-table-s1-cif-pv.csv"
)
DATASET_ID = "mg088fe010al014si090o3_bridgmanite_fu_2024_table_s1_cif_pv"
RECORD_IDS = (
    "mg088fe010al014si090o3_bridgmanite_fu_2024_bm2_1",
    "mg088fe010al014si090o3_bridgmanite_fu_2024_bm3_2",
)
CHECKSUM = "c735c80ede486b5ab85b755f1ab450a56e20864b41f7adeb2e2798e818090fff"


def load_source() -> tuple[dict, Material]:
    document = load_eosmat(MATERIAL_PATH)
    return document, Material.from_eosmat(document)


def test_fu_material_identity_structure_and_normalized_cell_contents():
    document, _ = load_source()

    assert document["identifier"] == "mg088fe010al014si090o3_bridgmanite"
    assert document["formula"] == "Mg0.88Fe0.10Al0.14Si0.90O3"
    assert document["phase"] == "Pbnm (Al,Fe)-bearing bridgmanite, synthesis run 5K2667"
    assert document["space_group"] == "Pbnm"
    assert document["space_group_number"] == 62
    assert document["formula_units_per_cell"] == 4
    assert document["lattice"] == {
        "a": 4.7564,
        "b": 4.9266,
        "c": 6.9005,
        "alpha": 90.0,
        "beta": 90.0,
        "gamma": 90.0,
    }

    multiplicities = {"4b": 4, "4c": 4, "8d": 8}
    contents: dict[str, float] = defaultdict(float)
    for site in document["atom_sites"]:
        contents[site["element"]] += (
            multiplicities[site["wyckoff"]]
            * site["occupancy"]
            / document["formula_units_per_cell"]
        )
    assert contents == pytest.approx(
        {"Mg": 0.88, "Fe": 0.09, "Al": 0.13, "Si": 0.90, "O": 3.0}
    )
    assert "measured Mg0.88Fe0.10Al0.14Si0.90O3 composition" in document["notes"]
    assert "normalized occupancies give Mg0.88Fe0.09Al0.13Si0.90O3" in document["notes"]


@pytest.mark.parametrize(
    (
        "index",
        "model",
        "parameters",
        "errors",
        "fixed",
        "benchmark_volume",
        "benchmark",
    ),
    (
        (
            0,
            "birch_murnaghan_2",
            {"V0": 163.85, "K0": 242.0},
            {"V0": 0.07, "K0": 3.0},
            [],
            136.29,
            64.4560123693,
        ),
        (
            1,
            "birch_murnaghan_3",
            {"V0": 164.64, "K0": 228.0, "K0_prime": 4.1},
            {"V0": 0.11, "K0": 5.0, "K0_prime": 0.2},
            [],
            161.70,
            4.2628529138,
        ),
    ),
)
def test_fu_records_load_execute_and_invert(
    index, model, parameters, errors, fixed, benchmark_volume, benchmark
):
    document, material = load_source()
    source = document["eos_records"][index]
    record = material.eos_records[index]

    assert source["identifier"] == RECORD_IDS[index]
    assert source["reference"]["doi"] == "10.2138/am-2023-8969"
    assert source["eos"]["model"] == model
    assert source["eos"]["parameters"] == parameters
    assert source["parameter_errors"] == errors
    assert source["parameter_error_confidence"] is None
    assert source["parameter_covariance"] is None
    assert source["fixed_parameters"] == fixed
    assert source["temperature_ref"] == 293.0
    assert source["experimental_pressure_range_gpa"] == [4.2, 64.6]
    assert source["fit_datasets"] == [DATASET_ID]
    assert source["author_fit_preference"].startswith("none;")
    assert source["scientific_validation"]["status"] == "primary_source_validated"

    assert record.pressure(parameters["V0"], 293.0) == pytest.approx(0.0, abs=1.0e-13)
    pressure = record.pressure(benchmark_volume, 293.0, check_validity=True)
    assert pressure == pytest.approx(benchmark, abs=1.0e-10)
    assert record.volume(pressure, 293.0, check_validity=True) == pytest.approx(
        benchmark_volume, abs=1.0e-9
    )


def test_fu_primary_data_transcription_checksum_and_source_disagreements():
    document = json.loads(MATERIAL_PATH.read_text(encoding="utf-8"))
    dataset = document["datasets"][0]

    assert dataset["identifier"] == DATASET_ID
    assert dataset["used_by_eos_records"] == list(RECORD_IDS)
    assert dataset["resource"]["sha256"] == CHECKSUM
    assert hashlib.sha256(DATA_PATH.read_bytes()).hexdigest() == CHECKSUM

    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 22
    assert [float(rows[0]["pressure_gpa"]), float(rows[-1]["pressure_gpa"])] == [
        4.2,
        64.6,
    ]
    volumes = [
        float(rows[0]["cif_volume_a3_conventional_cell"]),
        float(rows[-1]["cif_volume_a3_conventional_cell"]),
    ]
    assert volumes == [161.70, 136.29]
    assert {row["temperature_k"] for row in rows} == {"293"}
    assert rows[1]["table_s1_b_angstrom"] == "4.9098"
    assert rows[1]["cif_b_angstrom"] == "4.9097"
    assert rows[12]["table_s1_b_sd_angstrom"] == "0.0040"
    assert rows[12]["cif_b_sd_angstrom"] == "0.0010"
    assert rows[13]["table_s1_a_sd_angstrom"] == "0.0047"
    assert rows[13]["cif_a_sd_angstrom"] == "0.0015"
    assert rows[16]["table_s1_b_sd_angstrom"] == "0.0045"
    assert rows[16]["cif_b_sd_angstrom"] == "0.0012"


@pytest.mark.parametrize(
    ("model", "rmse", "max_residual", "benchmark_key", "benchmark_pressure"),
    (
        ("BM2", 0.7195012696, 1.4607777344, "high_pressure_state", 64.4560123693),
        ("BM3", 0.9290848204, 1.9358848929, "low_pressure_state", 4.2628529138),
    ),
)
def test_fu_independent_equation_and_diagnostic_refits(
    model, rmse, max_residual, benchmark_key, benchmark_pressure
):
    result = reproduce_model(model)

    assert result["rows"] == 22
    assert result["published_curve_pressure_rmse_gpa"] == pytest.approx(
        rmse, abs=1.0e-10
    )
    assert result["published_curve_max_abs_pressure_residual_gpa"] == pytest.approx(
        max_residual, abs=1.0e-10
    )
    assert result[benchmark_key]["calculated_pressure_gpa"] == pytest.approx(
        benchmark_pressure, abs=1.0e-10
    )
    assert result["unweighted_pressure_refit"]["rmse_gpa"] < rmse
    assert result["pressure_sigma_weighted_refit"]["rmse_gpa"] < rmse
    assert result["pressure_sigma_weighted_refit"]["reduced_chi_square"] > 0.0


def test_fu_pressure_calibration_resolves_to_bundled_gold_record():
    document, _ = load_source()
    gold = load_eosmat(MATERIAL_PATH.with_name("gold.eosmat"))
    gold_record_ids = {record["identifier"] for record in gold["eos_records"]}
    for source in document["eos_records"]:
        calibration = source["pressure_calibration"]
        method = calibration["methods"][0]
        assert calibration["status"] == "resolved"
        assert method["reference"]["doi"] == "10.1073/pnas.0609013104"
        assert method["reference_eos_record"] == "gold_fei_2007_vinet_2"
        assert method["reference_eos_record"] in gold_record_ids
        assert calibration["recalculation"]["status"] == (
            "missing_calibrant_observations"
        )
