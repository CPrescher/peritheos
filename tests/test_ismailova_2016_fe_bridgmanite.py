import csv
import hashlib
from pathlib import Path

import pytest

from peritheos import Material, load_eosmat
from peritheos.eosmat import validate_pressure_calibration_references
from scripts.reproduce_ismailova_2016_fe_bridgmanite import (
    diagnostic_refit,
    published_curve_diagnostic,
    reproduce,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIAL_PATH = (
    ROOT / "peritheos" / "data" / "materials" / "fe088sio3_bridgmanite.eosmat"
)
DATA_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "fe088sio3-bridgmanite-ismailova-2016-table-s2-selected-crystallography.csv"
)
RECORD_ID = "fe088sio3_bridgmanite_ismailova_2016_300k_bm2"
DATASET_ID = "fe088sio3_bridgmanite_ismailova_2016_table_s2_selected_crystallography"
CHECKSUM = "2fa03e31214c2bbddabcaa67aaff4921ab04d04431c124856fb228f7dad77591"


def _load():
    document = load_eosmat(MATERIAL_PATH)
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, document["eos_records"][0], document["datasets"][0], rows


def test_ismailova_record_has_correct_bm2_and_explicit_fit_limits():
    document, record, _, _ = _load()
    executable = Material.from_eosmat(document).eos_records[0]

    assert record["identifier"] == RECORD_ID
    assert record["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 178.98, "K0": 190.0},
    }
    assert record["implicit_parameters"] == {"K0_prime": 4.0}
    assert record["fit_provenance"]["fixed_parameters"] == ["K0_prime"]
    assert record["fit_provenance"]["objective"] == "not reported"
    assert record["fit_datasets"] == [DATASET_ID]
    assert record["experimental_pressure_range_gpa"] == [14.0, 132.0]
    assert executable.pressure(178.98, 300.0) == pytest.approx(0.0, abs=1.0e-13)


def test_ismailova_table_s2_checkpoint_transcription_is_lossless():
    _, _, dataset, rows = _load()

    assert hashlib.sha256(DATA_PATH.read_bytes()).hexdigest() == CHECKSUM
    assert dataset["resource"]["sha256"] == CHECKSUM
    assert dataset["identifier"] == DATASET_ID
    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert dataset["license"] == "CC BY-NC 4.0"
    assert dataset["provenance"]["official_supporting_pdf_md5"] == (
        "d772465eb207e08a511b987d73a729c9"
    )
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == 4
    assert [row["pressure_source_token"] for row in rows] == [
        "44.0(5)",
        "67.0(5)",
        "107.0(5)",
        "129.0(5)",
    ]
    assert [row["volume_source_token"] for row in rows] == [
        "148.44(18)",
        "142.36(5)",
        "132.17(11)",
        "126.93(13)",
    ]
    assert rows[0]["reported_high_temperature_condition_k"] == ""
    assert rows[-1]["reported_high_temperature_condition_source_token"] == ("1835(100)")


def test_ismailova_pressure_scale_lineage_is_resolved_but_not_recalculable():
    _, record, _, _ = _load()
    calibration = record["pressure_calibration"]
    artifact = record["scientific_validation"]["source_artifact_check"]

    assert calibration["status"] == "partially_resolved"
    assert calibration["methods"][0]["reference_eos_record"] == (
        "neon_fcc_fei_2007_vinet_2"
    )
    assert calibration["recalculation"]["status"] == ("missing_calibrant_observations")
    assert "no separate CSV" in artifact["repository_file_inventory"]
    assert "requested from the authors" in artifact["data_availability"]
    validate_pressure_calibration_references()


def test_ismailova_selected_rows_do_not_reproduce_the_source_fit():
    published = published_curve_diagnostic()
    unweighted = diagnostic_refit(errors_in_variables=False)
    weighted = diagnostic_refit(errors_in_variables=True)

    assert published["observations"] == 4
    assert published["pressure_rmse_gpa"] == pytest.approx(4.04371153455)
    assert published["max_abs_pressure_residual_gpa"] == pytest.approx(7.71306457768)
    assert unweighted["parameters"] == pytest.approx(
        {"V0": 172.092200829, "K0": 233.744397583, "K0_prime": 4.0}
    )
    assert weighted["parameters"] == pytest.approx(
        {"V0": 172.179953598, "K0": 234.492410265, "K0_prime": 4.0}
    )
    assert weighted["reduced_chi_square"] == pytest.approx(21.1758931382)

    outcome = reproduce()["source_fit"]
    assert outcome["status"] == "not_refittable"
    assert outcome["reported_exclusions"] == "none stated"
    assert len(outcome["blockers"]) == 4
