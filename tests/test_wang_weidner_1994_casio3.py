import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM2
from peritheos.fitting import fit_rt_eos

ROOT = Path(__file__).resolve().parents[1]
RECORD_ID = "ca_perovskite_wang_weidner_1994_bm2"
DATASET_ID = "ca_perovskite_wang_weidner_1994_figure3_digitized"
DATASET_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "ca-perovskite-wang-weidner-1994-figure3-digitized.csv"
)
DATASET_SHA256 = "d4e677cfa40f39ce67ca591fb24682e8b0848816798f3ba007ebefb3ab4fabea"


def source_record_and_dataset():
    document = get_material_document("ca_perovskite")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == RECORD_ID
    )
    dataset = next(
        row for row in document["datasets"] if row["identifier"] == DATASET_ID
    )
    return record, dataset


def test_wang_weidner_source_scope_and_calibration_are_explicit():
    record, dataset = source_record_and_dataset()

    assert record["eos"]["type"] == "BM2"
    assert record["eos"]["parameters"] == {"V0": 45.83, "K0": 280.0}
    assert record["parameter_errors"] == {"V0": 0.07, "K0": 23.0}
    assert "K0'=4" in record["parameter_provenance"]["equation"]
    assert "Exactly four" in record["parameter_provenance"]["fit_selection"]
    assert record["experimental_pressure_range_gpa"] == [2.0, 9.0]
    assert "11.7 GPa" in record["validity"]["notes"][0]

    calibration = record["pressure_calibration"]
    assert calibration["status"] == "resolved"
    assert calibration["methods"][0]["material"] == "NaCl"
    assert calibration["methods"][0]["reference"]["doi"] == "10.1063/1.1660714"
    assert calibration["recalculation"]["status"] == "missing_calibrant_observations"

    check = record["scientific_validation"]["primary_data_check"]
    assert check["status"] == "plot_only"
    assert check["digitized_dataset_identifiers"] == [DATASET_ID]
    assert dataset["used_by_eos_records"] == [RECORD_ID]


def test_wang_weidner_plot_dataset_is_complete_and_hashed():
    _, dataset = source_record_and_dataset()
    with DATASET_PATH.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    assert hashlib.sha256(DATASET_PATH.read_bytes()).hexdigest() == DATASET_SHA256
    assert dataset["resource"]["sha256"] == DATASET_SHA256
    assert (
        len(rows)
        == dataset["provenance"]["quality_control"]["source_stated_fit_observations"]
        == 4
    )
    assert all(row["used_in_published_fit"] == "1" for row in rows)
    assert [float(row["pressure_gpa"]) for row in rows] == pytest.approx(
        [2.74, 4.10, 6.14, 9.24]
    )
    assert [float(row["volume_a3_per_formula_unit"]) for row in rows] == pytest.approx(
        [45.341, 45.234, 44.887, 44.419]
    )


def test_wang_weidner_four_point_bm2_refit_recovers_published_coefficients():
    record, _ = source_record_and_dataset()
    with DATASET_PATH.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3_per_formula_unit"]) for row in rows])

    result = fit_rt_eos(
        BM2,
        volume,
        pressure,
        initial={"V0": 45.83, "K0": 280.0},
    )

    assert result.success
    assert result.parameters["V0"] == pytest.approx(45.832161459, abs=2.0e-8)
    assert result.parameters["K0"] == pytest.approx(279.27071645, abs=1.0e-6)
    assert abs(result.parameters["V0"] - 45.83) < record["parameter_errors"]["V0"]
    assert abs(result.parameters["K0"] - 280.0) < record["parameter_errors"]["K0"]
    published_residuals = BM2(45.83, 280.0).pressure(volume) - pressure
    assert np.sqrt(np.mean(published_residuals**2)) == pytest.approx(0.241748597567)
