import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document, validate_eosmat_document

SOURCE_ID = "mgo_b1_matsui_2000_bm3_300k"
REFIT_ID = "mgo_b1_matsui_2000_taylor_refit"
DATASET_ID = "mgo_matsui_2000_table3_pvt"


def _matsui_rows(document):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = Path("peritheos/data").joinpath(dataset["resource"]["path"])
    with path.open(encoding="utf-8", newline="") as stream:
        return dataset, list(csv.DictReader(stream))


def test_matsui_2000_records_and_table3_transcription():
    document = get_material_document("mgo")
    validate_eosmat_document(document)
    records = {record["identifier"]: record for record in document["eos_records"]}

    source = records[SOURCE_ID]
    refit = records[REFIT_ID]
    assert source["reference"]["doi"] == "10.2138/am-2000-2-308"
    assert source["equation_kind"] == "isothermal"
    assert source["eos"]["parameters"] == pytest.approx(
        {"V0": 74.6458805788525, "K0": 160.5, "K0_prime": 4.10}
    )
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": pytest.approx(0.2),
        "K0_prime": pytest.approx(0.02),
    }
    assert source["fixed_parameters"] == []

    assert refit["record_kind"] == "refit"
    assert refit["derived_from_record"] == SOURCE_ID
    assert refit["fit_provenance"]["dataset"] == DATASET_ID
    assert refit["thermal"]["model"] == "second_order_taylor_thermal_pressure"

    dataset, rows = _matsui_rows(document)
    resource = Path("peritheos/data").joinpath(dataset["resource"]["path"])
    assert (
        hashlib.sha256(resource.read_bytes()).hexdigest()
        == dataset["resource"]["sha256"]
    )
    assert len(rows) == 88
    assert sum(row["row_kind"] == "md" for row in rows) == 72
    assert sum(row["row_kind"] == "observed" for row in rows) == 16
    assert rows[0] == {
        "pressure_gpa": "0",
        "temperature_k": "300",
        "volume_ratio": "1.0000",
        "row_kind": "observed",
        "comparison_source": "Dubrovinsky_and_Saxena_1997",
    }
    assert rows[-1] == {
        "pressure_gpa": "100",
        "temperature_k": "3000",
        "volume_ratio": "0.7493",
        "row_kind": "md",
        "comparison_source": "present_study",
    }


def test_matsui_2000_published_bm3_and_thermal_refit_reproduce_table3():
    document = get_material_document("mgo")
    material = Material.from_eosmat(document, record_identifiers=[SOURCE_ID, REFIT_ID])
    source = material.get_eos_record(SOURCE_ID)
    refit = material.get_eos_record(REFIT_ID)
    _, rows = _matsui_rows(document)
    md_rows = [row for row in rows if row["row_kind"] == "md"]

    pressure = np.array([float(row["pressure_gpa"]) for row in md_rows])
    temperature = np.array([float(row["temperature_k"]) for row in md_rows])
    volume_ratio = np.array([float(row["volume_ratio"]) for row in md_rows])
    volume = source.eos.V0 * volume_ratio

    low_pressure_300k = (temperature == 300.0) & (pressure <= 20.0)
    source_residuals = (
        source.pressure(volume[low_pressure_300k]) - pressure[low_pressure_300k]
    )
    assert np.sqrt(np.mean(source_residuals**2)) == pytest.approx(0.008744268563302926)
    assert np.max(np.abs(source_residuals)) == pytest.approx(0.01829082080574196)

    residuals = refit.pressure(volume, temperature) - pressure
    assert np.sqrt(np.mean(residuals**2)) == pytest.approx(0.2106462806006488)
    assert np.max(np.abs(residuals)) == pytest.approx(0.5638802265021639)

    high_temperature_volume = source.eos.V0 * 0.7493
    fitted_pressure = refit.pressure(high_temperature_volume, 3000.0)
    assert fitted_pressure == pytest.approx(99.95974750699662)
    assert refit.volume(fitted_pressure, 3000.0) == pytest.approx(
        high_temperature_volume
    )
