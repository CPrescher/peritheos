import csv
import hashlib
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import get_material_document
from peritheos.materials import Material

DOI = "10.1029/2008GL035658"
RECORD_IDS = [
    "bridgmanite_katsura_2009_bm3",
    "bridgmanite_katsura_2009_bm3_k0_fixed",
]
DATASET_ID = "bridgmanite_katsura_2009_corrected_table1_pvt"


def _source_records():
    document = get_material_document("bridgmanite")
    records = [
        record
        for record in document["eos_records"]
        if record["reference"].get("doi", "").lower() == DOI.lower()
    ]
    assert [record["identifier"] for record in records] == RECORD_IDS
    return document, records


def _dataset_rows(document):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    payload = path.read_bytes()
    return dataset, payload, list(csv.DictReader(payload.decode().splitlines()))


def _bm3_pressure(volume_ratio, k0, k0_prime):
    eta = volume_ratio ** (-1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def test_katsura_2009_registers_two_source_generated_bm3_fits():
    document, records = _source_records()
    assert [record["eos"]["parameters"] for record in records] == [
        {"V0": 162.38559886687202, "K0": 256.0, "K0_prime": 3.8},
        {"V0": 162.38559886687202, "K0": 253.0, "K0_prime": 4.1},
    ]
    assert records[0]["fixed_parameters"] == ["V0"]
    assert records[1]["fixed_parameters"] == ["V0", "K0"]
    assert records[0]["parameter_errors"] == {
        "V0": pytest.approx(0.011812486067495595),
        "K0": 2.0,
        "K0_prime": 0.2,
    }
    assert records[1]["parameter_errors"] == {
        "V0": pytest.approx(0.011812486067495595),
        "K0": None,
        "K0_prime": 0.2,
    }
    for identifier in RECORD_IDS:
        loaded = Material.from_eosmat(document, record_identifiers=[identifier])
        assert loaded.eos_records[0].reference_temperature == 300.0


def test_katsura_2009_corrected_table_is_complete_and_checksummed():
    document, records = _source_records()
    dataset, payload, rows = _dataset_rows(document)
    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    assert len(rows) == 89
    assert sum(row["apparatus"] == "WC" for row in rows) == 36
    assert sum(row["apparatus"] == "SD" for row in rows) == 53
    assert sum(float(row["temperature_k"]) <= 308.0 for row in rows) == 9
    assert rows[0] == {
        "source_order": "1",
        "apparatus": "WC",
        "temperature_k": "1900",
        "mgo_v_v0": "0.9313",
        "mgo_v_v0_uncertainty": "0.0001",
        "pressure_gpa": "22.8",
        "pressure_uncertainty_gpa": "0.3",
        "bridgmanite_v_v0": "0.9614",
        "bridgmanite_v_v0_uncertainty": "0.0001",
    }
    assert rows[-1]["temperature_k"] == "305"
    assert rows[-1]["pressure_gpa"] == "28.9"
    assert rows[-1]["bridgmanite_v_v0"] == "0.9104"
    assert dataset["reference"]["doi"] == "10.1029/2009GL039318"
    assert dataset["used_by_eos_records"] == RECORD_IDS
    assert all(record["fit_datasets"] == [DATASET_ID] for record in records)


def test_katsura_2009_preferred_curve_and_conditional_refit_have_parity():
    document, records = _source_records()
    _, _, rows = _dataset_rows(document)
    ambient = [row for row in rows if float(row["temperature_k"]) <= 308.0]
    pressure = np.asarray([float(row["pressure_gpa"]) for row in ambient])
    volume_ratio = np.asarray([float(row["bridgmanite_v_v0"]) for row in ambient])
    published_pressure = _bm3_pressure(volume_ratio, 256.0, 3.8)
    fit = least_squares(
        lambda value: _bm3_pressure(volume_ratio, value[0], 3.8) - pressure,
        np.asarray([256.0]),
    )
    reproduction = records[0]["scientific_validation"]["numerical_reproduction"]
    assert np.sqrt(np.mean((published_pressure - pressure) ** 2)) == pytest.approx(
        0.3202084026, abs=5.0e-10
    )
    assert np.max(np.abs(published_pressure - pressure)) == pytest.approx(
        0.5699353837, abs=5.0e-10
    )
    assert fit.x[0] == pytest.approx(255.5075102, abs=5.0e-7)
    assert reproduction["conditional_refit"]["fitted_K0_gpa"] == pytest.approx(fit.x[0])
    assert abs(fit.x[0] - 256.0) < 2.0


def test_katsura_2009_k0_fixed_sensitivity_fit_has_parity():
    document, records = _source_records()
    _, _, rows = _dataset_rows(document)
    ambient = [row for row in rows if float(row["temperature_k"]) <= 308.0]
    pressure = np.asarray([float(row["pressure_gpa"]) for row in ambient])
    pressure_error = np.asarray(
        [float(row["pressure_uncertainty_gpa"]) for row in ambient]
    )
    volume_ratio = np.asarray([float(row["bridgmanite_v_v0"]) for row in ambient])
    fit = least_squares(
        lambda value: (
            (_bm3_pressure(volume_ratio, 253.0, value[0]) - pressure) / pressure_error
        ),
        np.asarray([4.1]),
    )
    reproduction = records[1]["scientific_validation"]["numerical_reproduction"]
    assert fit.x[0] == pytest.approx(4.12823286, abs=5.0e-8)
    assert reproduction["fitted_K0_prime"] == pytest.approx(fit.x[0])
    assert abs(fit.x[0] - 4.1) < 0.2


def test_katsura_2009_preserves_corrigendum_and_pressure_scale_scope():
    _, records = _source_records()
    preferred, sensitivity = records
    assert preferred["scientific_validation"]["status"] == "primary_source_validated"
    assert "formal corrigendum" in preferred["scientific_validation"]["note"]
    assert preferred["pressure_calibration"]["status"] == "partially_resolved"
    assert preferred["pressure_calibration"]["methods"][0]["kind"] == "other"
    assert (
        preferred["pressure_calibration"]["recalculation"]["status"]
        == "reference_eos_not_bundled"
    )
    reported = preferred["scientific_validation"]["reported_parameterizations"]
    assert {item["role"] for item in reported} == {
        "preferred_static_fit",
        "same_data_K0_fixed_sensitivity_fit",
        "full_thermal_fit_not_executable_here",
    }
    assert sensitivity["scientific_validation"]["numerical_reproduction"]["result"] == (
        "parity_accepted"
    )
