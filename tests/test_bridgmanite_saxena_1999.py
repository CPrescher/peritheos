import csv
import hashlib
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import brentq, least_squares

from peritheos import get_eos_record_document, get_material_document
from peritheos.materials import Material
from peritheos.units import molar_volume_to_cell_volume

DOI = "10.2138/am-1999-0303"
RECORD_ID = "bridgmanite_saxena_1999_bm3_4"
DATASET_ID = "bridgmanite_saxena_1999_table1_pvt"


def _document_record_and_loaded():
    document = get_material_document("bridgmanite")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    loaded = Material.from_eosmat(document, record_identifiers=[RECORD_ID])
    return document, record, loaded.eos_records[0]


def _dataset_rows(document):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    payload = path.read_bytes()
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    return dataset, payload, rows


def test_saxena_1999_registers_only_the_published_300_k_bm3_slice():
    _, record, loaded = _document_record_and_loaded()

    assert record["reference"]["doi"] == DOI
    assert record["eos"] == {
        "type": "BM3",
        "parameters": {"V0": 162.4007, "K0": 260.51, "K0_prime": 4.0},
        "model": "birch_murnaghan_3",
    }
    assert record["parameter_errors"] == {
        "V0": None,
        "K0": None,
        "K0_prime": None,
    }
    assert record["parameter_covariance"] is None
    assert record["fixed_parameters"] == ["K0_prime"]
    assert record["temperature_ref"] == 300.0
    assert record["validity"]["temperature_k"] == [300.0, 300.0]
    assert "thermal" not in record
    assert "cannot encode exactly" in record["notes"]
    assert get_eos_record_document(RECORD_ID)["identifier"] == RECORD_ID
    assert loaded.reference_temperature == 300.0


def test_saxena_1999_table1_is_complete_and_checksummed():
    document, record, _ = _document_record_and_loaded()
    dataset, payload, rows = _dataset_rows(document)

    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    assert len(rows) == 37
    assert sum(row["temperature_k"] == "300" for row in rows) == 10
    assert rows[0] == {
        "pressure_gpa": "36.2",
        "temperature_k": "300",
        "molar_volume_cm3_mol": "21.874",
        "molar_volume_uncertainty_cm3_mol": "0.011",
    }
    assert rows[23] == {
        "pressure_gpa": "82.8",
        "temperature_k": "1485",
        "molar_volume_cm3_mol": "20.09",
        "molar_volume_uncertainty_cm3_mol": "0.10",
    }
    assert rows[-1] == {
        "pressure_gpa": "109.1",
        "temperature_k": "1357",
        "molar_volume_cm3_mol": "19.172",
        "molar_volume_uncertainty_cm3_mol": "0.015",
    }
    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert record["fit_datasets"] == [DATASET_ID]


def test_saxena_1999_reproduces_the_highest_pressure_300_k_observation():
    _, record_document, record = _document_record_and_loaded()
    volume = molar_volume_to_cell_volume(19.846, 4, from_unit="cm^3/mol")

    # The rounded source volume predicts 0.0577 GPa above the printed upper
    # range, so evaluate the literature row without strict marginal checking.
    pressure = record.pressure(volume, 300.0)
    assert pressure == pytest.approx(82.5577, abs=5.0e-5)
    assert pressure == pytest.approx(82.5, abs=0.06)

    parameters = record_document["eos"]["parameters"]
    eta = (parameters["V0"] / volume) ** (1.0 / 3.0)
    independent_bm3 = (
        1.5
        * parameters["K0"]
        * (eta**7 - eta**5)
        * (1.0 + 0.75 * (parameters["K0_prime"] - 4.0) * (eta**2 - 1.0))
    )
    assert pressure == pytest.approx(independent_bm3, rel=2.0e-12)

    round_trip_volume = record.volume(82.5, 300.0, check_validity=True)
    assert record.pressure(
        round_trip_volume, 300.0, check_validity=True
    ) == pytest.approx(82.5, rel=1.0e-11)


def test_saxena_1999_independent_300_k_volume_weighted_refit_has_parity():
    document, record, _ = _document_record_and_loaded()
    _, _, rows = _dataset_rows(document)
    room_temperature = [row for row in rows if row["temperature_k"] == "300"]
    pressures = np.asarray([float(row["pressure_gpa"]) for row in room_temperature])
    volumes = np.asarray(
        [float(row["molar_volume_cm3_mol"]) for row in room_temperature]
    )
    volume_errors = np.asarray(
        [float(row["molar_volume_uncertainty_cm3_mol"]) for row in room_temperature]
    )

    def pressure_bm3(volume, v0, k0):
        eta = (v0 / volume) ** (1.0 / 3.0)
        return 1.5 * k0 * (eta**7 - eta**5)

    def volume_bm3(pressure, v0, k0):
        return brentq(
            lambda volume: pressure_bm3(volume, v0, k0) - pressure,
            0.5 * v0,
            v0,
        )

    def residuals(parameters):
        v0, k0 = parameters
        fitted_volumes = np.asarray(
            [volume_bm3(pressure, v0, k0) for pressure in pressures]
        )
        return (fitted_volumes - volumes) / volume_errors

    fit = least_squares(
        residuals, x0=(24.45, 260.51), bounds=((20.0, 100.0), (30.0, 500.0))
    )
    assert fit.success

    fitted_v0_molar, fitted_k0 = fit.x
    fitted_v0_cell = molar_volume_to_cell_volume(
        fitted_v0_molar, 4, from_unit="cm^3/mol"
    )
    covariance = np.linalg.inv(fit.jac.T @ fit.jac)
    standard_errors = np.sqrt(np.diag(covariance))
    fitted_v0_cell_error = molar_volume_to_cell_volume(
        standard_errors[0], 4, from_unit="cm^3/mol"
    )
    reduced_chi_square = np.sum(fit.fun**2) / (len(volumes) - len(fit.x))

    reproduction = record["scientific_validation"]["numerical_reproduction"][
        "independent_300_k_refit"
    ]
    result = reproduction["result"]
    parity = reproduction["published_parity"]

    assert fitted_v0_molar == pytest.approx(24.44757948, abs=5.0e-7)
    assert fitted_v0_cell == pytest.approx(162.38464327, abs=5.0e-7)
    assert fitted_k0 == pytest.approx(260.48766286, abs=5.0e-7)
    assert fitted_v0_cell_error == pytest.approx(0.2142853, abs=5.0e-7)
    assert standard_errors[1] == pytest.approx(2.62311623, abs=5.0e-7)
    assert reduced_chi_square == pytest.approx(0.00110423827, abs=5.0e-10)
    assert result["V0_molar_cm3_mol"] == pytest.approx(fitted_v0_molar)
    assert result["V0_cell_angstrom3"] == pytest.approx(fitted_v0_cell)
    assert result["K0_gpa"] == pytest.approx(fitted_k0)
    assert parity["delta_V0_cell_angstrom3"] == pytest.approx(
        fitted_v0_cell - record["eos"]["parameters"]["V0"], abs=5.0e-8
    )
    assert parity["delta_K0_gpa"] == pytest.approx(
        fitted_k0 - record["eos"]["parameters"]["K0"], abs=5.0e-8
    )
    assert parity["absolute_relative_delta_V0_percent"] < 0.01
    assert parity["absolute_relative_delta_K0_percent"] < 0.01
    assert reproduction["free_parameters"] == ["V0", "K0"]
    assert "K0_prime fixed" in reproduction["model"]
    assert "relative weights" in reproduction["uncertainty_treatment"]
    assert (
        "not a common pressure sigma" in reproduction["pressure_uncertainty_treatment"]
    )
    assert reproduction["status"] == "parity_accepted"


def test_saxena_1999_preserves_scope_calibration_and_missing_statistics():
    _, record, _ = _document_record_and_loaded()
    validation = record["scientific_validation"]
    calibration = record["pressure_calibration"]

    assert validation["status"] == "primary_source_validated"
    assert validation["primary_source_check"]["doi"] == DOI
    assert validation["primary_data_check"]["status"] == "bundled"
    assert validation["primary_data_check"]["dataset_identifiers"] == [DATASET_ID]
    assert {item["role"] for item in validation["reported_parameterizations"]} == {
        "model_1_not_executable",
        "model_2_not_executable",
    }
    assert (
        "full thermal/literature regression remains impossible"
        in validation["numerical_reproduction"]["limitation"]
    )
    assert calibration["status"] == "partially_resolved"
    assert calibration["methods"][0]["reference"]["authors"] == [
        "Jamieson",
        "Fritz",
        "Manghnani",
    ]
    assert calibration["recalculation"]["status"] == "missing_calibrant_observations"
