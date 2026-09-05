"""Primary-source checks for Xu et al. (2020) Al-bearing phase D."""

import csv
import hashlib
import math
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.materials import Material

DOI = "10.1029/2020GL088877"
RECORD_ID = "al_bearing_phase_d_xu_2020_bm3_300k"
DATASET_ID = "al_bearing_phase_d_xu_2020_earthchem_pvt"
DATASET_SHA256 = "4f462df434ca353b3b6d090980efc50ef2ad56f9ad99db8dcba1946dac0413c9"
WORKBOOK_SHA256 = "1349bc1260337c1e8a3de6830e089b4ad5c32c18b3aea69f0ff6986a9b9c4d3e"


def _document_record_dataset_rows():
    document = get_material_document("al_bearing_phase_d")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    payload = path.read_bytes()
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    return document, record, dataset, payload, rows


def _effective_variance_refit(rows):
    selected = [row for row in rows if row["used_in_room_temperature_fit"] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in selected])
    pressure_sigma = np.array(
        [
            0.0
            if row["pressure_uncertainty_gpa"] == ""
            else float(row["pressure_uncertainty_gpa"])
            for row in selected
        ]
    )
    volume = np.array([float(row["volume_a3"]) for row in selected])
    volume_sigma = np.array(
        [float(row["volume_uncertainty_a3"]) for row in selected]
    )
    parameters = np.array([143.0, 5.8])

    for _ in range(100):
        model = BM3(V0=86.71, K0=parameters[0], K0_prime=parameters[1])
        dp_dv = -np.asarray(model.bulk_modulus(volume), dtype=float) / volume
        sigma_effective = np.sqrt(
            pressure_sigma**2 + (dp_dv * volume_sigma) ** 2
        )

        def residual(candidate):
            candidate_model = BM3(
                V0=86.71, K0=candidate[0], K0_prime=candidate[1]
            )
            return (
                np.asarray(candidate_model.pressure(volume), dtype=float) - pressure
            ) / sigma_effective

        optimization = least_squares(
            residual,
            parameters,
            bounds=([100.0, 0.0], [200.0, 10.0]),
            x_scale="jac",
        )
        if np.allclose(optimization.x, parameters, rtol=1.0e-13, atol=1.0e-14):
            return optimization.x
        parameters = optimization.x

    raise AssertionError("effective-variance refit did not converge")


def test_xu_2020_material_identity_composition_and_volume_basis():
    document, _, _, _, _ = _document_record_dataset_rows()

    assert document["formula"] == "Mg0.90Al0.64Si1.29H3.10O6"
    assert document["phase"] == "Al-bearing trigonal phase D, P-31m"
    assert document["space_group"] == "P-31m"
    assert document["space_group_number"] == 162
    assert document["formula_units_per_cell"] == 1
    assert document["source"]["earthchem_record"]["sha256"] == WORKBOOK_SHA256

    lattice = document["lattice"]
    ambient_volume = (
        lattice["a"]
        * lattice["b"]
        * lattice["c"]
        * math.sin(math.radians(lattice["gamma"]))
    )
    assert ambient_volume == pytest.approx(86.51, abs=0.01)
    assert "86.71(2)" in document["notes"]
    assert "not a direct site-occupancy" in document["notes"]


def test_xu_2020_earthchem_pvt_is_complete_verbatim_and_checksummed():
    _, record, dataset, payload, rows = _document_record_dataset_rows()

    assert hashlib.sha256(payload).hexdigest() == DATASET_SHA256
    assert dataset["resource"]["sha256"] == DATASET_SHA256
    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == 28
    assert [int(row["source_order"]) for row in rows] == list(range(1, 29))
    assert sum(int(row["used_in_room_temperature_fit"]) for row in rows) == 8
    assert {float(row["temperature_k"]) for row in rows} == {
        300.0,
        500.0,
        700.0,
        900.0,
        1100.0,
        1300.0,
    }
    assert rows[0] == {
        "source_order": "1",
        "used_in_room_temperature_fit": "1",
        "pressure_gpa": "0.00",
        "pressure_uncertainty_gpa": "",
        "temperature_k": "300",
        "lattice_a_angstrom": "4.8078",
        "lattice_a_uncertainty_angstrom": "0.0006",
        "lattice_c_angstrom": "4.3316",
        "lattice_c_uncertainty_angstrom": "0.0002",
        "volume_a3": "86.71",
        "volume_uncertainty_a3": "0.02",
    }
    assert rows[1]["pressure_gpa"] == "17.96"
    assert rows[25]["pressure_gpa"] == "20.53"
    assert rows[-1]["temperature_k"] == "1300"
    assert rows[-1]["volume_a3"] == "79.63"
    assert dataset["provenance"]["source_file_official_sha1"] == (
        "cddfb0e859c36847ec75e776495ad8ef7f95b456"
    )
    assert record["fit_datasets"] == [DATASET_ID]

    for row in rows:
        calculated_volume = (
            float(row["lattice_a_angstrom"]) ** 2
            * float(row["lattice_c_angstrom"])
            * math.sin(math.radians(120.0))
        )
        assert calculated_volume == pytest.approx(float(row["volume_a3"]), abs=0.02)


def test_xu_2020_published_bm3_executes_and_reproduces_observations():
    document, source, _, _, _ = _document_record_dataset_rows()
    material = Material.from_eosmat(document, record_identifiers=[RECORD_ID])
    record = material.eos_records[0]

    assert source["reference"]["doi"] == DOI
    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 86.71, "K0": 143.0, "K0_prime": 5.8},
    }
    assert source["fixed_parameters"] == ["V0"]
    assert source["parameter_errors"] == {
        "V0": 0.02,
        "K0": 5.0,
        "K0_prime": 0.7,
    }
    assert source["parameter_covariance"] is None
    assert source["validity"]["pressure_gpa"] == [0.0, 20.53]
    assert source["pressure_calibration"]["methods"][0]["reference"]["doi"] == (
        "10.1029/2003JB002446"
    )
    assert source["pressure_calibration"]["recalculation"]["status"] == (
        "missing_calibrant_observations"
    )

    assert record.pressure(86.71, 300.0) == pytest.approx(0.0, abs=1.0e-12)
    assert record.eos.bulk_modulus(86.71) == pytest.approx(143.0, rel=1.0e-13)
    assert record.pressure(78.01, 300.0) == pytest.approx(20.53, abs=0.01)
    assert record.pressure(80.07, 300.0) == pytest.approx(14.18, abs=0.17)
    for pressure in (0.0, 5.0, 14.18, 20.5):
        volume = record.volume(pressure, 300.0, check_validity=True)
        assert record.pressure(
            volume, 300.0, check_validity=True
        ) == pytest.approx(pressure, rel=1.0e-11, abs=1.0e-12)


def test_xu_2020_effective_variance_refit_recovers_published_coefficients():
    _, source, _, _, rows = _document_record_dataset_rows()
    refit = _effective_variance_refit(rows)

    assert refit == pytest.approx([143.07336748, 5.86138397], abs=2.0e-7)
    published = source["eos"]["parameters"]
    assert refit[0] == pytest.approx(published["K0"], abs=0.08)
    assert refit[1] == pytest.approx(published["K0_prime"], abs=0.07)

    reproduction = source["scientific_validation"]["numerical_reproduction"]
    assert reproduction["published_curve_rmse_gpa"] == pytest.approx(
        0.12191813, abs=1.0e-8
    )
    assert reproduction["refit_curve_rmse_gpa"] == pytest.approx(
        0.11079143, abs=1.0e-8
    )


def test_xu_2020_does_not_silently_encode_ambiguous_thermal_or_acoustic_models():
    _, source, _, _, _ = _document_record_dataset_rows()
    validation = source["scientific_validation"]

    assert source["equation_kind"] == "isothermal"
    assert "thermal" not in source
    assert len(validation["unrepresented_parameterizations"]) == 2
    assert {item["role"] for item in validation["unrepresented_parameterizations"]} == {
        "high_temperature_bm3",
        "pressure_scale_free_acoustic_finite_strain",
    }
    assert len(validation["reported_inconsistencies"]) == 2
    assert validation["reported_inconsistencies"][0]["field"] == (
        "ambient conventional-cell volume"
    )
    assert validation["reported_inconsistencies"][1]["section_3_1"] == (
        "143(5) GPa"
    )
