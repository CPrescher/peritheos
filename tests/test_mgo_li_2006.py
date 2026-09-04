import csv
import hashlib
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import get_eos_record_document, get_material_document
from peritheos.materials import Material

DOI = "10.1029/2005JB004251"
RECORD_ID = "mgo_li_2006_bm3_absolute_acoustic"
DATASET_ID = "mgo_li_2006_table1_elasticity"


def _source_record():
    document = get_material_document("mgo")
    records = [
        record
        for record in document["eos_records"]
        if record["reference"].get("doi", "").lower() == DOI.lower()
    ]
    assert len(records) == 1
    return document, records[0]


def _dataset_rows(document):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    payload = path.read_bytes()
    return dataset, payload, list(csv.DictReader(payload.decode("utf-8").splitlines()))


def test_li_2006_absolute_acoustic_bm3_is_executable_and_distinct():
    document, source = _source_record()
    loaded = Material.from_eosmat(document, record_identifiers=[RECORD_ID])
    record = loaded.eos_records[0]

    assert get_eos_record_document(RECORD_ID)["identifier"] == RECORD_ID
    assert source["record_kind"] == "published"
    assert source["eos"] == {
        "type": "BM3",
        "parameters": {
            "V0": 74.702,
            "K0": 161.1767340842165,
            "K0_prime": 4.237911491402006,
        },
        "model": "birch_murnaghan_3",
    }
    assert source["fixed_parameters"] == ["V0"]
    assert source["parameter_errors"] == {
        "V0": 0.026,
        "K0": None,
        "K0_prime": None,
    }
    assert source["parameter_covariance"] is None
    assert source["fit_datasets"] == [DATASET_ID]
    assert source["pressure_calibration"]["methods"][0]["kind"] == "self_consistent"
    assert source["pressure_calibration"]["recalculation"]["status"] == "not_applicable"
    assert record.eos.bulk_modulus(record.reference_volume) == pytest.approx(
        source["eos"]["parameters"]["K0"], rel=1.0e-13
    )

    other_bm3_parameters = {
        (
            item["eos"]["parameters"]["V0"],
            item["eos"]["parameters"]["K0"],
            item["eos"]["parameters"]["K0_prime"],
        )
        for item in document["eos_records"]
        if item["identifier"] != RECORD_ID and item["eos"]["type"] == "BM3"
    }
    source_parameters = source["eos"]["parameters"]
    assert (
        source_parameters["V0"],
        source_parameters["K0"],
        source_parameters["K0_prime"],
    ) not in other_bm3_parameters


def test_li_2006_isothermal_conversion_reproduces_stored_parameters():
    _, source = _source_record()
    inputs = source["scientific_validation"]["isothermal_parameter_derivation"][
        "inputs"
    ]
    alpha_gamma_t = inputs["alpha_per_k"] * inputs["gamma0"] * inputs["temperature_k"]
    k0_t = inputs["K0S_gpa"] / (1.0 + alpha_gamma_t)
    k0_t_prime = (
        inputs["K0S_prime"]
        + inputs["q"] * alpha_gamma_t
        - inputs["gamma0"]
        * inputs["temperature_k"]
        * inputs["dK0T_dT_gpa_per_k"]
        / k0_t
    ) / (1.0 + alpha_gamma_t)

    assert source["eos"]["parameters"]["K0"] == pytest.approx(k0_t, rel=1.0e-14)
    assert source["eos"]["parameters"]["K0_prime"] == pytest.approx(
        k0_t_prime, rel=1.0e-14
    )
    assert "not a direct volumetric BM3 fit" in source["notes"]
    assert (
        source["scientific_validation"]["numerical_reproduction"]["refit_status"]
        == "not_directly_refittable_as_pressure_volume"
    )


def test_li_2006_reproduces_table1_endpoint_and_inverts_reported_range():
    document, source = _source_record()
    record = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]
    reference_volume = source["eos"]["parameters"]["V0"]

    assert record.pressure(reference_volume * 0.942, 300.0) == pytest.approx(
        10.95, abs=0.15
    )
    for pressure in (0.0, 1.7, 6.23, 10.9):
        volume = record.volume(pressure, 300.0, check_validity=True)
        assert record.pressure(volume, 300.0, check_validity=True) == pytest.approx(
            pressure, rel=1.0e-11, abs=1.0e-12
        )


def test_li_2006_independent_acoustic_finite_strain_refit():
    document, source = _source_record()
    _, _, rows = _dataset_rows(document)
    fit_rows = [
        row
        for row in rows
        if row["experimental_path"] in {"ambient", "decompression_after_annealing"}
    ]
    density = np.array([float(row["density_g_cm3"]) for row in fit_rows])
    p_velocity = np.array([float(row["p_wave_velocity_km_s"]) for row in fit_rows])
    s_velocity = np.array([float(row["s_wave_velocity_km_s"]) for row in fit_rows])
    density_0 = density[0]
    epsilon = (1.0 - (density / density_0) ** (2.0 / 3.0)) / 2.0

    def predicted_velocities(parameters):
        k_0s, k_0s_prime, g_0, g_0_prime = parameters
        l_1 = k_0s + 4.0 * g_0 / 3.0
        l_2 = 5.0 * l_1 - 3.0 * k_0s * (k_0s_prime + 4.0 * g_0_prime / 3.0)
        m_1 = g_0
        m_2 = 5.0 * g_0 - 3.0 * k_0s * g_0_prime
        strain_factor = (1.0 - 2.0 * epsilon) ** 2.5
        return (
            np.sqrt(strain_factor * (l_1 + l_2 * epsilon) / density),
            np.sqrt(strain_factor * (m_1 + m_2 * epsilon) / density),
        )

    def residuals(parameters):
        predicted_p, predicted_s = predicted_velocities(parameters)
        return np.concatenate((predicted_p - p_velocity, predicted_s - s_velocity))

    fit = least_squares(
        residuals,
        np.array([163.5, 4.20, 129.8, 2.42]),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
        max_nfev=5000,
    )
    assert fit.success
    assert fit.x == pytest.approx(
        [163.51148649, 4.15806624, 129.71748390, 2.39723610], abs=1.0e-7
    )

    published = source["scientific_validation"]["reported_parameterizations"][0]
    for actual, value_key, error_key in (
        (fit.x[0], "K0S_gpa", "K0S_standard_deviation_gpa"),
        (fit.x[1], "K0S_prime", "K0S_prime_standard_deviation"),
        (fit.x[2], "G0_gpa", "G0_standard_deviation_gpa"),
        (fit.x[3], "G0_prime", "G0_prime_standard_deviation"),
    ):
        assert abs(actual - published[value_key]) < published[error_key]

    conversion = source["scientific_validation"]["isothermal_parameter_derivation"]
    inputs = conversion["inputs"]
    alpha_gamma_t = inputs["alpha_per_k"] * inputs["gamma0"] * inputs["temperature_k"]
    refit_k_0t = fit.x[0] / (1.0 + alpha_gamma_t)
    refit_k_0t_prime = (
        fit.x[1]
        + inputs["q"] * alpha_gamma_t
        - inputs["gamma0"]
        * inputs["temperature_k"]
        * inputs["dK0T_dT_gpa_per_k"]
        / refit_k_0t
    ) / (1.0 + alpha_gamma_t)
    assert refit_k_0t == pytest.approx(161.18805735, abs=1.0e-7)
    assert refit_k_0t_prime == pytest.approx(4.19656803, abs=1.0e-7)
    assert abs(refit_k_0t - conversion["outputs"]["K0T_gpa"]) < 0.02
    assert abs(refit_k_0t_prime - conversion["outputs"]["K0T_prime"]) < 0.05

    predicted_p, predicted_s = predicted_velocities(fit.x)
    assert np.sqrt(np.mean((predicted_p[1:] - p_velocity[1:]) ** 2)) == pytest.approx(
        0.02164856, abs=1.0e-8
    )
    assert np.sqrt(np.mean((predicted_s[1:] - s_velocity[1:]) ** 2)) == pytest.approx(
        0.01473844, abs=1.0e-8
    )


def test_li_2006_table1_is_complete_verbatim_and_checksummed():
    document, source = _source_record()
    dataset, payload, rows = _dataset_rows(document)

    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    assert dataset["used_by_eos_records"] == [source["identifier"]]
    assert len(rows) == 18
    assert sum(row["experimental_path"] == "ambient" for row in rows) == 1
    assert sum(row["experimental_path"] == "initial_compression" for row in rows) == 7
    assert (
        sum(row["experimental_path"] == "decompression_after_annealing" for row in rows)
        == 10
    )
    assert rows[0]["density_g_cm3"] == "3.566"
    assert rows[0]["p_wave_velocity_km_s"] == "9.74"
    assert rows[0]["s_wave_velocity_km_s"] == "6.00"
    assert rows[0]["adiabatic_bulk_modulus_gpa"] == ""
    assert rows[0]["calculated_absolute_pressure_gpa"] == ""
    assert rows[1] == {
        "experimental_path": "initial_compression",
        "volume_ratio": "0.990",
        "density_g_cm3": "3.603",
        "density_standard_deviation_g_cm3": "0.006",
        "p_wave_velocity_km_s": "9.81",
        "p_wave_velocity_standard_deviation_km_s": "0.01",
        "s_wave_velocity_km_s": "6.08",
        "s_wave_velocity_standard_deviation_km_s": "0.01",
        "adiabatic_bulk_modulus_gpa": "169.2",
        "adiabatic_bulk_modulus_standard_deviation_gpa": "1.2",
        "shear_modulus_gpa": "133.2",
        "shear_modulus_standard_deviation_gpa": "0.5",
        "calculated_absolute_pressure_gpa": "1.70",
        "calculated_absolute_pressure_standard_deviation_gpa": "0.09",
    }
    assert rows[4]["volume_ratio"] == "0.962"
    assert rows[4]["density_g_cm3"] == "3.696"
    assert rows[4]["calculated_absolute_pressure_gpa"] == "6.23"
    assert rows[7]["calculated_absolute_pressure_gpa"] == "10.95"
    assert rows[-1]["calculated_absolute_pressure_gpa"] == "2.29"


def test_li_2006_primary_source_inconsistencies_are_not_silently_resolved():
    _, source = _source_record()
    inconsistencies = source["scientific_validation"]["reported_inconsistencies"]

    assert inconsistencies[0]["abstract_and_section_4"] == "163.5(11) GPa"
    assert inconsistencies[0]["table_2"] == "163.6(11) GPa"
    assert source["eos"]["parameters"]["K0"] < 163.5
    assert inconsistencies[1]["table_1"].startswith("V/V0=0.962")
