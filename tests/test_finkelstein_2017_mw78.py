"""Primary-source checks for Finkelstein et al. (2017) Mw78."""

import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document, validate_eosmat_document
from peritheos.eos.rt import BM3
from peritheos.pressure_calibrations import get_pressure_calibration

ROOT = Path(__file__).resolve().parents[1]
DATASET_SHA256 = "55ec130144ce8c0c7bd8fad6ebf45e3747ea631dfd901daf74993d7bd862f01d"
DATASET_ID = "magnesiowustite_mw78_finkelstein_2017_table2_compression"
HELIUM_ID = "mg0215fe0762vac0023o_finkelstein_2017_bm3_helium_1"
NEON_ID = "mg0215fe0762vac0023o_finkelstein_2017_bm3_neon_cubic_2"


def _document_dataset_rows_and_material():
    document = get_material_document("mg0215fe0762vac0023o_b1")
    dataset = document["datasets"][0]
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    material = Material.from_eosmat(document)
    return document, dataset, path, rows, material


def _effective_variance_refit(rows, *, v0_prior, k0_prime_prior=None):
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    pressure_sigma = np.array(
        [
            float(value) if value else 0.0
            for value in (row["pressure_standard_deviation_gpa"] for row in rows)
        ]
    )
    volume = np.array([float(row["cubic_volume_a3"]) for row in rows])
    volume_sigma = np.array(
        [float(row["cubic_volume_standard_deviation_a3"]) for row in rows]
    )

    def residual(parameters):
        model = BM3(V0=parameters[0], K0=parameters[1], K0_prime=parameters[2])
        predicted = np.asarray(model.pressure(volume), dtype=float)
        dp_dv = -np.asarray(model.bulk_modulus(volume), dtype=float) / volume
        sigma_effective = np.sqrt(pressure_sigma**2 + (dp_dv * volume_sigma) ** 2)
        result = list((predicted - pressure) / sigma_effective)
        result.append((parameters[0] - v0_prior[0]) / v0_prior[1])
        if k0_prime_prior is not None:
            result.append((parameters[2] - k0_prime_prior[0]) / k0_prime_prior[1])
        return np.asarray(result)

    return least_squares(
        residual,
        x0=[78.85, 155.0, 4.05],
        bounds=([70.0, 100.0, 0.0], [90.0, 250.0, 10.0]),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    ).x


def test_finkelstein_2017_document_is_valid_and_phase_specific():
    document, _, _, _, material = _document_dataset_rows_and_material()
    validate_eosmat_document(document)

    assert document["identifier"] == "mg0215fe0762vac0023o_b1"
    assert document["formula"] == "Mg0.215Fe0.762□0.023O"
    assert document["space_group"] == "Fm-3m"
    assert document["space_group_number"] == 225
    assert document["formula_units_per_cell"] == 4
    assert material.default_record().identifier == HELIUM_ID

    lattice = document["lattice"]
    assert lattice["a"] ** 3 == pytest.approx(78.94, abs=0.02)
    contents = {}
    for site in document["atom_sites"]:
        contents[site["element"]] = contents.get(site["element"], 0.0) + (
            site["multiplicity"] * site["occupancy"]
        )
    assert contents == pytest.approx({"Mg": 0.86, "Fe": 3.048, "O": 4.0})
    assert 4.0 - contents["Mg"] - contents["Fe"] == pytest.approx(4 * 0.023)
    assert "average diffraction proxy" in document["notes"]
    assert "high-spin Fe2+ and Fe3+ only at ambient pressure" in document["notes"]


def test_finkelstein_2017_table2_transcription_and_fit_selections():
    _, dataset, path, rows, _ = _document_dataset_rows_and_material()

    assert dataset["identifier"] == DATASET_ID
    assert hashlib.sha256(path.read_bytes()).hexdigest() == DATASET_SHA256
    assert dataset["resource"]["sha256"] == DATASET_SHA256
    assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
    assert len(rows) == 30
    assert sum(int(row["used_as_v0_prior"]) for row in rows) == 1
    assert sum(int(row["used_in_helium_cubic_fit"]) for row in rows) == 17
    assert sum(int(row["used_in_neon_cubic_fit"]) for row in rows) == 5
    assert sum(int(row["used_in_excluded_neon_hex_fit"]) for row in rows) == 7
    assert dataset["used_by_eos_records"] == [HELIUM_ID, NEON_ID]

    assert rows[0]["cubic_a_angstrom"] == "4.2898"
    assert rows[0]["cubic_volume_a3"] == "78.94"
    assert rows[5]["pressure_gpa"] == "16.199"
    assert rows[5]["pressure_standard_deviation_gpa"] == "0.065"
    estimated = next(
        row for row in rows if row["ruby1_observation_status"] == "estimated_from_ruby2"
    )
    assert estimated["pressure_gpa"] == "31.820"
    assert estimated["ruby1_wavelength_nm"] == "705.13"
    missing_pressure_error = next(
        row
        for row in rows
        if row["medium"] == "neon" and row["pressure_gpa"] == "19.244"
    )
    assert missing_pressure_error["pressure_standard_deviation_gpa"] == ""
    assert rows[-1]["hexagonal_volume_a3"] == "47.84"
    assert rows[-1]["hexagonal_c_over_a"] == "2.529"


def test_finkelstein_2017_cubic_records_preserve_parameters_and_priors():
    document, _, _, _, material = _document_dataset_rows_and_material()
    source = {record["identifier"]: record for record in document["eos_records"]}
    executable = {record.identifier: record for record in material.eos_records}

    assert set(source) == set(executable) == {HELIUM_ID, NEON_ID}
    assert source[HELIUM_ID]["eos"]["parameters"] == {
        "V0": 78.87,
        "K0": 148.0,
        "K0_prime": 4.09,
    }
    assert source[NEON_ID]["eos"]["parameters"] == {
        "V0": 78.742,
        "K0": 163.0,
        "K0_prime": 4.02,
    }
    assert source[HELIUM_ID]["parameter_errors"] == {
        "V0": 0.06,
        "K0": 3.0,
        "K0_prime": 0.12,
    }
    assert source[NEON_ID]["parameter_errors"] == {
        "V0": 0.014,
        "K0": 1.0,
        "K0_prime": 0.10,
    }
    assert source[HELIUM_ID]["fixed_parameters"] == []
    assert source[NEON_ID]["fixed_parameters"] == []
    assert source[HELIUM_ID]["parameter_covariance"] is None
    assert source[NEON_ID]["parameter_covariance"] is None
    assert (
        "V0 prior 78.94 +/- 0.1"
        in source[HELIUM_ID]["experimental_configuration"]["fit_constraints"]
    )
    assert (
        "K0_prime prior 4.0 +/- 0.1"
        in source[NEON_ID]["experimental_configuration"]["fit_constraints"]
    )

    for record_id, pressures in (
        (HELIUM_ID, (1.823, 30.0, 55.542)),
        (NEON_ID, (1.314, 10.0, 19.244)),
    ):
        record = executable[record_id]
        assert record.pressure(record.eos.V0) == pytest.approx(0.0, abs=1.0e-12)
        assert record.eos.bulk_modulus(record.eos.V0) == pytest.approx(record.eos.K0)
        for pressure in pressures:
            volume = record.volume(pressure, 300.0, check_validity=True)
            assert record.pressure(volume, 300.0, check_validity=True) == pytest.approx(
                pressure, rel=1.0e-11
            )


def test_finkelstein_2017_curves_reproduce_rows_and_source_density_benchmarks():
    document, _, _, rows, material = _document_dataset_rows_and_material()
    source = {record["identifier"]: record for record in document["eos_records"]}
    executable = {record.identifier: record for record in material.eos_records}

    selections = {
        HELIUM_ID: [row for row in rows if row["used_in_helium_cubic_fit"] == "1"],
        NEON_ID: [row for row in rows if row["used_in_neon_cubic_fit"] == "1"],
    }
    for record_id, selected in selections.items():
        record = executable[record_id]
        residuals = np.array(
            [
                record.pressure(float(row["cubic_volume_a3"]))
                - float(row["pressure_gpa"])
                for row in selected
            ]
        )
        expected = source[record_id]["scientific_validation"]["numerical_reproduction"][
            "published_curve_on_rounded_fit_rows"
        ]
        assert math.sqrt(np.mean(residuals**2)) == pytest.approx(
            expected["pressure_rmse_gpa"], abs=1.0e-10
        )
        assert np.max(np.abs(residuals)) == pytest.approx(
            expected["maximum_absolute_pressure_residual_gpa"], abs=1.0e-10
        )

    # Independent high-pressure benchmarks printed in the Implications section.
    nominal_molar_mass = 0.22 * 24.305 + 0.78 * 55.845 + 15.999
    avogadro = 6.02214076e23
    for record_id, published_density in ((HELIUM_ID, 8.17), (NEON_ID, 8.03)):
        volume = executable[record_id].volume(135.8)
        density = 4 * nominal_molar_mass / (avogadro * volume * 1.0e-24)
        assert density == pytest.approx(published_density, abs=0.008)
        stored = source[record_id]["scientific_validation"]["numerical_reproduction"][
            "published_135_8_gpa_density_benchmark"
        ]
        assert density == pytest.approx(stored["reproduced_value_g_cm3"], abs=1.0e-10)

    actual_molar_mass = 0.215 * 24.305 + 0.762 * 55.845 + 15.999
    actual_density = (
        4
        * actual_molar_mass
        / (avogadro * executable[HELIUM_ID].volume(135.8) * 1.0e-24)
    )
    assert actual_density == pytest.approx(8.0249703654)
    assert (
        "nominal Mg0.22Fe0.78O molar mass"
        in source[HELIUM_ID]["scientific_validation"]["unresolved_issues"][-1]
    )


def test_finkelstein_2017_diagnostic_refits_record_method_dependence():
    document, _, _, rows, _ = _document_dataset_rows_and_material()
    source = {record["identifier"]: record for record in document["eos_records"]}

    helium_rows = [row for row in rows if row["used_in_helium_cubic_fit"] == "1"]
    neon_rows = [row for row in rows if row["used_in_neon_cubic_fit"] == "1"]
    helium_refit = _effective_variance_refit(helium_rows, v0_prior=(78.94, 0.1))
    neon_refit = _effective_variance_refit(
        neon_rows, v0_prior=(78.94, 0.1), k0_prime_prior=(4.0, 0.1)
    )

    for record_id, refit in ((HELIUM_ID, helium_refit), (NEON_ID, neon_refit)):
        expected = source[record_id]["scientific_validation"]["independent_refit"]
        assert refit == pytest.approx(
            [
                expected["refit_parameters"]["V0"],
                expected["refit_parameters"]["K0"],
                expected["refit_parameters"]["K0_prime"],
            ],
            abs=1.0e-5,
        )

    assert (
        source[HELIUM_ID]["scientific_validation"]["independent_refit"]["result"]
        == "near_parity"
    )
    assert (
        source[NEON_ID]["scientific_validation"]["independent_refit"]["result"]
        == "method_dependent_near_parity"
    )


def test_finkelstein_2017_preserves_calibration_and_excludes_unresolved_hex_fit():
    document, _, _, rows, _ = _document_dataset_rows_and_material()
    source = {record["identifier"]: record for record in document["eos_records"]}

    for record in source.values():
        method = record["pressure_calibration"]["methods"][0]
        assert method["reference_calibration_record"] == "ruby_dewaele_2004"
        assert record["pressure_calibration"]["recalculation"]["status"] == (
            "missing_calibrant_observations"
        )

    calibration = get_pressure_calibration("ruby_dewaele_2004")
    first_helium = next(row for row in rows if row["pressure_gpa"] == "1.823")
    universal_reference_result = calibration.pressure_from_wavelength(
        float(first_helium["ruby1_wavelength_nm"])
    )
    assert universal_reference_result == pytest.approx(1.9557081021)
    assert universal_reference_result != pytest.approx(
        float(first_helium["pressure_gpa"]), abs=0.01
    )

    excluded = source[NEON_ID]["scientific_validation"][
        "excluded_hexagonal_parameterization"
    ]
    assert excluded["fit_range_gpa"] == [24.1, 53.3]
    assert excluded["V0_a3"] == 58.7
    assert excluded["K0_gpa"] == 176.8
    assert excluded["K0_prime"] == 4.0
    assert (
        "No space group"
        in source[HELIUM_ID]["scientific_validation"]["reported_parameterizations"][-1][
            "reason"
        ]
    )
    assert all("hex" not in identifier for identifier in source)
