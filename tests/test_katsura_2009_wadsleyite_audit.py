"""Source consistency and normalization checks for the scoped wadsleyite audit."""

import hashlib
import json

import numpy as np
import pytest
from scipy.integrate import quad

from scripts.audit_katsura_2009_wadsleyite import (
    DATA,
    ROOT,
    SOURCE,
    N,
    R,
    audit,
    cold_pressure,
    energy,
    load_data,
    pressure,
)


@pytest.fixture(scope="module")
def result():
    return audit()


def test_all_source_rows_errors_and_original_catalog_are_preserved():
    data = load_data()  # Verifies every value and error against PDF tokens.
    assert len(data["pressure_gpa"]) == 85
    assert hashlib.sha256(DATA.read_bytes()).hexdigest() == (
        "b0fba67929425c9f1c4c51137453b8a114d754cfd6ecd31280ebe89671538310"
    )
    document = json.loads(
        (ROOT / "peritheos/data/materials/wadsleyite.eosmat").read_text()
    )
    record = document["eos_records"][0]
    assert record["parameters"]["V0"]["value"] == 538.49
    assert record["thermal"]["parameters"]["n"]["value"] == 7
    assert record["thermal"]["parameters"]["gamma0"]["value"] == 1.64
    assert record["evidence"]["scientific_validation"]["status"] == "deferred"
    assert document["recommendation_hold"]["status"] == "on_hold"
    assert len(document["eos_records"]) == 1


def test_independent_debye_integral_and_classical_limit():
    theta, t = 814.0, 500.0
    y = theta / t
    integral = quad(lambda u: u**3 / np.expm1(u), 0, y, epsabs=1e-12)[0]
    assert energy(theta, t) == pytest.approx(9 * N * R * t * integral / y**3, rel=1e-12)
    assert energy(theta, 1e10) / 1e10 == pytest.approx(3 * N * R, rel=1e-7)


def test_reference_temperature_and_formula_cell_basis_are_consistent():
    ratios = np.array([0.915, 0.95, 1.0])
    assert pressure(ratios, 300) == pytest.approx(cold_pressure(ratios), abs=1e-12)
    assert pressure(ratios, 1700, n=56, z=1) == pytest.approx(
        pressure(ratios, 1700, n=7, z=8), abs=1e-12
    )


def test_cold_control_and_thermal_discrepancy_are_distinct(result):
    assert result["cold_stage"]["rows"] == 11
    # Source K'=4.1(1) is recovered when the near-ambient rows are treated as RT.
    assert result["cold_stage"]["k0_prime_when_assigned_300k"] == pytest.approx(
        4.1, abs=0.1
    )
    normal = result["variants"]["physical_n7_integrated"]
    assert normal["near_rt_rows"]["rms_gpa"] < 0.1
    assert normal["heated_rows"]["minimum_gpa"] > 0.32
    assert normal["all_rows"]["rms_gpa"] == pytest.approx(1.8944316914, abs=1e-8)
    for key in [
        "physical_n7_variable_exponent",
        "source_notation_cold_reference_gamma_n7",
    ]:
        assert result["variants"][key]["all_rows"]["rms_gpa"] > 1.5
    assert result["implementation_checks"]["native_available"]
    for key in [
        "native_maximum_pressure_difference_gpa",
        "python_maximum_pressure_difference_gpa",
    ]:
        assert result["implementation_checks"][key] < 1e-8


def test_source_figure_checkpoint_does_not_depend_on_recalibrating_observations(result):
    figure = result["figure2"]
    source = json.loads(SOURCE.read_text())["figure2"]
    assert source["kind"] == "digitized_source_calculation_not_observation"
    assert figure["source_ratios"][-1] == pytest.approx(1.0359056, abs=3e-7)
    assert figure["comparisons"]["n7"]["calculated_ratios"][-1] > 1.052
    assert (
        figure["comparisons"]["n7"]["maximum_absolute_ratio_difference"]
        > 50 * figure["graphical_bound"]
    )
    # The attractive normalization shortcut still fails the independent curve check.
    assert not figure["comparisons"]["n5"]["within_graphical_bound"]
    assert (
        result["variants"]["unphysical_n5_diagnostic_only"]["all_rows"]["rms_gpa"]
        < 0.21
    )
    assert not result["pressure_calibration"]["actual_recalibration_performed"]


def test_missing_joint_inputs_are_not_claimed_as_author_workflow_parity(result):
    assert any("Inoue" in value for value in result["missing_original_fit_inputs"])
    fit = result["conditional_high_pressure_fits"]["heated_74"]
    assert fit["success"] and fit["active_bounds"] == [0, 0]
    assert fit["residuals"]["rows"] == 74
    covariance = np.asarray(fit["conditional_covariance"])
    assert np.all(np.linalg.eigvalsh(covariance) > 0)
    assert "high-pressure data only" in fit["objective"]
    assert result["normalization_fingerprint"]["effective_n"] == pytest.approx(
        4.83421153, abs=1e-7
    )


def test_later_author_pressure_routine_preserves_historical_scope(result):
    check = result["later_author_code_crosscheck"]
    assert check["source_doi"] == "10.5281/zenodo.5903286"
    assert not check["original_2009_program_recovered"]
    assert check["source_2022_wadsleyite_n"] == 7
    assert check["source_2022_wadsleyite_z"] == 8
    assert check["source_2009_coefficients_residuals"]["rms_gpa"] == pytest.approx(
        1.8944462899, abs=1e-8
    )
    assert check["maximum_pressure_difference_from_exact_quadrature_gpa"] < 0.00012


def test_separate_htbm_expansivity_is_not_an_exact_mgd_benchmark(result):
    note = result["ambient_source_calculation"]["n7"]["note"]
    assert "separate HTBM" in note
    assert "not a decisive MGD benchmark" in note
