"""Compare the actual native evaluator and fitting API with archived EosFit data."""

import hashlib
import json

import numpy as np
import pytest

from scripts.compare_fes_peritheos_eosfit import OUTPUT, ROOT, reproduce
from scripts.run_fes_eosfit_console import data_lines


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_same_observations_and_explicit_volume_conversion(report):
    for path, digest in report["input_sha256"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == digest
    for selection in ["cold", "morard", "combined"]:
        archived = ROOT / f"docs/data/fes-eosfit7c/{selection}-unit/input.dat"
        actual = np.loadtxt(archived, skiprows=5)
        expected = np.array(
            [list(map(float, line.split())) for line in data_lines(selection)[5:]]
        )
        np.testing.assert_array_equal(actual, expected)
    assert report["catalog_adapter_max_pressure_difference_gpa"] < 1e-10


def test_native_pressure_agreement_and_identical_failed_bound_rows(report):
    assert len(report["cases"]) == 8
    for case in report["cases"].values():
        for stage in case.values():
            # Console coefficients/residuals are rounded and the implementations
            # use different gas constants. This is an empirical external-data
            # tolerance, not a demand for bitwise-identical calculations.
            assert stage["max_pressure_difference_gpa"] < 0.0012
            assert stage["max_native_python_difference_gpa"] < 1e-10
            assert (
                stage["peritheos"]["rows_exceeding_3_gpa"]
                == stage["eosfit"]["rows_exceeding_3_gpa"]
            )
            assert stage["peritheos"]["rows"] == stage["eosfit"]["rows"]
            assert (
                abs(stage["peritheos"]["rmse_gpa"] - stage["eosfit"]["rmse_gpa"])
                < 0.0003
            )


def test_joint_equal_weight_refit_recovers_coefficients_and_covariance(report):
    fit = report["independent_joint_equal_weight_refit"]
    assert fit["solver"]["success"]
    assert fit["degrees_of_freedom"] == 155
    assert fit["covariance_parameter_order"] == ["V0", "K0", "Kp", "Gamm0"]
    covariance = np.array(fit["covariance"])
    assert np.linalg.matrix_rank(covariance) == 4
    assert np.all(np.linalg.eigvalsh(covariance) > 0)
    for i, key in enumerate(fit["covariance_parameter_order"]):
        actual = fit["parameters_cm3_mol_gpa"][key]
        expected = fit["eosfit_parameters_cm3_mol_gpa"][key]
        assert actual["value"] == pytest.approx(expected["value"], rel=2e-4)
        assert actual["esd"] == pytest.approx(expected["esd"], rel=2e-4)
        assert actual["esd"] == pytest.approx(np.sqrt(covariance[i, i]), rel=1e-12)
    assert fit["all_rows"]["rmse_gpa"] == pytest.approx(1.31720292, abs=1e-5)
    assert len(fit["thermal_rows"]["rows_exceeding_3_gpa"]) == 4


def test_retained_report_matches_recomputed_scientific_results(report):
    retained = json.loads(OUTPUT.read_text())
    assert retained["input_sha256"] == report["input_sha256"]
    for name, case in report["cases"].items():
        for stage, result in case.items():
            previous = retained["cases"][name][stage]
            np.testing.assert_allclose(
                [r["peritheos_pressure_gpa"] for r in result["rows"]],
                [r["peritheos_pressure_gpa"] for r in previous["rows"]],
                rtol=0,
                atol=1e-10,
            )
    np.testing.assert_allclose(
        report["independent_joint_equal_weight_refit"]["covariance"],
        retained["independent_joint_equal_weight_refit"]["covariance"],
        rtol=1e-4,
        atol=1e-8,
    )
