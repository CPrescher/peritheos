"""Independent refit regression, original-coordinate residuals and limitations."""

import numpy as np
import pytest

from peritheos import get_material, get_material_document
from peritheos.eos.rt import BM3
from scripts.refit_chen_2010_argon import (
    DATASET,
    MASS_PER_CELL,
    RECORD,
    observations,
    pressure,
    reproduce,
    spacing_weights,
)


@pytest.fixture(scope="module")
def result():
    return reproduce()


def test_refit_matches_independent_solver_and_preserves_original_status(result):
    doc = get_material_document("argon_fcc")
    record = next(r for r in doc["eos_records"] if r["identifier"] == RECORD)
    assert record["record_kind"] == "refit"
    assert not record["default"]
    assert record["fit_datasets"] == [DATASET]
    assert record["parameter_covariance"] is None
    assert all(v is None for v in record["parameter_errors"].values())
    assert (
        record["scientific_validation"]["original_publication_reproduction_status"]
        == "not_reproduced"
    )
    for k, value in result["primary"]["parameters"].items():
        assert record["eos"]["parameters"][k] == pytest.approx(value, rel=1e-7)
    assert result["multistart_max_parameter_difference"] < 1e-5
    assert not result["published_curve_comparison"]["used_in_fit"]
    assert result["primary"]["pressure_rms_original_coordinates_gpa"] == pytest.approx(
        0.111198575, abs=1e-6
    )


def test_native_evaluation_inverse_and_positive_bulk_over_data_range(result):
    record = get_material("argon_fcc").get_eos_record(RECORD)
    p, rho, _, _ = observations()
    np.testing.assert_allclose(
        record.pressure(MASS_PER_CELL / rho, 290),
        pressure(rho, result["primary"]["density_parameters"]),
        atol=2e-12,
    )
    grid = np.linspace(min(p), max(p), 30)
    v = record.volume(grid, 290)
    np.testing.assert_allclose(record.pressure(v, 290), grid, atol=2e-8)
    assert np.all(np.diff(v) < 0)
    assert np.all(BM3(**result["primary"]["parameters"]).bulk_modulus(v) > 0)
    assert result["native_pressure_max_difference_gpa"] < 1e-11


def test_objective_uses_both_coordinates_and_reports_unadjusted_residual(result):
    p, rho, sp, sr = observations()
    f = result["primary"]
    adjusted = np.array(f["adjusted_density_g_cm3"])
    expected = np.sum(
        ((pressure(adjusted, f["density_parameters"]) - p) / sp) ** 2
        + ((adjusted - rho) / sr) ** 2
    )
    assert f["objective"] == pytest.approx(expected)
    assert (
        f["pressure_rms_original_coordinates_gpa"]
        > 5 * f["pressure_rms_adjusted_coordinates_gpa"]
    )
    assert len(f["row_weights"]) == 80
    assert all(w == 1 for w in f["row_weights"])
    for group in np.unique(np.floor(p)):
        assert sum(spacing_weights(p, 1)[np.floor(p) == group]) == pytest.approx(1)


def test_sensitivity_is_limited_in_range_but_not_claimed_as_confidence(result):
    assert (
        result["weighting_sensitivity"]["max_volume_difference_from_primary_percent"]
        < 0.7
    )
    assert (
        result["coordinate_shift_sensitivity"][
            "max_volume_difference_from_primary_percent"
        ]
        < 0.6
    )
    assert (
        result["deletion_sensitivity"]["max_volume_difference_from_primary_percent"]
        < 1.8
    )
    assert result["deletion_sensitivity"]["parameter_max"]["K0"] > 6
    assert (
        sum(
            f["withheld_positions"]
            for f in result["leave_one_pressure_bin_out"].values()
        )
        == 80
    )
    assert (
        max(
            f["withheld_pressure_rms_gpa"]
            for f in result["leave_one_pressure_bin_out"].values()
        )
        < 0.3
    )
