"""Distinguish atom-count scaling from an energy-offset signature."""

import numpy as np
import pytest

from scripts.ferrous_ringwoodite_brown_hypotheses import (
    investigate,
    series_energy,
    zero_point_pressure,
)
from scripts.ferrous_ringwoodite_brown_replay import (
    MOLAR_MASS,
    Brown1999,
    vibrational_energy,
)


@pytest.fixture(scope="module")
def hypotheses():
    return investigate()


def test_debye_series_independently_distinguishes_linear_term_from_atom_count():
    for temperature in (300, 600, 1200):
        for theta in (214, 375, 446):
            assert series_energy(temperature, theta) == pytest.approx(
                vibrational_energy(temperature, theta), abs=2e-7
            )
    model = Brown1999()
    volume = 0.3709
    _, gamma, theta = model.parameters(volume)
    for temperature in (300, 600, 1200):
        extra_energy = series_energy(temperature, theta, linear=-1 / 4) - series_energy(
            temperature, theta
        )
        assert extra_energy * gamma / (volume * MOLAR_MASS) / 1000 == pytest.approx(
            zero_point_pressure(model, volume) / 3, abs=1e-12
        )
    # Adding this term changes the zero-K reference, not the heat capacity.
    assert model.thermal(volume, 1) < 1e-7
    assert zero_point_pressure(model, volume) / 3 > 0.12


def test_brown_temperature_slopes_support_n2_and_require_an_offset(hypotheses):
    brown = hypotheses["brown_tables"]
    assert brown["n_from_within_volume_temperature_slopes"] == pytest.approx(
        2, abs=0.002
    )
    required = brown["required_n_by_temperature_at_one_volume"]
    assert required[0]["required_n_if_no_offset"] == pytest.approx(2.482322822106086)
    assert required[-1]["required_n_if_no_offset"] == pytest.approx(2.077135687746827)
    assert brown["wrong_n_only"]["1"]["rms_gpa"] > 1
    assert brown["wrong_n_only"]["4"]["rms_gpa"] > 1.8
    assert brown["n_only_fitted_with_fixed_table1_p0"]["rms_gpa"] > 0.04


def test_fixed_zpe_third_signature_and_temperature_holdout(hypotheses):
    brown = hypotheses["brown_tables"]
    assert brown["fixed_one_third_zero_point"]["count"] == 571
    assert brown["fixed_one_third_zero_point"]["rms_gpa"] < 0.011
    assert brown["fixed_one_third_zero_point"]["maximum_abs_gpa"] < 0.033
    assert brown["fitted_constant_pressure_offset"]["rms_gpa"] > 0.019
    assert brown["fitted_constant_molar_energy_offset"]["rms_gpa"] > 0.023
    assert brown["fitted_zero_point_fraction_at_n2"]["fraction"] == pytest.approx(
        1 / 3, abs=0.004
    )
    assert brown["joint_n_and_offset_diagnostic"]["n"] == pytest.approx(2, abs=0.002)
    assert brown["temperature_holdout"]["test"]["rms_gpa"] < 0.013
    assert (
        brown["independent_series_checks"]["normalized_series_minus_quadrature_max_gpa"]
        < 1e-9
    )
    assert (
        brown["independent_series_checks"][
            "altered_linear_term_minus_added_zpe_third_max_gpa"
        ]
        < 1e-9
    )
    assert hypotheses["brown_offset_cause_confirmed"] is False


def test_temperature_hypotheses_do_not_identify_the_original_workflow(hypotheses):
    result = hypotheses["table_s1"]
    summaries = result["accepted_statistics"]
    assert all(r["count"] == 14 for r in summaries.values())
    half = summaries["table1_n2_half_kelvin_heating"]["rms_gpa"]
    reference_error = summaries["table1_n2_midpoint_minus_E300_global"]["rms_gpa"]
    assert half == pytest.approx(0.10548218467206508)
    assert reference_error < 0.114
    assert abs(reference_error - half) < 0.009
    assert summaries["table1_n1_full_sample"]["rms_gpa"] > 0.23
    fits = result["conditional_unweighted_fits"]
    assert fits["temperature_fraction_without_intercept"]["parameters"][
        0
    ] == pytest.approx(0.49044752, abs=1e-6)
    assert fits["n_at_stated_midpoint"]["parameters"][0] == pytest.approx(
        1.63443858, abs=1e-6
    )
    assert hypotheses["historical_temperature_convention_confirmed"] is False
    assert hypotheses["source_pressures_or_eos_coefficients_changed"] is False


def test_all_hypothesis_rows_and_exclusions_preserved(hypotheses):
    rows = hypotheses["table_s1"]["rows"]
    assert len(rows) == 18 * 13
    for name in {r["hypothesis"] for r in rows}:
        assert len([r for r in rows if r["hypothesis"] == name]) == 18
    cold = [
        r
        for r in rows
        if r["measurement"] in ("Fay014quench", "Fay019preheat", "Fay029quench")
    ]
    bad_reference = [
        r for r in cold if r["hypothesis"] == "table1_n2_midpoint_minus_E300_global"
    ]
    assert all(r["included_by_table_s1"] == 0 for r in cold)
    assert all(r["predicted_minus_source_gpa"] < -0.53 for r in bad_reference)


def test_consistent_atom_count_and_volume_basis_are_equivalent():
    energy_per_atom_mole = vibrational_energy(1000, 375) / 2
    molar_volume_formula = 0.3709 * MOLAR_MASS
    formula = 1.19 * 2 * energy_per_atom_mole / molar_volume_formula
    cell = 1.19 * 8 * energy_per_atom_mole / (4 * molar_volume_formula)
    np.testing.assert_allclose(formula, cell, rtol=1e-15)


def test_published_ringwoodite_thermal_curves_support_n7_at_fixed_parameters(
    hypotheses,
):
    curves = hypotheses["ringwoodite_n_curve_check"]["by_n"]
    assert curves["7"]["thermal_checkpoints"] == 9
    assert curves["7"]["within_volume_ratio_tolerance_0_002"] == 9
    for n in ("1", "2", "4"):
        assert curves[n]["within_volume_ratio_tolerance_0_002"] == 0
