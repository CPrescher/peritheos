"""Check reconstruction agreement, phase holdouts, SI units and provenance."""

import json
from dataclasses import replace

import numpy as np
import pytest

from scripts.fit_yokoo_2009_gold_thermal import (
    OUTPUT,
    PUBLISHED,
    Parameters,
    check_reconstruction,
    independent_quad_pressure,
    pressure_gpa,
    reconstruct,
)
from scripts.reproduce_yokoo_2009_gold import V0, phonon_pressure, source_bm3


def test_primary_fit_recovers_near_published_coefficients_and_cold_volume():
    report = reconstruct()
    fit = report["fits"]["primary"]
    p = fit["parameters"]
    assert p["k0_gpa"] == 180
    assert p["theta0_k"] == 170
    assert p["k0_prime"] == pytest.approx(5.61, abs=0.001)
    assert p["a"] == pytest.approx(0.45, abs=0.003)
    assert p["b"] == pytest.approx(4.2, abs=0.02)
    assert p["gamma0"] == pytest.approx(2.96, rel=0.013)
    assert p["vc_over_v0"] == pytest.approx(0.99017812, abs=1e-7)
    assert fit["unmarked_states"]["rmse_gpa"] < 0.009
    assert fit["unmarked_states"]["max_abs_residual_gpa"] < 0.022
    baseline = report["fits"]["printed_coefficients_cold_volume_only"]
    assert baseline["unmarked_states"]["rmse_gpa"] > 0.08
    assert report["cold_volume_cell_a3"] == pytest.approx(p["vc_over_v0"] * V0)


def test_phase_annotations_and_missing_states_remain_outside_primary_fit():
    report = reconstruct()
    fit = report["fits"]["primary"]
    assert len(report["states"]) == 162
    assert fit["fit_states"]["states"] == 156
    assert fit["first_liquid_states"]["states"] == 6
    assert fit["first_liquid_states"]["max_abs_residual_gpa"] < 0.028
    for row in report["states"]:
        assert row["used_in_primary_fit"] == (
            row["source_phase_annotation"] == "unmarked"
        )
    assert report["fits"]["include_first_liquid_states"]["fit_states"]["states"] == 162


def test_volume_holdout_and_different_starts_agree():
    report = reconstruct()
    holdout = report["volume_interpolation_holdout"]
    assert holdout["training_states"]["states"] == 82
    assert holdout["held_out_states"]["states"] == 74
    assert holdout["held_out_states"]["rmse_gpa"] < 0.01
    assert holdout["held_out_states"]["max_abs_residual_gpa"] < 0.02
    primary = report["fits"]["primary"]["parameters"]
    for start in report["multistart"]:
        for name, value in primary.items():
            assert start["fitted"][name] == pytest.approx(value, abs=3e-7)


def test_pressure_normalization_absolute_phonon_and_adaptive_integral():
    report = reconstruct()
    assert report["adaptive_quad_max_difference_gpa"] < 1e-10
    p = Parameters(**report["fits"]["primary"]["parameters"])
    assert pressure_gpa(p.vc_over_v0, 0, p) == pytest.approx(0, abs=1e-10)
    for ratio in [0.6, 0.8, 1.0]:
        for temperature in [0, 50, 300, 1250, 3000]:
            assert float(pressure_gpa(ratio, temperature, p)) == pytest.approx(
                independent_quad_pressure(ratio, temperature, p), abs=1e-10
            )
        # Electronic pressure is zero at 300 K; this isolates the absolute
        # phonon energy and checks four-atom cell -> SI molar-volume conversion.
        cold = source_bm3(ratio, PUBLISHED.vc_over_v0, 180, 5.61)
        assert pressure_gpa(ratio, 300, PUBLISHED) - cold == pytest.approx(
            phonon_pressure(ratio, 300), abs=1e-10
        )
    assert pressure_gpa(1, 1e-12, p) == pytest.approx(pressure_gpa(1, 0, p), abs=1e-10)


def test_electronic_interpolation_is_explicit_and_bounded():
    # Eliminate phonon pressure to isolate corrected electronic table values.
    p = replace(PUBLISHED, gamma0=0)
    cold = source_bm3(1, p.vc_over_v0, p.k0_gpa, p.k0_prime)
    assert pressure_gpa(1, 1000, p) - cold == pytest.approx(0.01)
    midpoint = (pressure_gpa(1, 1200, p) + pressure_gpa(1, 1300, p)) / 2
    assert pressure_gpa(1, 1250, p) == pytest.approx(midpoint)
    for ratio, temperature in [
        (0.59, 1000),
        (1.01, 1000),
        (0.8, -1),
        (0.8, 3001),
        (np.nan, 300),
        (0.8, np.inf),
    ]:
        with pytest.raises(ValueError, match="Reconstruction requires"):
            pressure_gpa(ratio, temperature, p)


def test_saved_report_does_not_claim_original_fit_or_experimental_errors():
    report = reconstruct()
    check_reconstruction(json.loads(OUTPUT.read_text()), report)
    assert report["author_fit_reproduced"] is False
    assert report["experimental_observations"] == 0
    assert report["kind"] == "derived_table_output_pressure_reconstruction"
    conventions = report["conventions"]
    assert conventions["electronic_interpolation_author_verified"] is False
    assert conventions["domain"]["phase_stability_verified"] is False
    assert "No parameter errors or covariance" in conventions["uncertainty"]
    assert report["parameter_comparison"][0]["published"] is None
