"""Independent source checkpoints for the scoped follow-up, not DFT validation."""

import json
from pathlib import Path

import numpy as np
import pytest

from scripts.reproduce_dewaele_2008_followup import (
    REPORT,
    SOURCE,
    input_scale_check,
    reproduce,
    upgraded_pressure,
)
from scripts.reproduce_dewaele_2008_metals import inverse, vinet


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_ruby_mapping_preserves_same_wavelength_ratio():
    # Independent direct power laws at a specified wavelength ratio.
    ratio = 1.04
    p2004 = 1904 / 9.5 * (ratio**9.5 - 1)
    p2008 = 1920 / 9.61 * (ratio**9.61 - 1)
    assert float(upgraded_pressure(p2004)) == pytest.approx(p2008, abs=1e-10)
    assert float(upgraded_pressure(160) - 160) == pytest.approx(1.94745698)


def test_graph_axes_and_legend_recovery_match_table4_shock_comparisons():
    source = json.loads(SOURCE.read_text())
    metals = source["figure2b"]["metals"]
    assert set(metals) == {"Al", "Cu", "Ta", "Au", "Ni", "Ag", "Pt", "Mo", "W", "Co"}
    assert metals["Co"]["object_kind"] == "lines"
    # Independent Table IV constrained DAC / RSW pairs reproduce these paths.
    for metal, k0, dac_kp, shock_kp in [
        ("Ag", 101, 5.97, 6.19),
        ("Mo", 261, 4.06, 4.18),
    ]:
        xy = np.array(metals[metal]["pressure_difference_points_gpa"])
        graph_dp = np.interp(50.0, xy[:, 0], xy[:, 1])
        volume = inverse(50.0, [1, k0, dac_kp])
        expected = float(vinet(volume, [1, k0, shock_kp]) - 50.0)
        assert graph_dp == pytest.approx(expected, abs=0.03)


def test_hcp_fe_table_matches_thin_curve_but_semicore_identity_remains_unresolved(
    report,
):
    fe = report["fe_lda"]
    assert fe["catalog_changes"] is False
    for state, bound in [("static", 0.0025), ("298k", 0.0014)]:
        thin = fe["comparisons"][f"thin_without_semicore_{state}"]
        thick = fe["comparisons"][f"thick_with_semicore_{state}"]
        assert thin["source_curve_index"] == 24
        assert thick["source_curve_index"] == 25
        assert thin["max_absolute_volume_difference_a3_per_atom"] < bound
        assert thick["volume_difference_range_a3_per_atom"][0] > 0.10
        assert not thick["within_0_04_a3_per_atom_graph_reading_bound"]
        assert thin["within_0_04_a3_per_atom_graph_reading_bound"]


def test_printed_pressure_exception_is_preserved_and_helium_is_not_recalibrated():
    result = input_scale_check()
    assert result["non_helium_material_rows"] == 147
    assert result["within_joint_printed_rounding_bounds"] == 145
    assert {r["material"] for r in result["exceptions"]} == {"nickel", "zinc_hcp"}
    for row in result["exceptions"]:
        assert row["mao_pressure_token"] == "49.0"
        assert row["revised_pressure_token"] == "51.0"
        assert row["computed_revised_pressure_gpa"] == pytest.approx(50.1036671)


def test_crosscalibration_reports_assumptions_without_promoting_a_replay(report):
    ruby = report["ruby_crosscalibration"]
    assert ruby["catalog_changes"] is False
    common = ruby["common_range_mean"]
    assert common["pressure_range_gpa"][1] == pytest.approx(62.8474988)
    assert common["A_gpa"] == pytest.approx(1909.48, abs=0.02)
    assert common["B"] == pytest.approx(10.2912, abs=0.0002)
    full = ruby["full_range_graph_continuation"]
    assert full["mean_correction_at_160_gpa"] == pytest.approx(4.05063, abs=0.001)
    assert sum(r["extrapolated_to_160_gpa"] for r in full["models"].values()) == 6
    for row in full["models"].values():
        assert row["graph_fit_max_residual_gpa"] < 0.03
    for key in [
        "common_range_mean",
        "individual_ranges_pooled",
        "full_range_graph_continuation",
    ]:
        assert ruby[key]["solver_success"]
        assert ruby[key]["author_workflow_reproduced"] is False
        assert ruby[key]["coefficient_uncertainties"] is None


def test_saved_numerical_report_matches_executable_analysis(report):
    from scripts.validate_primary_eos_refits import check_ledger_values

    saved = json.loads(Path(REPORT).read_text())
    # These three diagnostic optimizations vary at optimizer/BLAS roundoff
    # across supported Python/SciPy builds. Keep source fields and assessments
    # exact; apply the ledger's calculated-value tolerance only to fit outputs.
    current = json.loads(json.dumps(report))
    for key in (
        "common_range_mean",
        "individual_ranges_pooled",
        "full_range_graph_continuation",
    ):
        check_ledger_values(
            saved["ruby_crosscalibration"].pop(key),
            current["ruby_crosscalibration"].pop(key),
            calculated=True,
        )
    check_ledger_values(saved, current)
