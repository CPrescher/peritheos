"""Preserve the scoped Sata/Speziale pressure-reduction findings."""

import csv
import json

import numpy as np
import pytest

from scripts.audit_sata_2010_fe3c import observations
from scripts.audit_sata_2010_mgo_pressure import POINTS, REPORT, reproduce


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_original_source_and_formulation_controls(report):
    source = report["source"]["preferred_BM3"]
    assert source["V0_a3"] == 74.71
    assert source["K0_fixed_gpa"] == 160.2
    assert source["Kprime"] == 3.99
    cases = report["cases"]
    assert cases["literal_BM3"]["max_abs_gpa"] == pytest.approx(0.352638989)
    assert cases["source_Vinet_Kprime_4_12_adopted_V0"]["max_abs_gpa"] == pytest.approx(
        4.985425126
    )
    assert cases["source_unconstrained_BM3"]["max_abs_gpa"] > 12
    assert report["thermal_controls_at_300K_max_difference_gpa"] < 1e-10
    archived = json.loads(REPORT.read_text())
    for name, case in cases.items():
        assert case["predicted_gpa"] == pytest.approx(
            archived["cases"][name]["predicted_gpa"], abs=1e-10
        )


def test_rounding_is_distinct_from_uncertainty_and_fitted_reducer(report):
    rounding = report["rounding"]
    assert rounding["fixed_coefficient_data_rounding_max_gpa"] < 0.014
    assert rounding["all_targets_in_separate_row_coefficient_boxes"] is True
    common = rounding["common_coefficient_set"]
    assert common["remaining_max_interval_gap_gpa"] == pytest.approx(
        0.011560308, abs=1e-7
    )
    assert common["multistart_gap_spread_gpa"] < 1e-7
    diagnostic = report["Kprime_only_diagnostic"]
    assert diagnostic["value"] == pytest.approx(3.99680594)
    assert diagnostic["inside_source_error_width"] is True
    assert diagnostic["inside_nearest_rounding_half_digit"] is False
    assert diagnostic["all_rows_overlap_P_V_display_rounding"] is True
    assert (
        np.max(np.abs(report["cases"]["Kprime_only_diagnostic"]["residual_gpa"])) < 0.01
    )


def test_recalculated_sample_coordinates_keep_source_pairs_and_pressure_bases(report):
    derived = report["recalculated_fe3c_coordinates"]
    rows = derived["rows"]
    original = observations()[1]
    assert derived["kind"] == "recalculated_pressure_coordinates"
    assert len(rows) == len(original) == 9
    for row, source in zip(rows, original):
        assert row["measurement_id"] == source["measurement_id"]
        assert row["run"] == source["run"]
        assert row["fe3c_volume_a3"] == source["volume_a3"]
        assert row["published_pressure_gpa"] == source["pressure_gpa"]
        assert (
            row["published_pressure_reported_width_gpa"] == source["error_pressure_gpa"]
        )
        # Independently write the printed BM3 reduction in finite strain.
        f = ((74.71 / source["mgo_volume_a3"]) ** (2 / 3) - 1) / 2
        expected = 3 * 160.2 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (3.99 - 4) * f)
        assert row["recalculated_pressure_gpa"] == pytest.approx(expected, abs=1e-10)
        assert row["recalculated_pressure_residual_gpa"] == pytest.approx(
            row["published_pressure_residual_gpa"]
            - row["recalculated_minus_published_pressure_gpa"],
            abs=1e-10,
        )
    assert derived["published_curve_on_original_nine"]["rmse_gpa"] == pytest.approx(
        1.255813063
    )
    assert derived["published_curve_on_recalculated_nine"]["rmse_gpa"] == pytest.approx(
        1.225288935
    )
    assert derived["published_curve_on_recalculated_nine"]["max_abs_gpa"] > 2.47
    with POINTS.open(newline="") as stream:
        archived = list(csv.DictReader(stream))
    assert len(archived) == 9
    assert [float(r["recalculated_pressure_gpa"]) for r in archived] == pytest.approx(
        [r["recalculated_pressure_gpa"] for r in rows], abs=1e-12
    )
