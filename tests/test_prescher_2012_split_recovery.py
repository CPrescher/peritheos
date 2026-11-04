"""Printed-digit and source-input boundaries for the split-fit investigation."""

import json

import numpy as np
import pytest

from scripts.recover_prescher_2012_fe3c import ROOT, TABLE, fit
from scripts.recover_prescher_2012_fe3c_splits import (
    STRAIN,
    TARGETS,
    legacy_pressure_replay,
    matches_printed,
    parameter_assessments,
    rounding_feasibility,
)


def test_error_bar_overlap_is_not_exact_printed_digit_recovery():
    assert matches_printed([147.18, 7.937, 3, 0.7], TARGETS["low"]) == [
        False,
        False,
        True,
        True,
    ]
    assert matches_printed([144.762, 8.494, 4.75, 0.88], TARGETS["low"]) == [
        True,
        True,
        False,
        False,
    ]
    assert all(matches_printed([145, 8.5, 3, 0.7], TARGETS["low"]))
    assert not matches_printed([145.5, 8.5, 3, 0.7], TARGETS["low"])[0]
    assert matches_printed([144.5, 8.5, 3, 0.7], TARGETS["low"])[0]


def test_strain_inset_has_19_rows_and_pressure_widths_near_point_two():
    evidence = json.loads(STRAIN.read_text())
    assert (
        evidence["kind"]
        == "digitized_graphic_coordinates_and_error_bars_not_original_inputs"
    )
    rows = evidence["rows"]
    assert [r["matching_table_s1_row"] for r in rows] == list(range(1, 20))
    for row in rows:
        factor = 3 * row["f"] * (1 + 2 * row["f"]) ** 2.5
        assert row["pressure_from_f_F_gpa"] == pytest.approx(factor * row["F_gpa"])
        assert row["pressure_width_if_F_bar_is_pressure_only_gpa"] == pytest.approx(
            0.2, abs=0.005
        )
    assert rows[2]["pressure_from_f_F_gpa"] == pytest.approx(6.8, abs=0.01)


def test_including_22_point_six_recovers_low_central_pair_but_not_errors():
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    previous = json.loads(
        (ROOT / "docs/data/prescher-2012-fe3c-recovery.json").read_text()
    )
    v0 = previous["source_drawn_curve_reconstruction"]["parameters"]["V0"]
    mask = table["pressure_gpa"] < 24
    assert table["source_row"][mask].tolist() == list(range(1, 11))
    assert table["pressure_gpa"][mask][-1] == 22.6
    result = fit(table["pressure_gpa"][mask], table["volume_a3"][mask], fixed_v0=v0)
    pair = [result["parameters"][k] for k in ("K0", "K0_prime")]
    assert pair == pytest.approx([144.7616, 8.494394], rel=1e-5)
    assert result["residual_scaled_errors"] == pytest.approx(
        [4.7494, 0.881408], rel=1e-4
    )
    assert matches_printed(pair + result["residual_scaled_errors"], TARGETS["low"]) == [
        True,
        True,
        False,
        False,
    ]
    p, v = table["pressure_gpa"][mask], table["volume_a3"][mask]
    absolute_case = fit(p, v, "pressure_widths", v0, np.full(len(p), 0.2))
    absolute_errors = np.sqrt(
        np.diag(absolute_case["unscaled_inverse_normal_matrix"])
    ).tolist()
    assert matches_printed(pair + absolute_errors, TARGETS["low"]) == [
        True,
        True,
        False,
        True,
    ]


def test_optimistic_input_rounding_does_not_resolve_low_error_in_this_protocol():
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    previous = json.loads(
        (ROOT / "docs/data/prescher-2012-fe3c-recovery.json").read_text()
    )
    v0 = previous["source_drawn_curve_reconstruction"]["parameters"]["V0"]
    mask = table["pressure_gpa"] < 24
    result = rounding_feasibility(
        table["pressure_gpa"][mask], table["volume_a3"][mask], v0, TARGETS["low"]
    )
    assert result["minimum_rectangle_distance_sum_squares_gpa2"] == pytest.approx(
        0.35546865, rel=1e-5
    )
    # This is explicitly a nominal-J sensitivity estimate, not a global theorem.
    assert result["optimistic_residual_scaled_error_floor_estimates"][0] > 3.5


def test_archived_high_pressure_candidate_matches_four_digits_conditionally():
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    cells = json.loads(
        (
            ROOT
            / "peritheos/data/research/cementite/prescher-2012-own-shelx-cells.json"
        ).read_text()
    )["rows"]
    assert [row["source_row"] for row in cells] == list(range(1, 20))
    # Retain the working-version discrepancy that matters to this candidate.
    assert cells[15]["volume_a3"] == 133.58
    assert table["volume_a3"][15] == 133.25
    volume = np.array([row["volume_a3"] for row in cells])
    previous = json.loads(
        (ROOT / "docs/data/prescher-2012-fe3c-recovery.json").read_text()
    )
    v0 = previous["source_drawn_curve_reconstruction"]["parameters"]["V0"]
    mask = table["pressure_gpa"] >= 22
    p, v = table["pressure_gpa"][mask], volume[mask]
    result = fit(p, v, "pressure_widths", v0, np.full(len(p), 0.1))
    central = [result["parameters"][key] for key in ("K0", "K0_prime")]
    absolute = np.sqrt(np.diag(result["unscaled_inverse_normal_matrix"])).tolist()
    assert all(matches_printed(central + absolute, TARGETS["high"]))
    assert not all(
        matches_printed(central + result["residual_scaled_errors"], TARGETS["high"])
    )
    legacy = legacy_pressure_replay(p, v, v0)
    assert legacy["errors_from_legacy_curve_fit"] == pytest.approx(
        result["residual_scaled_errors"], rel=1e-4
    )
    assert legacy["errors_from_unscaled_leastsq_matrix"] == pytest.approx(
        absolute, rel=1e-4
    )
    assert legacy["reduced_chi_square"] > 9


def test_parameter_parity_keeps_the_unreproduced_uncertainty_separate():
    report = json.loads(
        (ROOT / "docs/data/prescher-2012-fe3c-split-recovery.json").read_text()
    )
    assessments = parameter_assessments(report["cases"])
    for side in ("low", "high"):
        assessment = assessments[side]
        assert assessment["status"] == "parity"
        assert assessment["scope"] == "parameter_agreement"
        assert assessment["compared_parameters"] == ["K0", "K0_prime"]
        assert assessment["conditional"] is True
        assert assessment["qualifications"]
        assert assessment["author_workflow_reproduced"] is False
    low = assessments["low"]["uncertainty_comparisons"]
    assert low[0]["published"] == 3
    assert low[0]["reconstructed"] == pytest.approx(3.67049, rel=1e-5)
    assert low[0]["finding"] == "not_reproduced"
    assert low[1]["matches_printed_precision"] is True
    assert all(
        r["matches_printed_precision"]
        for r in assessments["high"]["uncertainty_comparisons"]
    )
    assert not report["cases"][assessments["low"]["case_index"]]["comparisons"][
        "input_widths_as_absolute_1sigma"
    ]["all_four_match"]
