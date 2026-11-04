"""Check the source replay and keep the two residual-coordinate meanings separate."""

import json
from collections import Counter

import numpy as np
import pytest

from scripts.reproduce_frost_2023_palladium import (
    OUTPUT,
    SOURCE_ROWS,
    observations,
    replay,
)


def test_complete_source_transcription_and_independent_fit_replay():
    p, v, sp, sv = observations()
    assert len(p) == 93
    assert (min(p), max(p)) == (2.99, 181.79)
    assert np.all(sv > 0) and np.all(sp == 0.1)
    rows = json.loads(SOURCE_ROWS.read_text())["rows"]
    assert Counter(r["loading"] for r in rows) == {"BX18": 26, "BX32": 20, "HP": 47}
    assert len(v) == len(rows)
    report = replay()  # Includes literal equations/native and independent ODR checks.
    saved = json.loads(OUTPUT.read_text())
    for kind, expected in saved["models"].items():
        current = report["models"][kind]
        assert current["orthogonal_refit"] == pytest.approx(
            expected["orthogonal_refit"], rel=2e-7
        )
        assert current["refit_pressure_rms_gpa"] > current["published_pressure_rms_gpa"]
        assert (
            current["orthogonal_objective"]
            < current["published_fixed_orthogonal_objective"]
        )
        assert current["within_combined_two_error_widths"] == [True, True, False]
        assert current["within_combined_two_error_widths_residual_scaled"] == [
            True,
            True,
            True,
        ]
        assert current["reduced_chi_square_refit"] > 19
        points = current["row_diagnostics"]
        measured_rms = np.sqrt(
            np.mean(
                [
                    point["refit_pressure_residual_at_measured_volume_gpa"] ** 2
                    for point in points
                ]
            )
        )
        joint_objective = sum(
            (point["refit_pressure_residual_at_adjusted_volume_gpa"] / 0.1) ** 2
            + (
                (
                    point["refit_adjusted_volume_cell_a3"]
                    - point["observed_volume_cell_a3"]
                )
                / point["volume_error_cell_a3"]
            )
            ** 2
            for point in points
        )
        assert measured_rms == pytest.approx(current["refit_pressure_rms_gpa"])
        assert joint_objective == pytest.approx(current["orthogonal_objective"])
