"""Check the source-increment audit independently of table-reconstruction fits."""

import json

import numpy as np
import pytest

from scripts.audit_yokoo_2009_pressure_conventions import OUTPUT, SOURCES, audit, phonon
from scripts.fit_yokoo_2009_gold_thermal import check_reconstruction


def test_archived_convention_audit_reproduces_and_preserves_source_states():
    report = audit()
    check_reconstruction(json.loads(OUTPUT.read_text()), report)
    assert report["published_analytical_pvt_reproduced"] is False
    for metal, count in [("gold", 162), ("platinum", 168)]:
        m = report["metals"][metal]
        assert len(m["states"]) == count
        assert m["missing_source_cells_added"] == 0
        assert m["all_source_temperatures_are_electronic_nodes"]
        assert m["adaptive_quad_max_difference_gpa"] < 1e-10


def test_ambient_rounding_intervals_reject_the_printed_pt_thermal_increment():
    rows = audit()["metals"]["platinum"]["ambient_interval_checks"]
    row = next(r for r in rows if r["temperature_k"] == 1500)
    # This comparison includes pressure-table, electronic-table, phonon and
    # density half-last-digit intervals. It does not use fitted corrections.
    assert row["interval_gap_gpa"] > 0.27
    assert (
        row["source_rounding_interval_gpa"][0]
        > row["printed_model_rounding_interval_gpa"][1]
    )
    assert row["source_electronic_pressure_gpa"] == 0.51


def test_cold_curve_and_constant_reference_offsets_cancel_in_thermal_increment():
    source = SOURCES["platinum"]
    v, t = np.array([0.6, 0.8, 1.0]), np.array([300, 1500, 3000])
    energy = phonon(v, t, source)
    cold = np.array([527.51, 108.48, -1.76])
    pel = np.array([0.04, 0.51, 1.64])
    for cold_offset in [0, 15, -35]:
        for electronic_reference in [0, 0.04]:
            total = cold + cold_offset + energy + pel - electronic_reference
            at_zero = cold + cold_offset - electronic_reference
            assert total - at_zero == pytest.approx(energy + pel, abs=1e-12)


def test_ambient_phonon_terms_do_not_depend_on_a_or_b():
    source = SOURCES["platinum"]
    t = np.array([300, 500, 1000, 1500, 2000, 2500, 3000])
    for a, b in [(0, 0.2), (0.1, 2), (0.8, 10)]:
        assert phonon(1, t, source, parameters=[2.63, a, b, 230]) == pytest.approx(
            phonon(1, t, source), abs=1e-12
        )


def test_electronic_polynomial_and_cold_normalization_do_not_close_pt_gap():
    pt = audit()["metals"]["platinum"]
    assert (
        pt["best_within_assumed_last_digit_rounding"]["residuals"][
            "max_abs_difference_gpa"
        ]
        > 0.29
    )
    assert pt["cold_normalized_phonon_alternative"]["max_abs_difference_gpa"] > 0.3
    assert (
        pt["pytheos_electronic_polynomial_increment_alternative"]["residuals"][
            "max_abs_difference_gpa"
        ]
        > 0.3
    )
