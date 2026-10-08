"""Check the recovered source polynomial and zero-pressure anchor convention."""

import pytest

from scripts.birch_1986_nacl import (
    WALKER_B2_ANCHORS,
    normalized_pressure,
    pressure_ratio,
    table5_check,
)


def test_final_polynomial_matches_all_printed_table5_reference_pressures():
    rows = table5_check()
    assert len(rows) == 8
    assert max(abs(row["difference_gpa"]) for row in rows) < 0.0005
    # Using rounded K0'=5.20 to reconstruct a=1.80 fails this tolerance.
    assert pressure_ratio(0.65) == pytest.approx(28.351, abs=0.0005)


def test_equation8_uses_birch_reference_and_additive_temperature_term():
    assert pressure_ratio(1.0, 298.15) == 0.0
    assert pressure_ratio(0.8, 773.15) - pressure_ratio(0.8) == pytest.approx(1.3585)


@pytest.mark.parametrize("anchor", WALKER_B2_ANCHORS)
def test_each_conditional_run_anchor_is_zero_at_its_actual_temperature(anchor):
    assert normalized_pressure(
        anchor["nacl_lattice_a_angstrom"], anchor["temperature_k"], anchor
    ) == pytest.approx(0.0, abs=1e-10)
