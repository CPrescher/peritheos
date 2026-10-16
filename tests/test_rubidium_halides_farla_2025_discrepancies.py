"""Source-data discrepancy localization without changing primary observations."""

import hashlib
import json

import numpy as np
import pytest

from scripts.audit_rubidium_halides_farla_2025_discrepancies import (
    BURNMAN_REPLAY,
    reproduce,
)
from scripts.investigate_rubidium_halides_farla_2025 import cscl_pressure
from scripts.reproduce_rubidium_halides_farla_2025 import (
    DATA,
    PREFIX,
    SOURCE,
    read_rows,
)


@pytest.fixture(scope="module")
def audit():
    files = list(DATA.glob(f"{PREFIX}-*-observations.csv")) + list(
        DATA.glob(f"{PREFIX}-*-benchmarks.csv")
    )
    before = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    result = reproduce()
    assert before == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return result


def test_executed_historical_code_agrees_with_independent_reference():
    replay = json.loads(BURNMAN_REPLAY.read_text())
    rows = [r for r in read_rows("rbcl", "observations") if r["source_table"] != "S1"]
    reference = cscl_pressure(
        np.array([float(r["cscl_volume_a3"]) for r in rows]),
        np.array([float(r["temperature_k"]) for r in rows]),
    )
    assert replay["burnman_version"] == "1.2.0"
    assert replay["rows"] == 109
    assert replay["source_sha256"] == replay["source"]["sha256"]
    np.testing.assert_allclose(replay["pressures_gpa"], reference, rtol=0, atol=1e-9)


def test_inverse_inputs_separate_rounding_control_and_uncertainty(audit):
    result = audit["cscl_input_diagnostics"]
    assert result["outliers_remaining_with_half_integer_T_rounding"] == 5
    assert result["outliers_remaining_with_target_control_tolerance"] == 4
    outliers = result["outliers"]
    np.testing.assert_allclose(
        [r["minimum_absolute_temperature_shift_k"] for r in outliers],
        [0.54631664, 9.33356896, 15.59285538, 28.55674457, 29.96316706],
        atol=1e-6,
    )
    assert all(abs(r["inverse_volume_shift_a3"]) > 0.0005 for r in outliers)
    assert "not recovered measurements" in result["qualification"]


def test_no_single_reference_lattice_or_temperature_shift_fixes_all_rows(audit):
    result = audit["cscl_input_diagnostics"]
    assert not result["common_a0_possible_at_printed_T"]
    assert not result["common_temperature_offset_possible"]
    low, high = result["common_a0_intersection_angstrom"]
    assert low > high
    assert low == pytest.approx(4.12215169285)
    assert high == pytest.approx(4.11897594830)


def test_only_the_5gpa_1500k_omission_restores_common_vector_consistency(audit):
    result = audit["rbcl_benchmark_diagnostics"]
    assert len(result["leave_one_out"]) == 54
    assert all(r["optimizer_success"] for r in result["leave_one_out"])
    assert result["single_omissions_with_compatible_remaining_table"] == [3]
    local = result["localized_candidate"]
    assert local["retained_rows"] == 53
    assert local["all_retained_within_ratio_rounding"]
    assert (
        max(abs(np.array(local["candidate_parameters"][:4]) - SOURCE["rbcl"][2]))
        <= 0.005 + 1e-10
    )
    missing = local["states"][3]
    assert missing["pressure_gpa"] == 5
    assert missing["temperature_k"] == 1500
    assert missing["source_ratio"] == 0.9282
    assert missing["calculated_ratio_rounded"] == "0.9283"


def test_reference_convention_sensitivities_do_not_restore_literal_table(audit):
    cases = audit["rbcl_benchmark_diagnostics"]["full_table_cases"]
    assert len(cases) == 5
    assert all(r["optimizer_success"] for r in cases.values())
    assert all(r["minimum_scaled_constraint"] > -1e-7 for r in cases.values())
    assert all(r["minimax_half_last_digits"] > 1.07 for r in cases.values())
    assert all(not r["all_retained_within_ratio_rounding"] for r in cases.values())


def test_omitted_ratio_extrema_both_round_to_the_same_alternative(audit):
    result = audit["rbcl_benchmark_diagnostics"]["omitted_prediction_range"]
    low, high = result["minimization"], result["maximization"]
    assert low["optimizer_success"] and high["optimizer_success"]
    assert 0.92825 < low["ratio"] < high["ratio"] < 0.92835
    assert low["minimum_scaled_constraint"] > -1e-7
    assert high["minimum_scaled_constraint"] > -1e-7
    assert "not certified global" in result["qualification"]


def test_one_cell_alternative_is_explicitly_a_diagnostic(audit):
    result = audit["rbcl_benchmark_diagnostics"]
    assert result["one_cell_replacement_diagnostic"][
        "all_retained_within_ratio_rounding"
    ]
    assert not result["source_table_modified"]
    assert not result["literal_published_table_parity"]
    assert not audit["author_fit_parity"]
    assert not audit["original_input_files_recovered"]
    assert audit["author_contact"] == "none"
