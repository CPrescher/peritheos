"""Historical calibration replay and explicit limits of fit parity."""

import json

import numpy as np
import pytest

from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import MieGruneisenDebye
from scripts.investigate_rubidium_halides_farla_2025 import (
    INPUTS,
    OLD_NA,
    calibration_replay,
    cscl_pressure,
    joint_benchmark_rounding,
    weighted_refit,
)


def test_historical_calibration_and_corrected_evaluator_are_distinct():
    facts = json.loads(INPUTS.read_text(encoding="utf-8"))
    sources = {s["id"]: s for s in facts["sources"]}
    assert sources["cscl_module"]["commit"] == (
        "9518e36112be23504b84c4f595eab0b0207a3539"
    )
    params = sources["cscl_module"]["parameter_dictionaries_si"]["_params_Dewaele_CsCl"]
    assert params["K_0"] == 16.74e9
    assert params["q_0"] == 1.18
    initial = sources["rubidium_initial"]["parameter_dictionaries_si"][
        "_params_Farla_RbCl"
    ]
    corrected = sources["rubidium_corrected"]["parameter_dictionaries_si"][
        "_params_Farla_RbCl"
    ]
    assert initial["K_0"] == 20.46743e9
    assert corrected["K_0"] == 19.89e9
    assert corrected["Debye_0"] == 177.04
    assert facts["cscl_gui_defaults"]["runtime_reference_volume_override"]
    assert not facts["file_recovery"]["original_eosfit_files_found"]
    assert facts["file_recovery"]["no_author_contact"]


def test_cscl_reference_matches_native_vinet_integrated_debye():
    volume = np.array([61.810, 46.959, 59.103])
    temperature = np.array([300, 500, 1800])
    molar_v0 = 4.1221**3 * OLD_NA * 1e-24
    # Native molar volumes are J/bar/mol, equal to 10 cm3/mol.
    eos = MieGruneisenDebye(Vinet(molar_v0 / 10, 16.74, 5.703), 300, 151, 1.99, 1.18, 2)
    np.testing.assert_allclose(
        cscl_pressure(volume, temperature),
        eos.pressure(volume * OLD_NA * 1e-25, temperature),
        atol=1e-12,
        rtol=0,
    )


def test_calibration_replay_preserves_five_discrepancies():
    result = calibration_replay()
    gui = result["gui_defaults_legacy_na"]
    assert gui["unique_rows"] == 109
    assert gui["repeated_source_rows"] == 327
    assert gui["within_pressure_and_volume_rounding"] == 104
    assert gui["within_reported_pressure_error_plus_rounding"] == 109
    assert gui["propagated_sigma_matches_printed_rounding"] == 109
    assert gui["max_abs_difference_gpa"] == pytest.approx(0.07090187135216)
    assert {
        (r["run"], r["temperature_k"]) for r in gui["outliers_beyond_rounding"]
    } == {
        ("BT792", 1100),
        ("BT792", 1200),
        ("BT792", 1300),
        ("BT792", 1400),
        ("BT795", 300),
    }
    assert result["module_legacy_na"]["within_pressure_and_volume_rounding"] == 100


@pytest.mark.parametrize("compound", ["rbcl", "rbbr", "rbi"])
def test_iterated_weights_improve_agreement_without_claiming_original_fit(compound):
    fitted = weighted_refit(compound, author_theta=True, eosfit_r=True)
    assert fitted["converged"]
    assert fitted["pressure_floor_gpa"] == 0
    assert fitted["max_coefficient_difference_in_reported_esd"] < 1.65
    assert not fitted["all_coefficients_round_to_source"]
    assert not fitted["author_input_files_recovered"]
    assert not fitted["author_workflow_parity"]


def test_zero_pressure_error_floor_is_a_material_analysis_choice():
    raw = weighted_refit("rbcl")
    floored = weighted_refit("rbcl", pressure_floor=0.005)
    assert raw["fit_parameters"][0] - floored["fit_parameters"][0] > 0.15
    assert raw["max_coefficient_difference_in_reported_esd"] < 0.5


@pytest.mark.parametrize("compound", ["rbbr", "rbi"])
def test_one_common_vector_can_match_all_rounded_calculated_ratios(compound):
    result = joint_benchmark_rounding(compound)
    assert result["rows"] == 54
    assert result["joint_rounded_table_solution_found"]
    assert result["max_absolute_ratio_difference"] <= 5e-5 + 1e-9
    assert "not original coefficients" in result["qualification"]
