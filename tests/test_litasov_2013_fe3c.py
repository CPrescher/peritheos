"""Scientific evidence and normalization boundaries for the Fe3C replay."""

import json
import shutil

import numpy as np
import pytest
from scipy.constants import Boltzmann
from scipy.integrate import quad

from scripts.reproduce_litasov_2013_fe3c import (
    CASES,
    DATA,
    PDF,
    TABLE3,
    TOKENS,
    checkpoints,
    comparison_inventory,
    debye_energy,
    fit,
    ke_state,
    library_check,
    selection,
    sha256,
    source_audit,
)


def test_preserved_observations_and_literal_tokens():
    assert (
        sha256(DATA)
        == "e0017d735a5f866053caca0e00fcb55beccb9d55a4e2fe2cbe30f25656bc5789"
    )
    assert (
        sha256(TOKENS)
        == "8171286e89547b6684b153a1f3f6a8f7f3f3a4624e2c12e66405ebb844643994"
    )
    tokens = json.loads(TOKENS.read_text())
    assert len(tokens) == 81
    assert [r["row"] for r in tokens if r["tokens"][1] is None] == [29, 39, 47]
    assert tokens[-1]["tokens"] == [
        "1073",
        "70.44(3)",
        "15.3(1)",
        "5.033",
        "6.654",
        "4.437",
        "148.59(2)",
    ]


@pytest.mark.skipif(
    not PDF.exists() or not shutil.which("pdftotext"),
    reason="Local primary PDF/Poppler not available",
)
def test_all_81_rows_against_final_pdf():
    audit = source_audit()
    assert audit["all_pdf_tokens_and_csv_values_match"]
    assert list(audit["run_counts"].values()) == [29, 9, 8, 11, 8, 11, 5]
    assert audit["volume_molar_cm3_per_formula"] == pytest.approx(23.26955, abs=1e-5)
    assert audit["max_abs_abc_minus_printed_volume_a3"] < 0.005


@pytest.mark.parametrize(
    "theta,temp", [(200, 300), (314, 873), (500, 1473), (800, 473), (1000, 300)]
)
def test_cell_debye_energy_against_adaptive_quadrature(theta, temp):
    y = theta / temp
    integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-12)[0]
    independent = 9 * 16 * Boltzmann * temp * integral / y**3
    assert debye_energy(theta, temp) == pytest.approx(independent, rel=2e-13, abs=1e-32)


def test_printed_alternatives_and_lost_minus_signs_stay_separate():
    assert CASES["fm_table2"][0] == 155.2
    assert CASES["fm_text"][0] == 155.27
    assert CASES["ht_table2"][0] == 154.56
    assert CASES["ht_text"][0] == 154.5
    assert CASES["ht_text"][-1] == -0.029
    assert CASES["mgd_abstract"][-2:] == [-0.1, 490]
    assert CASES["mgd_body"][-2:] == [-0.03, 314]
    assert CASES["mgd_fixed500_conditional"][-2:] == [-0.36, 500]
    assert len(TABLE3) == 25  # calculated benchmarks, separate from measured rows
    assert TABLE3[1, -1] == -26.775


def test_complete_ke_thermodynamic_checkpoints_and_printing_envelope():
    result = checkpoints()
    for case in ["ke_text", "ke_table2"]:
        assert result[case]["table3_within_linear_printing_envelope"]
        assert result[case]["max_absolute_difference"]["Cv_J_mol_K"] < 0.012
        assert result[case]["max_absolute_difference"]["gamma_oscillator"] < 0.0042
    assert result["table4"]["max_abs_x_difference"] < 0.00021
    # External printed entropy/Cv, including electrons; not a reference identity.
    state = ke_state(30, 1000)
    assert state[2] == pytest.approx(207.29, abs=0.04)
    assert state[4] == pytest.approx(107.23, abs=0.005)
    assert state[8] == pytest.approx(545.944, abs=0.03)


def test_direct_cell_equations_match_public_python_and_native_rust():
    for comparison in library_check().values():
        assert comparison["max_abs_difference_gpa"] < 1e-8
        if "max_abs_python_native_difference_gpa" in comparison:
            assert comparison["max_abs_python_native_difference_gpa"] < 1e-8


def test_conditional_refit_agreement_is_not_full_parameter_or_workflow_parity():
    assert len(selection("pm_conservative")) == 63
    assert len(selection("pm_473_included")) == 68
    for case in ["ht_table2", "ht_text", "mgd_body"]:
        result = fit(case, "pm_473_included")
        assert result["converged"]
        assert not result["author_workflow_reproduced"]
        for comparison in result["comparison"].values():
            if comparison["source_error_width"] is not None:
                assert comparison["within_source_error"]
        # V0 without a source error is not rescued by inventing an uncertainty.
        if case != "ht_text":
            assert result["comparison"]["V0"]["source_error_width"] is None
            assert abs(result["comparison"]["V0"]["delta"]) > 0.005
        covariance = np.array(result["covariance"])
        assert np.min(np.linalg.eigvalsh(covariance)) > -1e-9
        assert np.allclose(covariance, covariance.T)
    abstract = fit("mgd_abstract", "pm_473_included")
    assert not abstract["comparison"]["V0"]["within_source_error"]
    assert not abstract["comparison"]["K0_prime"]["within_source_error"]
    assert not abstract["comparison"]["theta0"]["within_source_error"]


def test_fm_fit_does_not_meet_printed_k0_width():
    result = fit("fm_text", "fm")
    assert result["n"] == 7
    assert result["parameters"]["K0"] == pytest.approx(178.66, abs=0.02)
    assert not result["comparison"]["K0"]["within_source_error"]
    assert result["parameters"]["K0_prime"] == 5


def test_comparison_rows_do_not_invent_bm_order_or_replace_originals():
    inventory = comparison_inventory()
    assert len(inventory) == 12
    assert sum(item["family"] == "Vinet" for item in inventory) == 5
    assert sum(item["family"] == "not_printed" for item in inventory) == 2
    assert all(not item["catalog_integration"] for item in inventory)
    for item in inventory:
        assert ("pressure_gpa" in item) == (item["family"] == "Vinet")
