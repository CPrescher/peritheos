"""Verify the author-source evidence and external-program audit conclusions."""

import json

import numpy as np
import pytest

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from scripts.audit_fes_author_inputs import (
    OUTPUT,
    compare_cold_sources,
    compare_table,
    digest,
    initial_rows,
    read_data,
)
from scripts.run_fes_eosfit_console import parse_console


@pytest.fixture(scope="module")
def report():
    return json.loads((OUTPUT / "manifest.json").read_text())


def test_source_bytes_and_actual_input_difference(report):
    for filename, source in report["source_files"].items():
        assert digest(OUTPUT / "sources" / filename) == source["sha256"]
    for filename, expected_hash in report["artifact_sha256"].items():
        assert digest(OUTPUT / filename) == expected_hash
    full = read_data(OUTPUT / "sources/FeS6EOSfittot-SataOhfuji.dat")
    short = read_data(OUTPUT / "sources/FeS6EOSfittot-SataOhfuji-300K.dat")
    assert len(full) == 178 and len(short) == 167
    keep = ~((full[:, 2] == 300) & (full[:, 1] == 0.0056))
    np.testing.assert_array_equal(full[keep], short)
    assert len(full[~keep]) == 11
    assert np.all(full[~keep, 3] == 5)
    assert np.sum(short[:, 2] == 300) == 21


def test_thermal_data_preserved_with_author_rounded_conversion(report):
    for filename in report["source_files"]:
        if not filename.endswith(".dat"):
            continue
        data = read_data(OUTPUT / "sources" / filename)
        comparison = compare_table(data)
        assert comparison["same_thermal_pressure_temperature_order"]
        assert (
            comparison["max_volume_difference_from_author_rounded_conversion_cm3_mol"]
            < 5.1e-9
        )
        hot = data[data[:, 2] > 300]
        np.testing.assert_array_equal(hot[6], hot[7])
        cold = compare_cold_sources(data)
        assert len(cold["sata_vi_volume_matches"]) == 13
        assert len(cold["additional_cold_rows"]) == 8
        literature = data[(data[:, 2] == 300) & (data[:, 1] != 0.0056)]
        np.testing.assert_allclose(
            literature[:, 1], 0.01 * literature[:, 0], atol=1e-12
        )
        np.testing.assert_allclose(
            literature[:, 5], 0.005 * literature[:, 4], atol=6e-10
        )


def test_native_replay_agrees_with_actual_external_program(report):
    eos = MieGruneisenDebye(BM3(1.537, 115.49, 4.99), 298, 417, 2.41671, 1, 2)
    assert hasattr(eos, "_native")
    for name, filename in [
        ("all-author-rows", "FeS6EOSfittot-SataOhfuji.dat"),
        ("without-own-cold", "FeS6EOSfittot-SataOhfuji-300K.dat"),
    ]:
        data = read_data(OUTPUT / "sources" / filename)
        case = report["cases"][name + "-replay"]
        external_rows = initial_rows(
            (OUTPUT / (name + "-replay") / "stdout.txt").read_text(), len(data)
        )
        residual = np.array([x["residual_gpa"] for x in external_rows])
        native = eos.pressure(data[:, 4] / 10, data[:, 2]) - data[:, 0]
        assert np.max(np.abs(native - residual)) < 0.0012
        assert case["eosfit"]["thermal"]["rmse_gpa"] == pytest.approx(
            1.6317545, abs=1e-7
        )
        assert len(case["eosfit"]["thermal"]["rows_exceeding_3_gpa"]) == 14
        assert np.max(np.abs(residual[data[:, 2] > 300])) == pytest.approx(5.0579)


def test_saved_model_is_full_mgd_with_conditional_gamma_variance(report):
    model = report["saved_author_model"]
    assert model["q_compromise_flag"] == 0
    assert model["q"] == 1 and model["Tref_k"] == 298
    assert model["covariance_nonzero_entries"] == [
        {"parameter_i": 18, "parameter_j": 18, "value": 0.0016139}
    ]
    assert model["gamma0_saved_esd"] == pytest.approx(np.sqrt(0.0016139))


def test_actual_refinements_converge_and_cold_fit_approaches_saved_model(report):
    assert len(report["cases"]) == 16
    for name, case in report["cases"].items():
        if name.endswith("replay"):
            continue
        assert case["converged"]
        # The input label itself contains 'cold'; the refined keys identify
        # whether a thermal parameter is present without relying on that label.
        keys = set(case["refined_parameters_final_cycle"])
        cold = keys == {"V0", "K0", "Kp"}
        dataset = "combined" if cold or "joint" in name else "morard"
        rows = len(case["fitted_replay"]["rows"])
        parsed = parse_console(
            (OUTPUT / name / "stdout.txt").read_text(),
            dataset,
            n_observations=rows,
            parameter_names=keys,
        )
        assert (
            parsed["refined_parameters_final_cycle"]
            == case["refined_parameters_final_cycle"]
        )
        cov = np.asarray(case["saved_refined_parameter_covariance"])
        np.testing.assert_array_equal(cov, cov.T)
        esds = [
            case["refined_parameters_final_cycle"][key]["esd"]
            for key in case["covariance_parameter_order"]
        ]
        np.testing.assert_allclose(np.sqrt(np.diag(cov)), esds, rtol=2e-4)
        # Fixed-width saved covariances have five significant digits. Highly
        # correlated cold fits can acquire a small negative eigenvalue on
        # output rounding; retain the source matrix without repairing it.
        assert np.min(np.linalg.eigvalsh(cov)) > -5e-5 * np.max(np.abs(cov))
    fitted = report["cases"]["without-own-cold-cold-errors"][
        "refined_parameters_final_cycle"
    ]
    assert fitted["V0"]["value"] == pytest.approx(15.37, abs=0.001)
    assert fitted["K0"]["value"] == pytest.approx(115.49, abs=0.051)
    assert fitted["Kp"]["value"] == pytest.approx(4.99, abs=0.002)
