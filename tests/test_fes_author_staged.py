"""Selected-data staged fits retain the cold curve and reproduce pressure."""

import json
import re

import numpy as np
import pytest

from scripts.audit_fes_author_inputs import digest, read_data
from scripts.run_fes_author_staged import AUTHOR, OUTPUT, PUBLISHED, verify_native
from scripts.run_fes_eosfit_console import parse_console


@pytest.fixture(scope="module")
def report():
    return json.loads((OUTPUT / "manifest.json").read_text())


def test_only_author_confirmed_selected_data_are_used(report):
    source = AUTHOR / report["source_file"]
    assert source.name == "FeS6EOSfittot-SataOhfuji-300K.dat"
    assert digest(source) == report["source_sha256"]
    assert report["selection"] == {"literature_cold": 21, "thermal": 146, "quenched": 0}
    selected = read_data(source)
    assert len(report["cases"]) == 6
    for name, case in report["cases"].items():
        for filename, checksum in case["artifact_sha256"].items():
            assert digest(OUTPUT / name / filename) == checksum
        folder = OUTPUT / (name + "-inputs")
        for filename, checksum in case["input_sha256"].items():
            assert digest(folder / filename) == checksum
        cold, hot = read_data(folder / "cold.dat"), read_data(folder / "thermal.dat")
        np.testing.assert_array_equal(cold, selected[selected[:, 2] == 300])
        np.testing.assert_array_equal(hot, selected[selected[:, 2] > 300])
        assert cold.shape == (21, 6) and hot.shape == (146, 6)
        assert np.all(cold[:, 1] != 0.0056)
        np.testing.assert_array_equal(hot[6], hot[7])


def test_cold_then_fixed_cold_thermal_with_only_gamma_refined(report):
    assert report["fixed_thermal_parameters"] == {
        "theta0_k": 417,
        "n": 2,
        "q": 1,
        "q_compromise": False,
    }
    for name, case in report["cases"].items():
        stdout = (OUTPUT / name / "stdout.txt").read_text()
        cold_text, thermal_text = stdout.split("EOSFIT-7.6>clear", 1)
        thermal = parse_console(thermal_text, "morard", n_observations=146)
        assert thermal["converged"]
        assert set(thermal["refined_parameters_final_cycle"]) == {"Gamm0"}
        assert thermal["fitted_replay"] == case["thermal_stage"]["fitted_replay"]
        final = thermal_text.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1]
        for key, value in case["fixed_cold_parameters"].items():
            found = re.search(rf"^\s*{key}\s+0\s+([-\d.]+)", final, re.M)
            assert found and float(found[1]) == pytest.approx(value, abs=1e-5)
        assert (
            "REFERENCE TEMPERATURE = " + str(case["reference_temperature_k"]) in final
        )
        stage = case["cold_stage"]
        if case["cold_source"] == "published":
            assert stage is None and case["fixed_cold_parameters"] == PUBLISHED
        else:
            assert stage["converged"]
            assert (
                parse_console(
                    cold_text,
                    "combined",
                    n_observations=21,
                    parameter_names=["V0", "K0", "Kp"],
                )["fitted_replay"]
                == stage["fitted_replay"]
            )
            cov = np.array(stage["saved_refined_parameter_covariance"])
            np.testing.assert_allclose(
                np.sqrt(np.diag(cov)),
                [
                    stage["refined_parameters_final_cycle"][k]["esd"]
                    for k in ("V0", "K0", "Kp")
                ],
                atol=6e-5,
                rtol=1e-5,
            )
        gamma = thermal["refined_parameters_final_cycle"]["Gamm0"]
        variance = case["thermal_stage"]["saved_refined_parameter_covariance"][0][0]
        assert np.sqrt(variance) == pytest.approx(gamma["esd"], abs=6e-6)
        assert len(thermal["fitted_replay"]["rows_exceeding_3_gpa"]) > 0


def test_native_model_agrees_with_all_staged_eosfit_results(report):
    verified = verify_native(json.loads(json.dumps(report)))
    for name, case in verified["cases"].items():
        actual = case["native_crosscheck"]
        expected = report["cases"][name]["native_crosscheck"]
        assert actual["implementation"] == expected["implementation"]
        assert actual["qualification"] == expected["qualification"]
        assert actual["max_residual_difference_gpa"] == pytest.approx(
            expected["max_residual_difference_gpa"], abs=1e-8
        )
        for key, value in actual["thermal_residual_summary"].items():
            if isinstance(value, float):
                assert value == pytest.approx(
                    expected["thermal_residual_summary"][key], abs=1e-8
                )
            else:
                assert value == expected["thermal_residual_summary"][key]
        assert case["native_crosscheck"]["max_residual_difference_gpa"] < 0.0012


def test_weighted_staged_comparison_with_published_gamma(report):
    published = report["cases"]["published-errors-tr300"]
    reproduced = report["cases"]["reproduced-errors-tr300"]
    assert published["thermal_stage"]["refined_parameters_final_cycle"]["Gamm0"][
        "value"
    ] == pytest.approx(2.29252)
    assert reproduced["thermal_stage"]["refined_parameters_final_cycle"]["Gamm0"][
        "value"
    ] == pytest.approx(2.46418)
    assert (
        reproduced["thermal_residual_summary"]["rmse_gpa"]
        < published["thermal_residual_summary"]["rmse_gpa"]
    )
    assert report["published_thermal_parameter"] == {"gamma0": 2.42, "esd": 0.03}
