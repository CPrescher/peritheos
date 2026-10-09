"""Verify q-compromise selection, unchanged observations and independent replay."""

import hashlib
import json
from pathlib import PurePosixPath

import numpy as np
import pytest

from scripts.compare_fes_q_compromise import (
    DIRECT,
    OUTPUT,
    PUBLISHED,
    ROOT,
    STAGED,
    model,
    reproduce,
)
from scripts.run_fes_eosfit_console import macro, parse_console
from scripts.run_fes_eosfit_staged import staged_macro


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_q_compromise_has_volume_independent_thermal_pressure():
    eos = model(PUBLISHED)
    volumes = np.array([0.9, 1.2, 1.54])
    assert not hasattr(eos, "q")
    assert np.all(eos.thermal_pressure(volumes, 300) == 0)
    np.testing.assert_array_equal(
        eos.thermal_pressure(volumes, 2000),
        np.full(3, eos.thermal_pressure(1.54, 2000)),
    )
    np.testing.assert_allclose(
        eos.pressure(volumes, 2000) - eos.rt_eos.pressure(volumes),
        eos.thermal_pressure(volumes, 2000),
        rtol=0,
        atol=1e-13,
    )


def test_original_full_model_macros_are_preserved():
    def archived_macro(case):
        text = (case / "run.mcr").read_text(encoding="utf-8")
        # The immutable evidence retains paths from the original Mac checkout.
        # Compare commands on that recorded basis, including on Windows or CI.
        log = text.splitlines()[0]
        assert log.startswith("log ")
        folder = PurePosixPath(log[4:]).parent
        assert folder.name == case.name
        return folder, text

    folder = ROOT / "docs/data/fes-eosfit7c"
    for dataset in ["cold", "morard", "combined"]:
        for weighted in [False, True]:
            case = folder / (dataset + ("-errors" if weighted else "-unit"))
            original_folder, text = archived_macro(case)
            assert macro(dataset, weighted, original_folder) == text
    folder = ROOT / "docs/data/fes-eosfit7c-staged"
    for weighted in [False, True]:
        case = folder / ("errors" if weighted else "unit")
        original_folder, text = archived_macro(case)
        assert staged_macro(original_folder, weighted) == text


def test_actual_console_selected_q_compromise_and_kept_inputs():
    direct = json.loads((DIRECT / "manifest.json").read_text(encoding="utf-8"))
    full = json.loads(
        (ROOT / "docs/data/fes-eosfit7c/manifest.json").read_text(encoding="utf-8")
    )
    assert direct["executable_sha256"] == full["executable_sha256"]
    assert direct["thermal_source_sha256"] == full["thermal_source_sha256"]
    assert direct["cold_source_sha256"] == full["cold_source_sha256"]
    assert direct["version"] == "7.60"
    for name, case in direct["cases"].items():
        for filename, digest in case["files_sha256"].items():
            assert (
                hashlib.sha256((DIRECT / name / filename).read_bytes()).hexdigest()
                == digest
            )
        assert (DIRECT / name / "input.dat").read_bytes() == (
            ROOT / "docs/data/fes-eosfit7c" / name / "input.dat"
        ).read_bytes()
        text = (DIRECT / name / "stdout.txt").read_text(encoding="utf-8")
        parsed = parse_console(text, name.split("-")[0])
        assert parsed["converged"]
        assert parsed["fitted_replay"] == case["fitted_replay"]
        if name.startswith("cold"):
            assert (
                case["refined_parameters_final_cycle"]
                == full["cases"][name]["refined_parameters_final_cycle"]
            )
            continue
        assert "with q-compromise" in text
        assert "q is not defined for q-compromise" in text
        saved = (DIRECT / name / "fitted.eos").read_text(encoding="utf-8")
        assert "Param =14 1.00000" in saved  # Explicit q-compromise flag.
        assert "Param =19 0.00000" in saved  # Not a free or fixed physical q.


def test_staged_q_compromise_preserves_cold_solution():
    manifest = json.loads((STAGED / "manifest.json").read_text(encoding="utf-8"))
    full = json.loads(
        (ROOT / "docs/data/fes-eosfit7c-staged/manifest.json").read_text(
            encoding="utf-8"
        )
    )
    for name, case in manifest["cases"].items():
        for filename, digest in case["files_sha256"].items():
            assert (
                hashlib.sha256((STAGED / name / filename).read_bytes()).hexdigest()
                == digest
            )
        cold, hot = case["cold_stage"], case["thermal_stage"]
        assert cold["converged"] and hot["converged"]
        assert (
            cold["refined_parameters_final_cycle"]
            == full["cases"][name]["cold_stage"]["refined_parameters_final_cycle"]
        )
        assert hot["fixed_cold_parameters"] == cold["refined_parameters_final_cycle"]
        assert set(hot["refined_parameters_final_cycle"]) == {"Gamm0"}
        np.testing.assert_array_equal(
            [r["residual_gpa"] for r in hot["fitted_replay"]["rows"][146:]],
            [r["residual_gpa"] for r in cold["fitted_replay"]["rows"]],
        )


def test_peritheos_component_evaluation_matches_actual_console(report):
    assert len(report["cases"]) == 6
    retained = json.loads(OUTPUT.read_text(encoding="utf-8"))
    for name, case in report["cases"].items():
        for stage, result in case.items():
            assert result["max_pressure_difference_gpa"] < 0.0012
            assert result["max_native_python_energy_difference_in_pressure_gpa"] < 1e-10
            assert (
                result["eosfit_all"]["rows_exceeding_3_gpa"]
                == result["peritheos_all"]["rows_exceeding_3_gpa"]
            )
            np.testing.assert_allclose(
                [r["peritheos_pressure_gpa"] for r in result["rows"]],
                [
                    r["peritheos_pressure_gpa"]
                    for r in retained["cases"][name][stage]["rows"]
                ],
                rtol=0,
                atol=1e-10,
            )


def test_independent_joint_refit_and_failed_published_bound(report):
    fit = report["independent_joint_equal_weight_refit"]
    console = json.loads((DIRECT / "manifest.json").read_text(encoding="utf-8"))[
        "cases"
    ]["combined-unit"]
    assert fit["solver"]["success"]
    assert fit["degrees_of_freedom"] == 155
    cov = np.array(fit["covariance"])
    assert np.linalg.matrix_rank(cov) == 4
    assert np.all(np.linalg.eigvalsh(cov) > 0)
    for i, key in enumerate(fit["covariance_parameter_order"]):
        actual = fit["parameters_cm3_mol_gpa"][key]
        expected = console["refined_parameters_final_cycle"][key]
        assert actual["value"] == pytest.approx(expected["value"], rel=2e-4)
        assert actual["esd"] == pytest.approx(expected["esd"], rel=2e-4)
        assert actual["esd"] == pytest.approx(np.sqrt(cov[i, i]), rel=1e-12)
    assert fit["all_rows"]["rmse_gpa"] == pytest.approx(1.3373041, abs=1e-5)
    published = report["cases"]["morard-unit"]["initial"]["eosfit_thermal"]
    assert published["rmse_gpa"] == pytest.approx(2.17538184)
    assert len(published["rows_exceeding_3_gpa"]) == 21
    for case in report["cases"].values():
        assert case["fitted"]["eosfit_thermal"]["max_abs_gpa"] > 3
