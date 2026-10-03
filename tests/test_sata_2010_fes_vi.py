"""Check the recovered cold observations, reference equation and EosFit evidence."""

import hashlib
import json

import numpy as np
import pytest

from peritheos.eos.rt import BM3
from scripts.audit_sata_2010_fes_vi import (
    DATA,
    OUTPUT,
    PARAMETERS,
    ROOT,
    VR,
    elastic_properties,
    linear_design,
    observations,
    parse_token,
    pressure,
    reproduce,
)
from scripts.run_fes_eosfit_console import data_lines, parse_console, saved_covariance


def test_source_precision_phase_selection_and_unknown_errors():
    rows = observations()
    assert len(rows) == 19
    included = [r for r in rows if r["included"]]
    assert len(included) == 13
    assert included[-3]["run"] == "4-100"
    assert included[-1]["pressure_gpa"] == 214.5
    assert {r["z"] for r in included} == {4}
    assert rows[13]["error_volume_angstrom3"] is None
    assert rows[-1]["error_volume_angstrom3"] == pytest.approx(0.0067)
    assert parse_token("71.30(3)") == (71.3, 0.03)
    assert parse_token("214.5(6)") == pytest.approx((214.5, 0.6))
    assert rows[0]["volume_raw"] == "76.71(3)"
    manifest = json.loads(
        (DATA.parent / "fes-sata-2010-source-manifest.json").read_text()
    )
    assert manifest["table_sha256"] == hashlib.sha256(DATA.read_bytes()).hexdigest()


def test_finite_reference_has_declared_pressure_modulus_and_derivative():
    vr, values = 83.79080132959231, [36.0, 306.0, 3.81]
    assert pressure(vr, values, vr) == pytest.approx(36)
    assert elastic_properties(vr, values, vr) == pytest.approx((306, 3.81))
    step = vr * 1e-5
    for v in [vr, 70, 61.31]:
        numerical_k = (
            -v
            * (pressure(v + step, values, vr) - pressure(v - step, values, vr))
            / (2 * step)
        )
        k, kp = elastic_properties(v, values, vr)
        assert k == pytest.approx(numerical_k, rel=1e-8)
        k1, _ = elastic_properties(v - step, values, vr)
        k2, _ = elastic_properties(v + step, values, vr)
        numerical_kp = (k1 - k2) / (
            pressure(v - step, values, vr) - pressure(v + step, values, vr)
        )
        assert kp == pytest.approx(numerical_kp, rel=1e-8)


def test_zero_reference_reduces_to_native_bm3_and_linear_benchmark():
    v = np.array([VR, 82.92, 70, 61.31])
    native = BM3(VR, 148, 4.53)
    np.testing.assert_allclose(pressure(v), [native.pressure(x) for x in v], atol=1e-11)
    beta = [PARAMETERS[0], PARAMETERS[1], PARAMETERS[1] * PARAMETERS[2]]
    np.testing.assert_allclose(linear_design(v) @ beta, pressure(v), atol=1e-11)


def test_fit_errors_covariance_and_source_replay_remain_qualified():
    report = reproduce()
    saved = json.loads(OUTPUT.read_text())
    assert report["source_rows"] == saved["source_rows"]
    assert report["published_replay"]["rmse_gpa"] == pytest.approx(1.9860371659)
    fit = report["conditional_refits"]["unweighted_pressure"]
    assert fit["parameters"] == pytest.approx([-0.04423032, 148.3688273, 4.52511843])
    assert fit["standard_errors_residual_scaled"] == pytest.approx(
        [4.213544, 15.694117, 0.3340795]
    )
    for fit in report["conditional_refits"].values():
        assert fit["rank"] == 3 and fit["degrees_of_freedom"] == 10
        assert fit["converged"] and not fit["at_bound"]
        covariance = np.asarray(fit["covariance_residual_scaled"])
        assert np.min(np.linalg.eigvalsh(covariance)) > 0
        np.testing.assert_allclose(covariance, covariance.T, atol=1e-12)
    check = report["finite_reference_check"]
    assert check["max_reference_transform_pressure_difference_gpa"] < 1e-10
    assert report["published_reference"]["parameter_covariance"] is None


def test_direct_eosfit_inputs_preserve_all_rows_units_and_derived_conversions():
    for dataset, n in [("cold", 13), ("morard", 146), ("combined", 159)]:
        lines = data_lines(dataset)
        assert len(lines) - 5 == n
        if dataset != "cold":
            assert lines[2] == "Vscale cm^3/mol"
        if dataset == "combined":
            cold = [r for r in observations() if r["included"]]
            for line, row in zip(lines[-13:], cold):
                p, dp, v, dv, t, dt = map(float, line.split())
                assert (p, dp, t, dt) == (
                    row["pressure_gpa"],
                    row["error_pressure_gpa"],
                    300,
                    0,
                )
                assert v == pytest.approx(
                    row["molar_volume_cm3_per_mol_fes"], rel=1e-14
                )


def test_direct_console_outputs_are_complete_hashed_and_actually_converged():
    folder = ROOT / "docs/data/fes-eosfit7c"
    manifest = json.loads((folder / "manifest.json").read_text())
    assert manifest["version"] == "7.60"
    assert len(manifest["cases"]) == 6
    for name, case in manifest["cases"].items():
        for file, digest in case["files_sha256"].items():
            assert (
                hashlib.sha256((folder / name / file).read_bytes()).hexdigest()
                == digest
            )
        parsed = parse_console(
            (folder / name / "stdout.txt").read_text(), name.split("-")[0]
        )
        assert parsed["converged"]
        assert parsed["fitted_replay"] == case["fitted_replay"]
        assert len(parsed["fitted_replay"]["rows"]) == case["dataset_rows"]
        covariance = saved_covariance(folder / name / "fitted.eos", name.split("-")[0])
        np.testing.assert_array_equal(
            covariance, case["saved_refined_parameter_covariance"]
        )
        assert np.min(np.linalg.eigvalsh(covariance)) > 0
        model = (folder / name / "fitted.eos").read_text()
        assert "Model =  2,  (Birch-Murnaghan)" in model
        if name.startswith("cold"):
            assert "Thermal =  0" in model
        else:
            assert "Thermal =  7,  (Mie-Gruneisen-Debye)" in model
            assert "Param =11 417.000" in model
            assert "Param =13 2.00000" in model
            assert "Param =19 1.00000" in model
    gamma = manifest["cases"]["morard-unit"]["refined_parameters_final_cycle"]["Gamm0"]
    assert gamma == {"value": 2.28211, "esd": 0.03231}
    source_frame = manifest["cases"]["cold-unit"]["sata_fixed_reference_representation"]
    independent = reproduce()["conditional_refits"]["unweighted_pressure"]
    np.testing.assert_allclose(
        source_frame["parameters"], independent["parameters"], atol=1e-3, rtol=0
    )
    np.testing.assert_allclose(
        source_frame["esds"], independent["standard_errors_residual_scaled"], rtol=1e-3
    )
    assert (
        len(
            manifest["cases"]["morard-unit"]["initial_source_replay"][
                "rows_exceeding_3_gpa"
            ]
        )
        == 15
    )
    assert any("esds=0" in w for w in manifest["cases"]["combined-errors"]["warnings"])


def test_manual_defaults_are_not_mistaken_for_author_model_files():
    report = json.loads((ROOT / "docs/data/eosfit-mgd-manual-audit.json").read_text())
    assert report["full_mgd"]["q_compromise_prompt_default"] == "N"
    assert "authors exact" in report["qualification"]


def test_same_process_staged_refits_keep_the_cold_curve_fixed():
    folder = ROOT / "docs/data/fes-eosfit7c-staged"
    manifest = json.loads((folder / "manifest.json").read_text())
    assert set(manifest["cases"]) == {"unit", "errors"}
    for name, case in manifest["cases"].items():
        for filename, digest in case["files_sha256"].items():
            assert (
                hashlib.sha256((folder / name / filename).read_bytes()).hexdigest()
                == digest
            )
        cold, thermal = case["cold_stage"], case["thermal_stage"]
        assert cold["converged"] and thermal["converged"]
        assert (
            thermal["fixed_cold_parameters"] == cold["refined_parameters_final_cycle"]
        )
        assert set(thermal["refined_parameters_final_cycle"]) == {"Gamm0"}
        assert len(thermal["fitted_replay"]["rows"]) == 159
        np.testing.assert_array_equal(
            [r["residual_gpa"] for r in cold["fitted_replay"]["rows"]],
            [r["residual_gpa"] for r in thermal["fitted_replay"]["rows"][146:]],
        )
        final = (
            (folder / name / "stdout.txt")
            .read_text()
            .split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1]
        )
        for parameter, values in cold["refined_parameters_final_cycle"].items():
            assert f"{values['value']:.5f}  [NOT REFINED]" in final
        covariance = np.array(cold["saved_refined_parameter_covariance"])
        assert np.min(np.linalg.eigvalsh(covariance)) > 0
    gamma = manifest["cases"]["unit"]["thermal_stage"][
        "refined_parameters_final_cycle"
    ]["Gamm0"]
    assert gamma == {"value": 2.44764, "esd": 0.04153}
