"""Evidence from the real EosFit replay of the user-relayed author protocol."""

import hashlib
import json
import re

import numpy as np
import pytest
from scipy.optimize import least_squares

from scripts.reconstruct_miozzi_2020_iron import ROOT, read_data
from scripts.refit_miozzi_2020_tange import reproduce
from scripts.replay_miozzi_2020_staged import FACTOR, model_pressure

REPORT = ROOT / "docs/data/miozzi-2020-eosfit-staged/report.json"


def test_console_evidence_preserves_inputs_and_stage_constraints():
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    _, data, hashes = read_data()
    assert report["source_csv_sha256"] == hashes
    assert "User explicitly clarified" in report["protocol_provenance"]
    for name, case in report["cases"].items():
        assert "failure" not in case
        assert case["returncode"] == 0
        assert case["completed_stages"] == case["converged_stages"] == 3
        assert case["stages_saved"]
        folder = REPORT.parent / name
        for filename, sha in case["files_sha256"].items():
            assert hashlib.sha256((folder / filename).read_bytes()).hexdigest() == sha
        cold, thermal, final = case["stages"]
        assert [s["rows"] for s in case["stages"]] == [36, 131, 131]
        assert thermal["free_parameters"] == ["Gamm0", "q"]
        assert final["free_parameters"] == ["V0", "K0", "Kp", "Gamm0", "q"]
        for key in ["V0", "K0", "Kp"]:
            assert thermal["parameters"][key] == cold["parameters"][key]
        assert thermal["parameters"]["ThMGD"] == final["parameters"]["ThMGD"] == 420
        for stage in case["stages"]:
            cov = np.asarray(stage["covariance"])
            errors = np.array(
                list(stage["conditional_eosfit_standard_errors"].values())
            )
            assert np.all(np.linalg.eigvalsh(cov) > 0)
            np.testing.assert_allclose(np.diag(cov), errors**2, rtol=3e-4)
        inputs = np.loadtxt(folder / "all.dat", skiprows=5)
        np.testing.assert_array_equal(inputs[:, 4], data["temperature_k"])
        np.testing.assert_array_equal(inputs[:15, 0], data["pressure_gpa"][:15])
        assert len(inputs) == 131


def test_tange_staged_console_matches_registered_joint_optimum():
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    final = report["cases"]["tange_vinet-unit-theta_fixed"]["stages"][-1]
    stored = reproduce()["primary"]
    mapping = {"V0": "V0", "K0": "K0", "Kp": "K0_prime", "Gamm0": "gamma0", "q": "q"}
    for console, name in mapping.items():
        # R rounding and tabulated coefficients limit the console comparison.
        assert final["parameters"][console] == pytest.approx(
            stored["parameters"][name], rel=6e-5
        )
        assert final["conditional_eosfit_standard_errors"][console] == pytest.approx(
            stored["standard_errors"][name], rel=3e-4
        )
    assert final["rmse_gpa_from_printed_coefficients"] == pytest.approx(
        stored["fit_row_metrics"]["rmse_gpa"], abs=1e-5
    )
    np.testing.assert_allclose(final["covariance"], stored["covariance"], rtol=2e-3)


def test_speziale_cold_recovery_does_not_establish_source_thermal_parity():
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    stages = report["cases"]["speziale_variable_q_debye-unit-theta_fixed"]["stages"]
    assert stages[0]["parameters"]["K0"] == pytest.approx(129, abs=1)
    assert stages[-1]["parameters"]["Gamm0"] > 1.8
    assert stages[-1]["parameters"]["Gamm0"] != pytest.approx(1.11, abs=0.01)


@pytest.mark.parametrize("switch_to_full", [False, True])
def test_q_compromise_console_evidence_and_independent_pressure(switch_to_full):
    suffix = "-start" if switch_to_full else ""
    path = ROOT / f"docs/data/miozzi-2020-eosfit-q-compromise{suffix}/report.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    original = json.loads(REPORT.read_text(encoding="utf-8"))
    _, data, hashes = read_data()
    assert report["source_csv_sha256"] == hashes
    for name, case in report["cases"].items():
        assert "failure" not in case
        assert case["returncode"] == 0
        assert case["completed_stages"] == case["converged_stages"] == 3
        folder = path.parent / name
        for filename, sha in case["files_sha256"].items():
            assert hashlib.sha256((folder / filename).read_bytes()).hexdigest() == sha
        assert case["stages"][0] == original["cases"][name]["stages"][0]
        assert case["stages"][1]["free_parameters"] == ["Gamm0"]
        assert "q" not in case["stages"][1]["parameters"]
        inputs = np.loadtxt(folder / "all.dat", skiprows=5)
        np.testing.assert_array_equal(inputs[:, 4], data["temperature_k"])
        np.testing.assert_array_equal(inputs[:15, 0], data["pressure_gpa"][:15])
        blocks = (
            (folder / "stdout.txt")
            .read_text(encoding="utf-8")
            .split("RESULTS AFTER FINAL REFINEMENT CYCLE")[1:]
        )
        for i, (stage, block) in enumerate(zip(case["stages"], blocks)):
            q_compromise = i > 0 and not (switch_to_full and i == 2)
            saved = (folder / f"stage{i + 1}.eos").read_text(encoding="utf-8")
            flag = re.search(r"Param =14\s+([\d.]+)", saved)
            assert float(flag[1]) == int(q_compromise)
            if i:
                assert stage["parameters"]["ThMGD"] == 420
                assert stage["parameters"]["Gamm0"] > 1.8
            cov = np.asarray(stage["covariance"])
            errors = np.array(
                list(stage["conditional_eosfit_standard_errors"].values())
            )
            assert np.all(np.linalg.eigvalsh(cov) > 0)
            # Console ESDs have five decimals; saved covariances have five
            # significant figures. Compare their rounding intervals.
            rounding = np.array(
                [
                    1e-5 / FACTOR if key == "V0" else 1e-5
                    for key in stage["conditional_eosfit_standard_errors"]
                ]
            )
            assert np.all(
                np.abs(np.diag(cov) - errors**2)
                <= errors * rounding + np.abs(np.diag(cov)) * 5e-5
            )
            # The console's actual calculated pressures independently confirm
            # the different approximation. R and printed precision differ.
            rows = [
                line.split()
                for line in block.splitlines()
                if re.match(r"^\s*\d+\s+1\s+[-\d.]", line)
            ][: stage["rows"]]
            select = data["temperature_k"] <= 300 if i == 0 else np.ones(131, bool)
            p = model_pressure(
                data["volume_a3"][select],
                data["temperature_k"][select],
                stage["parameters"],
                q_compromise,
            )
            np.testing.assert_allclose(
                p,
                [float(row[8 if i else 7]) for row in rows],
                rtol=0,
                atol=0.0025,
            )
            residual = p - inputs[select, 0]
            assert stage["rmse_gpa_from_printed_coefficients"] == pytest.approx(
                np.sqrt(np.mean(residual**2)),
                abs=1e-12,
            )
        final = case["stages"][-1]
        if switch_to_full:
            assert final["free_parameters"] == ["V0", "K0", "Kp", "Gamm0", "q"]
            full = original["cases"][name]["stages"][-1]
            for key, value in full["parameters"].items():
                assert final["parameters"][key] == pytest.approx(value, rel=3e-4)
            assert final["rmse_gpa_from_printed_coefficients"] == pytest.approx(
                full["rmse_gpa_from_printed_coefficients"],
                abs=1e-6,
            )
        else:
            assert final["free_parameters"] == ["V0", "K0", "Kp", "Gamm0"]
            assert "q" not in final["parameters"]


@pytest.mark.parametrize(
    "calibration", ["printed", "tange_vinet", "speziale_variable_q_debye"]
)
def test_q_compromise_matches_independent_least_squares(calibration):
    folder = ROOT / "docs/data/miozzi-2020-eosfit-q-compromise"
    name = f"{calibration}-unit-theta_fixed"
    report = json.loads((folder / "report.json").read_text(encoding="utf-8"))
    final = report["cases"][name]["stages"][-1]
    _, data, _ = read_data()
    pressures = np.loadtxt(folder / name / "all.dat", skiprows=5)[:, 0]

    def residual(x):
        values = dict(zip(["V0", "K0", "Kp", "Gamm0"], x))
        values["ThMGD"] = 420
        return (
            model_pressure(data["volume_a3"], data["temperature_k"], values, True)
            - pressures
        )

    fit = least_squares(
        residual,
        [22.81, 129, 6.24, 1.11],
        bounds=([20, 50, 2, 0.1], [25, 300, 10, 4]),
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    assert fit.success
    for key, value in zip(["V0", "K0", "Kp", "Gamm0"], fit.x):
        assert final["parameters"][key] == pytest.approx(value, rel=1e-4)
    assert np.sqrt(np.mean(fit.fun**2)) == pytest.approx(
        final["rmse_gpa_from_printed_coefficients"],
        abs=1e-6,
    )
    covariance = np.linalg.inv(fit.jac.T @ fit.jac) * np.sum(fit.fun**2) / (131 - 4)
    np.testing.assert_allclose(final["covariance"], covariance, rtol=2e-3)


@pytest.mark.parametrize("multiplier", [1, 2])
def test_full_mgd_n2_author_replay_and_volume_control(multiplier):
    suffix = "n2-author-volume" if multiplier == 1 else "n2-double-volume-control"
    folder = ROOT / f"docs/data/miozzi-2020-eosfit-{suffix}"
    report = json.loads((folder / "report.json").read_text(encoding="utf-8"))
    old = json.loads(REPORT.read_text(encoding="utf-8"))
    _, data, hashes = read_data()
    assert report["source_csv_sha256"] == hashes
    assert report["thermal_model"] == "full_mgd"
    assert report["normalization"]["eosfit_atoms_per_formula_unit"] == 2
    assert (
        report["normalization"]["molar_volume_multiplier_relative_to_moles_of_Fe"]
        == multiplier
    )
    assert "n=2 and V0 approximately 6.87" in report["protocol_provenance"]
    for name, case in report["cases"].items():
        assert "failure" not in case
        assert case["returncode"] == 0
        assert case["completed_stages"] == case["converged_stages"] == 3
        for filename, sha in case["files_sha256"].items():
            assert (
                hashlib.sha256((folder / name / filename).read_bytes()).hexdigest()
                == sha
            )
        inputs = np.loadtxt(folder / name / "all.dat", skiprows=5)
        prior_inputs = np.loadtxt(REPORT.parent / name / "all.dat", skiprows=5)
        np.testing.assert_array_equal(
            inputs[:, [0, 1, 4, 5]], prior_inputs[:, [0, 1, 4, 5]]
        )
        np.testing.assert_allclose(
            inputs[:, [2, 3]], prior_inputs[:, [2, 3]] * multiplier, rtol=1e-15
        )
        blocks = (
            (folder / name / "stdout.txt")
            .read_text(encoding="utf-8")
            .split("RESULTS AFTER FINAL REFINEMENT CYCLE")[1:]
        )
        for i, (stage, block) in enumerate(zip(case["stages"], blocks)):
            saved = (folder / name / f"stage{i + 1}.eos").read_text(encoding="utf-8")
            assert float(re.search(r"Param =14\s+([\d.]+)", saved)[1]) == 0
            if i:
                assert float(re.search(r"Param =13\s+([\d.]+)", saved)[1]) == 2
                assert stage["parameters"]["ThMGD"] == 420
                assert "q" in stage["parameters"]
            assert stage["rows"] == (36 if i == 0 else 131)
            rows = [
                line.split()
                for line in block.splitlines()
                if re.match(r"^\s*\d+\s+1\s+[-\d.]", line)
            ][: stage["rows"]]
            select = data["temperature_k"] <= 300 if i == 0 else np.ones(131, bool)
            p = model_pressure(
                data["volume_a3"][select],
                data["temperature_k"][select],
                stage["parameters"],
                False,
                2,
                multiplier,
            )
            np.testing.assert_allclose(
                p, [float(row[8 if i else 7]) for row in rows], rtol=0, atol=0.0025
            )
            if multiplier == 2:
                prior = old["cases"][name]["stages"][i]
                for key, value in prior["parameters"].items():
                    assert stage["parameters"][key] == pytest.approx(value, rel=3e-4)
                assert stage["rmse_gpa_from_printed_coefficients"] == pytest.approx(
                    prior["rmse_gpa_from_printed_coefficients"], abs=3e-5
                )
        final = case["stages"][-1]
        assert final["free_parameters"] == ["V0", "K0", "Kp", "Gamm0", "q"]
        if multiplier == 1:
            assert 0.9 < final["parameters"]["Gamm0"] < 1.2


def test_full_mgd_normalization_equivalence_including_volume_dependent_theta():
    _, data, _ = read_data()
    values = {"V0": 22.81, "K0": 129, "Kp": 6.24, "ThMGD": 420, "Gamm0": 1.11, "q": 0.3}
    p1 = model_pressure(data["volume_a3"], data["temperature_k"], values, False, 1, 1)
    p2 = model_pressure(data["volume_a3"], data["temperature_k"], values, False, 2, 2)
    np.testing.assert_allclose(p1, p2, rtol=1e-15)


@pytest.mark.parametrize(
    "calibration", ["printed", "tange_vinet", "speziale_variable_q_debye"]
)
def test_n2_full_mgd_console_matches_independent_fit(calibration):
    folder = ROOT / "docs/data/miozzi-2020-eosfit-n2-author-volume"
    name = f"{calibration}-unit-theta_fixed"
    final = json.loads((folder / "report.json").read_text(encoding="utf-8"))["cases"][
        name
    ]["stages"][-1]
    _, data, _ = read_data()
    pressure = np.loadtxt(folder / name / "all.dat", skiprows=5)[:, 0]

    def residual(x):
        values = dict(zip(["V0", "K0", "Kp", "Gamm0", "q"], x))
        values["ThMGD"] = 420
        return (
            model_pressure(
                data["volume_a3"], data["temperature_k"], values, False, 2, 1
            )
            - pressure
        )

    fit = least_squares(
        residual,
        [22.81, 129, 6.24, 1.11, 0.3],
        bounds=([20, 50, 2, 0.1, -3], [25, 300, 10, 4, 5]),
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    assert fit.success
    for key, value in zip(["V0", "K0", "Kp", "Gamm0", "q"], fit.x):
        assert final["parameters"][key] == pytest.approx(value, rel=4e-4)
    assert np.sqrt(np.mean(fit.fun**2)) == pytest.approx(
        final["rmse_gpa_from_printed_coefficients"], abs=2e-6
    )
    covariance = np.linalg.inv(fit.jac.T @ fit.jac) * np.sum(fit.fun**2) / (131 - 5)
    np.testing.assert_allclose(final["covariance"], covariance, rtol=5e-3)
