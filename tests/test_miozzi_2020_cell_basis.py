"""Whole-cell refit evidence, molar normalization and independent cell physics."""

import hashlib
import json

import numpy as np
import pytest

from scripts.reconstruct_miozzi_2020_iron import NA, ROOT, read_data
from scripts.refit_miozzi_2020_cell_basis import OUTPUT, direct_cell_pressure
from scripts.refit_miozzi_2020_tange import target_pressures


def test_cell_basis_evidence_preserves_input_and_code_fingerprints():
    report = json.loads((OUTPUT / "report.json").read_text(encoding="utf-8"))
    _, data, hashes = read_data()
    assert report["source_csv_sha256"] == hashes
    path = ROOT / "scripts/refit_miozzi_2020_cell_basis.py"
    assert report["script_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    for path, fingerprint in report["dependency_script_sha256"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == fingerprint
    np.testing.assert_array_equal(
        report["observations"]["volume_cell_a3"], data["volume_a3"]
    )
    np.testing.assert_array_equal(
        report["observations"]["temperature_k"], data["temperature_k"]
    )
    np.testing.assert_array_equal(
        report["cases"]["printed"]["observed_pressure_gpa"], data["pressure_gpa"]
    )
    np.testing.assert_allclose(
        report["cases"]["tange_vinet"]["observed_pressure_gpa"],
        target_pressures(data, "vinet"),
        rtol=1e-14,
        atol=0,
    )


@pytest.mark.parametrize(
    "case", ["printed", "tange_vinet", "speziale_variable_q_debye"]
)
def test_cell_basis_fit_matches_per_fe_and_direct_cell_physics(case):
    report = json.loads((OUTPUT / "report.json").read_text(encoding="utf-8"))
    result = report["cases"][case]
    n1, n2 = result["fits"]["n1"], result["fits"]["n2"]
    v = np.array(report["observations"]["volume_cell_a3"])
    t = np.array(report["observations"]["temperature_k"])
    observed = np.array(result["observed_pressure_gpa"])
    for n, fit in ((1, n1), (2, n2)):
        normalization = fit["normalization"]
        assert normalization["n"] == n
        scale = NA * 1e-25 / (2 if n == 1 else 1)
        assert normalization["cell_a3_to_model_volume_scale"] == scale
        rt, thermal, final = fit["stages"]
        assert [len(stage["residuals"]) for stage in fit["stages"]] == [36, 131, 131]
        assert all(stage["solver"]["success"] for stage in fit["stages"])
        assert final["diagnostics"]["degrees_of_freedom"] == 126
        assert thermal["free_parameters"] == ["gamma0", "q"]
        for stage in (thermal, final):
            assert stage["parameters"]["n"] == n
            assert stage["parameters"]["theta0"] == 420
            assert stage["parameters"]["Tr"] == 300
        for key in ("V0", "K0", "K0_prime"):
            assert (
                thermal["model"]["parameters"]["rt_eos." + key] == rt["parameters"][key]
            )
        np.testing.assert_allclose(
            final["adjusted_volume"], v * scale, rtol=1e-14, atol=0
        )
        independent = direct_cell_pressure(v, t, fit["parameters"])
        np.testing.assert_allclose(
            independent, fit["predicted_pressure_gpa"], rtol=0, atol=1e-8
        )
        assert fit["rmse_gpa"] == pytest.approx(
            np.sqrt(np.mean((independent - observed) ** 2)), abs=1e-9
        )
    # Optimizer stopping differs slightly between the two internal volume scales.
    np.testing.assert_allclose(
        n1["predicted_pressure_gpa"], n2["predicted_pressure_gpa"], rtol=0, atol=1e-4
    )
    np.testing.assert_allclose(
        n1["covariance_in_cell_volume_basis"],
        n2["covariance_in_cell_volume_basis"],
        rtol=1e-4,
        atol=1e-7,
    )
    assert n1["rmse_gpa"] == pytest.approx(n2["rmse_gpa"], abs=1e-9)
    assert (
        n2["normalization"]["initial_V0_model"]
        == 2 * n1["normalization"]["initial_V0_model"]
    )
