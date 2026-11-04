"""Scientific checks for the isolated Brown NaCl reconstruction."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import interpolate
from scipy.integrate import quad

from scripts.brown_nacl_1999_reconstruct import (
    PREF,
    VREF,
    M,
    Model,
    N,
    R,
    Settings,
    Thermal,
    comparison,
    cv,
    energy,
    expansion_volume,
    length_integral,
    read,
)


def test_recovered_resource_bytes_match_the_source_manifest():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(
        (
            root / "output/material-candidates/brown-nacl-1999/recovery-manifest.json"
        ).read_text(encoding="utf-8")
    )
    for name, resource in manifest["datasets"].items():
        # Source fingerprints include the supplied CSV line endings.
        payload = (root / "peritheos/data/datasets" / name).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == resource["sha256"], name


@pytest.fixture(scope="module")
def model():
    if not hasattr(interpolate, "make_smoothing_spline"):
        pytest.skip(
            "Optional natural-smoothing research diagnostic requires SciPy >=1.10"
        )
    return Model(
        Settings(smoothing=1e-7, spline_kind="natural_smoothing", join_width=0.025)
    )


def test_optional_smoothing_requirement_does_not_change_the_requested_method(
    monkeypatch,
):
    monkeypatch.delattr(interpolate, "make_smoothing_spline", raising=False)
    with pytest.raises(RuntimeError, match="requires SciPy >=1.10"):
        Model(Settings(spline_kind="natural_smoothing"))


@pytest.mark.parametrize("T,theta", [(100, 300), (300, 375), (1200, 375), (5000, 450)])
def test_debye_against_independent_adaptive_quadrature(T, theta):
    y = theta / T
    expected = (
        9
        * N
        * R
        * T
        / M
        / y**3
        * quad(lambda x: x**3 / np.expm1(x), 0, y, epsabs=1e-12)[0]
    )
    assert float(energy(T, theta)) == pytest.approx(expected, rel=1e-12)


def test_heat_capacity_is_energy_derivative():
    assert float(cv(300, 375)) == pytest.approx(
        float((energy(300.001, 375) - energy(299.999, 375)) / 0.002), rel=1e-9
    )
    assert float(cv(1e6, 300)) == pytest.approx(3 * N * R / M, rel=1e-7)


def test_gamma_theta_identity():
    t = Thermal(Settings(join_width=0.025))
    for v in [0.33, 0.45, 0.462, 0.49, 0.515]:
        h = 1e-6
        derivative = (np.log(t.theta(v + h)) - np.log(t.theta(v - h))) / (2 * h) * v
        assert derivative == pytest.approx(-float(t.gamma(v)), abs=1e-7)


def test_expansion_uses_room_length_convention():
    v = expansion_volume([773.15])[0]
    assert v / VREF == pytest.approx(
        ((1 + length_integral(500)) / (1 + length_integral(293 - 273.15))) ** 3
    )
    assert abs(v - expansion_volume([773.15], "instantaneous")[0]) > 1e-4


def test_source_rows_and_printed_footnotes():
    bk = read("boehler-kennedy-1980-table1")
    shock = read("fritz-1971-table1")
    assert len(bk) == 82 and len(set(r["run"] for r in bk)) == 19
    assert len(shock) == 43
    assert shock[16]["orientation"] == "110"
    assert shock[17]["shock_velocity_km_s"] == "5.59"
    assert shock[20]["shock_velocity_km_s"] == "5.94"
    assert sum(r["standard"] == "Cu" for r in shock) == 5


def test_hugoniot_reference_and_iterative_fixed_point(model):
    assert model.history[-1] < 1e-10
    assert model.stats()["rh_residual_equals_fit_residual_max_gpa"] < 1e-9
    assert np.all(model.e_thermal_h > 0)
    # Numerical work integral matches the analytical spline antiderivative.
    integ = model.spline.antiderivative()
    expected = -1000 * VREF * (integ(0.35 / VREF - 1) - integ(0))
    assert expected == pytest.approx(
        -1000 * quad(lambda v: float(model.p0(v)), VREF, 0.35, epsabs=1e-8)[0], abs=2e-5
    )


def test_thermodynamic_identities(model):
    p = model.properties(PREF, 600)
    assert p["ks"] / p["kt"] == pytest.approx(p["cp"] / p["cv"])
    assert p["gamma"] == pytest.approx(p["alpha"] * p["ks"] * p["v"] * 1000 / p["cp"])
    assert p["cp"] > p["cv"] > 0 and p["kt"] > 0
    # dV/dT at fixed pressure independently verifies alpha.
    dV = (model.volume(PREF, 600.01) - model.volume(PREF, 599.99)) / 0.02
    assert dV / p["v"] == pytest.approx(p["alpha"], rel=1e-6)


def test_primary_fit_does_not_read_brown_tables(monkeypatch):
    import scripts.brown_nacl_1999_reconstruct as module

    original = module.read

    def guarded(name):
        assert not name.startswith("brown-table")
        return original(name)

    monkeypatch.setattr(module, "read", guarded)
    result = Model(Settings())
    assert result.stats()["reduced_pressure_rms_gpa"]["BK"] < 0.03


def test_published_table_arithmetic_signature(model):
    result, _ = comparison(model)
    assert result["table_arithmetic_cells"] == 571
    assert result["table1_standard_to_table2_rms_gpa"] > 0.12
    assert result["table1_zp3_signature_to_table2_rms_gpa"] < 0.012


def test_source_density_pressure_are_retained_separately():
    rows = read("fritz-1971-table1")
    r = rows[0]
    p = 2.1645 * float(r["shock_velocity_km_s"]) * float(r["particle_velocity_km_s"])
    assert p != float(r["pressure_kbar"]) / 10
    # Published reduced/rounded columns are preserved, never overwritten by RH.
    assert r["pressure_kbar"] == "29" and r["density_g_cm3"] == "2.358"


def test_mode_factor_is_not_a_uniform_thermal_energy_rescaling():
    from scripts.brown_nacl_1999_reconstruct import thermal_mode_factor_check

    result = thermal_mode_factor_check()
    assert result["increments"] == 508
    assert result["best_uniform_energy_scale"] == pytest.approx(1, abs=0.002)
    rms = result["increment_rms_gpa_by_scale"]
    assert rms["standard"] < 0.006
    assert rms["one_third"] > 0.9 and rms["three"] > 2.9
