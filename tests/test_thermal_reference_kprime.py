"""Linear Kprime(T): independent equations, native/fallback parity and fits."""

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import EOSUncertainty, Material, eosmat_schema
from peritheos.eos.rt import BM2, BM3, BM4, Vinet
from peritheos.eos.thermal import ThermalReferenceStateEOS
from peritheos.errors import EosmatError, EosValidationError
from peritheos.fitting import fit_thermal_eos


def conditional_continuous_pressure(volume, temperature):
    """Fukui (2003), Eqs. 1-2/p. 513 slopes, assuming 300 K intercepts.

    This standalone equation deliberately uses no Peritheos model or staged
    material resources. It tests the mechanism, not author-fit reproduction.
    """
    delta = np.asarray(temperature) - 300.0
    v0 = 40.917 * np.exp(6.69e-5 * delta)
    k0 = 41.8 - 0.0319 * delta
    kp = 6.6 + 0.00571 * delta
    eta = (v0 / np.asarray(volume)) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def model(**updates):
    kwargs = dict(
        rt_eos=BM3(40.917, 41.8, 6.6),
        Tr=300.0,
        alpha0=6.69e-5,
        dK_dT=-0.0319,
        dKprime_dT=0.00571,
    )
    kwargs.update(updates)
    return ThermalReferenceStateEOS(**kwargs)


def document():
    """A conditional audit model; no new claim about the author's intercepts."""
    doc = {
        "format": "peritheos.material",
        "format_version": 3,
        "identifier": "conditional_brucite_test",
        "name": "Conditional equation test",
        "formula": "Mg(OH)2",
        "formula_units_per_cell": 1,
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        },
        "eos_records": [
            {
                "identifier": "conditional_brucite_linear_kprime",
                "label": "Conditional equation test",
                "default": True,
                "reference": "Fukui (2003) slopes with assumed 300 K intercepts; not author-fit reproduction",
                "eos": {
                    "type": "BM3",
                    "model": "birch_murnaghan_3",
                    "parameters": {"V0": 40.917, "K0": 41.8, "K0_prime": 6.6},
                },
                "parameter_errors": {},
                "fixed_parameters": [],
                "scientific_validation": {"status": "pending_primary_source_check"},
            }
        ],
    }
    record = doc["eos_records"][0]
    record["record_kind"] = "diagnostic"
    record["equation_kind"] = "thermal"
    record["label"] = "Conditional 300 K anchored linear Kprime(T) audit model"
    record["thermal"] = {
        "type": "AlphaKT",
        "model": "thermal_reference_state",
        "parameters": {
            "Tr": 300.0,
            "alpha0": 6.69e-5,
            "dK_dT": -0.0319,
            "dKprime_dT": 0.00571,
        },
        "parameter_errors": {"dKprime_dT": 0.00045},
    }
    return doc


def test_linear_kprime_matches_independent_brucite_equation_and_native_batches():
    eos = model()
    volumes = np.array([[40.0], [37.0], [34.0]])
    temperatures = np.array([300.0, 473.0, 673.0, 873.0])
    expected = conditional_continuous_pressure(volumes, temperatures)
    actual = eos.pressure(volumes, temperatures)
    assert actual.shape == (3, 4)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    # Explicit separate equation from the source audit, away from Tr and V0.
    assert eos.pressure(34.0, 873.0) == pytest.approx(14.051454515344124)
    assert isinstance(eos.pressure(34.0, 873.0), float)
    native = eos._native
    del eos._native
    np.testing.assert_allclose(eos.pressure(volumes, temperatures), actual)
    np.testing.assert_allclose(
        eos.bulk_modulus(volumes, temperatures),
        model().bulk_modulus(volumes, temperatures),
    )
    # Exercise the native ordered batch independently of the Python wrapper.
    v, t = np.broadcast_arrays(volumes, temperatures)
    np.testing.assert_allclose(
        native.evaluate_array("pressure", v.ravel(), t.ravel()),
        actual.ravel(),
    )


def test_linear_kprime_reference_state_derivatives_and_inverse():
    eos = model()
    temperatures = np.array([300.0, 473.0, 673.0, 873.0])
    delta = temperatures - 300.0
    v0 = 40.917 * np.exp(6.69e-5 * delta)
    k0 = 41.8 - 0.0319 * delta
    kp = 6.6 + 0.00571 * delta
    np.testing.assert_allclose(eos.pressure(v0, temperatures), 0, atol=1e-12)
    np.testing.assert_allclose(eos.bulk_modulus(v0, temperatures), k0)
    volumes = np.array([39.0, 37.0, 35.0, 34.0])
    h = 1e-4
    bulk = (
        -volumes
        * (
            eos.pressure(volumes + h, temperatures)
            - eos.pressure(volumes - h, temperatures)
        )
        / (2 * h)
    )
    np.testing.assert_allclose(eos.bulk_modulus(volumes, temperatures), bulk, rtol=1e-8)
    # dK/dP at each zero-pressure reference is K0_prime(T), not its 300 K value.
    slope = (
        eos.bulk_modulus(v0 - h, temperatures) - eos.bulk_modulus(v0 + h, temperatures)
    ) / (eos.pressure(v0 - h, temperatures) - eos.pressure(v0 + h, temperatures))
    np.testing.assert_allclose(slope, kp, rtol=1e-8)
    pressure = eos.pressure(volumes, temperatures)
    np.testing.assert_allclose(eos.volume(pressure, temperatures), volumes, rtol=1e-10)
    np.testing.assert_allclose(eos.thermal_pressure(volumes, 300), 0, atol=1e-12)
    reconstructed = eos.with_parameters(dKprime_dT=0.004)
    assert reconstructed.dKprime_dT == 0.004
    assert (
        reconstructed.parameter_values(include_reference=False)["dKprime_dT"] == 0.004
    )


@pytest.mark.parametrize(
    "reference", [BM3(40.917, 41.8, 6.6), Vinet(40.917, 41.8, 6.6)]
)
def test_linear_and_logarithmic_kprime_terms_compose(reference):
    eos = model(rt_eos=reference, kprime_log_coefficient=0.0003)
    t = 700.0
    d = t - 300.0
    expected = reference.with_parameters(
        V0=40.917 * np.exp(6.69e-5 * d),
        K0=41.8 - 0.0319 * d,
        K0_prime=6.6 + 0.00571 * d + 0.0003 * d * np.log(t / 300.0),
    )
    assert eos.pressure(35, t) == pytest.approx(expected.pressure(35))
    del eos._native
    assert eos.pressure(35, t) == pytest.approx(expected.pressure(35))


@pytest.mark.parametrize("reference", [BM2(40.917, 41.8), BM4(40.917, 41.8, 6.6, -0.1)])
def test_unsupported_reference_rejects_nonzero_slope_but_accepts_zero(reference):
    with pytest.raises(EosValidationError, match="Unsupported reference"):
        model(rt_eos=reference)
    zero = model(rt_eos=reference, dKprime_dT=0)
    assert np.isfinite(zero.pressure(39, 300))


@pytest.mark.parametrize("coefficient", [np.nan, np.inf, -np.inf])
def test_nonfinite_linear_slope_is_rejected(coefficient):
    with pytest.raises(EosValidationError):
        model(dKprime_dT=coefficient)


@pytest.mark.parametrize(
    "volume,temperature", [(0, 300), (34, 0), (np.nan, 300), (34, np.inf), (34, 2000)]
)
def test_invalid_inputs_and_nonpositive_hot_modulus(volume, temperature):
    eos = model()
    with pytest.raises(ValueError):
        eos.pressure(volume, temperature)
    del eos._native
    with pytest.raises(ValueError):
        eos.pressure(volume, temperature)


def test_optional_slope_defaults_to_existing_constant_kprime_behavior():
    reference = BM3(40.917, 41.8, 6.6)
    eos = ThermalReferenceStateEOS(reference, 300, 6.69e-5, -0.0319)
    assert eos.dKprime_dT == 0
    t = 673
    expected = BM3(40.917 * np.exp(6.69e-5 * (t - 300)), 41.8 - 0.0319 * (t - 300), 6.6)
    assert eos.pressure(35, t) == pytest.approx(expected.pressure(35))


def test_linear_kprime_schema_loading_metadata_and_roundtrip():
    doc = document()
    Draft202012Validator(eosmat_schema()).validate(doc)
    material = Material.from_eosmat(doc, require_primary_validation=False)
    record = material.eos_records[0]
    assert record.eos.dKprime_dT == 0.00571
    assert record.pressure(34, 873) == pytest.approx(14.051454515344124)
    exported = material.to_eosmat()
    restored = Material.from_eosmat(exported, require_primary_validation=False)
    assert restored.to_eosmat() == exported
    assert (
        exported["eos_records"][0]["thermal"]["parameter_errors"]["dKprime_dT"]
        == 0.00045
    )
    # Schema and Python reject using this parameter on a different thermal type.
    malformed = document()
    malformed["eos_records"][0]["thermal"]["type"] = "LinearThermalPressure"
    malformed["eos_records"][0]["thermal"]["model"] = "linear_thermal_pressure"
    assert list(Draft202012Validator(eosmat_schema()).iter_errors(malformed))
    with pytest.raises(EosmatError, match="requires AlphaKT"):
        Material.from_eosmat(malformed, require_primary_validation=False)


def test_native_fitting_and_covariance_include_linear_kprime_slope():
    eos = model()
    volume, temperature = np.meshgrid(
        [33.0, 35.0, 37.0, 39.0], [300.0, 473.0, 673.0, 873.0]
    )
    volume, temperature = volume.ravel(), temperature.ravel()
    # The expected pressure comes from the independent equation, not the fit model.
    pressure = conditional_continuous_pressure(volume, temperature)
    fit = fit_thermal_eos(
        ThermalReferenceStateEOS,
        eos.rt_eos,
        volume=volume,
        temperature=temperature,
        pressure=pressure,
        initial={"dKprime_dT": 0.003, "dK_dT": -0.02},
        fixed={
            key: value
            for key, value in eos.parameter_values(include_reference=False).items()
            if key not in {"dKprime_dT", "dK_dT"}
        },
        configuration=eos.configuration_values(),
    )
    assert fit.parameters["dKprime_dT"] == pytest.approx(0.00571, rel=1e-7)
    assert fit.parameters["dK_dT"] == pytest.approx(-0.0319, rel=1e-7)
    assert fit.covariance.shape == (2, 2)


def test_linear_slope_uncertainty_matches_analytic_bm3_sensitivity():
    eos = model()
    v, t, sigma = 34.0, 873.0, 0.00045
    eta = (40.917 * np.exp(6.69e-5 * (t - 300)) / v) ** (1 / 3)
    k = 41.8 - 0.0319 * (t - 300)
    derivative = 1.5 * k * (eta**7 - eta**5) * 0.75 * (eta**2 - 1) * (t - 300)
    uncertainty = EOSUncertainty(eos, parameter_errors={"dKprime_dT": sigma})
    prediction = uncertainty.pressure(v, t)
    assert prediction.standard_error == pytest.approx(abs(derivative) * sigma, rel=1e-7)
    assert uncertainty.pressure(v, 300).standard_error == pytest.approx(0, abs=1e-12)
