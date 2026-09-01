"""Direct tests for the private native extension during the additive phase."""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, _rust, get_material_document, list_material_documents
from peritheos.eos.rt import (
    BM2,
    BM3,
    BM4,
    Holzapfel,
    ModifiedTait,
    Murnaghan,
    NaturalStrain2,
    NaturalStrain3,
    NaturalStrain4,
    Vinet,
)
from peritheos.eos.thermal import (
    MieGruneisenDebye,
    MieGruneisenEinstein,
    Sokolova2016,
    ThermalModifiedTait,
)


@pytest.mark.parametrize(
    ("python_model", "native_model"),
    [
        (BM2(10.0, 120.0), _rust.RtEos.bm2(10.0, 120.0)),
        (BM3(10.0, 120.0, 4.3), _rust.RtEos.bm3(10.0, 120.0, 4.3)),
        (
            BM4(10.0, 120.0, 4.3, -0.02),
            _rust.RtEos.bm4(10.0, 120.0, 4.3, -0.02),
        ),
        (
            Murnaghan(10.0, 120.0, 4.3),
            _rust.RtEos.murnaghan(10.0, 120.0, 4.3),
        ),
        (
            ModifiedTait(10.0, 120.0, 4.3, -0.02),
            _rust.RtEos.modified_tait(10.0, 120.0, 4.3, -0.02),
        ),
        (
            NaturalStrain2(10.0, 120.0),
            _rust.RtEos.natural_strain2(10.0, 120.0),
        ),
        (
            NaturalStrain3(10.0, 120.0, 4.3),
            _rust.RtEos.natural_strain3(10.0, 120.0, 4.3),
        ),
        (
            NaturalStrain4(10.0, 120.0, 4.3, -0.02),
            _rust.RtEos.natural_strain4(10.0, 120.0, 4.3, -0.02),
        ),
        (Vinet(10.0, 120.0, 4.3), _rust.RtEos.vinet(10.0, 120.0, 4.3)),
        (
            Holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
            _rust.RtEos.holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
        ),
    ],
)
def test_native_isothermal_binding_matches_python(python_model, native_model):
    fractions = np.array([[0.55, 0.7], [0.85, 1.0]], dtype=float)
    volumes = python_model.V0 * fractions
    expected_pressures = python_model.pressure(volumes)

    assert native_model.reference_volume == python_model.V0
    assert np.allclose(
        native_model.pressure_array(volumes), expected_pressures, rtol=1.0e-12
    )
    assert np.allclose(
        native_model.bulk_modulus_array(volumes),
        python_model.bulk_modulus(volumes),
        rtol=1.0e-12,
    )
    assert np.allclose(
        native_model.volume_array(expected_pressures), volumes, rtol=1.0e-10
    )
    assert isinstance(native_model.pressure_scalar(float(volumes[0, 0])), float)
    assert native_model.pressure_array(volumes).shape == volumes.shape


def test_native_binding_preserves_error_categories():
    with pytest.raises(ValueError):
        _rust.RtEos.bm2(0.0, 100.0)

    model = _rust.RtEos.bm2(10.0, 100.0)
    with pytest.raises(ValueError):
        model.pressure_scalar(0.0)
    with pytest.raises(ValueError):
        model.volume_scalar(-100.0)


def test_large_parallel_isothermal_batches_preserve_order_shape_and_strides():
    model = _rust.RtEos.bm3(10.0, 120.0, 4.3)
    source = np.linspace(6.0, 12.0, 80_000, dtype=float).reshape(400, 200)
    volumes = source[:, ::2]
    strain = 0.5 * ((10.0 / volumes) ** (2.0 / 3.0) - 1.0)
    expected = (
        3.0
        * 120.0
        * strain
        * (1.0 + 2.0 * strain) ** 2.5
        * (1.0 + 1.5 * (4.3 - 4.0) * strain)
    )

    actual = model.pressure_array(volumes)

    assert actual.shape == volumes.shape
    assert np.allclose(actual, expected, rtol=2.0e-11, atol=1.0e-13)


def test_large_parallel_volume_batch_round_trips_in_input_order():
    model = _rust.RtEos.bm3(10.0, 120.0, 4.3)
    volumes = np.linspace(7.0, 10.5, 2_000, dtype=float).reshape(40, 50)
    pressures = model.pressure_array(volumes)

    recovered = model.volume_array(pressures)

    assert recovered.shape == volumes.shape
    assert np.allclose(recovered, volumes, rtol=1.0e-10)


def test_concurrent_parallel_batches_are_deterministic_and_thread_safe():
    model = _rust.RtEos.bm3(10.0, 120.0, 4.3)
    volumes = np.linspace(6.5, 11.0, 100_000, dtype=float)
    expected_pressure = model.pressure_array(volumes)
    expected_modulus = model.bulk_modulus_array(volumes)

    def evaluate(quantity):
        return getattr(model, f"{quantity}_array")(volumes)

    quantities = ["pressure", "bulk_modulus"] * 4
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(evaluate, quantities))

    for quantity, result in zip(quantities, results):
        expected = expected_pressure if quantity == "pressure" else expected_modulus
        assert np.array_equal(result, expected)


@pytest.mark.parametrize(
    ("python_model", "native_model", "caloric"),
    [
        (
            MieGruneisenDebye(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0),
            _rust.ThermalEos.mie_gruneisen_debye(
                _rust.RtEos.bm3(1.0, 160.0, 4.0),
                300.0,
                800.0,
                1.5,
                1.0,
                2.0,
            ),
            True,
        ),
        (
            MieGruneisenEinstein(BM3(1.0, 160.0, 4.0), 300.0, 800.0, 1.5, 1.0, 2.0),
            _rust.ThermalEos.mie_gruneisen_einstein(
                _rust.RtEos.bm3(1.0, 160.0, 4.0),
                300.0,
                800.0,
                1.5,
                1.0,
                2.0,
            ),
            True,
        ),
        (
            ThermalModifiedTait(
                ModifiedTait(1.0, 160.0, 4.0, -0.01),
                298.15,
                700.0,
                2.5e-5,
                2.0,
            ),
            _rust.ThermalEos.thermal_modified_tait(
                _rust.RtEos.modified_tait(1.0, 160.0, 4.0, -0.01),
                298.15,
                700.0,
                2.5e-5,
                2.0,
            ),
            True,
        ),
        (
            Sokolova2016(
                Holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
                298.15,
                684.0,
                0.564,
                1561.0,
                2.436,
                -0.506,
                1.085,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            _rust.ThermalEos.multi_oscillator_gruneisen(
                _rust.RtEos.holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
                298.15,
                684.0,
                0.564,
                1561.0,
                2.436,
                -0.506,
                1.085,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            False,
        ),
    ],
)
def test_native_thermal_binding_matches_python(python_model, native_model, caloric):
    volumes, temperatures = np.broadcast_arrays(
        python_model.rt_eos.V0 * np.array([[0.8], [0.9]]),
        np.array([[800.0, 1800.0]]),
    )
    expected_pressure = python_model.pressure(volumes, temperatures)

    assert np.allclose(
        native_model.evaluate_array("thermal_pressure", volumes, temperatures),
        python_model.thermal_pressure(volumes, temperatures),
        rtol=3.0e-5,
    )
    assert np.allclose(
        native_model.evaluate_array("pressure", volumes, temperatures),
        expected_pressure,
        rtol=3.0e-5,
    )
    recovered = native_model.evaluate_array("volume", expected_pressure, temperatures)
    assert np.allclose(recovered, volumes, rtol=1.0e-9)

    pressure = float(expected_pressure[0, 1])
    volume = float(volumes[0, 1])
    temperature = float(temperatures[0, 1])
    assert np.isclose(
        native_model.evaluate_scalar("temperature", pressure, volume),
        temperature,
        rtol=1.0e-9,
    )
    if caloric:
        assert np.isclose(
            native_model.evaluate_scalar("molar_heat_capacity_v", volume, temperature),
            python_model.molar_heat_capacity_v(volume, temperature),
            rtol=1.0e-9,
        )
    else:
        with pytest.raises(NotImplementedError):
            native_model.evaluate_scalar("molar_heat_capacity_v", volume, temperature)


def test_large_parallel_thermal_batch_matches_scalar_evaluation():
    model = _rust.ThermalEos.mie_gruneisen_debye(
        _rust.RtEos.bm3(1.0, 160.0, 4.0),
        300.0,
        800.0,
        1.5,
        1.0,
        2.0,
    )
    volumes = np.linspace(0.7, 1.0, 3_000, dtype=float).reshape(60, 50)
    temperatures = np.linspace(300.0, 3_000.0, 3_000, dtype=float).reshape(60, 50)

    actual = model.evaluate_array("pressure", volumes, temperatures)
    expected = np.array(
        [
            model.evaluate_scalar("pressure", float(volume), float(temperature))
            for volume, temperature in zip(volumes.flat, temperatures.flat)
        ]
    ).reshape(volumes.shape)

    assert actual.shape == volumes.shape
    assert np.array_equal(actual, expected)


def test_native_thermal_binding_enforces_reference_model_types():
    with pytest.raises(TypeError, match="ModifiedTait"):
        _rust.ThermalEos.thermal_modified_tait(
            _rust.RtEos.bm3(1.0, 160.0, 4.0), 298.15, 700.0, 2.5e-5, 2.0
        )
    with pytest.raises(ValueError, match="n is required"):
        _rust.ThermalEos.multi_oscillator_gruneisen(
            _rust.RtEos.bm3(1.0, 160.0, 4.0),
            298.15,
            684.0,
            0.564,
            1561.0,
            2.436,
            -0.506,
            1.085,
            0.0,
            0.0,
            0.0,
            0.0,
        )


def test_native_multi_oscillator_accepts_a_generic_reference_isotherm():
    model = _rust.ThermalEos.multi_oscillator_gruneisen(
        _rust.RtEos.bm3(1.0, 160.0, 4.0),
        298.15,
        684.0,
        0.564,
        1561.0,
        2.436,
        -0.506,
        1.085,
        0.0,
        0.0,
        0.0,
        0.0,
        n=2.0,
    )

    assert model.evaluate_scalar("thermal_pressure", 0.9, 298.15) == pytest.approx(0.0)


def test_every_bundled_material_record_has_an_exact_native_model():
    for identifier in list_material_documents():
        material = Material.from_eosmat(get_material_document(identifier))
        for record in material.eos_records:
            assert hasattr(record.eos, "_native"), (
                identifier,
                record.identifier,
                type(record.eos).__name__,
            )


@pytest.mark.parametrize("loss", ["linear", "soft_l1", "huber", "cauchy", "arctan"])
def test_native_least_squares_matches_scipy_robust_regression(loss):
    coordinates = np.linspace(-1.0, 1.0, 21)
    observations = 2.5 - 1.2 * coordinates + 0.05 * np.sin(np.arange(21))
    observations[10] += 3.0

    def residuals(parameters):
        return parameters[0] + parameters[1] * coordinates - observations

    initial = np.array([1.0, 0.0])
    lower = np.array([-5.0, -5.0])
    upper = np.array([5.0, 5.0])
    scipy_result = least_squares(
        residuals,
        initial,
        bounds=(lower, upper),
        x_scale="jac",
        loss=loss,
        f_scale=0.2,
        max_nfev=1000,
    )
    native_result = _rust.fit_least_squares(
        residuals,
        initial,
        lower,
        upper,
        loss=loss,
        f_scale=0.2,
        max_nfev=1000,
    )

    assert native_result.success
    assert native_result.x == pytest.approx(scipy_result.x, rel=2.0e-5, abs=2.0e-6)
    assert native_result.cost == pytest.approx(scipy_result.cost, rel=1.0e-9)
    assert native_result.fun == pytest.approx(scipy_result.fun, rel=2.0e-5, abs=3.0e-6)
    assert native_result.jac.shape == scipy_result.jac.shape


def test_native_least_squares_reports_evaluation_limit():
    result = _rust.fit_least_squares(
        lambda parameters: np.array([parameters[0] - 2.0]),
        np.array([0.0]),
        np.array([-5.0]),
        np.array([5.0]),
        max_nfev=1,
    )

    assert not result.success
    assert result.status == 0
    assert result.nfev == 1
    assert result.jac.shape == (1, 1)

    with pytest.raises(ValueError, match="dimensions"):
        _rust.fit_least_squares(
            lambda parameters: np.array([parameters[0]]),
            np.array([0.0]),
            np.array([-1.0]),
            np.array([1.0]),
            global_parameter_count=1,
        )


def test_native_least_squares_handles_differently_scaled_columns():
    def residuals(parameters):
        return np.array([1.0e-9 * (parameters[0] - 2.0), 1.0e9 * (parameters[1] - 3.0)])

    result = _rust.fit_least_squares(
        residuals,
        np.array([0.0, 0.0]),
        np.array([-10.0, -10.0]),
        np.array([10.0, 10.0]),
        max_nfev=1000,
    )

    assert result.success
    assert result.x == pytest.approx([2.0, 3.0], abs=1.0e-8)
    assert result.cost < 1.0e-10


def test_native_covariance_is_a_rank_aware_pseudoinverse():
    result = _rust.fit_least_squares(
        lambda parameters: np.full(3, parameters[0] + 2.0 * parameters[1] - 3.0),
        np.array([1.0, 1.0]),
        np.array([-10.0, -10.0]),
        np.array([10.0, 10.0]),
    )
    expected = np.linalg.pinv(result.jac.T @ result.jac, hermitian=True)

    assert result.parameter_covariance == pytest.approx(expected)


def test_native_linear_uncertainty_matches_dense_reference():
    jacobian = np.array([[1.0, 2.0], [-1.0, 0.5]])
    parameter_covariance = np.array([[4.0, 0.5], [0.5, 1.0]])
    state_variance = np.array([0.25, 0.0])

    result = _rust.linear_uncertainty(
        jacobian,
        parameter_covariance,
        state_variance,
        full_covariance=True,
    )
    expected = jacobian @ parameter_covariance @ jacobian.T
    expected += np.diag(state_variance)

    assert result.variance == pytest.approx(np.diag(expected))
    assert result.covariance == pytest.approx(expected)


@pytest.mark.parametrize(
    ("covariance", "state_variance"),
    [
        (np.array([[-1.0]]), np.array([0.0])),
        (np.array([[1.0]]), np.array([-1.0])),
    ],
)
def test_native_linear_uncertainty_rejects_invalid_variances(
    covariance, state_variance
):
    with pytest.raises(ValueError):
        _rust.linear_uncertainty(
            np.array([[1.0]]), covariance, state_variance, full_covariance=True
        )


def test_native_monte_carlo_summary_matches_numpy():
    samples = np.array(
        [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]], dtype=float
    )
    result = _rust.monte_carlo_summary(samples, confidence=0.5, full_covariance=True)

    assert result.standard_error == pytest.approx(samples.std(axis=0, ddof=1))
    assert result.lower == pytest.approx(np.quantile(samples, 0.25, axis=0))
    assert result.upper == pytest.approx(np.quantile(samples, 0.75, axis=0))
    assert result.covariance == pytest.approx(np.cov(samples, rowvar=False, ddof=1))
