"""Tests for equation-of-state diagnostic transforms."""

import numpy as np
import pytest

from peritheos.diagnostics import (
    BirchMurnaghanFiniteStrainDiagnostic,
    birch_murnaghan_finite_strain_diagnostic,
)
from peritheos.eos.rt import BM2, BM3
from peritheos.errors import EosValidationError
from peritheos.fitting import fit_rt_eos


def test_bm2_has_constant_normalized_stress_and_omits_reference_point():
    model = BM2(V0=10.0, K0=120.0)
    volumes = np.array([10.0, 9.8, 9.2, 8.5])
    diagnostic = birch_murnaghan_finite_strain_diagnostic(
        volumes, model.pressure(volumes), model=model
    )

    assert diagnostic.reference_volume == 10.0
    assert isinstance(diagnostic, BirchMurnaghanFiniteStrainDiagnostic)
    assert diagnostic.omitted_indices.tolist() == [0]
    assert diagnostic.strain.shape == (3,)
    assert diagnostic.normalized_stress == pytest.approx(120.0)


def test_bm3_normalized_stress_slope_encodes_k0_prime():
    model = BM3(V0=10.0, K0=160.0, K0_prime=5.2)
    volumes = np.linspace(9.8, 7.5, 12)
    diagnostic = birch_murnaghan_finite_strain_diagnostic(
        volumes, model.pressure(volumes), model=model
    )

    slope, intercept = np.polyfit(diagnostic.strain, diagnostic.normalized_stress, 1)
    assert intercept == pytest.approx(model.K0)
    assert slope == pytest.approx(1.5 * model.K0 * (model.K0_prime - 4.0))

    curve_strain, curve_stress = diagnostic.model_curve(points=41)
    assert curve_strain[0] == 0.0
    assert curve_stress[0] == pytest.approx(model.K0)
    assert curve_stress == pytest.approx(
        model.K0 + 1.5 * model.K0 * (model.K0_prime - 4.0) * curve_strain
    )


def test_measurement_errors_are_propagated_to_both_axes():
    model = BM3(V0=10.0, K0=150.0, K0_prime=4.5)
    volume = np.array([9.5])
    pressure = np.asarray(model.pressure(volume))
    volume_sigma = np.array([0.01])
    pressure_sigma = np.array([0.2])
    diagnostic = birch_murnaghan_finite_strain_diagnostic(
        volume,
        pressure,
        model=model,
        volume_sigma=volume_sigma,
        pressure_sigma=pressure_sigma,
    )

    f = diagnostic.strain[0]
    F = diagnostic.normalized_stress[0]
    factor = 3.0 * f * (1.0 + 2.0 * f) ** 2.5
    df_dv = -(1.0 + 2.0 * f) / (3.0 * volume[0])
    dF_df = -F * (1.0 / f + 5.0 / (1.0 + 2.0 * f))
    expected_f_sigma = abs(df_dv) * volume_sigma[0]
    expected_F_sigma = np.hypot(
        pressure_sigma[0] / factor, dF_df * df_dv * volume_sigma[0]
    )

    assert diagnostic.strain_standard_error == pytest.approx(expected_f_sigma)
    assert diagnostic.normalized_stress_standard_error == pytest.approx(
        expected_F_sigma
    )


def test_fitted_model_can_be_attached_to_diagnostic():
    expected = BM3(V0=10.0, K0=140.0, K0_prime=4.6)
    volumes = np.linspace(9.8, 8.0, 16)
    pressures = expected.pressure(volumes)
    result = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        initial={"V0": 10.1, "K0": 130.0, "K0_prime": 4.0},
    )

    diagnostic = birch_murnaghan_finite_strain_diagnostic(
        volumes, pressures, model=result.model
    )

    assert diagnostic.model is result.model
    assert diagnostic.reference_volume == pytest.approx(result.model.V0)
    assert diagnostic.normalized_stress.size == volumes.size


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({}, "reference_volume is required"),
        ({"reference_volume": 0.0}, "greater than zero"),
        ({"reference_volume": 10.0, "pressure_sigma": 0.0}, "greater than zero"),
    ],
)
def test_invalid_diagnostic_inputs_raise(kwargs, message):
    with pytest.raises(EosValidationError, match=message):
        birch_murnaghan_finite_strain_diagnostic([9.0], [10.0], **kwargs)


def test_all_reference_volume_observations_raise():
    with pytest.raises(EosValidationError, match="every volume equals V0"):
        birch_murnaghan_finite_strain_diagnostic(
            [10.0, 10.0], [0.0, 0.0], reference_volume=10.0
        )


def test_model_curve_requires_model_and_at_least_two_points():
    diagnostic = birch_murnaghan_finite_strain_diagnostic(
        [9.0], [10.0], reference_volume=10.0
    )
    with pytest.raises(EosValidationError, match="model is required"):
        diagnostic.model_curve()

    model_diagnostic = birch_murnaghan_finite_strain_diagnostic(
        [9.0], [10.0], model=BM2(10.0, 100.0)
    )
    with pytest.raises(EosValidationError, match="greater than one"):
        model_diagnostic.model_curve(points=1)
