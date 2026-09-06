"""Tests for the generalized Rydberg--Stacey equation of state."""

import numpy as np
import pytest

from peritheos.eos.rt import RydbergStacey, Vinet
from peritheos.utils import derivative


def test_reference_state_and_array_evaluation():
    eos = RydbergStacey(47.218, 133.6, 5.2035, 2.0)
    volumes = np.array([47.218, 0.9 * 47.218, 0.7 * 47.218])

    assert eos.pressure(volumes)[0] == pytest.approx(0.0)
    assert eos.bulk_modulus(volumes)[0] == pytest.approx(133.6)
    assert np.all(eos.pressure(volumes)[1:] > 0.0)
    assert np.all(eos.bulk_modulus(volumes)[1:] > 133.6)


def test_analytic_bulk_modulus_matches_pressure_derivative():
    eos = RydbergStacey(29.468, 352.6, 4.411, 5.0 / 3.0)
    volumes = np.linspace(0.55, 1.1, 12) * eos.V0
    expected = -volumes * derivative(eos.pressure, volumes)

    assert np.allclose(eos.bulk_modulus(volumes), expected, rtol=2.0e-8)
    assert eos.bulk_modulus_derivative(eos.V0) == pytest.approx(
        eos.K0_prime, rel=1.0e-8
    )


def test_two_thirds_limit_is_vinet():
    generalized = RydbergStacey(10.0, 160.0, 4.3, 2.0 / 3.0)
    vinet = Vinet(10.0, 160.0, 4.3)
    volumes = np.linspace(0.45, 1.2, 20) * 10.0

    assert np.allclose(generalized.pressure(volumes), vinet.pressure(volumes))
    assert np.allclose(generalized.bulk_modulus(volumes), vinet.bulk_modulus(volumes))


@pytest.mark.parametrize("value", [0.0, -1.0, np.nan, np.inf])
def test_invalid_reference_parameters_are_rejected(value):
    with pytest.raises(ValueError):
        RydbergStacey(value, 160.0, 4.0, 5.0 / 3.0)
    with pytest.raises(ValueError):
        RydbergStacey(10.0, value, 4.0, 5.0 / 3.0)


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_nonfinite_derivatives_are_rejected(value):
    with pytest.raises(ValueError):
        RydbergStacey(10.0, 160.0, value, 5.0 / 3.0)
    with pytest.raises(ValueError):
        RydbergStacey(10.0, 160.0, 4.0, value)
