"""Independent equation and API checks for Supplemental S4 Eq. (2)."""

import numpy as np
import pytest

from peritheos.eos.rt import Vinet, Vinet3


def test_density_equation_and_analytic_bulk_modulus():
    eos = Vinet3(47.218085048721, 133.6, 6.29, 2.06, 1.65)
    density = np.array([[8.939, 10.0], [20.0, 30.0]])
    ratio = 8.939 / density
    y = 1.0 - ratio ** (1.0 / 3.0)
    # Independent transcription of the printed pressure-density equation.
    expected = (
        3
        * 133.6
        * y
        / ratio ** (2.0 / 3.0)
        * np.exp(6.29 * y + 2.06 * y**2 + 1.65 * y**3)
    )
    volume = eos.V0 * ratio
    np.testing.assert_allclose(eos.pressure(volume), expected, atol=1e-12)
    step = volume * 1e-5
    numerical_k = (
        -volume
        * (eos.pressure(volume + step) - eos.pressure(volume - step))
        / (2 * step)
    )
    np.testing.assert_allclose(eos.bulk_modulus(volume), numerical_k, rtol=2e-9)
    assert eos.pressure(eos.V0) == 0.0
    assert eos.bulk_modulus(eos.V0) == pytest.approx(eos.K0)
    assert eos.bulk_modulus_derivative(eos.V0) == pytest.approx(
        1 + 2 * eos.eta / 3, rel=1e-5
    )
    np.testing.assert_allclose(eos.volume(expected), volume, rtol=1e-9)
    assert isinstance(eos.pressure(eos.V0), float)
    assert eos.pressure(np.empty((0, 2))).shape == (0, 2)


def test_vinet_limit_and_parameter_reconstruction():
    eos = Vinet3(10.0, 160.0, 4.95, 0.0, 0.0)
    ordinary = Vinet(10.0, 160.0, 4.3)
    volumes = np.linspace(3, 11, 15)
    np.testing.assert_allclose(
        eos.pressure(volumes), ordinary.pressure(volumes), atol=1e-12
    )
    np.testing.assert_allclose(
        eos.bulk_modulus(volumes), ordinary.bulk_modulus(volumes)
    )
    assert eos.parameter_values() == dict(V0=10, K0=160, eta=4.95, beta=0, psi=0)
    updated = eos.with_parameters(beta=2.06, psi=1.65)
    assert isinstance(updated, Vinet3)
    assert updated.pressure(5) != pytest.approx(eos.pressure(5))


@pytest.mark.parametrize("name", ["V0", "K0", "eta", "beta", "psi"])
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_nonfinite_parameters(name, value):
    parameters = dict(V0=10, K0=160, eta=6.29, beta=2.06, psi=1.65)
    parameters[name] = value
    with pytest.raises(ValueError):
        Vinet3(**parameters)


@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf])
def test_invalid_volumes(value):
    eos = Vinet3(10, 160, 6.29, 2.06, 1.65)
    for method in (eos.pressure, eos.bulk_modulus):
        with pytest.raises(ValueError):
            method(np.array([10, value]))


@pytest.mark.parametrize("name", ["V0", "K0"])
@pytest.mark.parametrize("value", [0, -1])
def test_nonpositive_reference_parameters(name, value):
    parameters = dict(V0=10, K0=160, eta=6.29, beta=2.06, psi=1.65)
    parameters[name] = value
    with pytest.raises(ValueError):
        Vinet3(**parameters)
