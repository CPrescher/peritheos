import numpy as np
import pytest
from scipy.constants import Avogadro

from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import Tange2009Debye

V0 = 74.698 * Avogadro * 1.0e-25 / 4.0


def make_eos(**updates):
    parameters = {
        "Tr": 300.0,
        "theta0": 761.0,
        "gamma0": 1.442,
        "a": 0.138,
        "b": 5.4,
        "n": 2.0,
    }
    parameters.update(updates)
    return Tange2009Debye(Vinet(V0, 160.63, 4.367), **parameters)


def test_tange_gamma_and_theta_reference_identities():
    eos = make_eos()

    assert eos.gruneisen_parameter(V0) == pytest.approx(eos.gamma0)
    assert eos.characteristic_temperature(V0) == pytest.approx(eos.theta0)
    volumes = V0 * np.array([1.0, 0.9, 0.7])
    assert np.asarray(eos.gruneisen_parameter(volumes)).shape == (3,)
    assert np.asarray(eos.characteristic_temperature(volumes)).shape == (3,)


def test_tange_characteristic_temperature_satisfies_gamma_derivative():
    eos = make_eos()
    volume = 0.8 * V0
    step = 1.0e-6 * volume
    derivative = (
        np.log(eos.characteristic_temperature(volume + step))
        - np.log(eos.characteristic_temperature(volume - step))
    ) / (np.log(volume + step) - np.log(volume - step))

    assert -derivative == pytest.approx(eos.gruneisen_parameter(volume), rel=1e-9)


def test_tange_b_zero_reduces_to_constant_gamma_model():
    eos = make_eos(b=0.0)
    volume = 0.8 * V0

    assert eos.gruneisen_parameter(volume) == pytest.approx(eos.gamma0)
    assert eos.characteristic_temperature(volume) == pytest.approx(
        eos.theta0 * (volume / V0) ** (-eos.gamma0)
    )


@pytest.mark.parametrize("a", [-0.1, 1.1, np.nan])
def test_tange_rejects_invalid_a(a):
    with pytest.raises(ValueError):
        make_eos(a=a)


def test_tange_model_reconstructs_with_parameter_updates():
    eos = make_eos()
    modified = eos.with_parameters(gamma0=1.5, **{"rt_eos.K0": 161.0})

    assert isinstance(modified, Tange2009Debye)
    assert modified.gamma0 == 1.5
    assert modified.rt_eos.K0 == 161.0
    assert "q" not in eos.parameter_values()
