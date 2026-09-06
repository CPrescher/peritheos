import numpy as np
import pytest

from peritheos.eos.rt import Morse3, SunMorse3, SunMorse4


@pytest.mark.parametrize("model", [Morse3, SunMorse3, SunMorse4])
def test_morse_models_recover_reference_state_and_vectorize(model):
    eos = model(10.0, 120.0, 4.3)

    assert eos.pressure(10.0) == pytest.approx(0.0)
    assert eos.bulk_modulus(10.0) == pytest.approx(120.0)
    pressures = eos.pressure(np.array([10.0, 9.0, 8.0]))
    assert pressures.shape == (3,)
    assert np.all(np.diff(pressures) > 0.0)


@pytest.mark.parametrize("model", [Morse3, SunMorse3, SunMorse4])
def test_morse_analytical_bulk_modulus_matches_pressure_derivative(model):
    eos = model(10.0, 120.0, 4.3)
    volume = 8.0
    step = 1.0e-5
    numerical = (
        -volume
        * (eos.pressure(volume + step) - eos.pressure(volume - step))
        / (2.0 * step)
    )
    assert eos.bulk_modulus(volume) == pytest.approx(numerical, rel=1.0e-9)


def test_morse3_rejects_singular_pressure_derivative():
    with pytest.raises(ValueError, match="K0_prime"):
        Morse3(10.0, 120.0, 1.0)


def test_sun_morse_reference_equations_match_direct_evaluation():
    volume = 8.465 * 0.8
    for model, n, k0, k0_prime in [
        (SunMorse3, 3.0, 147.24, 5.9298),
        (SunMorse4, 4.0, 147.3, 5.8813),
    ]:
        x = (volume / 8.465) ** (1.0 / n)
        alpha = (n * k0_prime + 1.0) / 3.0
        expected = (
            n
            * k0
            / alpha
            * (np.exp(2.0 * alpha * (1.0 - x)) - np.exp(alpha * (1.0 - x)))
        )
        assert model(8.465, k0, k0_prime).pressure(volume) == pytest.approx(expected)
