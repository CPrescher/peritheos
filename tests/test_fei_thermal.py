import numpy as np
import pytest

from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import MieGruneisenDebye


def _gold_debye(**overrides):
    parameters = dict(
        Tr=300.0,
        theta0=170.0,
        gamma0=2.97,
        q=0.6,
        n=1.0,
    )
    parameters.update(overrides)
    return MieGruneisenDebye(Vinet(1.0, 167.0, 6.0), **parameters)


def test_integrated_gruneisen_is_the_default_debye_temperature_law():
    implicit = _gold_debye()
    explicit = _gold_debye(debye_temperature_law="integrated_gruneisen")

    assert implicit.debye_temperature_law == "integrated_gruneisen"
    assert implicit.characteristic_temperature(0.7) == pytest.approx(
        explicit.characteristic_temperature(0.7)
    )


def test_variable_exponent_law_follows_fei_published_convention():
    eos = _gold_debye(debye_temperature_law="variable_exponent")
    ratio = np.array([1.0, 0.9, 0.7])
    gamma = 2.97 * ratio**0.6
    expected = 170.0 * ratio ** (-gamma)

    assert np.allclose(eos.characteristic_temperature(ratio), expected)


def test_variable_exponent_is_not_silently_replaced_by_integrated_law():
    variable = _gold_debye(debye_temperature_law="variable_exponent")
    integrated = _gold_debye()

    assert variable.characteristic_temperature(0.7) != pytest.approx(
        integrated.characteristic_temperature(0.7)
    )


def test_debye_temperature_law_validation_and_parameter_reconstruction():
    for invalid in ("unknown", None, []):
        with pytest.raises(ValueError, match="debye_temperature_law"):
            _gold_debye(debye_temperature_law=invalid)

    eos = _gold_debye(debye_temperature_law="variable_exponent")
    reconstructed = eos.with_parameters(gamma0=3.0)

    assert "debye_temperature_law" not in eos.parameter_values()
    assert eos.configuration_values() == {"debye_temperature_law": "variable_exponent"}
    assert reconstructed.gamma0 == 3.0
    assert reconstructed.debye_temperature_law == "variable_exponent"


def test_variable_exponent_reference_thermal_pressure_and_arrays():
    eos = _gold_debye(debye_temperature_law="variable_exponent")

    assert eos.thermal_pressure(0.8, 300.0) == pytest.approx(0.0, abs=1.0e-14)
    result = eos.thermal_pressure(np.array([0.8, 0.7]), np.array([1000.0, 2000.0]))
    assert result.shape == (2,)
    assert np.all(result > 0.0)
