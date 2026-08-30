"""Direct tests for the private native extension during the additive phase."""

import numpy as np
import pytest

from peritheos import _rust
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
            _rust.ThermalEos.sokolova2016(
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


def test_native_thermal_binding_enforces_reference_model_types():
    with pytest.raises(TypeError, match="ModifiedTait"):
        _rust.ThermalEos.thermal_modified_tait(
            _rust.RtEos.bm3(1.0, 160.0, 4.0), 298.15, 700.0, 2.5e-5, 2.0
        )
    with pytest.raises(TypeError, match="Holzapfel"):
        _rust.ThermalEos.sokolova2016(
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
