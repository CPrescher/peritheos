import pytest
import numpy as np
from peritheos.eos.thermal.sokolova2016 import (
    f_gamV,
    I_gamV,
)

from peritheos.eos.rt.holzapfel import Holzapfel
from peritheos.eos.thermal.sokolova2016 import Sokolova2016

# diamond parameters
V0 = 0.3414  # in JBar^-1 (same as [cm^3/mol]/10)
K0 = 441.5  # in GPa
K0_kbar = (
    K0 * 10
)  # in kbar (original code from Sokolova et al. 2016 works in kbar for K0)
K0_prime = 3.9
Theta_1 = 684  # in K
m1 = 0.564
Theta_2 = 1561  # in K
m2 = 2.436
delta = -0.506
t = 1.085
a_0 = 0
m = 0
e_0 = 0
g = 0

n = 1
z = 6
Tr = 298.15  # in K - Reference temperature

V = V0 * 0.9
T = 3000  # in K
x = V / V0


def test_diamond_thermal_pressure_regression():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    assert np.isclose(
        sokolova_eos.thermal_pressure(V, T), 14.860490047369, rtol=1e-4
    )


def test_f_gamV():
    V = V0 * 0.8
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    Px = rt_eos.pressure(V)
    KT = rt_eos.bulk_modulus(V)
    kkx = rt_eos.bulk_modulus_derivative(V)

    x = V / V0

    g0 = delta
    gb = t

    assert np.isclose(f_gamV(x, Px, KT, kkx, g0, gb), 0.879820464953846)


def test_I_gamV_compression():
    V = V0 * 0.9
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)

    x = V / V0

    g0 = delta
    gb = t

    integral = I_gamV(x, g0, gb, rt_eos=rt_eos)

    assert np.isclose(np.exp(integral), 1.09427641040855)


def test_I_gamV_expansion():
    V = V0 * 1.1
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)

    x = V / V0

    g0 = delta
    gb = t

    integral = I_gamV(x, g0, gb, rt_eos=rt_eos)

    assert np.isclose(np.exp(integral), 0.910574451622214)


def test_I_gamV_multiple_volumes():
    V = np.linspace(V0 * 0.8, V0 * 0.6, 21)
    x = V / V0

    g0 = delta
    gb = t
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)

    I_gamV_1 = I_gamV(x, g0, gb, rt_eos=rt_eos)

    assert len(I_gamV_1) == len(V)


def test_complete_pressure_terms_regression():
    """Exercise pressure terms that were previously omitted from the class API."""
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    parameters = {
        "beta": 0.35,
        "QBo": 480.0,
        "d": 2.4,
        "mb": 0.75,
        "QB1o": 1120.0,
        "d1": 1.6,
        "mb1": 0.4,
        "a_0": 5.2,
        "m": 1.3,
        "g": 0.8,
        "e_0": 2.7,
    }
    eos = Sokolova2016(
        rt_eos,
        Tr,
        Theta_1,
        m1,
        Theta_2,
        m2,
        delta,
        t,
        parameters["a_0"],
        parameters["m"],
        parameters["g"],
        parameters["e_0"],
        beta=parameters["beta"],
        QBo=parameters["QBo"],
        d=parameters["d"],
        mb=parameters["mb"],
        QB1o=parameters["QB1o"],
        d1=parameters["d1"],
        mb1=parameters["mb1"],
    )

    expected_pressures = (
        (0.72, 900.0, 3.097984574947913),
        (0.88, 2400.0, 17.277569642979014),
        (1.05, 700.0, 2.6221802405808816),
    )
    for volume_ratio, temperature, expected_pressure in expected_pressures:
        volume = V0 * volume_ratio

        assert np.isclose(
            eos.thermal_pressure(volume, temperature),
            expected_pressure,
            rtol=1e-4,
            atol=1e-7,
        )


@pytest.mark.parametrize(
    ("keyword", "value"),
    [
        ("QBo", 0),
        ("d", -1),
        ("mb", -0.1),
        ("QB1o", np.nan),
        ("d1", 0),
        ("mb1", -1),
        ("beta", np.inf),
    ],
)
def test_complete_pressure_parameters_are_validated(keyword, value):
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    kwargs = {keyword: value}

    with pytest.raises(ValueError):
        Sokolova2016(
            rt_eos,
            Tr,
            Theta_1,
            m1,
            Theta_2,
            m2,
            delta,
            t,
            a_0,
            m,
            g,
            e_0,
            **kwargs,
        )


def test_thermal_pressure_multiple_volumes_single_temperature():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    V = np.linspace(V0 * 0.8, V0 * 0.7, 11)
    T = 3000

    thermal_pressures = sokolova_eos.thermal_pressure(V, T)
    assert len(thermal_pressures) == len(V)

    expected_thermal_pressures = np.array(
        [
            143576,
            143214,
            142876,
            142564,
            142277,
            142015,
            141778,
            141566,
            141380,
            141219,
            141084,
        ]
    ) / 10000
    assert np.allclose(thermal_pressures, expected_thermal_pressures, rtol=1e-4)


def test_thermal_pressure_single_volume_multiple_temperatures():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    V = V0 * 0.8
    T = np.linspace(2000, 4000, 11)

    thermal_pressures = sokolova_eos.thermal_pressure(V, T)
    assert len(thermal_pressures) == len(T)

    expected_thermal_pressures = np.array(
        [
            81856,
            94020,
            106295,
            118658,
            131089,
            143576,
            156108,
            168678,
            181279,
            193906,
            206556,
        ]
    ) / 10000
    assert np.allclose(thermal_pressures, expected_thermal_pressures, rtol=1e-4)


def test_thermal_pressure_multiple_volumes_multiple_temperature():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    V = np.linspace(V0 * 0.8, V0 * 0.7, 11)
    T = np.linspace(2000, 4000, 11)

    thermal_pressures = sokolova_eos.thermal_pressure(V, T)

    assert len(thermal_pressures) == len(V)

    expected_thermal_pressures = np.array(
        [
            81856,
            93690,
            105634,
            117674,
            129802,
            142015,
            154312,
            166694,
            179164,
            191726,
            204387,
        ]
    ) / 10000
    assert np.allclose(thermal_pressures, expected_thermal_pressures, rtol=1e-4)


def test_thermal_pressure_multiple_volumes_multiple_temperatures_different_number():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    V = np.linspace(V0 * 0.8, V0 * 0.7, 11)
    T = np.linspace(2000, 4000, 21)

    with pytest.raises(ValueError):
        sokolova_eos.thermal_pressure(V, T)


def test_total_pressure_uses_gpa_consistently():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    thermal_pressure = sokolova_eos.thermal_pressure(V, T)
    expected = rt_eos.pressure(V) + thermal_pressure

    assert np.isclose(sokolova_eos.pressure(V, T), expected)
    assert np.isclose(expected, 71.8816, rtol=1e-5)


def test_thermal_volume_pressure_round_trip():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )
    pressure = sokolova_eos.pressure(V, T)

    calculated_volume = sokolova_eos.calculate_volume(pressure, T)

    assert np.isclose(calculated_volume, V, rtol=1e-10)


def test_thermal_volume_pressure_round_trip_arrays():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )
    volumes = np.array([0.95, 0.85, 0.75]) * V0
    temperatures = np.array([500.0, 1500.0, 3000.0])
    pressures = sokolova_eos.pressure(volumes, temperatures)

    calculated_volumes = sokolova_eos.volume(pressures, temperatures)

    assert np.allclose(calculated_volumes, volumes, rtol=1e-10)


def test_zero_pressure_volume_expands_above_reference_temperature():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    expanded_volume = sokolova_eos.volume(0.0, 3000.0)

    assert expanded_volume > V0
    assert np.isclose(sokolova_eos.pressure(expanded_volume, 3000.0), 0.0, atol=1e-8)


@pytest.mark.parametrize("temperature", [0, -1, np.nan, np.inf])
def test_thermal_pressure_rejects_invalid_temperature(temperature):
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    with pytest.raises(ValueError):
        sokolova_eos.thermal_pressure(V, temperature)


def test_thermal_pressure_is_stable_at_low_temperature():
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )

    assert np.isfinite(sokolova_eos.thermal_pressure(V, 1.0))
