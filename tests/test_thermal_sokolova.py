import pytest
import numpy as np
from peritheos.eos.thermal.sokolova2016 import (
    f_gamV,
    I_gamV,
)

import peritheos.eos._reference.sokolova2016 as sokolova2016_original
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


def test_original_thermal_pressure_calculation():
    expp = 1.09427641040855
    Pth = sokolova2016_original.P_thermal(
        n,
        z,
        V0,
        K0_kbar,
        K0_prime,
        Tr,
        x,
        expp,
        V,
        1,
        1,
        0,
        1,
        1,
        0,
        Theta_1,
        m1,
        Theta_2,
        m2,
        T,
        delta,
        0,
        t,
        a_0,
        m,
        g,
        e_0,
    )
    assert np.isclose(Pth, 148604.90047369, rtol=1e-4)


def test_compare_original_with_modified_thermal_pressure_calculation():
    # original implementation
    expp = 1.09427641040855
    P1 = sokolova2016_original.P_thermal(
        n,
        z,
        V0,
        K0_kbar,
        K0_prime,
        Tr,
        x,
        expp,
        V,
        1,
        1,
        0,
        1,
        1,
        0,
        Theta_1,
        m1,
        Theta_2,
        m2,
        T,
        delta,
        0,
        t,
        a_0,
        m,
        g,
        e_0,
    )

    # class based implementation
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    sokolova_eos = Sokolova2016(
        rt_eos, Tr, Theta_1, m1, Theta_2, m2, delta, t, a_0, m, g, e_0
    )
    P2 = sokolova_eos.thermal_pressure(V, T)
    assert np.isclose(P1, P2)


def test_f_gamV():
    V = V0 * 0.8
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)
    Px = rt_eos.pressure(V)
    KT = rt_eos.bulk_modulus(V)
    kkx = rt_eos.bulk_modulus_derivative(V)

    x = V / V0

    g0 = delta
    gb = t

    gamV0_1 = sokolova2016_original.f_gamV(x, n, z, V0, K0_kbar, K0_prime, g0, gb, 0)
    gamV0_2 = f_gamV(x, Px, KT, kkx, g0, gb)

    assert np.isclose(gamV0_1, gamV0_2)


def test_I_gamV_compression():
    V = V0 * 0.9
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)

    x = V / V0

    g0 = delta
    gb = t

    I_gamV_1 = sokolova2016_original.I_gamV(z, n, x, V0, K0_kbar, K0_prime, g0, gb, 0)
    I_gamV_2 = I_gamV(x, g0, gb, rt_eos=rt_eos)

    assert np.isclose(I_gamV_1, I_gamV_2)
    assert np.isclose(np.exp(I_gamV_2), 1.09427641040855)


def test_I_gamV_expansion():
    V = V0 * 1.1
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)

    x = V / V0

    g0 = delta
    gb = t

    I_gamV_1 = sokolova2016_original.I_gamV(z, n, x, V0, K0_kbar, K0_prime, g0, gb, 0)
    I_gamV_2 = I_gamV(x, g0, gb, rt_eos=rt_eos)

    assert np.isclose(I_gamV_1, I_gamV_2)
    assert np.isclose(np.exp(I_gamV_2), 0.910574451622214)


def test_I_gamV_multiple_volumes():
    V = np.linspace(V0 * 0.8, V0 * 0.6, 21)
    x = V / V0

    g0 = delta
    gb = t
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)

    I_gamV_1 = I_gamV(x, g0, gb, rt_eos=rt_eos)

    assert len(I_gamV_1) == len(V)


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
    )
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
    )
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
    )
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
