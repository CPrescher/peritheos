import numpy as np
from peritheos.eos.thermal import (
    Pth_original,
    Pth_modified,
    f_gamV_original,
    f_gamV,
    I_gamV_original,
    I_gamV,
)
from peritheos.eos.holzapfel import Holzapfel


# diamond parameters
V0 = 0.3414  # in JBar^-1 (same as [cm^3/mol]/10)
K0 = 4415  # in kbar
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
    print(V)
    Pth = Pth_original(
        n,
        z,
        V0,
        K0,
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
    expp = 1.09427641040855
    P1 = Pth_original(
        n,
        z,
        V0,
        K0,
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

    P2 = Pth_modified(
        n,
        z,
        V0,
        K0,
        K0_prime,
        Tr,
        V,
        Theta_1,
        m1,
        Theta_2,
        m2,
        T,
        delta,
        t,
        a_0,
        m,
        g,
        e_0,
    )
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

    gamV0_1 = f_gamV_original(x, n, z, V0, K0, K0_prime, g0, gb, 0)
    gamV0_2 = f_gamV(x, Px, KT, kkx, g0, gb)

    assert np.isclose(gamV0_1, gamV0_2)


def test_I_gamV_compression():
    V = V0 * 0.9
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)

    x = V / V0

    g0 = delta
    gb = t

    I_gamV_1 = I_gamV_original(z, n, x, V0, K0, K0_prime, g0, gb, 0)
    I_gamV_2 = I_gamV(x, g0, gb, rt_eos=rt_eos)

    assert np.isclose(I_gamV_1, I_gamV_2)
    assert np.isclose(np.exp(I_gamV_2), 1.09427641040855)


def test_I_gamV_expansion():
    V = V0 * 1.1
    rt_eos = Holzapfel(V0, K0, K0_prime, n, z)

    x = V / V0

    g0 = delta
    gb = t

    I_gamV_1 = I_gamV_original(z, n, x, V0, K0, K0_prime, g0, gb, 0)
    I_gamV_2 = I_gamV(x, g0, gb, rt_eos=rt_eos)

    assert np.isclose(I_gamV_1, I_gamV_2)
    assert np.isclose(np.exp(I_gamV_2), 0.910574451622214)
