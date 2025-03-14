import numpy as np
from peritheos.eos.thermal import Pth_original, Pth_modified


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


def test_original():
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


def test_compare_original_with_modified():
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
        x,
        expp,
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