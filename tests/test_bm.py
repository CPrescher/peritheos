import numpy as np
from theospy.eos.bm import BM2, BM3


def test_bm2_P():
    V0 = 100
    K0 = 100
    eos = BM2(V0, K0)
    assert eos.pressure(V0) == 0
    assert eos.bulk_modulus() == K0

    V1 = V0 * 0.9
    P_test = 3 / 2 * K0 * ((V0 / V1) ** (7 / 3) - (V0 / V1) ** (5 / 3))
    assert np.isclose(eos.pressure(V1), P_test)


def test_bm2_K():
    V0 = 100
    K0 = 100
    eos = BM2(V0, K0)
    assert eos.bulk_modulus(V0) == K0

    V1 = V0 * 0.9
    K1_test = K0 / 2 * (7 * (V0 / V1) ** (7 / 3) - 5 * (V0 / V1) ** (5 / 3))
    assert np.isclose(eos.bulk_modulus(V1), K1_test)


def test_bm3_P():
    V0 = 100
    K0 = 100
    K0_prime = 5
    eos = BM3(V0, K0, K0_prime)
    assert eos.pressure(V0) == 0
    assert eos.bulk_modulus(V0) == K0

    # different formulation in anderson et al.
    V1 = V0 * 0.9
    n = V1 / V0
    first_term = -3 / 2 * K0 * (n ** (-5 / 3) - n ** (-7 / 3))
    second_term = 1 - 3 / 4 * (K0_prime - 4) * (1 - n ** (-2 / 3))
    P_test = first_term * second_term

    assert np.isclose(eos.pressure(V1), P_test)


def test_bm3_K():
    V0 = 100
    K0 = 100
    K0_prime = 5
    eos = BM3(V0, K0, K0_prime)
    assert eos.bulk_modulus(V0) == K0

    # Different formulation in anderson et al. page 168
    V1 = V0 * 0.9
    n = V1 / V0
    first_term = 0.5 * (1 - n ** (-2 / 3))
    second_term = 5 - 3 * K0_prime - 27 / 4 * (4 - K0_prime) * (1 - n ** (-2 / 3))
    K_test = K0 * n ** (-5 / 3) * (1 + first_term * second_term)
    assert np.isclose(eos.bulk_modulus(V1), K_test)
