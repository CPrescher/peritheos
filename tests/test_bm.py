"""
Tests for the Birch-Murnaghan equation of state implementations.

This module contains tests for the 2nd and 3rd order Birch-Murnaghan
equations of state, verifying both pressure and bulk modulus calculations
against analytical solutions.

References:
-----------
1. Birch, F. (1947). Finite elastic strain of cubic crystals. Physical Review, 71(11), 809.
2. Anderson, O.L. (2000). Equations of State of Solids for Geophysics and Ceramic Science,
   Oxford University Press, Oxford, UK.
"""

import numpy as np
from theospy.eos.bm import BM2, BM3


def test_bm2_P():
    """
    Test the pressure calculation for the 2nd-order Birch-Murnaghan EOS.
    
    This test verifies:
    1. Pressure is zero at the reference volume
    2. Pressure calculation matches the analytical formula at a compressed volume
    
    The 2nd-order BM EOS pressure is given by:
    P(V) = (3/2) * K0 * [(V0/V)^(7/3) - (V0/V)^(5/3)]
    """
    V0 = 100  # Reference volume
    K0 = 100  # Bulk modulus at reference volume
    eos = BM2(V0, K0)
    
    # Test pressure at reference volume (should be zero)
    assert eos.pressure(V0) == 0

    # Test pressure at a compressed volume (90% of V0)
    V1 = V0 * 0.9
    P_test = 3 / 2 * K0 * ((V0 / V1) ** (7 / 3) - (V0 / V1) ** (5 / 3))
    assert np.isclose(eos.pressure(V1), P_test)


def test_bm2_K():
    """
    Test the bulk modulus calculation for the 2nd-order Birch-Murnaghan EOS.
    
    This test verifies:
    1. Bulk modulus equals K0 at the reference volume
    2. Bulk modulus calculation matches the analytical formula at a compressed volume
    
    The 2nd-order BM EOS bulk modulus is given by:
    K(V) = K0/2 * [7(V0/V)^(7/3) - 5(V0/V)^(5/3)]
    """
    V0 = 100  # Reference volume
    K0 = 100  # Bulk modulus at reference volume
    eos = BM2(V0, K0)
    
    # Test bulk modulus at reference volume (should equal K0)
    assert eos.bulk_modulus(V0) == K0

    # Test bulk modulus at a compressed volume (90% of V0)
    V1 = V0 * 0.9
    K1_test = K0 / 2 * (7 * (V0 / V1) ** (7 / 3) - 5 * (V0 / V1) ** (5 / 3))
    assert np.isclose(eos.bulk_modulus(V1), K1_test)


def test_bm3_P():
    """
    Test the pressure calculation for the 3rd-order Birch-Murnaghan EOS.
    
    This test verifies:
    1. Pressure is zero at the reference volume
    2. Pressure calculation matches the analytical formula at a compressed volume
    
    The 3rd-order BM EOS pressure is given by:
    P(V) = (3/2) * K0 * [(V0/V)^(7/3) - (V0/V)^(5/3)] * 
           {1 + (3/4) * (K0' - 4) * [(V0/V)^(2/3) - 1]}
    
    This test uses an alternative formulation from Anderson (2000):
    P(V) = -(3/2) * K0 * (n^(-5/3) - n^(-7/3)) * 
           [1 - (3/4) * (K0' - 4) * (1 - n^(-2/3))]
    
    where n = V/V0 is the relative volume.
    """
    V0 = 100  # Reference volume
    K0 = 100  # Bulk modulus at reference volume
    K0_prime = 5  # Pressure derivative of bulk modulus
    eos = BM3(V0, K0, K0_prime)
    
    # Test pressure at reference volume (should be zero)
    assert eos.pressure(V0) == 0

    # Test pressure at a compressed volume (90% of V0)
    # Using the formulation from Anderson et al.
    V1 = V0 * 0.9
    n = V1 / V0  # Relative volume
    first_term = -3 / 2 * K0 * (n ** (-5 / 3) - n ** (-7 / 3))
    second_term = 1 - 3 / 4 * (K0_prime - 4) * (1 - n ** (-2 / 3))
    P_test = first_term * second_term

    assert np.isclose(eos.pressure(V1), P_test)


def test_bm3_K():
    """
    Test the bulk modulus calculation for the 3rd-order Birch-Murnaghan EOS.
    
    This test verifies:
    1. Bulk modulus equals K0 at the reference volume
    2. Bulk modulus calculation matches the analytical formula at a compressed volume
    
    The test uses the formulation from Anderson (2000), page 168:
    K(V) = K0 * n^(-5/3) * [1 + 0.5 * (1 - n^(-2/3)) * 
           (5 - 3K0' - (27/4) * (4 - K0') * (1 - n^(-2/3)))]
    
    where n = V/V0 is the relative volume.
    """
    V0 = 100  # Reference volume
    K0 = 100  # Bulk modulus at reference volume
    K0_prime = 5  # Pressure derivative of bulk modulus
    eos = BM3(V0, K0, K0_prime)
    
    # Test bulk modulus at reference volume (should equal K0)
    assert eos.bulk_modulus(V0) == K0

    # Test bulk modulus at a compressed volume (90% of V0)
    # Using the formulation from Anderson et al. page 168
    V1 = V0 * 0.9
    n = V1 / V0  # Relative volume
    first_term = 0.5 * (1 - n ** (-2 / 3))
    second_term = 5 - 3 * K0_prime - 27 / 4 * (4 - K0_prime) * (1 - n ** (-2 / 3))
    K_test = K0 * n ** (-5 / 3) * (1 + first_term * second_term)
    assert np.isclose(eos.bulk_modulus(V1), K_test)
