import numpy as np
from scipy.integrate import quad

from .holzapfel import Holzapfel


def Pth_modified(
    n,
    z,
    V0,
    K0,
    K_prime,
    Tr,
    V,
    QE1o,
    mE1,
    QE2o,
    mE2,
    TK,
    gamVo,
    gb,
    ao,
    m,
    mm,
    ae,
):
    """
    Original thermal pressure equation from sokolova et al. 2016.

    Parameters
    ----------
    n : float
        Number of atoms in the chemical formula
    z : float
        Atomic number of the chemical formula unit
    V0 : float
        Reference volume in [JBar^-1] (same as [cm^3/mol]/10)
    K0 : float
        Bulk modulus at reference volume [kbar]
    K_prime : float
        Bulk modulus derivative at reference volume
    Tr : float
        Reference temperature in [K] for the EOS (typically 298.15 K)
    V : float
        Volume in [JBar^-1] (same as [cm^3/mol]/10)
    QE1o : float
        Einstein characteristic temperature, Theta_1 in [K]
    mE1 : float
        The first Einstein number
    QE2o : float
        Einstein characteristic temperature, Theta_02 in [K]
    mE2 : float
        The second Einstein number
    TK : float
        Temperature in [K]
    gamVo : float
        Additive normalizing constant for the Gruneisen parameter, (given as delta in Sokolova et al. 2016)
    gb : float
        Generalized Gruneisen parameter, given as t in Sokolova et al. 2016
    ao : float
        Intrinsic anharmonicity parameter (10e-6 [K]), given as a_0 in Sokolova et al. 2016
    m : float
        Anharmonic analogue of the Grüneisen parameter
    mm : float
        Electronic analogue of the Grüneisen parameter, given as g in Sokolova et al. 2016
    ae : float
        Free electrons parameter (10e-6 [K]), given as e_0 in Sokolova et al. 2016

    Returns
    -------
    Pth : float
        Thermal pressure in [bar]
    """
    R = 8.31451

    x = V / V0  # fractional volume

    rt_eos = Holzapfel(V0=V0, K0=K0, K0_prime=K_prime, n=n, Z=z)
    Px = rt_eos.pressure(V)
    KT = rt_eos.bulk_modulus(V) * 1000  # conver kbar to bar
    kkx = rt_eos.bulk_modulus_derivative(V)

    # Equation (10) - seems not the same as in the original paper
    gamV = (-3 * KT + 2 * Px * gb + 9 * KT * kkx - 6 * gb * KT) / 6 / (
        3 * KT - 2 * Px * gb
    ) + gamVo

    # Exponent part in equation (9) - However this is not the correct value - the Excel spreadsheet
    # calculation has a more elaborate calculation included

    # expp_test = np.exp(0.5 * ao * TK * 1e6 * (V / V0) ** (1 / 3))
    expp = np.exp(I_gamV(x, gamVo, gb, rt_eos))

    # Equation (9)
    QE1 = QE1o * expp
    QE2 = QE2o * expp

    # Equation (12) for the different Einstein contributions at the temperature TK
    e1 = np.exp(QE1 / TK)
    PE1 = mE1 * R * ((QE1 / 2 + QE1 / (e1 - 1)) * gamV / V)

    e2 = np.exp(QE2 / TK)
    PE2 = mE2 * R * ((QE2 / 2 + QE2 / (e2 - 1)) * gamV / V)

    # Equation (12) for the different Einstein contributions at the reference temperature
    e1r = np.exp(QE1 / Tr)
    PE1r = mE1 * R * ((QE1 / 2 + QE1 / (e1r - 1)) * gamV / V)

    e2r = np.exp(QE2 / Tr)
    PE2r = mE2 * R * ((QE2 / 2 + QE2 / (e2r - 1)) * gamV / V)

    # Equation (12) second additive term
    Pea = 3 / 2 * n * R * ao / 1000000 * x ** (m) * (m) / V * (TK**2 - Tr**2)
    Pee = 3 / 2 * n * R * ae / 1000000 * x ** (mm) * (mm) / V * (TK**2 - Tr**2)

    Pth = PE1 + PE2 - PE1r - PE2r + Pee + Pea
    return Pth


def I_gamV(x, gamVo, gb, rt_eos):
    """
    Integral of the Gruneisen parameter over the volume ratio (from x to 1).

    Parameters
    ----------
    x : float
        Fractional volume (V/Vo)
    gamVo : float
        Additive normalizing constant for the Gruneisen parameter, (given as delta in Sokolova et al. 2016)
    gb : float
        Generalized Gruneisen parameter, given as t in Sokolova et al. 2016
    rt_eos : RT_EOS
        room temperature equation of state object used for the calculation

    Returns
    -------
    I_gamV : float
        Integral of the Gruneisen parameter over the volume ratio (from x to 1)
    """

    V0 = rt_eos.V0

    def f_gamV_x(x):
        Px_x = rt_eos.pressure(x * V0)
        KT_x = rt_eos.bulk_modulus(x * V0)
        kkx_x = rt_eos.bulk_modulus_derivative(x * V0)
        return f_gamV(x, Px_x, KT_x, kkx_x, gamVo, gb)

    I_gamV = quad(f_gamV_x, x, 1)

    return I_gamV[0]


def f_gamV(x, Px, KT, kkx, gamVo, gb):
    """
    Helper function to calculate the Gruneisen parameter at a given temperature and pressure.

    Parameters
    ----------
    x : float
        Fractional volume (V/Vo)
    Px : float
        Pressure in [bar]
    KT: float
        Bulk modulus at temperature in [kbar]
    kkx: float
        Bulk modulus derivative at temperature
    gamVo: float
        Additive normalizing constant for the Gruneisen parameter, (given as delta in Sokolova et al. 2016)
    gb: float
        Generalized Gruneisen parameter, given as t in Sokolova et al. 2016
    """
    KT = KT * 1000  # convert kbar to bar

    f_gamV_value = (
        gamVo
        + (-3 * KT + 2 * Px * gb + 9 * KT * kkx - 6 * gb * KT)
        / (6 * (3 * KT - 2 * Px * gb))
    ) / x

    return f_gamV_value
