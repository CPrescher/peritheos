import numpy as np
from scipy.integrate import quad

from .holzapfel import Holzapfel


def Pth_original(
    n,
    z,
    Vo,
    Ko,
    kk,
    Tr,
    x,
    expp,
    V,
    QBo,
    d,
    mb,
    QB1o,
    d1,
    mb1,
    QE1o,
    mE1,
    QE2o,
    mE2,
    TK,
    gamVo,
    beta,
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
    Vo : float
        Reference volume in [JBar^-1] (same as [cm^3/mol]/10)
    Ko : float
        Bulk modulus at reference volume [kbar]
    kk : float
        Bulk modulus derivative at reference volume
    Tr : float
        Reference temperature in [K] for the EOS (typically 298.15 K)
    x : float
        Fractional volume (V/Vo)
    expp : float
        unused in function
    V : float
        Volume in [JBar^-1] (same as [cm^3/mol]/10)
    QBo : float
        should be 1
    d : float
        should be 1
    mb : float
        should be 0
    QB1o : float
        should be 1
    d1 : float
        should be 1
    mb1 : float
        should be 0
    QE1o : float
        Einstein characteristic temperature, Theta_1 in [K]
    mE1 : float
        The first Einstein number
    QE2o : float
        Einstain characteristic temperature, Theta_02 in [K]
    mE2 : float
        The second Einstein number
    TK : float
        Temperature in [K]
    gamVo : float
        Additive normalizing constant for the Gruneisen parameter, (given as delta in Sokolova et al. 2016)
    beta : float
        Should be 0
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

    fw = -np.log(3 * Ko / 10 / (1003.6 * (z * n / (Vo * 10)) ** (5 / 3)))
    ff = x ** (1 / 3)
    aa = 1.5 * (kk - 3) - fw
    Px = (
        3
        * Ko
        * 1000
        * np.exp(fw * (1 - ff))
        * (1 / ff**5 - 1 / ff**4)
        * (1 + aa * ff - aa * ff**2)
    )
    KT = (
        Ko
        * 1000
        / ff**6
        * np.exp(fw * (1 - ff))
        * (
            (-5 / ff**2 + 4 / ff) * (1 + aa * ff - aa * ff**2)
            + (1 / ff - 1) * (1 + aa * ff - aa * ff**2) * (-fw)
            + (1 / ff - 1) * (aa - 2 * aa * ff)
        )
        * (-x)
    )
    ex = (
        3
        / ff**4
        * Ko
        * 1000
        * np.exp(fw * (1 - ff))
        * (
            (-5 / ff**2 + 4 / ff) * (1 + aa * ff - aa * ff**2)
            - (1 / ff - 1) * (1 + aa * ff - aa * ff**2) * fw
            + (1 / ff - 1) * (aa - 2 * aa * ff)
        )
    )
    ex1 = (
        1
        / ff**3
        * Ko
        * 1000
        * np.exp(fw * (1 - ff))
        * fw
        * (
            (-5 / ff**2 + 4 / ff) * (1 + aa * ff - aa * ff**2)
            - (1 / ff - 1) * (1 + aa * ff - aa * ff**2) * fw
            + (1 / ff - 1) * (aa - 2 * aa * ff)
        )
    )
    ex2 = (
        Ko
        * 1000
        * np.exp(fw * (1 - ff))
        * (
            (10 / ff**3 - 4 / ff**2) * (1 + aa * ff - aa * ff**2)
            + (-5 / ff**2 + 4 / ff) * (aa - 2 * aa * ff)
            + fw / ff**2 * (1 + aa * ff - aa * ff**2)
            - (1 / ff - 1) * (aa - 2 * aa * ff) * fw
            - (aa - 2 * aa * ff) / ff**2
            - 2 * aa * (1 / ff - 1)
        )
        / ff**3
    )
    kkx = (ex + ex1 - ex2) / (-KT / ff) / 3

    gt = gb - beta * x ** (1 / 3)
    gtx = -beta / 3 * x ** (-2 / 3)
    gamV = (-3 * KT + 2 * Px * gt + 9 * KT * kkx - 6 * gt * KT) / 6 / (
        3 * KT - 2 * Px * gt
    ) + gamVo
    ff2 = (x + 0.00001) ** (1 / 3)
    KT2 = (
        Ko
        * 1000
        / ff2**6
        * np.exp(fw * (1 - ff2))
        * (
            (-5 / ff2**2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2**2)
            + (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2**2) * (-fw)
            + (1 / ff2 - 1) * (aa - 2 * aa * ff2)
        )
        * (-(x + 0.00001))
    )
    ex = (
        3
        / ff2**4
        * Ko
        * 1000
        * np.exp(fw * (1 - ff2))
        * (
            (-5 / ff2**2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2**2)
            - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2**2) * fw
            + (1 / ff2 - 1) * (aa - 2 * aa * ff2)
        )
    )
    ex1 = (
        1
        / ff2**3
        * Ko
        * 1000
        * np.exp(fw * (1 - ff2))
        * fw
        * (
            (-5 / ff2**2 + 4 / ff2) * (1 + aa * ff2 - aa * ff2**2)
            - (1 / ff2 - 1) * (1 + aa * ff2 - aa * ff2**2) * fw
            + (1 / ff2 - 1) * (aa - 2 * aa * ff2)
        )
    )
    ex2 = (
        Ko
        * 1000
        * np.exp(fw * (1 - ff2))
        * (
            (10 / ff2**3 - 4 / ff2**2) * (1 + aa * ff2 - aa * ff2**2)
            + (-5 / ff2**2 + 4 / ff2) * (aa - 2 * aa * ff2)
            + fw / ff2**2 * (1 + aa * ff2 - aa * ff2**2)
            - (1 / ff2 - 1) * (aa - 2 * aa * ff2) * fw
            - (aa - 2 * aa * ff2) / ff2**2
            - 2 * aa * (1 / ff2 - 1)
        )
        / ff2**3
    )
    kkx2 = (ex + ex1 - ex2) / (-KT2 / ff2) / 3
    ff1 = (x - 0.00001) ** (1 / 3)
    KT1 = (
        Ko
        * 1000
        / ff1**6
        * np.exp(fw * (1 - ff1))
        * (
            (-5 / ff1**2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1**2)
            + (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1**2) * (-fw)
            + (1 / ff1 - 1) * (aa - 2 * aa * ff1)
        )
        * (-(x - 0.00001))
    )
    ex = (
        3
        / ff1**4
        * Ko
        * 1000
        * np.exp(fw * (1 - ff1))
        * (
            (-5 / ff1**2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1**2)
            - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1**2) * fw
            + (1 / ff1 - 1) * (aa - 2 * aa * ff1)
        )
    )
    ex1 = (
        1
        / ff1**3
        * Ko
        * 1000
        * np.exp(fw * (1 - ff1))
        * fw
        * (
            (-5 / ff1**2 + 4 / ff1) * (1 + aa * ff1 - aa * ff1**2)
            - (1 / ff1 - 1) * (1 + aa * ff1 - aa * ff1**2) * fw
            + (1 / ff1 - 1) * (aa - 2 * aa * ff1)
        )
    )
    ex2 = (
        Ko
        * 1000
        * np.exp(fw * (1 - ff1))
        * (
            (10 / ff1**3 - 4 / ff1**2) * (1 + aa * ff1 - aa * ff1**2)
            + (-5 / ff1**2 + 4 / ff1) * (aa - 2 * aa * ff1)
            + fw / ff1**2 * (1 + aa * ff1 - aa * ff1**2)
            - (1 / ff1 - 1) * (aa - 2 * aa * ff1) * fw
            - (aa - 2 * aa * ff1) / ff1**2
            - 2 * aa * (1 / ff1 - 1)
        )
        / ff1**3
    )
    kkx1 = (ex + ex1 - ex2) / (-KT1 / ff1) / 3

    dkkdx = (kkx2 - kkx1) / (0.00002)

    qV = (
        -1
        / 2
        * (
            -6 * KT * kkx**2 * Px * gt / x
            + 6 * KT * dkkdx * Px * gt
            + 6 * KT * kkx * KT * gt / x
            - 6 * KT * kkx * Px * gtx
            + 4 * gt**2 * KT * kkx * Px / x
            - 4 * gt**2 * KT**2 / x
            - 9 * KT**2 * dkkdx
            + 6 * gtx * KT**2
        )
        / (3 * KT - 2 * Px * gt) ** 2
        * x
        / gamV
    )
    a = ao / 1000000 * x**m

    QB = QBo * expp
    QB1 = QB1o * expp
    QE1 = QE1o * expp
    QE2 = QE2o * expp

    ggr = d * np.log(1 + QB / Tr / d)
    br = 1 / (np.exp(ggr) - 1)
    ex1r = np.exp(QB / Tr)
    PBr = (
        mb * R * ((QB * (d - 1) / 2 / d + Tr * QB * d * br / (Tr * d + QB)) * gamV / V)
    )

    gg = d * np.log(1 + QB / TK / d)
    b = 1 / (np.exp(gg) - 1)
    ex1 = np.exp(QB / TK)
    PB = mb * R * ((QB * (d - 1) / 2 / d + TK * QB * d * b / (TK * d + QB)) * gamV / V)

    gg1r = d1 * np.log(1 + QB1 / Tr / d1)
    b1r = 1 / (np.exp(gg1r) - 1)
    ex2r = np.exp(QB1 / Tr)
    PB1r = (
        mb1
        * R
        * ((QB1 * (d1 - 1) / 2 / d1 + Tr * QB1 * d1 * b1r / (Tr * d1 + QB1)) * gamV / V)
    )

    gg1 = d1 * np.log(1 + QB1 / TK / d1)
    b1 = 1 / (np.exp(gg1) - 1)
    ex2 = np.exp(QB1 / TK)
    PB1 = (
        mb1
        * R
        * ((QB1 * (d1 - 1) / 2 / d1 + TK * QB1 * d1 * b1 / (TK * d1 + QB1)) * gamV / V)
    )

    e1 = np.exp(QE1 / TK)
    PE1 = mE1 * R * ((QE1 / 2 + QE1 / (e1 - 1)) * gamV / V)

    e2 = np.exp(QE2 / TK)
    PE2 = mE2 * R * ((QE2 / 2 + QE2 / (e2 - 1)) * gamV / V)

    e1r = np.exp(QE1 / Tr)
    PE1r = mE1 * R * ((QE1 / 2 + QE1 / (e1r - 1)) * gamV / V)

    e2r = np.exp(QE2 / Tr)
    PE2r = mE2 * R * ((QE2 / 2 + QE2 / (e2r - 1)) * gamV / V)
    Pea = 3 / 2 * n * R * ao / 1000000 * x ** (m) * (m) / V * (TK**2 - Tr**2)
    Pee = 3 / 2 * n * R * ae / 1000000 * x ** (mm) * (mm) / V * (TK**2 - Tr**2)

    Pth = PB + PB1 + PE1 + PE2 - PBr - PB1r - PE1r - PE2r + Pee + Pea
    return Pth


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


def I_gamV(x, gamVo, gb, rt_eos, e=0.001):
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
    e : float
        Tolerance for the integration

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


def I_gamV_original(z, n, b, Vo, Ko, kk, gamVo, gb, beta, e=0.001):
    """
    Integral of the Gruneisen parameter over the volume.

    Parameters
    ----------
    z : float
        Atomic number of the chemical formula unit
    n : float
        Number of atoms in the chemical formula
    b : float
        Fractional volume (V/V0)
    Vo : float
        Reference volume in [JBar^-1] (same as [cm^3/mol]/10)
    Ko : float
        Bulk modulus at reference volume [kbar]
    kk : float
        Bulk modulus derivative at reference volume
    gamVo : float
        Additive normalizing constant for the Gruneisen parameter, (given as delta in Sokolova et al. 2016)
    gb : float
        Generalized Gruneisen parameter, given as t in Sokolova et al. 2016
    beta : float
        Anharmonic analogue of the Grüneisen parameter
    e : float
        Tolerance for the integration

    Returns
    -------
    I_gamV : float
        Integral of the Gruneisen parameter over the volume

    """
    # check whether b is greater than 1
    a = min(b, 1)
    b = max(b, 1)

    # calculate the integral using the trapezoidal rule
    s2 = 1
    h = b - a
    S = f_gamV_original(a, n, z, Vo, Ko, kk, gamVo, gb, beta) + f_gamV_original(
        b, n, z, Vo, Ko, kk, gamVo, gb, beta
    )

    while True:
        s3 = s2
        h /= 2
        s1 = 0
        x = a + h

        while x < b:
            s1 += 2 * f_gamV_original(x, n, z, Vo, Ko, kk, gamVo, gb, beta)
            x += 2 * h

        S += s1
        s2 = (S + s1) * h / 3
        error = abs(s3 - s2) / 15

        if error < e:
            break

    return -s2 if b > 1 else s2


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


def f_gamV_original(x, n, z, Vo, Ko, kk, gamVo, gb, beta):
    """
    Helper function to calculate the Gruneisen parameter at a given temperature and pressure.

    Parameters
    ----------
    x : float
        Fractional volume (V/Vo)
    n : float
        Number of atoms in the chemical formula
    z : float
        Atomic number of the chemical formula unit
    Vo : float
        Reference volume in [JBar^-1] (same as [cm^3/mol]/10)
    Ko : float
        Bulk modulus at reference volume [kbar]
    kk : float
        Bulk modulus derivative at reference volume
    gamVo : float
        Additive normalizing constant for the Gruneisen parameter, (given as delta in Sokolova et al. 2016)
    gb : float
        Generalized Gruneisen parameter, given as t in Sokolova et al. 2016
    beta : float
        Anharmonic analogue of the Grüneisen parameter

    Returns
    -------
    f_gamV_value : float
        Gruneisen parameter at a given temperature and pressure
    """
    fw = -np.log(3 * Ko / 10 / (1003.6 * (z * n / (Vo * 10)) ** (5 / 3)))
    ff = x ** (1 / 3)
    aa = 1.5 * (kk - 3) - fw

    KT = (
        Ko
        * 1000
        / ff**6
        * np.exp(fw * (1 - ff))
        * (
            (-5 / ff**2 + 4 / ff) * (1 + aa * ff - aa * ff**2)
            + (1 / ff - 1) * (1 + aa * ff - aa * ff**2) * (-fw)
            + (1 / ff - 1) * (aa - 2 * aa * ff)
        )
        * (-x)
    )

    Px = (
        3
        * Ko
        * 1000
        * np.exp(fw * (1 - ff))
        * (1 / ff**5 - 1 / ff**4)
        * (1 + aa * ff - aa * ff**2)
    )

    ex = (
        3
        / ff**4
        * Ko
        * 1000
        * np.exp(fw * (1 - ff))
        * (
            (-5 / ff**2 + 4 / ff) * (1 + aa * ff - aa * ff**2)
            - (1 / ff - 1) * (1 + aa * ff - aa * ff**2) * fw
            + (1 / ff - 1) * (aa - 2 * aa * ff)
        )
    )

    ex1 = (
        1
        / ff**3
        * Ko
        * 1000
        * np.exp(fw * (1 - ff))
        * fw
        * (
            (-5 / ff**2 + 4 / ff) * (1 + aa * ff - aa * ff**2)
            - (1 / ff - 1) * (1 + aa * ff - aa * ff**2) * fw
            + (1 / ff - 1) * (aa - 2 * aa * ff)
        )
    )

    ex2 = (
        Ko
        * 1000
        * np.exp(fw * (1 - ff))
        * (
            (10 / ff**3 - 4 / ff**2) * (1 + aa * ff - aa * ff**2)
            + (-5 / ff**2 + 4 / ff) * (aa - 2 * aa * ff)
            + fw / ff**2 * (1 + aa * ff - aa * ff**2)
            - (1 / ff - 1) * (aa - 2 * aa * ff) * fw
            - (aa - 2 * aa * ff) / ff**2
            - 2 * aa * (1 / ff - 1)
        )
        / ff**3
    )

    kkx = (ex + ex1 - ex2) / (-KT / ff) / 3

    gt = gb - beta * x ** (1 / 3)
    gtx = -beta / 3 * x ** (-2 / 3)

    f_gamV_value = (
        gamVo
        + (-3 * KT + 2 * Px * gt + 9 * KT * kkx - 6 * gt * KT)
        / (6 * (3 * KT - 2 * Px * gt))
    ) / x

    return f_gamV_value
