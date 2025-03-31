import numpy as np
from scipy.integrate import quad
from scipy.constants import R

from ..rt.holzapfel import Holzapfel
from .. import ThermalEOS
from .. import NumericType


class Sokolova2016(ThermalEOS):
    def __init__(
        self,
        rt_eos: Holzapfel,
        Tr: float,
        QE1o: float,
        mE1: float,
        QE2o: float,
        mE2: float,
        delta: float,
        t: float,
        a_0: float,
        m: float,
        g: float,
        e_0: float,
    ):
        """
        Original thermal pressure equation from sokolova et al. 2016.

        Parameters
        ----------
        rt_eos : Holzapfel
            Room temperature equation of state (Holzapfel EOS)
        Tr : float
            Reference temperature in [K] for the EOS (typically 298.15 K)
        QE1o : float
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
        delta : float
            Additive normalizing constant for the Gruneisen parameter
        t : float
            Generalized Gruneisen parameter
        a_0 : float
            Intrinsic anharmonicity parameter (10e-6 [K])
        m : float
            Anharmonic analogue of the Grüneisen parameter
        g : float
            Electronic analogue of the Grüneisen parameter
        e_0 : float
            Free electrons parameter (10e-6 [K])
        """
        super().__init__(rt_eos)
        self.Tr = Tr
        self.QE1o = QE1o
        self.mE1 = mE1
        self.QE2o = QE2o
        self.mE2 = mE2
        self.delta = delta
        self.t = t
        self.a_0 = a_0
        self.m = m
        self.g = g
        self.e_0 = e_0

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """
        Calculate the thermal pressure using the Sokolova et al. 2016 model.

        Parameters
        ----------
        V : NumericType
            Volume in [cm^3/mol]
        T : NumericType
            Temperature in [K]

        Returns
        -------
        thermal_pressure : NumericType
            Thermal pressure in [bar]
        """
        if isinstance(V, np.ndarray) and isinstance(T, np.ndarray):
            if len(V) != len(T):
                raise ValueError(
                    "V and T either must have the same length or one must be a scalar"
                )

        x = V / self.rt_eos.V0  # fractional volume
        Px = self.rt_eos.pressure(V) * 10000  # convert GPa to bar
        KT = self.rt_eos.bulk_modulus(V) * 10000  # convert GPa to bar
        kkx = self.rt_eos.bulk_modulus_derivative(V)

        # Equation (10) - seems not the same as in the original paper
        gamV = (-3 * KT + 2 * Px * self.t + 9 * KT * kkx - 6 * self.t * KT) / 6 / (
            3 * KT - 2 * Px * self.t
        ) + self.delta

        # Exponent part in equation (9) - However this is not the correct value - the Excel spreadsheet
        # calculation has a more elaborate calculation included

        # expp_test = np.exp(0.5 * ao * TK * 1e6 * (V / V0) ** (1 / 3))
        expp = np.exp(I_gamV(x, self.delta, self.t, self.rt_eos))

        # Equation (9)
        QE1 = self.QE1o * expp
        QE2 = self.QE2o * expp

        # Equation (12) for the different Einstein contributions at the temperature TK
        e1 = np.exp(QE1 / T)
        PE1 = self.mE1 * R * ((QE1 / 2 + QE1 / (e1 - 1)) * gamV / V)

        e2 = np.exp(QE2 / T)
        PE2 = self.mE2 * R * ((QE2 / 2 + QE2 / (e2 - 1)) * gamV / V)

        # Equation (12) for the different Einstein contributions at the reference temperature
        e1r = np.exp(QE1 / self.Tr)
        PE1r = self.mE1 * R * ((QE1 / 2 + QE1 / (e1r - 1)) * gamV / V)

        e2r = np.exp(QE2 / self.Tr)
        PE2r = self.mE2 * R * ((QE2 / 2 + QE2 / (e2r - 1)) * gamV / V)

        # Equation (12) second additive term
        Pea = (
            3
            / 2
            * self.rt_eos.n
            * R
            * self.a_0
            / 1000000
            * x ** (self.m)
            * (self.m)
            / V
            * (T**2 - self.Tr**2)
        )
        Pee = (
            3
            / 2
            * self.rt_eos.n
            * R
            * self.e_0
            / 1000000
            * x ** (self.g)
            * (self.g)
            / V
            * (T**2 - self.Tr**2)
        )

        Pth = PE1 + PE2 - PE1r - PE2r + Pee + Pea
        return Pth


def I_gamV(x, delta, t, rt_eos):
    """
    Integral of the Gruneisen parameter over the volume ratio (from x to 1).

    Parameters
    ----------
    x : float
        Fractional volume (V/Vo)
    delta : float
        Additive normalizing constant for the Gruneisen parameter
    t : float
        Generalized Gruneisen parameter
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
        return f_gamV(x, Px_x, KT_x, kkx_x, delta, t)

    if isinstance(x, np.ndarray) or isinstance(x, list):
        return np.array([quad(f_gamV_x, x_i, 1)[0] for x_i in x])
    else:
        return quad(f_gamV_x, x, 1)[0]


def f_gamV(x, Px, KT, kkx, delta, t):
    """
    Helper function to calculate the Gruneisen parameter at a given temperature and pressure.

    Parameters
    ----------
    x : float
        Fractional volume (V/Vo)
    Px : float
        Pressure in [GPa]
    KT: float
        Bulk modulus at temperature in [GPa]
    kkx: float
        Bulk modulus derivative at temperature
    delta: float
        Additive normalizing constant for the Gruneisen parameter
    gb: float
        Generalized Gruneisen parameter
    """
    KT = KT * 10000  # convert GPa to bar
    Px = Px * 10000  # convert GPa to bar

    f_gamV_value = (
        delta
        + (-3 * KT + 2 * Px * t + 9 * KT * kkx - 6 * t * KT)
        / (6 * (3 * KT - 2 * Px * t))
    ) / x

    return f_gamV_value
