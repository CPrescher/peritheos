import numpy as np
from scipy.integrate import quad
from scipy.constants import R

from ..rt.holzapfel import Holzapfel
from .. import (
    NumericType,
    ThermalEOS,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


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
        if not isinstance(rt_eos, Holzapfel):
            raise TypeError("Sokolova2016 requires a Holzapfel room-temperature EOS")
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.QE1o = validate_positive_scalar(QE1o, "QE1o")
        self.mE1 = validate_finite_scalar(mE1, "mE1")
        self.QE2o = validate_positive_scalar(QE2o, "QE2o")
        self.mE2 = validate_finite_scalar(mE2, "mE2")
        if self.mE1 < 0 or self.mE2 < 0:
            raise ValueError("Einstein multiplicities must not be negative")
        self.delta = validate_finite_scalar(delta, "delta")
        self.t = validate_finite_scalar(t, "t")
        self.a_0 = validate_finite_scalar(a_0, "a_0")
        self.m = validate_finite_scalar(m, "m")
        self.g = validate_finite_scalar(g, "g")
        self.e_0 = validate_finite_scalar(e_0, "e_0")

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """
        Calculate the thermal pressure using the Sokolova et al. 2016 model.

        Parameters
        ----------
        V : NumericType
            Molar volume in [J bar^-1], equal to [cm^3/mol] / 10
        T : NumericType
            Temperature in [K]

        Returns
        -------
        thermal_pressure : NumericType
            Thermal pressure in [GPa]
        """
        V = validate_volume(V)
        temperatures = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise ValueError("Temperature must be finite and greater than zero")
        try:
            V, T = np.broadcast_arrays(np.asarray(V, dtype=float), temperatures)
        except ValueError as error:
            raise ValueError("V and T must have broadcast-compatible shapes") from error

        x = V / self.rt_eos.V0  # fractional volume
        Px = self.rt_eos.pressure(V)
        KT = self.rt_eos.bulk_modulus(V)
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
        PE1 = self.mE1 * R * (_einstein_energy(QE1, T) * gamV / V)

        PE2 = self.mE2 * R * (_einstein_energy(QE2, T) * gamV / V)

        # Equation (12) for the different Einstein contributions at the reference temperature
        PE1r = self.mE1 * R * (_einstein_energy(QE1, self.Tr) * gamV / V)

        PE2r = self.mE2 * R * (_einstein_energy(QE2, self.Tr) * gamV / V)

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

        # R [J mol^-1 K^-1] divided by V [J bar^-1 mol^-1] produces bar.
        Pth_bar = PE1 + PE2 - PE1r - PE2r + Pee + Pea
        Pth_gpa = Pth_bar / 10000
        if Pth_gpa.ndim == 0:
            return float(Pth_gpa)
        return Pth_gpa


def _einstein_energy(theta, temperature):
    """Return Einstein oscillator energy in kelvin without exponential overflow."""
    ratio = theta / temperature
    decay = np.exp(-ratio)
    thermal_part = theta * decay / (-np.expm1(-ratio))
    return theta / 2 + thermal_part


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

    x_values = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(x_values)) or np.any(x_values <= 0):
        raise ValueError("Volume ratio must be finite and greater than zero")
    integrals = np.array(
        [quad(f_gamV_x, float(x_i), 1)[0] for x_i in x_values.flat]
    ).reshape(x_values.shape)
    if integrals.ndim == 0:
        return float(integrals)
    return integrals


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
    f_gamV_value = (
        delta
        + (-3 * KT + 2 * Px * t + 9 * KT * kkx - 6 * t * KT)
        / (6 * (3 * KT - 2 * Px * t))
    ) / x

    return f_gamV_value
