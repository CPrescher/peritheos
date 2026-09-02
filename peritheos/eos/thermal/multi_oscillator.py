"""Composable multi-oscillator Gruneisen thermal-pressure equation."""

from __future__ import annotations

import numpy as np
from scipy.constants import R
from scipy.integrate import quad

from peritheos.errors import EosValidationError

from .. import (
    EosBase,
    NumericType,
    ThermalEOS,
    _native_for_exact_model,
    _native_thermal_evaluate,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class MultiOscillatorGruneisenThermalEOS(ThermalEOS):
    def __init__(
        self,
        rt_eos: EosBase,
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
        beta: float = 0.0,
        QBo: float = 1.0,
        d: float = 1.0,
        mb: float = 0.0,
        QB1o: float = 1.0,
        d1: float = 1.0,
        mb1: float = 0.0,
        n: float | None = None,
    ):
        """
        Multi-oscillator Gruneisen thermal-pressure correction.

        The reference isotherm is an independent component and may be any
        :class:`~peritheos.eos.EosBase` that provides pressure and bulk
        modulus. The generic EOS interface supplies ``dK/dP`` numerically when
        a reference EOS does not implement a specialized form.

        The optional ``beta`` and Bose-mode parameters reproduce the complete
        pressure expression in the original calculation model. Their defaults
        preserve the reduced form used by existing Peritheos releases.

        Parameters
        ----------
        rt_eos : EosBase
            Freely chosen reference-temperature equation of state
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
        delta : float
            Additive normalizing constant for the Gruneisen parameter
        t : float
            Generalized Gruneisen parameter
        a_0 : float
            Intrinsic anharmonicity parameter in 10^-6 K^-1
        m : float
            Anharmonic analogue of the Grüneisen parameter
        g : float
            Electronic analogue of the Grüneisen parameter
        e_0 : float
            Free-electron parameter in 10^-6 K^-1
        beta : float
            Volume-dependent correction to ``t``. The model uses
            ``t - beta * (V/V0)^(1/3)``.
        QBo, QB1o : float
            Reference characteristic temperatures of the two generalized
            Bose-Einstein modes in [K].
        d, d1 : float
            Positive dispersion parameters for the Bose-Einstein modes.
        mb, mb1 : float
            Non-negative multiplicities of the Bose-Einstein modes.
        n : float, optional
            Number of atoms per chemical formula. This should be supplied for
            new model compositions. If omitted, the value is read from a
            reference EOS that exposes ``n`` for compatibility with earlier
            Holzapfel-based usage.
        """
        super().__init__(rt_eos)
        self.Tr = validate_positive_scalar(Tr, "Tr")
        self.QE1o = validate_positive_scalar(QE1o, "QE1o")
        self.mE1 = validate_finite_scalar(mE1, "mE1")
        self.QE2o = validate_positive_scalar(QE2o, "QE2o")
        self.mE2 = validate_finite_scalar(mE2, "mE2")
        if self.mE1 < 0 or self.mE2 < 0:
            raise EosValidationError("Einstein multiplicities must not be negative")
        self.delta = validate_finite_scalar(delta, "delta")
        self.t = validate_finite_scalar(t, "t")
        self.a_0 = validate_finite_scalar(a_0, "a_0")
        self.m = validate_finite_scalar(m, "m")
        self.g = validate_finite_scalar(g, "g")
        self.e_0 = validate_finite_scalar(e_0, "e_0")
        self.beta = validate_finite_scalar(beta, "beta")
        self.QBo = validate_positive_scalar(QBo, "QBo")
        self.d = validate_positive_scalar(d, "d")
        self.mb = validate_finite_scalar(mb, "mb")
        self.QB1o = validate_positive_scalar(QB1o, "QB1o")
        self.d1 = validate_positive_scalar(d1, "d1")
        self.mb1 = validate_finite_scalar(mb1, "mb1")
        if self.mb < 0 or self.mb1 < 0:
            raise EosValidationError(
                "Bose-Einstein multiplicities must not be negative"
            )
        if n is None:
            n = getattr(rt_eos, "n", None)
        if n is None:
            raise EosValidationError(
                "n must be supplied when the reference isotherm does not expose it"
            )
        self.n = validate_positive_scalar(n, "n")
        reference_native = _native_for_exact_model(rt_eos)
        if (
            reference_native is not None
            and type(self) is MultiOscillatorGruneisenThermalEOS
        ):
            from peritheos import _rust

            self._native = _rust.ThermalEos.multi_oscillator_gruneisen(
                reference_native,
                self.Tr,
                self.QE1o,
                self.mE1,
                self.QE2o,
                self.mE2,
                self.delta,
                self.t,
                self.a_0,
                self.m,
                self.g,
                self.e_0,
                self.beta,
                self.QBo,
                self.d,
                self.mb,
                self.QB1o,
                self.d1,
                self.mb1,
                self.n,
            )

    def _volume_terms(self, V: NumericType) -> tuple[np.ndarray, ...]:
        """Precompute the volume-dependent terms in the pressure expression."""
        volume = np.asarray(validate_volume(V), dtype=float)
        x = volume / self.rt_eos.V0  # fractional volume
        Px = self.rt_eos.pressure(volume)
        KT = self.rt_eos.bulk_modulus(volume)
        kkx = self.rt_eos.bulk_modulus_derivative(volume)

        # Equation (10), following the original calculation model.
        generalized_t = self.t - self.beta * np.cbrt(x)
        gamV = (
            -3 * KT + 2 * Px * generalized_t + 9 * KT * kkx - 6 * generalized_t * KT
        ) / 6 / (3 * KT - 2 * Px * generalized_t) + self.delta

        # Exponent in equation (9), obtained by integrating gamma(V) / V.
        expp = np.exp(I_gamV(x, self.delta, self.t, self.rt_eos, self.beta))

        QB = self.QBo * expp
        QB1 = self.QB1o * expp
        QE1 = self.QE1o * expp
        QE2 = self.QE2o * expp
        reference_oscillator_pressure = (
            self.mb * R * (_bose_energy(QB, self.Tr, self.d) * gamV / volume)
            + self.mb1 * R * (_bose_energy(QB1, self.Tr, self.d1) * gamV / volume)
            + self.mE1 * R * (_einstein_energy(QE1, self.Tr) * gamV / volume)
            + self.mE2 * R * (_einstein_energy(QE2, self.Tr) * gamV / volume)
        )
        squared_temperature_coefficient = (
            3
            / 2
            * self.n
            * R
            / 1000000
            / volume
            * (self.a_0 * x**self.m * self.m + self.e_0 * x**self.g * self.g)
        )

        return (
            volume,
            np.asarray(gamV),
            QB,
            QB1,
            QE1,
            QE2,
            np.asarray(reference_oscillator_pressure),
            np.asarray(squared_temperature_coefficient),
        )

    def _thermal_pressure_from_volume_terms(
        self, volume_terms: tuple[np.ndarray, ...], T: NumericType
    ) -> NumericType:
        """Evaluate thermal pressure using prepared fixed-volume terms."""
        temperatures = np.asarray(T, dtype=float)
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise EosValidationError("Temperature must be finite and greater than zero")
        try:
            V, gamV, QB, QB1, QE1, QE2, reference_pressure, t2_coefficient, T = (
                np.broadcast_arrays(*volume_terms, temperatures)
            )
        except ValueError as error:
            raise EosValidationError(
                "V and T must have broadcast-compatible shapes"
            ) from error

        # Equation (12) for the different oscillator contributions at T and Tr.
        PB = self.mb * R * (_bose_energy(QB, T, self.d) * gamV / V)
        PB1 = self.mb1 * R * (_bose_energy(QB1, T, self.d1) * gamV / V)
        PE1 = self.mE1 * R * (_einstein_energy(QE1, T) * gamV / V)
        PE2 = self.mE2 * R * (_einstein_energy(QE2, T) * gamV / V)

        # R [J mol^-1 K^-1] divided by V [J bar^-1 mol^-1] produces bar.
        Pth_bar = (
            PB
            + PB1
            + PE1
            + PE2
            - reference_pressure
            + t2_coefficient * (T**2 - self.Tr**2)
        )
        return self._scalar_or_array(np.asarray(Pth_bar / 10000))

    def _thermal_pressure_function(self, V: float):
        """Prepare the costly volume integral once for temperature inversion."""
        if hasattr(self, "_native"):
            return lambda temperature: _native_thermal_evaluate(
                self._native, "thermal_pressure", V, temperature
            )
        volume_terms = self._volume_terms(V)
        return lambda temperature: self._thermal_pressure_from_volume_terms(
            volume_terms, temperature
        )

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        """Calculate the multi-oscillator thermal pressure.

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
        if hasattr(self, "_native"):
            volumes, temperatures = self._broadcast_state(V, T)
            return _native_thermal_evaluate(
                self._native, "thermal_pressure", volumes, temperatures
            )
        return self._thermal_pressure_from_volume_terms(self._volume_terms(V), T)


def _einstein_energy(theta, temperature):
    """Return Einstein oscillator energy in kelvin without exponential overflow."""
    ratio = theta / temperature
    decay = np.exp(-ratio)
    thermal_part = theta * decay / (-np.expm1(-ratio))
    return theta / 2 + thermal_part


def _bose_energy(theta, temperature, dispersion):
    """Return generalized Bose-mode energy in kelvin without overflow."""
    exponent = dispersion * np.log1p(theta / (temperature * dispersion))
    decay = np.exp(-exponent)
    occupation = decay / (-np.expm1(-exponent))
    zero_point = theta * (dispersion - 1) / (2 * dispersion)
    thermal_part = (
        temperature
        * theta
        * dispersion
        * occupation
        / (temperature * dispersion + theta)
    )
    return zero_point + thermal_part


def I_gamV(x, delta, t, rt_eos, beta=0.0):
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
    beta : float
        Volume-dependent correction to ``t``

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
        return f_gamV(x, Px_x, KT_x, kkx_x, delta, t, beta)

    x_values = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(x_values)) or np.any(x_values <= 0):
        raise EosValidationError("Volume ratio must be finite and greater than zero")
    integrals = np.array(
        [
            0.0
            if np.isclose(x_i, 1.0, rtol=0.0, atol=8 * np.finfo(float).eps)
            else quad(f_gamV_x, float(x_i), 1)[0]
            for x_i in x_values.flat
        ]
    ).reshape(x_values.shape)
    if integrals.ndim == 0:
        return float(integrals)
    return integrals


def f_gamV(x, Px, KT, kkx, delta, t, beta=0.0):
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
    t: float
        Generalized Gruneisen parameter
    beta: float
        Volume-dependent correction to ``t``
    """
    generalized_t = t - beta * np.cbrt(x)
    f_gamV_value = (
        delta
        + (-3 * KT + 2 * Px * generalized_t + 9 * KT * kkx - 6 * generalized_t * KT)
        / (6 * (3 * KT - 2 * Px * generalized_t))
    ) / x

    return f_gamV_value
