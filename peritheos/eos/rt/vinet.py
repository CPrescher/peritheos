"""
Vinet equation of state
"""

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class Vinet(EosBase):
    """
    Vinet equation of state. It has been presented multiple times in the literature. But has been
    attributed finally to Vinet 1986, 1987. An earlier description can be found in the work of
    Stacey et al. 1981. The Vinet EOS is based on a variant of the Morse potential, the Rydberg
    potential which was introduced by Rydberg in 1932.

    The Vinet EOS is defined as:
    P(V) = 3 * K0 * (1 - f)/ f^2 * exp(3/2 * (K0_prime - 1) * (1 - f))
    where f = (V / V0)^(1/3)

    The bulk modulus is given by:
    K(V) = K0 * f^(-2) * (1 + (eta * f + 1) * (1 - f) * exp(eta * (1 - f)))
    where f = (V / V0)^(1/3)
    and eta = 3/2 * (K0_prime - 1)
    """

    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
        """
        Initialize the Vinet equation of state.

        Parameters
        ----------
        V0 : float
            Reference volume (in cubic angstroms or any consistent unit)
        K0 : float
            Bulk modulus at reference volume (in GPa or any consistent unit)
        K0_prime : float
            Pressure derivative of bulk modulus at reference volume
        """
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")
        self._native = _rust.RtEos.vinet(self.V0, self.K0, self.K0_prime)

    def pressure(self, V: NumericType) -> NumericType:
        """
        Calculate pressure using the Vinet equation of state.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in cubic angstroms or any consistent unit)

        Returns
        -------
        float or numpy.ndarray
            Pressure (in the same units as K0)
        """
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "pressure", V)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """
        Calculate the bulk modulus using the Vinet equation of state.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in cubic angstroms or any consistent unit)

        Returns
        -------
        float or numpy.ndarray
            Bulk modulus (in the same units as K0)
        """
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "bulk_modulus", V)
