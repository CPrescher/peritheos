"""
Equations of state module for Peritheos
"""

from typing import Union
import numpy as np
from numpy.typing import NDArray
from scipy import optimize

# Type alias for numeric values (scalar or array)
NumericType = Union[float, NDArray[np.float64]]


class EosBase:
    """
    Base class for equation of state implementations.

    This abstract class defines the interface that all equation of state
    implementations should follow.
    """

    def pressure(self, V: NumericType) -> NumericType:
        """
        Calculate pressure at a given volume.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in cubic angstroms or any consistent unit)

        Returns
        -------
        float or numpy.ndarray
            Pressure (in the same units as K0)
        """
        raise NotImplementedError("Subclasses must implement the pressure method.")

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """
        Calculate the bulk modulus.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in cubic angstroms or any consistent unit)

        Returns
        -------
        float or numpy.ndarray
            Bulk modulus (in the same units as K0)
        """
        raise NotImplementedError("Subclasses must implement the bulk_modulus method.")

    def calculate_volume(self, P: NumericType) -> NumericType:
        """
        Calculate the volume at a given pressure, uses numerical optimization to find the volume..

        Parameters
        ----------
        P : float or numpy.ndarray
            Pressure (in the same units as K0)

        Returns
        -------
        float or numpy.ndarray
            Volume (in the same units as V0)
        """
        start_volume = self.V0 * 0.8
        # minimize the difference between the pressure and the pressure calculated from the EOS

        if isinstance(P, np.ndarray) or isinstance(P, list):
            return np.array([self.calculate_volume(P_i) for P_i in P])
        else:
            return optimize.minimize(
                lambda V: (self.pressure(V) - P) ** 2,
                start_volume,
                method="Nelder-Mead",
            ).x[0]


class ThermalEOS(EosBase):
    def __init__(self, rt_eos: EosBase):
        self.rt_eos = rt_eos

    def thermal_pressure(self, V: float, T: float) -> float:
        raise NotImplementedError("This method should be implemented by the subclass")

    def pressure(self, V: float, T: float) -> float:
        return self.thermal_pressure(V, T) + self.rt_eos.pressure(V)
