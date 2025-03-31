"""
Equations of state module for Peritheos
""" 

from typing import Union
import numpy as np
from numpy.typing import NDArray

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

class ThermalEOS(EosBase):
    def __init__(self, rt_eos: EosBase):
        self.rt_eos = rt_eos

    def thermal_pressure(self, V: float, T: float) -> float:
        raise NotImplementedError("This method should be implemented by the subclass")

    def pressure(self, V: float, T: float) -> float:
        return self.thermal_pressure(V, T) + self.rt_eos.pressure(V)