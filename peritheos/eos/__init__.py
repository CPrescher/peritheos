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