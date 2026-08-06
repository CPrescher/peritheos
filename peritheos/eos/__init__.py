"""
Equations of state module for Peritheos
"""

from typing import Callable, Union
import numpy as np
from numpy.typing import NDArray
from scipy import optimize

# Type alias for numeric values (scalar or array)
NumericType = Union[float, NDArray[np.float64]]


def validate_finite_scalar(value: float, name: str) -> float:
    """Return *value* as a finite float or raise a descriptive error."""
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def validate_positive_scalar(value: float, name: str) -> float:
    """Return *value* as a positive finite float."""
    value = validate_finite_scalar(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be greater than zero")
    return value


def validate_volume(V: NumericType) -> NumericType:
    """Validate volume input while preserving scalar or array behaviour."""
    values = np.asarray(V, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("Volume must be finite")
    if np.any(values <= 0):
        raise ValueError("Volume must be greater than zero")
    if values.ndim == 0:
        return float(values)
    return values


def _scalar_pressure(pressure_function: Callable[[float], NumericType], V: float) -> float:
    """Evaluate a pressure callable and require a finite scalar result."""
    pressure = np.asarray(pressure_function(V), dtype=float)
    if pressure.ndim != 0:
        raise TypeError("The pressure function must return a scalar for scalar volume")
    result = float(pressure)
    if not np.isfinite(result):
        raise ArithmeticError(f"EOS returned a non-finite pressure at V={V}")
    return result


def solve_volume(
    pressure_function: Callable[[float], NumericType],
    pressure: float,
    reference_volume: float,
) -> float:
    """Solve ``pressure_function(V) == pressure`` on the branch nearest V0.

    The EOS is assumed to be locally monotonic around its reference volume:
    pressure increases on compression and decreases on expansion. Starting at
    ``V0``, the search geometrically decreases volume for a higher target
    pressure or increases it for a lower target pressure until the pressure
    residual changes sign. ``scipy.optimize.brentq`` then solves inside that
    positive-volume bracket. If no sign change is found, the requested state is
    reported as outside the model's invertible range.
    """
    target = validate_finite_scalar(pressure, "Pressure")
    V0 = validate_positive_scalar(reference_volume, "Reference volume")
    p0 = _scalar_pressure(pressure_function, V0)
    f0 = p0 - target
    pressure_tolerance = 1e-10 * max(1.0, abs(target), abs(p0))
    if abs(f0) <= pressure_tolerance:
        return V0

    if f0 < 0:
        # The target is above P(V0), so search towards compression.
        upper, f_upper = V0, f0
        lower = V0
        f_lower = f0
        for _ in range(160):
            lower *= 0.8
            if lower <= V0 * 1e-14:
                break
            f_lower = _scalar_pressure(pressure_function, lower) - target
            if f_lower >= 0:
                break
        else:
            f_lower = np.nan
        if not np.isfinite(f_lower) or f_lower < 0:
            raise ValueError(
                f"Could not bracket a positive volume for pressure {target}"
            )
    else:
        # The target is below P(V0), so search along the first expansion branch.
        lower, f_lower = V0, f0
        upper = V0
        f_upper = f0
        for _ in range(160):
            upper *= 1.05
            if upper >= V0 * 1e4:
                break
            f_upper = _scalar_pressure(pressure_function, upper) - target
            if f_upper <= 0:
                break
        else:
            f_upper = np.nan
        if not np.isfinite(f_upper) or f_upper > 0:
            raise ValueError(
                f"Pressure {target} is outside the invertible expansion range"
            )

    result = optimize.brentq(
        lambda volume: _scalar_pressure(pressure_function, volume) - target,
        lower,
        upper,
        xtol=max(np.finfo(float).eps * V0, 1e-14),
        rtol=1e-12,
    )
    residual = abs(_scalar_pressure(pressure_function, result) - target)
    if residual > 1e-8 * max(1.0, abs(target)):
        raise ArithmeticError(
            f"Volume inversion did not converge to the requested pressure; residual={residual}"
        )
    return float(result)


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
        Calculate volume at a given pressure using a bracketed root solver.

        Parameters
        ----------
        P : float or numpy.ndarray
            Pressure (in the same units as K0)

        Returns
        -------
        float or numpy.ndarray
            Volume (in the same units as V0)
        """
        pressures = np.asarray(P, dtype=float)
        if not np.all(np.isfinite(pressures)):
            raise ValueError("Pressure must be finite")
        if pressures.ndim == 0:
            return solve_volume(self.pressure, float(pressures), self.V0)
        return np.array(
            [solve_volume(self.pressure, value, self.V0) for value in pressures.flat]
        ).reshape(pressures.shape)

    def volume(self, P: NumericType) -> NumericType:
        """Alias for :meth:`calculate_volume` using the P-to-V terminology."""
        return self.calculate_volume(P)


class ThermalEOS(EosBase):
    def __init__(self, rt_eos: EosBase):
        self.rt_eos = rt_eos

    def thermal_pressure(self, V: NumericType, T: NumericType) -> NumericType:
        raise NotImplementedError("This method should be implemented by the subclass")

    def pressure(self, V: NumericType, T: NumericType) -> NumericType:
        return self.thermal_pressure(V, T) + self.rt_eos.pressure(V)

    def calculate_volume(self, P: NumericType, T: NumericType) -> NumericType:
        """Calculate volume at pressure and temperature using bracketed roots."""
        pressures, temperatures = np.broadcast_arrays(
            np.asarray(P, dtype=float), np.asarray(T, dtype=float)
        )
        if not np.all(np.isfinite(pressures)):
            raise ValueError("Pressure must be finite")
        if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
            raise ValueError("Temperature must be finite and greater than zero")

        volumes = np.array(
            [
                solve_volume(
                    lambda volume, temperature=float(temperature): self.pressure(
                        volume, temperature
                    ),
                    float(pressure),
                    self.rt_eos.V0,
                )
                for pressure, temperature in zip(pressures.flat, temperatures.flat)
            ]
        ).reshape(pressures.shape)
        if volumes.ndim == 0:
            return float(volumes)
        return volumes

    def volume(self, P: NumericType, T: NumericType) -> NumericType:
        """Alias for :meth:`calculate_volume` for a thermal EOS."""
        return self.calculate_volume(P, T)
