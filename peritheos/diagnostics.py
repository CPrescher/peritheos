"""Numerical diagnostic transforms for equation-of-state data."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from peritheos.eos import EosBase, validate_positive_scalar
from peritheos.errors import EosValidationError


def _optional_standard_error(
    value: Any | None, shape: tuple[int, ...], name: str
) -> NDArray[np.float64] | None:
    if value is None:
        return None
    try:
        result = np.broadcast_to(np.asarray(value, dtype=float), shape).copy()
    except ValueError as error:
        raise EosValidationError(
            f"{name} must broadcast to the volume and pressure shape"
        ) from error
    if not np.all(np.isfinite(result)) or np.any(result <= 0.0):
        raise EosValidationError(f"{name} must be finite and greater than zero")
    return result


@dataclass(frozen=True)
class BirchMurnaghanFiniteStrainDiagnostic:
    """Eulerian finite strain and normalized stress for P-V observations.

    Arrays are one-dimensional and contain only observations for which the
    normalized stress is defined. Points at exactly ``V0`` are recorded in
    ``omitted_indices`` because both pressure and the normalizing factor vanish
    there. Standard errors use first-order propagation and treat supplied
    pressure and volume errors as independent.
    """

    strain: NDArray[np.float64]
    normalized_stress: NDArray[np.float64]
    reference_volume: float
    omitted_indices: NDArray[np.int64]
    model: EosBase | None = None
    strain_standard_error: NDArray[np.float64] | None = None
    normalized_stress_standard_error: NDArray[np.float64] | None = None

    def model_curve(
        self, *, points: int = 200
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return the attached EOS as an F-f curve over the observed range."""
        if self.model is None:
            raise EosValidationError("A model is required to construct an F-f curve")
        if isinstance(points, bool) or not isinstance(points, (int, np.integer)):
            raise EosValidationError("points must be an integer greater than one")
        if points < 2:
            raise EosValidationError("points must be an integer greater than one")

        lower = min(0.0, float(np.min(self.strain)))
        upper = max(0.0, float(np.max(self.strain)))
        if lower == upper:
            padding = max(abs(lower) * 0.05, 1.0e-6)
            lower -= padding
            upper += padding
        strain = np.linspace(lower, upper, int(points))
        volume = self.reference_volume / (1.0 + 2.0 * strain) ** 1.5
        pressure = np.asarray(self.model.pressure(volume), dtype=float)
        factor = 3.0 * strain * (1.0 + 2.0 * strain) ** 2.5

        with np.errstate(divide="ignore", invalid="ignore"):
            stress = pressure / factor
        at_reference = strain == 0.0
        if np.any(at_reference):
            bulk_modulus = getattr(self.model, "K0", None)
            if bulk_modulus is None:
                stress[at_reference] = np.nan
            else:
                stress[at_reference] = float(bulk_modulus)
        return strain, stress


def birch_murnaghan_finite_strain_diagnostic(
    volume: Any,
    pressure: Any,
    *,
    reference_volume: float | None = None,
    model: EosBase | None = None,
    pressure_sigma: Any | None = None,
    volume_sigma: Any | None = None,
) -> BirchMurnaghanFiniteStrainDiagnostic:
    r"""Transform P-V observations into a Birch-Murnaghan F-f diagnostic.

    The Eulerian strain and normalized stress are

    .. math::

       f = \frac{1}{2}\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right],
       \qquad
       F = \frac{P}{3f(1+2f)^{5/2}}.

    For BM3, ``F = K0 + 3*K0*(K0_prime - 4)*f/2``. A horizontal trend is
    therefore the graphical signature of BM2, where ``K0_prime`` is fixed at
    four. ``reference_volume`` may be omitted when ``model`` exposes ``V0``.
    """
    if reference_volume is None:
        if model is None or not hasattr(model, "V0"):
            raise EosValidationError(
                "reference_volume is required when the model does not expose V0"
            )
        reference_volume = float(model.V0)
    reference_volume = validate_positive_scalar(reference_volume, "Reference volume")

    try:
        volumes, pressures = np.broadcast_arrays(
            np.asarray(volume, dtype=float), np.asarray(pressure, dtype=float)
        )
    except ValueError as error:
        raise EosValidationError(
            "Volume and pressure must have broadcast-compatible shapes"
        ) from error
    if volumes.size == 0:
        raise EosValidationError("Volume and pressure must not be empty")
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise EosValidationError("Volume must be finite and greater than zero")
    if not np.all(np.isfinite(pressures)):
        raise EosValidationError("Pressure must be finite")

    pressure_errors = _optional_standard_error(
        pressure_sigma, pressures.shape, "pressure_sigma"
    )
    volume_errors = _optional_standard_error(
        volume_sigma, volumes.shape, "volume_sigma"
    )

    flat_volume = volumes.ravel()
    flat_pressure = pressures.ravel()
    strain = 0.5 * ((reference_volume / flat_volume) ** (2.0 / 3.0) - 1.0)
    valid = strain != 0.0
    if not np.any(valid):
        raise EosValidationError(
            "Normalized stress is undefined because every volume equals V0"
        )

    factor = 3.0 * strain[valid] * (1.0 + 2.0 * strain[valid]) ** 2.5
    stress = flat_pressure[valid] / factor
    if not np.all(np.isfinite(stress)):
        raise EosValidationError("The F-f transform produced non-finite values")

    strain_error = None
    stress_error = None
    if volume_errors is not None:
        selected_volume_error = volume_errors.ravel()[valid]
        df_dv = -(1.0 + 2.0 * strain[valid]) / (3.0 * flat_volume[valid])
        strain_error = np.abs(df_dv) * selected_volume_error

    if pressure_errors is not None or volume_errors is not None:
        variance = np.zeros_like(stress)
        if pressure_errors is not None:
            d_stress_dp = 1.0 / factor
            variance += (d_stress_dp * pressure_errors.ravel()[valid]) ** 2
        if volume_errors is not None:
            df_dv = -(1.0 + 2.0 * strain[valid]) / (3.0 * flat_volume[valid])
            d_stress_df = -stress * (
                1.0 / strain[valid] + 5.0 / (1.0 + 2.0 * strain[valid])
            )
            variance += (d_stress_df * df_dv * volume_errors.ravel()[valid]) ** 2
        stress_error = np.sqrt(variance)

    return BirchMurnaghanFiniteStrainDiagnostic(
        strain=np.asarray(strain[valid], dtype=float),
        normalized_stress=np.asarray(stress, dtype=float),
        reference_volume=reference_volume,
        omitted_indices=np.flatnonzero(~valid),
        model=model,
        strain_standard_error=strain_error,
        normalized_stress_standard_error=stress_error,
    )


__all__ = [
    "BirchMurnaghanFiniteStrainDiagnostic",
    "birch_murnaghan_finite_strain_diagnostic",
]
