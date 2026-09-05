"""Bounded errors-in-variables fitting for EOS models."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, cast

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares
from scipy.sparse import csr_matrix, issparse, lil_matrix
from scipy.sparse.linalg import spsolve

from peritheos import _rust
from peritheos.eos import EosBase, ThermalEOS
from peritheos.errors import FitValidationError
from peritheos.hugoniot import HugoniotBase, LinearUsUpHugoniot

if TYPE_CHECKING:
    from peritheos.uncertainty import EOSUncertainty


def _json_safe(value: Any) -> Any:
    """Return nested built-in values that strict JSON can represent."""
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


@dataclass(frozen=True)
class FitResult:
    """Result of an EOS fit.

    ``adjusted_volume`` and ``adjusted_temperature`` are the fitted latent
    values of measurements whose uncertainties were supplied.  They equal the
    observations when the corresponding uncertainty was omitted.
    """

    model: EosBase
    parameters: dict[str, float]
    standard_errors: dict[str, float]
    covariance: NDArray[np.float64]
    correlation: NDArray[np.float64]
    free_parameters: tuple[str, ...]
    residuals: NDArray[np.float64]
    weighted_residuals: NDArray[np.float64]
    adjusted_volume: NDArray[np.float64]
    adjusted_temperature: NDArray[np.float64] | None
    volume_corrections: NDArray[np.float64]
    temperature_corrections: NDArray[np.float64] | None
    chi_square: float
    reduced_chi_square: float
    degrees_of_freedom: int
    aic: float
    bic: float
    success: bool
    message: str
    loss: str = "linear"
    f_scale: float = 1.0
    max_nfev: int | None = None
    status: int = 0
    nfev: int = 0
    njev: int | None = None
    cost: float = np.nan
    optimality: float = np.nan

    def eos_uncertainty(
        self,
        *,
        additional: EOSUncertainty | None = None,
        assume_blocks_independent: bool = False,
    ) -> EOSUncertainty:
        """Return an uncertainty-aware EOS using this fit's covariance."""
        from peritheos.uncertainty import EOSUncertainty

        return EOSUncertainty.from_fit(
            self,
            additional=additional,
            assume_blocks_independent=assume_blocks_independent,
        )

    def summary(self, *, precision: int = 6) -> str:
        """Return a compact, human-readable fit report."""
        if isinstance(precision, bool) or not isinstance(precision, (int, np.integer)):
            raise FitValidationError("precision must be a positive integer")
        if precision <= 0:
            raise FitValidationError("precision must be a positive integer")
        precision = int(precision)

        def formatted(value: float) -> str:
            return f"{float(value):.{precision}g}"

        free = set(self.free_parameters)
        parameter_order = self.free_parameters + tuple(
            name for name in self.parameters if name not in free
        )
        name_width = max(9, *(len(name) for name in parameter_order))
        lines = [
            f"FitResult ({type(self.model).__name__})",
            f"Success: {self.success} (status {self.status})",
            f"Message: {self.message}",
            "",
            "Parameters:",
            f"  {'Name':<{name_width}}  {'Value':>14}  {'Std. error':>14}  Status",
        ]
        for name in parameter_order:
            lines.append(
                f"  {name:<{name_width}}  "
                f"{formatted(self.parameters[name]):>14}  "
                f"{formatted(self.standard_errors[name]):>14}  "
                f"{'free' if name in free else 'fixed'}"
            )
        lines.extend(
            [
                "",
                "Diagnostics:",
                f"  chi-square: {formatted(self.chi_square)}",
                f"  reduced chi-square: {formatted(self.reduced_chi_square)}",
                f"  degrees of freedom: {self.degrees_of_freedom}",
                f"  AIC: {formatted(self.aic)}",
                f"  BIC: {formatted(self.bic)}",
                "",
                "Solver:",
                f"  loss: {self.loss}",
                f"  f_scale: {formatted(self.f_scale)}",
                f"  function evaluations: {self.nfev}",
                f"  cost: {formatted(self.cost)}",
                f"  optimality: {formatted(self.optimality)}",
            ]
        )
        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Return a versioned, JSON-safe representation of the fit."""
        result = {
            "schema_version": 1,
            "model": {
                "module": type(self.model).__module__,
                "class": type(self.model).__name__,
                "parameters": self.model.parameter_values(include_reference=True),
                "configuration": self.model.configuration_values(),
            },
            "parameters": self.parameters,
            "standard_errors": self.standard_errors,
            "free_parameters": self.free_parameters,
            "covariance": self.covariance,
            "correlation": self.correlation,
            "residuals": self.residuals,
            "weighted_residuals": self.weighted_residuals,
            "adjusted_volume": self.adjusted_volume,
            "adjusted_temperature": self.adjusted_temperature,
            "volume_corrections": self.volume_corrections,
            "temperature_corrections": self.temperature_corrections,
            "diagnostics": {
                "chi_square": self.chi_square,
                "reduced_chi_square": self.reduced_chi_square,
                "degrees_of_freedom": self.degrees_of_freedom,
                "aic": self.aic,
                "bic": self.bic,
            },
            "solver": {
                "success": self.success,
                "status": self.status,
                "message": self.message,
                "loss": self.loss,
                "f_scale": self.f_scale,
                "max_nfev": self.max_nfev,
                "nfev": self.nfev,
                "njev": self.njev,
                "cost": self.cost,
                "optimality": self.optimality,
            },
        }
        return _json_safe(result)

    def to_json(self, path: str | Path | None = None, *, indent: int | None = 2) -> str:
        """Return strict JSON and optionally write it to *path*."""
        serialized = json.dumps(
            self.to_dict(), indent=indent, sort_keys=True, allow_nan=False
        )
        if path is not None:
            Path(path).write_text(serialized + "\n", encoding="utf-8")
        return serialized


@dataclass(frozen=True)
class HugoniotFitResult:
    """Result of fitting a linear ``Us = c0 + s * up`` relation."""

    model: LinearUsUpHugoniot
    parameters: dict[str, float]
    standard_errors: dict[str, float]
    covariance: NDArray[np.float64]
    correlation: NDArray[np.float64]
    free_parameters: tuple[str, ...]
    residuals: NDArray[np.float64]
    weighted_residuals: NDArray[np.float64]
    adjusted_particle_velocity: NDArray[np.float64]
    particle_velocity_corrections: NDArray[np.float64]
    chi_square: float
    reduced_chi_square: float
    degrees_of_freedom: int
    aic: float
    bic: float
    success: bool
    message: str
    loss: str = "linear"
    f_scale: float = 1.0
    max_nfev: int | None = None
    status: int = 0
    nfev: int = 0
    njev: int | None = None
    cost: float = np.nan
    optimality: float = np.nan

    def eos_uncertainty(self) -> EOSUncertainty:
        """Return the fitted Hugoniot with its parameter covariance."""
        from peritheos.uncertainty import EOSUncertainty

        return EOSUncertainty.from_fit(self)

    def summary(self, *, precision: int = 6) -> str:
        """Return a compact human-readable fit report."""
        if isinstance(precision, bool) or not isinstance(precision, (int, np.integer)):
            raise FitValidationError("precision must be a positive integer")
        if precision <= 0:
            raise FitValidationError("precision must be a positive integer")

        def formatted(value: float) -> str:
            return f"{float(value):.{int(precision)}g}"

        lines = [
            "HugoniotFitResult (LinearUsUpHugoniot)",
            f"Success: {self.success} (status {self.status})",
            f"Message: {self.message}",
            "",
            "Parameters:",
        ]
        free = set(self.free_parameters)
        for name, value in self.parameters.items():
            lines.append(
                f"  {name}: {formatted(value)} +/- "
                f"{formatted(self.standard_errors[name])} "
                f"({'free' if name in free else 'fixed'})"
            )
        lines.extend(
            [
                "",
                "Diagnostics:",
                f"  chi-square: {formatted(self.chi_square)}",
                f"  reduced chi-square: {formatted(self.reduced_chi_square)}",
                f"  degrees of freedom: {self.degrees_of_freedom}",
            ]
        )
        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Return a versioned, JSON-safe representation of the Hugoniot fit."""
        return _json_safe(
            {
                "schema_version": 1,
                "model": {
                    "module": type(self.model).__module__,
                    "class": type(self.model).__name__,
                    "parameters": self.model.parameter_values(),
                    "configuration": {},
                },
                "parameters": self.parameters,
                "standard_errors": self.standard_errors,
                "free_parameters": self.free_parameters,
                "covariance": self.covariance,
                "correlation": self.correlation,
                "residuals": self.residuals,
                "weighted_residuals": self.weighted_residuals,
                "adjusted_particle_velocity": self.adjusted_particle_velocity,
                "particle_velocity_corrections": self.particle_velocity_corrections,
                "diagnostics": {
                    "chi_square": self.chi_square,
                    "reduced_chi_square": self.reduced_chi_square,
                    "degrees_of_freedom": self.degrees_of_freedom,
                    "aic": self.aic,
                    "bic": self.bic,
                },
                "solver": {
                    "success": self.success,
                    "status": self.status,
                    "message": self.message,
                    "loss": self.loss,
                    "f_scale": self.f_scale,
                    "max_nfev": self.max_nfev,
                    "nfev": self.nfev,
                    "njev": self.njev,
                    "cost": self.cost,
                    "optimality": self.optimality,
                },
            }
        )

    def to_json(self, path: str | Path | None = None, *, indent: int | None = 2) -> str:
        """Return strict JSON and optionally write it to *path*."""
        serialized = json.dumps(
            self.to_dict(), indent=indent, sort_keys=True, allow_nan=False
        )
        if path is not None:
            Path(path).write_text(serialized + "\n", encoding="utf-8")
        return serialized


def _validated_observations(pressure: Any) -> NDArray[np.float64]:
    observed = np.asarray(pressure, dtype=float)
    if observed.size == 0 or not np.all(np.isfinite(observed)):
        raise FitValidationError("Observed pressure must contain finite values")
    return observed


def _validated_uncertainty(
    uncertainty: Any | None,
    shape: tuple[int, ...],
    name: str,
    *,
    default: float | None = None,
) -> NDArray[np.float64] | None:
    if uncertainty is None:
        if default is None:
            return None
        return np.full(shape, default, dtype=float)
    try:
        values = np.broadcast_to(np.asarray(uncertainty, dtype=float), shape).copy()
    except ValueError as error:
        raise FitValidationError(
            f"{name} must broadcast to the pressure shape"
        ) from error
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise FitValidationError(f"{name} must be finite and greater than zero")
    return values


def _pressure_uncertainty(
    pressure_sigma: Any | None,
    sigma: Any | None,
    shape: tuple[int, ...],
) -> tuple[NDArray[np.float64], bool]:
    if pressure_sigma is not None and sigma is not None:
        raise FitValidationError("Use either pressure_sigma or sigma, not both")
    supplied = pressure_sigma is not None or sigma is not None
    values = _validated_uncertainty(
        pressure_sigma if pressure_sigma is not None else sigma,
        shape,
        "pressure_sigma",
        default=1.0,
    )
    assert values is not None
    return values, supplied


def _validated_observation_covariance(
    covariance: Any,
    shape: tuple[int, ...],
    component_names: Sequence[str],
) -> NDArray[np.float64]:
    """Return Cholesky factors for per-observation covariance matrices."""
    component_count = len(component_names)
    target_shape = shape + (component_count, component_count)
    try:
        matrices = np.broadcast_to(
            np.asarray(covariance, dtype=float), target_shape
        ).copy()
    except ValueError as error:
        raise FitValidationError(
            "observation_covariance must broadcast to the pressure shape plus "
            f"({component_count}, {component_count}) for "
            f"{tuple(component_names)}"
        ) from error
    if not np.all(np.isfinite(matrices)):
        raise FitValidationError("observation_covariance must contain finite values")
    if not np.allclose(
        matrices, np.swapaxes(matrices, -1, -2), rtol=1.0e-12, atol=1.0e-15
    ):
        raise FitValidationError("observation_covariance must be symmetric")
    try:
        factors = np.linalg.cholesky(matrices)
    except np.linalg.LinAlgError as error:
        raise FitValidationError(
            "observation_covariance must be positive definite"
        ) from error
    return factors.reshape(-1, component_count, component_count)


def _validated_solver_options(
    loss: str | Callable[..., Any],
    f_scale: float,
    max_nfev: int | None,
) -> tuple[str | Callable[..., Any], float, int | None]:
    allowed_losses = {"linear", "soft_l1", "huber", "cauchy", "arctan"}
    if not callable(loss) and loss not in allowed_losses:
        raise FitValidationError(
            f"loss must be one of {sorted(allowed_losses)} or callable"
        )
    f_scale = float(f_scale)
    if not np.isfinite(f_scale) or f_scale <= 0.0:
        raise FitValidationError("f_scale must be finite and greater than zero")
    if max_nfev is not None:
        if isinstance(max_nfev, bool) or not isinstance(max_nfev, (int, np.integer)):
            raise FitValidationError("max_nfev must be a positive integer or None")
        if max_nfev <= 0:
            raise FitValidationError("max_nfev must be a positive integer or None")
        max_nfev = int(max_nfev)
    return loss, f_scale, max_nfev


def _diagonal_uncertainties(
    cholesky: NDArray[np.float64], shape: tuple[int, ...]
) -> tuple[NDArray[np.float64], ...] | None:
    """Return standard deviations when a covariance factor is diagonal."""
    off_diagonal = cholesky.copy()
    indices = np.arange(cholesky.shape[-1])
    off_diagonal[:, indices, indices] = 0.0
    if not np.allclose(off_diagonal, 0.0, rtol=0.0, atol=1.0e-15):
        return None
    return tuple(cholesky[:, index, index].reshape(shape) for index in indices)


def _observation_error_model(
    shape: tuple[int, ...],
    component_names: Sequence[str],
    component_sigmas: Sequence[Any | None],
    sigma_alias: Any | None,
    observation_covariance: Any | None,
) -> tuple[tuple[NDArray[np.float64] | None, ...], NDArray[np.float64] | None, bool]:
    """Resolve independent sigmas or correlated observation covariance."""
    if len(component_names) != len(component_sigmas):
        raise AssertionError("Each observation component requires a sigma entry")
    if component_names[0] != "pressure":
        raise AssertionError("Pressure must be the first observation component")

    if observation_covariance is not None:
        if sigma_alias is not None or any(
            value is not None for value in component_sigmas
        ):
            raise FitValidationError(
                "observation_covariance cannot be combined with individual "
                "sigma arguments"
            )
        cholesky = _validated_observation_covariance(
            observation_covariance, shape, component_names
        )
        diagonal = _diagonal_uncertainties(cholesky, shape)
        if diagonal is not None:
            return diagonal, None, True
        return tuple(np.ones(shape) for _ in component_names), cholesky, True

    pressure_uncertainty, pressure_supplied = _pressure_uncertainty(
        component_sigmas[0], sigma_alias, shape
    )
    coordinate_uncertainties = tuple(
        _validated_uncertainty(raw_sigma, shape, f"{name}_sigma")
        for name, raw_sigma in zip(component_names[1:], component_sigmas[1:])
    )
    supplied = pressure_supplied or any(
        uncertainty is not None for uncertainty in coordinate_uncertainties
    )
    return (pressure_uncertainty, *coordinate_uncertainties), None, supplied


def _jacobian_sparsity(
    point_count: int,
    parameter_count: int,
    coordinate_slices: Mapping[str, slice],
    *,
    coupled_residuals: bool = False,
):
    """Describe the block-local dependence of latent observation residuals."""
    adjusted_count = len(coordinate_slices)
    if adjusted_count == 0:
        return None
    matrix = lil_matrix(
        (
            point_count * (1 + adjusted_count),
            parameter_count + point_count * adjusted_count,
        ),
        dtype=int,
    )
    rows = np.arange(point_count)
    if coupled_residuals:
        for block in range(1 + adjusted_count):
            block_rows = point_count * block + rows
            matrix[block_rows, :parameter_count] = 1
            for coordinate_slice in coordinate_slices.values():
                matrix[block_rows, coordinate_slice.start + rows] = 1
        return matrix.tocsr()

    matrix[:point_count, :parameter_count] = 1
    for block, coordinate_slice in enumerate(coordinate_slices.values()):
        columns = coordinate_slice.start + rows
        matrix[rows, columns] = 1
        correction_rows = point_count * (block + 1) + rows
        matrix[correction_rows, columns] = 1
    return matrix.tocsr()


def _parameter_covariance(jacobian, parameter_count: int) -> NDArray[np.float64]:
    """Return the parameter block after profiling out latent observations."""
    if not issparse(jacobian):
        _, singular_values, right_vectors = np.linalg.svd(
            np.asarray(jacobian), full_matrices=False
        )
        tolerance = np.finfo(float).eps * max(jacobian.shape) * singular_values[0]
        retained = singular_values > tolerance
        scaled_vectors = (
            right_vectors[retained, :parameter_count].T / singular_values[retained]
        )
        return scaled_vectors @ scaled_vectors.T

    parameter_jacobian = jacobian[:, :parameter_count].toarray()
    latent_jacobian = jacobian[:, parameter_count:].tocsc()
    parameter_information = np.einsum(
        "ij,ik->jk", parameter_jacobian, parameter_jacobian
    )
    latent_information = latent_jacobian.T @ latent_jacobian
    cross_information = np.asarray(latent_jacobian.T @ parameter_jacobian)
    solved = np.asarray(spsolve(latent_information, cross_information)).reshape(
        latent_information.shape[0], parameter_count
    )
    profiled_information = parameter_information - np.einsum(
        "ij,ik->jk", cross_information, solved
    )
    return np.linalg.pinv(profiled_information, hermitian=True)


@cache
def _native_fitting_types() -> tuple[type[EosBase], ...]:
    """Load exact built-in classes lazily to avoid package import cycles."""
    from peritheos.eos.rt import (
        BM2,
        BM3,
        BM4,
        Holzapfel,
        ModifiedTait,
        Murnaghan,
        NaturalStrain2,
        NaturalStrain3,
        NaturalStrain4,
        Vinet,
    )
    from peritheos.eos.thermal import (
        LinearThermalPressure,
        LogVolumeThermalPressure,
        MieGruneisenDebye,
        MieGruneisenEinstein,
        SecondOrderTaylorThermalPressure,
        Sokolova2016,
        Tange2009Debye,
        ThermalModifiedTait,
        ThermalReferenceStateEOS,
    )

    return (
        BM2,
        BM3,
        BM4,
        Murnaghan,
        ModifiedTait,
        NaturalStrain2,
        NaturalStrain3,
        NaturalStrain4,
        Vinet,
        Holzapfel,
        LinearThermalPressure,
        LogVolumeThermalPressure,
        MieGruneisenDebye,
        MieGruneisenEinstein,
        Tange2009Debye,
        SecondOrderTaylorThermalPressure,
        ThermalModifiedTait,
        ThermalReferenceStateEOS,
        Sokolova2016,
    )


def _native_fitting_model(model: EosBase):
    """Return the native model only for an exact built-in Peritheos class.

    User subclasses may inherit ``_native`` while overriding pressure
    evaluation, so merely checking for that attribute would silently bypass
    their Python behavior during fitting.
    """
    native = getattr(model, "_native", None)
    if native is None or type(model) not in _native_fitting_types():
        return None
    return native


def _fit_model(
    factory: Callable[[Mapping[str, float]], EosBase],
    evaluator: Callable[
        [EosBase, Mapping[str, NDArray[np.float64]]], NDArray[np.float64]
    ],
    observed: NDArray[np.float64],
    pressure_sigma: NDArray[np.float64],
    coordinates: Mapping[str, NDArray[np.float64]],
    coordinate_sigmas: Mapping[str, NDArray[np.float64] | None],
    initial: Mapping[str, float],
    fixed: Mapping[str, float] | None,
    bounds: Mapping[str, Sequence[float]] | None,
    scale_covariance: bool,
    observation_cholesky: NDArray[np.float64] | None,
    loss: str | Callable[..., Any],
    f_scale: float,
    max_nfev: int | None,
    coordinate_lower_bounds: Mapping[str, float] | None = None,
) -> FitResult:
    loss, f_scale, max_nfev = _validated_solver_options(loss, f_scale, max_nfev)
    fixed_values = {name: float(value) for name, value in (fixed or {}).items()}
    overlap = set(initial) & set(fixed_values)
    if overlap:
        raise FitValidationError(
            f"Parameters cannot be both initial and fixed: {sorted(overlap)}"
        )
    names = tuple(initial)
    if not names:
        raise FitValidationError("At least one free parameter is required")
    parameter_x0 = np.array([float(initial[name]) for name in names], dtype=float)
    if not np.all(np.isfinite(parameter_x0)):
        raise FitValidationError("Initial parameters must be finite")

    configured_bounds = bounds or {}
    parameter_lower = []
    parameter_upper = []
    for name in names:
        interval = configured_bounds.get(name, (-np.inf, np.inf))
        if len(interval) != 2 or interval[0] >= interval[1]:
            raise FitValidationError(f"Invalid bounds for {name}")
        parameter_lower.append(float(interval[0]))
        parameter_upper.append(float(interval[1]))

    adjusted_names = tuple(
        name
        for name, uncertainty in coordinate_sigmas.items()
        if uncertainty is not None
    )
    coordinate_slices: dict[str, slice] = {}
    x0_parts = [parameter_x0]
    lower_parts = [np.asarray(parameter_lower)]
    upper_parts = [np.asarray(parameter_upper)]
    offset = len(names)
    for name in adjusted_names:
        values = coordinates[name]
        size = values.size
        coordinate_slices[name] = slice(offset, offset + size)
        x0_parts.append(values.ravel())
        coordinate_lower = (coordinate_lower_bounds or {}).get(
            name, np.finfo(float).tiny
        )
        lower_parts.append(np.full(size, coordinate_lower))
        upper_parts.append(np.full(size, np.inf))
        offset += size

    x0 = np.concatenate(x0_parts)
    lower = np.concatenate(lower_parts)
    upper = np.concatenate(upper_parts)
    jacobian_sparsity = _jacobian_sparsity(
        observed.size,
        len(names),
        coordinate_slices,
        coupled_residuals=observation_cholesky is not None,
    )

    def parameter_mapping(values):
        return {**fixed_values, **dict(zip(names, map(float, values[: len(names)])))}

    def adjusted_coordinates(values):
        adjusted = dict(coordinates)
        for name in adjusted_names:
            adjusted[name] = values[coordinate_slices[name]].reshape(observed.shape)
        return adjusted

    def residual_function(values):
        model = factory(parameter_mapping(values))
        adjusted = adjusted_coordinates(values)
        predicted = np.asarray(evaluator(model, adjusted), dtype=float)
        if predicted.shape != observed.shape:
            raise FitValidationError("Model predictions must match the pressure shape")
        pressure_residual = predicted - observed
        if observation_cholesky is not None:
            raw_components = [pressure_residual.ravel()]
            raw_components.extend(
                (adjusted[name] - coordinates[name]).ravel() for name in adjusted_names
            )
            raw_residuals = np.column_stack(raw_components)
            whitened = np.linalg.solve(
                observation_cholesky, raw_residuals[..., np.newaxis]
            )[..., 0]
            return whitened.T.ravel()

        residual_parts = [(pressure_residual / pressure_sigma).ravel()]
        for name in adjusted_names:
            uncertainty = coordinate_sigmas[name]
            assert uncertainty is not None
            residual_parts.append(
                ((adjusted[name] - coordinates[name]) / uncertainty).ravel()
            )
        return np.concatenate(residual_parts)

    prototype = factory(parameter_mapping(x0))
    native_model = _native_fitting_model(prototype)
    if isinstance(loss, str):
        if native_model is not None and "temperature" not in coordinates:
            optimization = _rust.fit_rt_eos_native(
                native_model,
                names,
                parameter_x0,
                np.asarray(parameter_lower),
                np.asarray(parameter_upper),
                observed.ravel(),
                coordinates["volume"].ravel(),
                pressure_sigma.ravel(),
                None
                if coordinate_sigmas["volume"] is None
                else coordinate_sigmas["volume"].ravel(),
                observation_cholesky,
                loss=loss,
                f_scale=f_scale,
                max_nfev=max_nfev,
            )
        elif native_model is not None:
            optimization = _rust.fit_thermal_eos_native(
                native_model,
                names,
                parameter_x0,
                np.asarray(parameter_lower),
                np.asarray(parameter_upper),
                observed.ravel(),
                coordinates["volume"].ravel(),
                coordinates["temperature"].ravel(),
                pressure_sigma.ravel(),
                None
                if coordinate_sigmas["volume"] is None
                else coordinate_sigmas["volume"].ravel(),
                None
                if coordinate_sigmas["temperature"] is None
                else coordinate_sigmas["temperature"].ravel(),
                observation_cholesky,
                loss=loss,
                f_scale=f_scale,
                max_nfev=max_nfev,
            )
        else:
            native_options = {}
            if adjusted_names:
                native_options = {
                    "global_parameter_count": len(names),
                    "point_count": observed.size,
                    "latent_coordinate_count": len(adjusted_names),
                }
            optimization = _rust.fit_least_squares(
                residual_function,
                x0,
                lower,
                upper,
                loss=loss,
                f_scale=f_scale,
                max_nfev=max_nfev,
                **native_options,
            )
    else:
        # Callable robust losses are an intentional compatibility fallback:
        # arbitrary Python callables cannot be represented by the native enum.
        optimization = least_squares(
            residual_function,
            x0,
            bounds=(lower, upper),
            jac_sparsity=jacobian_sparsity,
            x_scale="jac",
            loss=loss,
            f_scale=f_scale,
            max_nfev=max_nfev,
        )
    parameters = parameter_mapping(optimization.x)
    model = factory(parameters)
    adjusted = adjusted_coordinates(optimization.x)
    native_prediction = getattr(optimization, "predicted_pressure", None)
    predicted = (
        np.asarray(evaluator(model, adjusted), dtype=float)
        if native_prediction is None
        else np.asarray(native_prediction, dtype=float).reshape(observed.shape)
    )
    residuals = predicted - observed
    weighted_residuals = np.asarray(optimization.fun, dtype=float)
    count = int(weighted_residuals.size)
    degrees_of_freedom = count - optimization.x.size
    chi_square = float(np.sum(weighted_residuals**2))
    reduced_chi_square = (
        chi_square / degrees_of_freedom if degrees_of_freedom > 0 else np.nan
    )

    native_covariance = getattr(optimization, "parameter_covariance", None)
    if native_covariance is not None:
        covariance = np.asarray(native_covariance, dtype=float)
    else:
        covariance_jacobian = optimization.jac
        if adjusted_names and isinstance(loss, str):
            covariance_jacobian = csr_matrix(covariance_jacobian)
        covariance = _parameter_covariance(covariance_jacobian, len(names))
    if scale_covariance and degrees_of_freedom > 0:
        covariance *= reduced_chi_square
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(errors, errors)
    correlation = np.divide(
        covariance,
        denominator,
        out=np.zeros_like(covariance),
        where=denominator > 0.0,
    )
    standard_errors = {name: float(error) for name, error in zip(names, errors)}
    standard_errors.update({name: 0.0 for name in fixed_values})

    log_variance = math.log(max(chi_square / count, float(np.finfo(float).tiny)))
    fitted_count = optimization.x.size
    aic = count * log_variance + 2.0 * fitted_count
    bic = count * log_variance + fitted_count * math.log(count)
    adjusted_volume = adjusted["volume"]
    adjusted_temperature = adjusted.get("temperature")
    return FitResult(
        model=model,
        parameters=parameters,
        standard_errors=standard_errors,
        covariance=covariance,
        correlation=correlation,
        free_parameters=names,
        residuals=residuals,
        weighted_residuals=weighted_residuals,
        adjusted_volume=adjusted_volume,
        adjusted_temperature=adjusted_temperature,
        volume_corrections=adjusted_volume - coordinates["volume"],
        temperature_corrections=(
            None
            if adjusted_temperature is None
            else adjusted_temperature - coordinates["temperature"]
        ),
        chi_square=chi_square,
        reduced_chi_square=float(reduced_chi_square),
        degrees_of_freedom=degrees_of_freedom,
        aic=float(aic),
        bic=float(bic),
        success=bool(optimization.success),
        message=str(optimization.message),
        loss=(
            loss
            if isinstance(loss, str)
            else getattr(loss, "__qualname__", type(loss).__name__)
        ),
        f_scale=f_scale,
        max_nfev=max_nfev,
        status=int(optimization.status),
        nfev=int(optimization.nfev),
        njev=(None if optimization.njev is None else int(optimization.njev)),
        cost=float(optimization.cost),
        optimality=float(optimization.optimality),
    )


def fit_rt_eos(
    eos_class: type[EosBase],
    volume: Any,
    pressure: Any,
    initial: Mapping[str, float],
    *,
    fixed: Mapping[str, float] | None = None,
    bounds: Mapping[str, Sequence[float]] | None = None,
    pressure_sigma: Any | None = None,
    volume_sigma: Any | None = None,
    sigma: Any | None = None,
    absolute_sigma: bool = False,
    observation_covariance: Any | None = None,
    loss: str | Callable[..., Any] = "linear",
    f_scale: float = 1.0,
    max_nfev: int | None = None,
) -> FitResult:
    """Fit an isothermal EOS with optional errors in pressure and volume.

    ``sigma`` is retained as a compatibility alias for ``pressure_sigma``.
    When ``volume_sigma`` is supplied, the true volumes are fitted as latent
    values and their normalized corrections form part of the objective.
    ``observation_covariance`` accepts one or per-point 2-by-2 covariance
    matrices ordered as pressure, volume and cannot be combined with sigmas.
    """
    volumes, observed = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(pressure, dtype=float)
    )
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise FitValidationError("Volume must be finite and greater than zero")
    observed = _validated_observations(observed)
    uncertainties, observation_cholesky, uncertainty_supplied = (
        _observation_error_model(
            observed.shape,
            ("pressure", "volume"),
            (pressure_sigma, volume_sigma),
            sigma,
            observation_covariance,
        )
    )
    pressure_uncertainties, volume_uncertainties = uncertainties
    assert pressure_uncertainties is not None
    return _fit_model(
        lambda parameters: eos_class(**parameters),
        lambda model, coordinates: np.asarray(
            model.pressure(coordinates["volume"]), dtype=float
        ),
        observed,
        pressure_uncertainties,
        {"volume": volumes},
        {"volume": volume_uncertainties},
        initial,
        fixed,
        bounds,
        scale_covariance=not absolute_sigma or not uncertainty_supplied,
        observation_cholesky=observation_cholesky,
        loss=loss,
        f_scale=f_scale,
        max_nfev=max_nfev,
    )


def fit_linear_us_up(
    particle_velocity: Any,
    shock_velocity: Any,
    *,
    V0: float,
    rho0: float,
    P0: float = 0.0,
    initial: Mapping[str, float] | None = None,
    fixed: Mapping[str, float] | None = None,
    bounds: Mapping[str, Sequence[float]] | None = None,
    shock_velocity_sigma: Any | None = None,
    particle_velocity_sigma: Any | None = None,
    sigma: Any | None = None,
    absolute_sigma: bool = False,
    observation_covariance: Any | None = None,
    loss: str | Callable[..., Any] = "linear",
    f_scale: float = 1.0,
    max_nfev: int | None = None,
) -> HugoniotFitResult:
    """Fit ``Us = c0 + s * up`` by OLS, WLS, or errors in variables.

    With no uncertainty arguments, minimizing the vertical ``Us`` residuals is
    ordinary least squares. ``shock_velocity_sigma`` produces weighted least
    squares. Supplying ``particle_velocity_sigma`` additionally fits latent
    particle velocities, and ``observation_covariance`` accepts one or
    per-point 2-by-2 covariance matrices ordered as ``(Us, up)``.

    ``V0``, ``rho0``, and ``P0`` define the initial state used to turn the
    fitted velocity relation into an executable pressure-volume Hugoniot; they
    are not inferred from velocity data.
    """
    particle, shock = np.broadcast_arrays(
        np.asarray(particle_velocity, dtype=float),
        np.asarray(shock_velocity, dtype=float),
    )
    if particle.size < 2 or not np.all(np.isfinite(particle)):
        raise FitValidationError(
            "Particle velocity must contain at least two finite values"
        )
    if np.any(particle < 0.0):
        raise FitValidationError("Particle velocity must be non-negative")
    if not np.all(np.isfinite(shock)) or np.any(shock <= 0.0):
        raise FitValidationError("Shock velocity must be finite and positive")
    if np.ptp(particle) == 0.0:
        raise FitValidationError("Particle velocities must not all be equal")

    allowed = {"c0", "s"}
    supplied_initial = dict(initial or {})
    supplied_fixed = dict(fixed or {})
    configured_bounds = dict(bounds or {})
    unknown = (
        set(supplied_initial) | set(supplied_fixed) | set(configured_bounds)
    ) - allowed
    if unknown:
        raise FitValidationError(
            f"Only c0 and s can be fit from Us-up data; found {sorted(unknown)}"
        )
    overlap = set(supplied_initial) & set(supplied_fixed)
    if overlap:
        raise FitValidationError(
            f"Parameters cannot be both initial and fixed: {sorted(overlap)}"
        )

    slope, intercept = np.polyfit(particle.ravel(), shock.ravel(), 1)
    guesses = {"c0": float(intercept), "s": float(slope)}
    free_initial = {
        name: float(supplied_initial.get(name, guesses[name]))
        for name in ("c0", "s")
        if name not in supplied_fixed
    }
    if not free_initial:
        raise FitValidationError("At least c0 or s must be a free parameter")
    if any(value <= 0.0 or not np.isfinite(value) for value in free_initial.values()):
        raise FitValidationError(
            "Initial c0 and s estimates must be finite and positive; provide initial"
        )

    uncertainties, observation_cholesky, uncertainty_supplied = (
        _observation_error_model(
            shock.shape,
            ("pressure", "volume"),
            (shock_velocity_sigma, particle_velocity_sigma),
            sigma,
            observation_covariance,
        )
    )
    shock_uncertainty, particle_uncertainty = uncertainties
    assert shock_uncertainty is not None
    fixed_parameters = {
        "V0": float(V0),
        "rho0": float(rho0),
        "P0": float(P0),
        **{name: float(value) for name, value in supplied_fixed.items()},
    }
    positive_bounds: dict[str, Sequence[float]] = {
        "c0": (float(np.finfo(float).tiny), float(np.inf)),
        "s": (float(np.finfo(float).tiny), float(np.inf)),
        **configured_bounds,
    }

    def factory(parameters: Mapping[str, float]) -> LinearUsUpHugoniot:
        return LinearUsUpHugoniot(**parameters)

    def evaluator(
        model: HugoniotBase,
        coordinates: Mapping[str, NDArray[np.float64]],
    ) -> NDArray[np.float64]:
        return np.asarray(
            model.shock_velocity_from_particle_velocity(coordinates["volume"]),
            dtype=float,
        )

    inner = _fit_model(
        cast(Callable[[Mapping[str, float]], EosBase], factory),
        cast(
            Callable[
                [EosBase, Mapping[str, NDArray[np.float64]]],
                NDArray[np.float64],
            ],
            evaluator,
        ),
        shock,
        shock_uncertainty,
        {"volume": particle},
        {"volume": particle_uncertainty},
        free_initial,
        fixed_parameters,
        positive_bounds,
        scale_covariance=not absolute_sigma or not uncertainty_supplied,
        observation_cholesky=observation_cholesky,
        loss=loss,
        f_scale=f_scale,
        max_nfev=max_nfev,
        coordinate_lower_bounds={"volume": 0.0},
    )
    model = cast(LinearUsUpHugoniot, inner.model)
    return HugoniotFitResult(
        model=model,
        parameters=inner.parameters,
        standard_errors=inner.standard_errors,
        covariance=inner.covariance,
        correlation=inner.correlation,
        free_parameters=inner.free_parameters,
        residuals=inner.residuals,
        weighted_residuals=inner.weighted_residuals,
        adjusted_particle_velocity=inner.adjusted_volume,
        particle_velocity_corrections=inner.volume_corrections,
        chi_square=inner.chi_square,
        reduced_chi_square=inner.reduced_chi_square,
        degrees_of_freedom=inner.degrees_of_freedom,
        aic=inner.aic,
        bic=inner.bic,
        success=inner.success,
        message=inner.message,
        loss=inner.loss,
        f_scale=inner.f_scale,
        max_nfev=inner.max_nfev,
        status=inner.status,
        nfev=inner.nfev,
        njev=inner.njev,
        cost=inner.cost,
        optimality=inner.optimality,
    )


def fit_thermal_eos(
    eos_class: type[ThermalEOS],
    rt_eos: EosBase,
    volume: Any,
    temperature: Any,
    pressure: Any,
    initial: Mapping[str, float],
    *,
    configuration: Mapping[str, Any] | None = None,
    fixed: Mapping[str, float] | None = None,
    bounds: Mapping[str, Sequence[float]] | None = None,
    pressure_sigma: Any | None = None,
    volume_sigma: Any | None = None,
    temperature_sigma: Any | None = None,
    sigma: Any | None = None,
    absolute_sigma: bool = False,
    observation_covariance: Any | None = None,
    loss: str | Callable[..., Any] = "linear",
    f_scale: float = 1.0,
    max_nfev: int | None = None,
) -> FitResult:
    """Fit a thermal EOS with optional errors in pressure, volume, and temperature.

    The reference EOS remains fixed. ``configuration`` supplies fixed,
    non-numeric constructor choices such as ``debye_temperature_law``;
    configuration values are not fitted. Supplied volume and temperature errors
    turn their corresponding true values into latent fit variables.
    ``sigma`` is retained as a compatibility alias for ``pressure_sigma``.
    ``observation_covariance`` accepts one or per-point 3-by-3 covariance
    matrices ordered as pressure, volume, temperature.
    """
    configuration = dict(configuration or {})
    overlap = set(configuration) & (set(initial) | set(fixed or {}))
    if overlap:
        raise FitValidationError(
            "Configuration choices must not also be supplied as fit parameters: "
            f"{sorted(overlap)}"
        )
    volumes, temperatures, observed = np.broadcast_arrays(
        np.asarray(volume, dtype=float),
        np.asarray(temperature, dtype=float),
        np.asarray(pressure, dtype=float),
    )
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise FitValidationError("Volume must be finite and greater than zero")
    if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0.0):
        raise FitValidationError("Temperature must be finite and greater than zero")
    observed = _validated_observations(observed)
    uncertainties, observation_cholesky, uncertainty_supplied = (
        _observation_error_model(
            observed.shape,
            ("pressure", "volume", "temperature"),
            (pressure_sigma, volume_sigma, temperature_sigma),
            sigma,
            observation_covariance,
        )
    )
    (
        pressure_uncertainties,
        volume_uncertainties,
        temperature_uncertainties,
    ) = uncertainties
    assert pressure_uncertainties is not None
    return _fit_model(
        lambda parameters: eos_class(rt_eos=rt_eos, **parameters, **configuration),
        lambda model, coordinates: np.asarray(
            cast(ThermalEOS, model).pressure(
                coordinates["volume"], coordinates["temperature"]
            ),
            dtype=float,
        ),
        observed,
        pressure_uncertainties,
        {"volume": volumes, "temperature": temperatures},
        {
            "volume": volume_uncertainties,
            "temperature": temperature_uncertainties,
        },
        initial,
        fixed,
        bounds,
        scale_covariance=not absolute_sigma or not uncertainty_supplied,
        observation_cholesky=observation_cholesky,
        loss=loss,
        f_scale=f_scale,
        max_nfev=max_nfev,
    )


def fit_joint_eos(
    eos_class: type[ThermalEOS],
    rt_eos_class: type[EosBase],
    volume: Any,
    temperature: Any,
    pressure: Any,
    initial: Mapping[str, float],
    *,
    configuration: Mapping[str, Any] | None = None,
    fixed: Mapping[str, float] | None = None,
    bounds: Mapping[str, Sequence[float]] | None = None,
    pressure_sigma: Any | None = None,
    volume_sigma: Any | None = None,
    temperature_sigma: Any | None = None,
    sigma: Any | None = None,
    absolute_sigma: bool = False,
    observation_covariance: Any | None = None,
    loss: str | Callable[..., Any] = "linear",
    f_scale: float = 1.0,
    max_nfev: int | None = None,
) -> FitResult:
    """Jointly fit a reference isotherm and thermal EOS to P-V-T data.

    Reference-EOS parameters use the same dotted names exposed by thermal
    models, for example ``rt_eos.V0`` and ``rt_eos.K0``. Thermal parameters
    retain their constructor names. The returned covariance therefore includes
    cross-correlations between reference and thermal parameters and can be
    passed directly to :meth:`FitResult.eos_uncertainty`. ``configuration``
    supplies fixed, non-numeric thermal constructor choices and is not fitted.
    """

    configuration = dict(configuration or {})
    overlap = set(configuration) & (set(initial) | set(fixed or {}))
    if overlap:
        raise FitValidationError(
            "Configuration choices must not also be supplied as fit parameters: "
            f"{sorted(overlap)}"
        )

    def factory(parameters: Mapping[str, float]) -> ThermalEOS:
        reference_parameters = {
            name.removeprefix("rt_eos."): value
            for name, value in parameters.items()
            if name.startswith("rt_eos.")
        }
        thermal_parameters = {
            name: value
            for name, value in parameters.items()
            if not name.startswith("rt_eos.")
        }
        return eos_class(
            rt_eos=rt_eos_class(**reference_parameters),
            **thermal_parameters,
            **configuration,
        )

    volumes, temperatures, observed = np.broadcast_arrays(
        np.asarray(volume, dtype=float),
        np.asarray(temperature, dtype=float),
        np.asarray(pressure, dtype=float),
    )
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise FitValidationError("Volume must be finite and greater than zero")
    if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0.0):
        raise FitValidationError("Temperature must be finite and greater than zero")
    observed = _validated_observations(observed)

    all_parameter_names = set(initial) | set(fixed or {})
    if not any(name.startswith("rt_eos.") for name in all_parameter_names):
        raise FitValidationError("Joint fitting requires rt_eos.* reference parameters")

    uncertainties, observation_cholesky, uncertainty_supplied = (
        _observation_error_model(
            observed.shape,
            ("pressure", "volume", "temperature"),
            (pressure_sigma, volume_sigma, temperature_sigma),
            sigma,
            observation_covariance,
        )
    )
    (
        pressure_uncertainties,
        volume_uncertainties,
        temperature_uncertainties,
    ) = uncertainties
    assert pressure_uncertainties is not None

    return _fit_model(
        factory,
        lambda model, coordinates: np.asarray(
            cast(ThermalEOS, model).pressure(
                coordinates["volume"], coordinates["temperature"]
            ),
            dtype=float,
        ),
        observed,
        pressure_uncertainties,
        {"volume": volumes, "temperature": temperatures},
        {
            "volume": volume_uncertainties,
            "temperature": temperature_uncertainties,
        },
        initial,
        fixed,
        bounds,
        scale_covariance=not absolute_sigma or not uncertainty_supplied,
        observation_cholesky=observation_cholesky,
        loss=loss,
        f_scale=f_scale,
        max_nfev=max_nfev,
    )
