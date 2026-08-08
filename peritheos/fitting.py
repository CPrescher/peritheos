"""Bounded errors-in-variables fitting for EOS models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares
from scipy.sparse import issparse, lil_matrix
from scipy.sparse.linalg import spsolve

from peritheos.eos import EosBase, ThermalEOS


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


def _validated_observations(pressure: Any) -> NDArray[np.float64]:
    observed = np.asarray(pressure, dtype=float)
    if observed.size == 0 or not np.all(np.isfinite(observed)):
        raise ValueError("Observed pressure must contain finite values")
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
        raise ValueError(f"{name} must broadcast to the pressure shape") from error
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError(f"{name} must be finite and greater than zero")
    return values


def _pressure_uncertainty(
    pressure_sigma: Any | None,
    sigma: Any | None,
    shape: tuple[int, ...],
) -> tuple[NDArray[np.float64], bool]:
    if pressure_sigma is not None and sigma is not None:
        raise ValueError("Use either pressure_sigma or sigma, not both")
    supplied = pressure_sigma is not None or sigma is not None
    values = _validated_uncertainty(
        pressure_sigma if pressure_sigma is not None else sigma,
        shape,
        "pressure_sigma",
        default=1.0,
    )
    assert values is not None
    return values, supplied


def _jacobian_sparsity(
    point_count: int, parameter_count: int, coordinate_slices: Mapping[str, slice]
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
    matrix[:point_count, :parameter_count] = 1
    rows = np.arange(point_count)
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
        tolerance = (
            np.finfo(float).eps
            * max(jacobian.shape)
            * singular_values[0]
        )
        retained = singular_values > tolerance
        scaled_vectors = (
            right_vectors[retained, :parameter_count].T
            / singular_values[retained]
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
) -> FitResult:
    fixed_values = {name: float(value) for name, value in (fixed or {}).items()}
    overlap = set(initial) & set(fixed_values)
    if overlap:
        raise ValueError(
            f"Parameters cannot be both initial and fixed: {sorted(overlap)}"
        )
    names = tuple(initial)
    if not names:
        raise ValueError("At least one free parameter is required")
    parameter_x0 = np.array([float(initial[name]) for name in names], dtype=float)
    if not np.all(np.isfinite(parameter_x0)):
        raise ValueError("Initial parameters must be finite")

    configured_bounds = bounds or {}
    parameter_lower = []
    parameter_upper = []
    for name in names:
        interval = configured_bounds.get(name, (-np.inf, np.inf))
        if len(interval) != 2 or interval[0] >= interval[1]:
            raise ValueError(f"Invalid bounds for {name}")
        parameter_lower.append(float(interval[0]))
        parameter_upper.append(float(interval[1]))

    adjusted_names = tuple(
        name for name, uncertainty in coordinate_sigmas.items() if uncertainty is not None
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
        lower_parts.append(np.full(size, np.finfo(float).tiny))
        upper_parts.append(np.full(size, np.inf))
        offset += size

    x0 = np.concatenate(x0_parts)
    lower = np.concatenate(lower_parts)
    upper = np.concatenate(upper_parts)
    jacobian_sparsity = _jacobian_sparsity(
        observed.size, len(names), coordinate_slices
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
            raise ValueError("Model predictions must match the pressure shape")
        residual_parts = [((predicted - observed) / pressure_sigma).ravel()]
        for name in adjusted_names:
            uncertainty = coordinate_sigmas[name]
            assert uncertainty is not None
            residual_parts.append(
                ((adjusted[name] - coordinates[name]) / uncertainty).ravel()
            )
        return np.concatenate(residual_parts)

    optimization = least_squares(
        residual_function,
        x0,
        bounds=(lower, upper),
        jac_sparsity=jacobian_sparsity,
        x_scale="jac",
    )
    parameters = parameter_mapping(optimization.x)
    model = factory(parameters)
    adjusted = adjusted_coordinates(optimization.x)
    predicted = np.asarray(evaluator(model, adjusted), dtype=float)
    residuals = predicted - observed
    weighted_residuals = residual_function(optimization.x)
    count = weighted_residuals.size
    degrees_of_freedom = count - optimization.x.size
    chi_square = float(np.sum(weighted_residuals**2))
    reduced_chi_square = (
        chi_square / degrees_of_freedom if degrees_of_freedom > 0 else np.nan
    )

    covariance = _parameter_covariance(optimization.jac, len(names))
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

    log_variance = np.log(max(chi_square / count, np.finfo(float).tiny))
    fitted_count = optimization.x.size
    aic = count * log_variance + 2.0 * fitted_count
    bic = count * log_variance + fitted_count * np.log(count)
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
) -> FitResult:
    """Fit an isothermal EOS with optional errors in pressure and volume.

    ``sigma`` is retained as a compatibility alias for ``pressure_sigma``.
    When ``volume_sigma`` is supplied, the true volumes are fitted as latent
    values and their normalized corrections form part of the objective.
    """
    volumes, observed = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(pressure, dtype=float)
    )
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("Volume must be finite and greater than zero")
    observed = _validated_observations(observed)
    pressure_uncertainties, pressure_sigma_supplied = _pressure_uncertainty(
        pressure_sigma, sigma, observed.shape
    )
    volume_uncertainties = _validated_uncertainty(
        volume_sigma, observed.shape, "volume_sigma"
    )
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
        scale_covariance=(
            not absolute_sigma
            or not (pressure_sigma_supplied or volume_uncertainties is not None)
        ),
    )


def fit_thermal_eos(
    eos_class: type[ThermalEOS],
    rt_eos: EosBase,
    volume: Any,
    temperature: Any,
    pressure: Any,
    initial: Mapping[str, float],
    *,
    fixed: Mapping[str, float] | None = None,
    bounds: Mapping[str, Sequence[float]] | None = None,
    pressure_sigma: Any | None = None,
    volume_sigma: Any | None = None,
    temperature_sigma: Any | None = None,
    sigma: Any | None = None,
    absolute_sigma: bool = False,
) -> FitResult:
    """Fit a thermal EOS with optional errors in pressure, volume, and temperature.

    The reference EOS remains fixed.  Supplied volume and temperature errors
    turn their corresponding true values into latent fit variables.
    ``sigma`` is retained as a compatibility alias for ``pressure_sigma``.
    """
    volumes, temperatures, observed = np.broadcast_arrays(
        np.asarray(volume, dtype=float),
        np.asarray(temperature, dtype=float),
        np.asarray(pressure, dtype=float),
    )
    if not np.all(np.isfinite(volumes)) or np.any(volumes <= 0.0):
        raise ValueError("Volume must be finite and greater than zero")
    if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0.0):
        raise ValueError("Temperature must be finite and greater than zero")
    observed = _validated_observations(observed)
    pressure_uncertainties, pressure_sigma_supplied = _pressure_uncertainty(
        pressure_sigma, sigma, observed.shape
    )
    volume_uncertainties = _validated_uncertainty(
        volume_sigma, observed.shape, "volume_sigma"
    )
    temperature_uncertainties = _validated_uncertainty(
        temperature_sigma, observed.shape, "temperature_sigma"
    )
    return _fit_model(
        lambda parameters: eos_class(rt_eos=rt_eos, **parameters),
        lambda model, coordinates: np.asarray(
            model.pressure(coordinates["volume"], coordinates["temperature"]),
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
        scale_covariance=(
            not absolute_sigma
            or not (
                pressure_sigma_supplied
                or volume_uncertainties is not None
                or temperature_uncertainties is not None
            )
        ),
    )
