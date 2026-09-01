"""Uncertainty propagation for deterministic equation-of-state models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.stats import norm

from peritheos import _rust
from peritheos.eos import EosBase, NumericType, ThermalEOS


def _scalar_or_array(values: NDArray[np.float64]) -> NumericType:
    if values.ndim == 0:
        return float(values)
    return values


def _validate_covariance(
    covariance: Any, size: int, name: str = "covariance"
) -> NDArray[np.float64]:
    matrix = np.asarray(covariance, dtype=float)
    if matrix.shape != (size, size):
        raise ValueError(f"{name} must have shape ({size}, {size})")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must contain only finite values")
    scale = max(1.0, float(np.max(np.abs(matrix))))
    tolerance = 1.0e-12 * scale
    if not np.allclose(matrix, matrix.T, rtol=1.0e-10, atol=tolerance):
        raise ValueError(f"{name} must be symmetric")
    matrix = (matrix + matrix.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(matrix)
    if eigenvalues[0] < -tolerance:
        raise ValueError(f"{name} must be positive semidefinite")
    return matrix


def _correlation_from_covariance(
    covariance: NDArray[np.float64], errors: NDArray[np.float64]
) -> NDArray[np.float64]:
    denominator = np.outer(errors, errors)
    correlation = np.divide(
        covariance,
        denominator,
        out=np.zeros_like(covariance),
        where=denominator > 0.0,
    )
    positive = errors > 0.0
    correlation[np.diag_indices_from(correlation)] = positive.astype(float)
    return correlation


class ParameterUncertainty:
    """Validated parameter covariance and its named representation."""

    def __init__(
        self,
        *,
        parameter_errors: Mapping[str, float] | None = None,
        covariance: Any | None = None,
        correlation: Any | None = None,
        parameter_names: Sequence[str] | None = None,
        assumptions: Sequence[str] = (),
    ) -> None:
        if covariance is not None and (
            parameter_errors is not None or correlation is not None
        ):
            raise ValueError(
                "covariance cannot be combined with parameter_errors or correlation"
            )
        if correlation is not None and parameter_errors is None:
            raise ValueError("correlation requires parameter_errors")

        inferred_assumptions = list(assumptions)
        if covariance is not None:
            if parameter_names is None:
                raise ValueError("parameter_names are required with covariance")
            names = tuple(parameter_names)
            if not names:
                raise ValueError("At least one uncertain parameter is required")
            matrix = _validate_covariance(covariance, len(names))
        elif parameter_errors is not None:
            names = tuple(parameter_errors)
            if not names:
                raise ValueError("At least one uncertain parameter is required")
            if parameter_names is not None and tuple(parameter_names) != names:
                raise ValueError(
                    "parameter_names must match the order of parameter_errors"
                )
            errors = np.array([float(parameter_errors[name]) for name in names])
            if not np.all(np.isfinite(errors)) or np.any(errors <= 0.0):
                raise ValueError(
                    "parameter_errors must be finite and greater than zero"
                )
            if correlation is None:
                correlation_matrix = np.eye(len(names))
                inferred_assumptions.append(
                    "parameter errors treated as mutually independent"
                )
            else:
                correlation_matrix = _validate_covariance(
                    correlation, len(names), "correlation"
                )
                if not np.allclose(
                    np.diag(correlation_matrix), 1.0, rtol=0.0, atol=1.0e-10
                ):
                    raise ValueError("correlation must have ones on its diagonal")
                if np.any(np.abs(correlation_matrix) > 1.0 + 1.0e-10):
                    raise ValueError("correlation coefficients must lie in [-1, 1]")
            matrix = errors[:, None] * correlation_matrix * errors[None, :]
        else:
            raise ValueError("Provide covariance or parameter_errors")

        if len(set(names)) != len(names):
            raise ValueError("parameter names must be unique")
        errors = np.sqrt(np.maximum(np.diag(matrix), 0.0))
        self.parameter_names = names
        self.covariance = matrix
        self.standard_errors = {
            name: float(error) for name, error in zip(names, errors)
        }
        self.correlation = _correlation_from_covariance(matrix, errors)
        self.assumptions = tuple(dict.fromkeys(inferred_assumptions))

    @classmethod
    def _from_matrix(
        cls,
        parameter_names: Sequence[str],
        covariance: Any,
        assumptions: Sequence[str],
    ) -> ParameterUncertainty:
        return cls(
            parameter_names=parameter_names,
            covariance=covariance,
            assumptions=assumptions,
        )


@dataclass(frozen=True)
class PredictionUncertainty:
    """Nominal EOS result and its propagated uncertainty."""

    value: NumericType
    standard_error: NumericType
    lower: NumericType
    upper: NumericType
    covariance: NDArray[np.float64] | None
    method: str
    confidence: float
    assumptions: tuple[str, ...]
    rejected_fraction: float = 0.0


class EOSUncertainty:
    """Attach parameter uncertainty propagation to a deterministic EOS."""

    def __init__(
        self,
        eos: EosBase,
        *,
        parameter_errors: Mapping[str, float] | None = None,
        covariance: Any | None = None,
        correlation: Any | None = None,
        parameter_names: Sequence[str] | None = None,
    ) -> None:
        if not isinstance(eos, EosBase):
            raise TypeError("eos must be an equation of state")
        if covariance is not None and parameter_names is None:
            all_names = tuple(eos.parameter_values(include_reference=True))
            covariance_shape = np.asarray(covariance).shape
            if covariance_shape == (len(all_names), len(all_names)):
                parameter_names = all_names
            else:
                raise ValueError(
                    "parameter_names are required when covariance does not cover "
                    "all EOS parameters"
                )
        uncertainty = ParameterUncertainty(
            parameter_errors=parameter_errors,
            covariance=covariance,
            correlation=correlation,
            parameter_names=parameter_names,
        )
        self._initialize(eos, uncertainty)

    def _initialize(
        self, eos: EosBase, parameter_uncertainty: ParameterUncertainty
    ) -> None:
        available = eos.parameter_values(include_reference=True)
        unknown = set(parameter_uncertainty.parameter_names) - set(available)
        if unknown:
            raise ValueError(
                f"Unknown parameters for {type(eos).__name__}: {sorted(unknown)}"
            )
        self.eos = eos
        self.parameter_uncertainty = parameter_uncertainty

    @classmethod
    def _from_parameter_uncertainty(
        cls, eos: EosBase, parameter_uncertainty: ParameterUncertainty
    ) -> EOSUncertainty:
        result = cls.__new__(cls)
        result._initialize(eos, parameter_uncertainty)
        return result

    @classmethod
    def from_fit(
        cls,
        fit_result: Any,
        *,
        additional: EOSUncertainty | None = None,
        assume_blocks_independent: bool = False,
    ) -> EOSUncertainty:
        """Create uncertainty from a fit covariance and optional reference block."""
        uncertainty = ParameterUncertainty._from_matrix(
            fit_result.free_parameters,
            fit_result.covariance,
            ("parameter covariance obtained from EOS fit",),
        )
        if additional is None:
            return cls._from_parameter_uncertainty(fit_result.model, uncertainty)
        if not assume_blocks_independent:
            raise ValueError(
                "Set assume_blocks_independent=True to combine separate fits"
            )
        if not isinstance(fit_result.model, ThermalEOS):
            raise TypeError("additional uncertainty requires a thermal fitted EOS")
        reference_values = fit_result.model.rt_eos.parameter_values(
            include_reference=True
        )
        if type(additional.eos) is not type(fit_result.model.rt_eos) or any(
            not np.isclose(reference_values[name], value)
            for name, value in additional.eos.parameter_values(
                include_reference=True
            ).items()
        ):
            raise ValueError("additional uncertainty must describe the reference EOS")

        first_size = len(uncertainty.parameter_names)
        second_size = len(additional.parameter_uncertainty.parameter_names)
        covariance = np.zeros((first_size + second_size,) * 2)
        covariance[:first_size, :first_size] = uncertainty.covariance
        covariance[first_size:, first_size:] = (
            additional.parameter_uncertainty.covariance
        )
        names = uncertainty.parameter_names + tuple(
            f"rt_eos.{name}"
            for name in additional.parameter_uncertainty.parameter_names
        )
        assumptions = (
            uncertainty.assumptions
            + additional.parameter_uncertainty.assumptions
            + ("thermal and reference-EOS covariance blocks treated as independent",)
        )
        combined = ParameterUncertainty._from_matrix(names, covariance, assumptions)
        return cls._from_parameter_uncertainty(fit_result.model, combined)

    @property
    def parameter_names(self) -> tuple[str, ...]:
        return self.parameter_uncertainty.parameter_names

    @property
    def covariance(self) -> NDArray[np.float64]:
        return self.parameter_uncertainty.covariance

    @property
    def standard_errors(self) -> dict[str, float]:
        return dict(self.parameter_uncertainty.standard_errors)

    @property
    def correlation(self) -> NDArray[np.float64]:
        return self.parameter_uncertainty.correlation

    def _evaluate_model(
        self,
        eos: EosBase,
        quantity: str,
        arguments: tuple[Any, ...],
        quantity_kwargs: Mapping[str, Any],
    ) -> NDArray[np.float64]:
        if not hasattr(eos, quantity) or quantity.startswith("_"):
            raise ValueError(f"Unknown public EOS quantity: {quantity}")
        values = np.asarray(
            getattr(eos, quantity)(*arguments, **quantity_kwargs), dtype=float
        )
        if not np.all(np.isfinite(values)):
            raise ArithmeticError(f"{quantity} returned non-finite values")
        return values

    def _parameter_jacobian(
        self,
        quantity: str,
        arguments: tuple[Any, ...],
        quantity_kwargs: Mapping[str, Any],
        nominal: NDArray[np.float64],
        relative_step: float,
    ) -> NDArray[np.float64]:
        values = self.eos.parameter_values(include_reference=True)
        derivatives = []
        for name in self.parameter_names:
            value = values[name]
            scale = max(abs(value), self.standard_errors[name], 1.0)
            step = relative_step * scale
            plus = minus = None
            try:
                plus = self._evaluate_model(
                    self.eos.with_parameters(**{name: value + step}),
                    quantity,
                    arguments,
                    quantity_kwargs,
                )
            except (ArithmeticError, TypeError, ValueError):
                pass
            try:
                minus = self._evaluate_model(
                    self.eos.with_parameters(**{name: value - step}),
                    quantity,
                    arguments,
                    quantity_kwargs,
                )
            except (ArithmeticError, TypeError, ValueError):
                pass
            if plus is not None and minus is not None:
                derivative = (plus - minus) / (2.0 * step)
            elif plus is not None:
                derivative = (plus - nominal) / step
            elif minus is not None:
                derivative = (nominal - minus) / step
            else:
                raise ArithmeticError(
                    f"Could not perturb parameter {name!r} for {quantity}"
                )
            if derivative.shape != nominal.shape:
                raise ArithmeticError("Perturbed EOS output shape changed")
            derivatives.append(derivative.ravel())
        return np.column_stack(derivatives)

    def _state_variance(
        self,
        quantity: str,
        arguments: tuple[Any, ...],
        quantity_kwargs: Mapping[str, Any],
        argument_sigmas: Mapping[int, Any],
        nominal: NDArray[np.float64],
        relative_step: float,
    ) -> NDArray[np.float64]:
        variance = np.zeros(nominal.size)
        for index, raw_sigma in argument_sigmas.items():
            argument = np.asarray(arguments[index], dtype=float)
            try:
                sigma = np.broadcast_to(
                    np.asarray(raw_sigma, dtype=float), nominal.shape
                )
            except ValueError as error:
                raise ValueError(
                    "state uncertainty must broadcast to the calculated shape"
                ) from error
            if not np.all(np.isfinite(sigma)) or np.any(sigma <= 0.0):
                raise ValueError(
                    "state uncertainties must be finite and greater than zero"
                )
            step = relative_step * np.maximum(np.abs(argument), 1.0)
            plus_arguments = list(arguments)
            minus_arguments = list(arguments)
            plus_arguments[index] = argument + step
            minus_arguments[index] = argument - step
            plus = minus = None
            try:
                plus = self._evaluate_model(
                    self.eos, quantity, tuple(plus_arguments), quantity_kwargs
                )
            except (ArithmeticError, TypeError, ValueError):
                pass
            try:
                minus = self._evaluate_model(
                    self.eos, quantity, tuple(minus_arguments), quantity_kwargs
                )
            except (ArithmeticError, TypeError, ValueError):
                pass
            if plus is not None and minus is not None:
                derivative = (plus - minus) / (2.0 * step)
            elif plus is not None:
                derivative = (plus - nominal) / step
            elif minus is not None:
                derivative = (nominal - minus) / step
            else:
                raise ArithmeticError(
                    f"Could not perturb argument {index} for {quantity}"
                )
            variance += (
                np.broadcast_to(derivative, nominal.shape) * sigma
            ).ravel() ** 2
        return variance

    @staticmethod
    def _validate_options(
        confidence: float, relative_step: float
    ) -> tuple[float, float]:
        confidence = float(confidence)
        relative_step = float(relative_step)
        if not np.isfinite(confidence) or not 0.0 < confidence < 1.0:
            raise ValueError("confidence must lie between zero and one")
        if not np.isfinite(relative_step) or relative_step <= 0.0:
            raise ValueError("relative_step must be finite and greater than zero")
        return confidence, relative_step

    def _linear_evaluate(
        self,
        quantity: str,
        arguments: tuple[Any, ...],
        quantity_kwargs: Mapping[str, Any],
        argument_sigmas: Mapping[int, Any],
        confidence: float,
        full_covariance: bool,
        relative_step: float,
    ) -> PredictionUncertainty:
        nominal = self._evaluate_model(self.eos, quantity, arguments, quantity_kwargs)
        jacobian = self._parameter_jacobian(
            quantity, arguments, quantity_kwargs, nominal, relative_step
        )
        state_variance = self._state_variance(
            quantity,
            arguments,
            quantity_kwargs,
            argument_sigmas,
            nominal,
            relative_step,
        )
        propagation = _rust.linear_uncertainty(
            jacobian,
            self.covariance,
            state_variance,
            full_covariance=full_covariance,
        )
        variance = np.asarray(propagation.variance, dtype=float)
        standard_error = np.sqrt(variance).reshape(nominal.shape)
        output_covariance = propagation.covariance
        quantile = float(norm.ppf((1.0 + confidence) / 2.0))
        assumptions = self.parameter_uncertainty.assumptions + (
            "local linear (delta-method) uncertainty propagation",
        )
        if argument_sigmas:
            assumptions += ("state-variable errors treated as independent",)
        return PredictionUncertainty(
            value=_scalar_or_array(nominal),
            standard_error=_scalar_or_array(standard_error),
            lower=_scalar_or_array(nominal - quantile * standard_error),
            upper=_scalar_or_array(nominal + quantile * standard_error),
            covariance=output_covariance,
            method="linear",
            confidence=confidence,
            assumptions=tuple(dict.fromkeys(assumptions)),
        )

    def _monte_carlo_evaluate(
        self,
        quantity: str,
        arguments: tuple[Any, ...],
        quantity_kwargs: Mapping[str, Any],
        argument_sigmas: Mapping[int, Any],
        confidence: float,
        full_covariance: bool,
        sample_count: int,
        random_state: Any,
    ) -> PredictionUncertainty:
        if isinstance(sample_count, bool) or int(sample_count) != sample_count:
            raise ValueError("sample_count must be an integer")
        sample_count = int(sample_count)
        if sample_count < 2:
            raise ValueError("sample_count must be at least two")
        nominal = self._evaluate_model(self.eos, quantity, arguments, quantity_kwargs)
        rng = np.random.default_rng(random_state)
        parameter_values = self.eos.parameter_values(include_reference=True)
        means = np.array([parameter_values[name] for name in self.parameter_names])
        prepared_sigmas = {}
        for index, raw_sigma in argument_sigmas.items():
            argument = np.asarray(arguments[index], dtype=float)
            try:
                sigma = np.broadcast_to(
                    np.asarray(raw_sigma, dtype=float), argument.shape
                )
            except ValueError as error:
                raise ValueError(
                    "state uncertainty must broadcast to its state argument"
                ) from error
            if not np.all(np.isfinite(sigma)) or np.any(sigma <= 0.0):
                raise ValueError(
                    "state uncertainties must be finite and greater than zero"
                )
            prepared_sigmas[index] = sigma

        accepted = []
        attempted = 0
        maximum_attempts = 20 * sample_count
        while len(accepted) < sample_count and attempted < maximum_attempts:
            remaining = sample_count - len(accepted)
            batch_size = min(max(remaining * 2, 64), maximum_attempts - attempted)
            draws = rng.multivariate_normal(
                means, self.covariance, size=batch_size, check_valid="ignore"
            )
            for draw in draws:
                attempted += 1
                sampled_arguments = list(arguments)
                for index, sigma in prepared_sigmas.items():
                    sampled_arguments[index] = rng.normal(arguments[index], sigma)
                try:
                    sampled_eos = self.eos.with_parameters(
                        **dict(zip(self.parameter_names, map(float, draw)))
                    )
                    result = self._evaluate_model(
                        sampled_eos,
                        quantity,
                        tuple(sampled_arguments),
                        quantity_kwargs,
                    )
                    if result.shape != nominal.shape:
                        raise ArithmeticError("Sampled EOS output shape changed")
                except (ArithmeticError, FloatingPointError, TypeError, ValueError):
                    if attempted >= maximum_attempts:
                        break
                    continue
                accepted.append(result.ravel())
                if len(accepted) == sample_count:
                    break
        if len(accepted) < sample_count:
            raise ArithmeticError(
                "Could not obtain enough valid Monte Carlo EOS samples"
            )

        samples = np.asarray(accepted)
        summary = _rust.monte_carlo_summary(
            samples, confidence, full_covariance=full_covariance
        )
        standard_error = np.asarray(summary.standard_error).reshape(nominal.shape)
        lower = np.asarray(summary.lower).reshape(nominal.shape)
        upper = np.asarray(summary.upper).reshape(nominal.shape)
        output_covariance = summary.covariance
        assumptions = self.parameter_uncertainty.assumptions + (
            "parameter uncertainty sampled as a multivariate normal distribution",
        )
        if argument_sigmas:
            assumptions += (
                "state-variable errors sampled as independent normal distributions",
            )
        return PredictionUncertainty(
            value=_scalar_or_array(nominal),
            standard_error=_scalar_or_array(standard_error),
            lower=_scalar_or_array(lower),
            upper=_scalar_or_array(upper),
            covariance=output_covariance,
            method="monte_carlo",
            confidence=confidence,
            assumptions=tuple(dict.fromkeys(assumptions)),
            rejected_fraction=(attempted - sample_count) / attempted,
        )

    def evaluate(
        self,
        quantity: str,
        *arguments: Any,
        argument_sigmas: Mapping[int, Any] | None = None,
        method: str = "linear",
        confidence: float = 0.95,
        full_covariance: bool = False,
        relative_step: float = 1.0e-6,
        sample_count: int = 5000,
        random_state: Any = None,
        **quantity_kwargs: Any,
    ) -> PredictionUncertainty:
        """Evaluate an EOS method with propagated parameter uncertainty."""
        confidence, relative_step = self._validate_options(confidence, relative_step)
        state_sigmas = dict(argument_sigmas or {})
        invalid_indices = set(state_sigmas) - set(range(len(arguments)))
        if invalid_indices:
            raise ValueError(f"Invalid argument uncertainty indices: {invalid_indices}")
        if method == "linear":
            return self._linear_evaluate(
                quantity,
                tuple(arguments),
                quantity_kwargs,
                state_sigmas,
                confidence,
                bool(full_covariance),
                relative_step,
            )
        if method == "monte_carlo":
            return self._monte_carlo_evaluate(
                quantity,
                tuple(arguments),
                quantity_kwargs,
                state_sigmas,
                confidence,
                bool(full_covariance),
                sample_count,
                random_state,
            )
        raise ValueError("method must be 'linear' or 'monte_carlo'")

    def pressure(
        self,
        V: Any,
        T: Any | None = None,
        *,
        volume_sigma: Any | None = None,
        temperature_sigma: Any | None = None,
        **options: Any,
    ) -> PredictionUncertainty:
        """Return pressure and propagated uncertainty."""
        if isinstance(self.eos, ThermalEOS):
            if T is None:
                raise ValueError("Temperature is required for a thermal EOS")
            volumes, temperatures = np.broadcast_arrays(
                np.asarray(V, dtype=float), np.asarray(T, dtype=float)
            )
            sigmas = {
                index: sigma
                for index, sigma in enumerate((volume_sigma, temperature_sigma))
                if sigma is not None
            }
            return self.evaluate(
                "pressure", volumes, temperatures, argument_sigmas=sigmas, **options
            )
        if T is not None or temperature_sigma is not None:
            raise ValueError("Temperature is only valid for a thermal EOS")
        sigmas = {} if volume_sigma is None else {0: volume_sigma}
        return self.evaluate("pressure", V, argument_sigmas=sigmas, **options)

    def volume(
        self,
        P: Any,
        T: Any | None = None,
        *,
        pressure_sigma: Any | None = None,
        temperature_sigma: Any | None = None,
        **options: Any,
    ) -> PredictionUncertainty:
        """Return volume and propagated uncertainty."""
        if isinstance(self.eos, ThermalEOS):
            if T is None:
                raise ValueError("Temperature is required for a thermal EOS")
            pressures, temperatures = np.broadcast_arrays(
                np.asarray(P, dtype=float), np.asarray(T, dtype=float)
            )
            sigmas = {
                index: sigma
                for index, sigma in enumerate((pressure_sigma, temperature_sigma))
                if sigma is not None
            }
            return self.evaluate(
                "volume", pressures, temperatures, argument_sigmas=sigmas, **options
            )
        if T is not None or temperature_sigma is not None:
            raise ValueError("Temperature is only valid for a thermal EOS")
        sigmas = {} if pressure_sigma is None else {0: pressure_sigma}
        return self.evaluate("volume", P, argument_sigmas=sigmas, **options)

    def bulk_modulus(
        self,
        V: Any,
        T: Any | None = None,
        *,
        volume_sigma: Any | None = None,
        temperature_sigma: Any | None = None,
        **options: Any,
    ) -> PredictionUncertainty:
        """Return bulk modulus and propagated uncertainty."""
        if isinstance(self.eos, ThermalEOS):
            if T is None:
                raise ValueError("Temperature is required for a thermal EOS")
            volumes, temperatures = np.broadcast_arrays(
                np.asarray(V, dtype=float), np.asarray(T, dtype=float)
            )
            sigmas = {
                index: sigma
                for index, sigma in enumerate((volume_sigma, temperature_sigma))
                if sigma is not None
            }
            return self.evaluate(
                "bulk_modulus",
                volumes,
                temperatures,
                argument_sigmas=sigmas,
                **options,
            )
        if T is not None or temperature_sigma is not None:
            raise ValueError("Temperature is only valid for a thermal EOS")
        sigmas = {} if volume_sigma is None else {0: volume_sigma}
        return self.evaluate("bulk_modulus", V, argument_sigmas=sigmas, **options)


__all__ = ["EOSUncertainty", "ParameterUncertainty", "PredictionUncertainty"]
