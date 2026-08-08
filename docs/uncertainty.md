# Uncertainty in EOS calculations

Peritheos keeps its EOS classes deterministic and attaches uncertainty through
`EOSUncertainty`. This avoids changing the return type of established methods
such as `pressure()` and makes every statistical assumption explicit.

## From individual literature errors

Published EOS tables often report parameter standard errors without their
correlations:

```python
from peritheos import EOSUncertainty
from peritheos.eos.rt import BM3

eos = BM3(V0=10.0, K0=120.0, K0_prime=4.3)
eos_uncertainty = EOSUncertainty(
    eos,
    parameter_errors={
        "V0": 0.01,
        "K0": 2.0,
        "K0_prime": 0.1,
    },
)

prediction = eos_uncertainty.pressure([9.5, 9.0, 8.5])
print(prediction.value)
print(prediction.standard_error)
print(prediction.lower, prediction.upper)
```

This constructs a diagonal covariance matrix and records the assumption
`parameter errors treated as mutually independent` in `prediction.assumptions`.
EOS parameters are frequently correlated, so this fallback should not be
presented as equivalent to propagation from a fitted covariance matrix.

## With correlation or covariance

Individual errors and a correlation matrix define a covariance matrix:

```python
eos_uncertainty = EOSUncertainty(
    eos,
    parameter_errors={"V0": 0.01, "K0": 2.0, "K0_prime": 0.1},
    correlation=correlation_matrix,
)
```

A covariance matrix already contains both scale and correlation, so individual
errors are unnecessary. Matrix ordering must be stated when the matrix covers
only some EOS parameters:

```python
eos_uncertainty = EOSUncertainty(
    eos,
    parameter_names=("V0", "K0", "K0_prime"),
    covariance=covariance_matrix,
)
```

Correlation alone is insufficient because its diagonal is one and contains no
parameter-error scale. Peritheos rejects a correlation matrix supplied without
`parameter_errors`.

## When only some parameters have errors

Uncertainties do not need to be supplied for every EOS parameter. Only the
parameters named in `parameter_errors` or `parameter_names` are treated as
uncertain; every omitted parameter is held exact during propagation:

```python
eos_uncertainty = EOSUncertainty(
    eos,
    parameter_errors={
        "K0": 2.0,
        "K0_prime": 0.1,
    },
)
```

In this example, uncertainty in `K0` and `K0_prime` is propagated while `V0`
is treated as exact. When covariance is supplied for a subset, its dimensions
and ordering must match that subset exactly:

```python
eos_uncertainty = EOSUncertainty(
    eos,
    parameter_names=("K0", "K0_prime"),
    covariance=covariance_2_by_2,
)
```

Omitting a parameter is a statistical assumption that its uncertainty is zero.
It does not mean that Peritheos estimates an unpublished or unknown error. If an
omitted parameter actually has uncertainty, the calculated prediction error
will generally be underestimated. Reports should distinguish parameters known
to be fixed from parameters whose uncertainties were simply unavailable.

The same rule applies to state variables. `volume_sigma`, `temperature_sigma`,
and `pressure_sigma` are independently optional; an omitted state uncertainty
is treated as zero.

## From an EOS fit

`FitResult.eos_uncertainty()` uses the fitted covariance in the exact order of
`free_parameters`:

```python
fit = fit_rt_eos(...)

eos_uncertainty = fit.eos_uncertainty()
prediction = eos_uncertainty.bulk_modulus(volumes)
```

Parameters fixed during fitting remain deterministic. The fitted covariance,
correlation, and standard errors are available through:

```python
eos_uncertainty.covariance
eos_uncertainty.correlation
eos_uncertainty.standard_errors
```

## Linear propagation

Linear propagation is the default. Peritheos evaluates numerical parameter
derivatives and applies the delta method:

```text
output_covariance = J @ parameter_covariance @ J.T
```

It is fast and normally appropriate near a well-constrained fit. A complete
covariance between calculated curve points can be requested explicitly:

```python
prediction = eos_uncertainty.pressure(
    volumes,
    full_covariance=True,
    confidence=0.95,
)
print(prediction.covariance)
```

The reported lower and upper values are pointwise normal intervals, not a
simultaneous confidence band for the entire curve.

## Monte Carlo propagation

For nonlinear inversion, broad uncertainties, or asymmetric intervals, sample
the parameter distribution:

```python
prediction = eos_uncertainty.volume(
    pressures,
    method="monte_carlo",
    sample_count=10_000,
    random_state=42,
)
```

Parameters are sampled from a multivariate normal distribution. Samples that
violate EOS constructor constraints or calculation domains are rejected and
resampled. `prediction.rejected_fraction` reports the fraction rejected.

## Errors in the requested state

Measurement errors in the state being evaluated can be combined with parameter
uncertainty:

```python
pressure = eos_uncertainty.pressure(
    volume,
    temperature,
    volume_sigma=volume_errors,
    temperature_sigma=temperature_errors,
)

volume = eos_uncertainty.volume(
    pressure,
    temperature,
    pressure_sigma=pressure_errors,
    temperature_sigma=temperature_errors,
)
```

State errors are treated as independent between variables and evaluation
points. The generic `evaluate()` method accepts `argument_sigmas` for other EOS
quantities:

```python
heat_capacity = eos_uncertainty.evaluate(
    "molar_heat_capacity_p",
    volume,
    temperature,
    argument_sigmas={0: volume_errors, 1: temperature_errors},
)
```

## Thermal and reference-EOS uncertainty

A thermal fit currently holds its reference EOS fixed. Its fit covariance
therefore covers only the fitted thermal parameters. A separately quantified
reference EOS can be included explicitly:

```python
reference_uncertainty = EOSUncertainty(
    reference_eos,
    parameter_errors={"V0": 0.001, "K0": 1.0, "K0_prime": 0.05},
)

combined = thermal_fit.eos_uncertainty(
    additional=reference_uncertainty,
    assume_blocks_independent=True,
)
```

This forms a block-diagonal covariance and records the independence assumption.
A joint P-V-T fit is required when cross-correlations between reference and
thermal parameters matter.

## Interpretation limits

Propagated parameter uncertainty describes uncertainty within the selected EOS
and its parameter distribution. It does not include model-form error, systematic
calibration error not represented in the inputs, or disagreement between EOS
families. Linear covariance also relies on a local Gaussian approximation;
inspect Monte Carlo intervals and parameter correlations when extrapolating.
