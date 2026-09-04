# Uncertainty in EOS calculations

`EOSUncertainty` accepts equilibrium EOSs and shock Hugoniots. For a Hugoniot,
use `evaluate()` with any public path quantity, for example
`uncertainty.evaluate("shock_velocity", volume)`; temperature is not an
independent state variable.

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

\[
\Sigma_y = J_\theta\,\Sigma_\theta\,J_\theta^{\mathsf T}.
\]

Here $y=g(\boldsymbol\theta,\boldsymbol x)$, $J_\theta$ is the Jacobian
with respect to uncertain EOS parameters, and $\Sigma_\theta$ is their
covariance. When independent state-variable errors are supplied, the complete
implemented approximation is

\[
\Sigma_y =
J_\theta\Sigma_\theta J_\theta^{\mathsf T}
+\operatorname{diag}\left[
\sum_k\left(
\frac{\partial g}{\partial x_k}\sigma_{x_k}
\right)^2
\right].
\]

Parameter covariance can therefore correlate different calculated curve
points. State errors contribute only to the diagonal because they are treated
as independent between variables and evaluation points.

For confidence level $c$, the pointwise normal interval is

\[
y\pm z_{(1+c)/2}\,\sigma_y.
\]

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
The reported confidence limits are empirical quantiles,

\[
\left[Q_{(1-c)/2},\;Q_{(1+c)/2}\right],
\]

so they can be asymmetric even though the sampled parameter distribution is
multivariate normal.

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

For `temperature_from_volumes`, both measured volumes contribute uncertainty,
and `f_dac` is an experimental assumption rather than an EOS parameter. The
reduced equation divides the reference-isotherm pressure difference by
`1 - f_dac`, so the result becomes especially sensitive as `f_dac` approaches
one. Propagate the two
volume errors with
`evaluate("temperature_from_volumes", V_ambient, V_heated, f_dac=f_dac,
argument_sigmas={0: ambient_sigma, 1: heated_sigma})` and report a sensitivity
sweep over the plausible `f_dac` interval; keyword-only `f_dac` uncertainty is
not automatically propagated as a model parameter.

The sensitivity sweep represents uncertainty in the experimental boundary
condition, not ordinary uncertainty in an EOS coefficient. Do not combine a
literature percentage with the volume errors until its denominator is known:
fractions of EOS thermal pressure and fractions of cold pressure are different
models. Report the `f_dac=0` isobaric result and each confinement scenario
separately.

Uncertainty in `V_ambient` should include both measurement error and uncertainty
from choosing the stable reference interval. If reference volumes before and
after heating differ systematically, a single pooled standard deviation is not
an adequate description of baseline drift. Analyze the branches separately or
use a documented time-dependent baseline. Near `f_dac=1`, or when
`V_heated - V_ambient` is comparable with its uncertainty, the inversion can be
strongly nonlinear; prefer Monte Carlo propagation and inspect rejected or
nonphysical samples.

## Thermal and reference-EOS uncertainty

A `fit_thermal_eos` fit holds its reference EOS fixed. Its fit covariance
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
When cross-correlations between reference and thermal parameters matter, use a
joint P-V-T fit instead:

```python
from peritheos.fitting import fit_joint_eos

joint_fit = fit_joint_eos(...)
combined = joint_fit.eos_uncertainty()
```

The dotted `rt_eos.*` parameter names in the joint covariance are understood by
the thermal model's reconstruction API, so both linear and Monte Carlo
propagation preserve the fitted cross-correlations.

## Interpretation limits

Propagated parameter uncertainty describes uncertainty within the selected EOS
and its parameter distribution. It does not include model-form error, systematic
calibration error not represented in the inputs, or disagreement between EOS
families. Linear covariance also relies on a local Gaussian approximation;
inspect Monte Carlo intervals and parameter correlations when extrapolating.
