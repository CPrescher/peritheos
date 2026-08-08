# Fitting P-V and P-V-T data

Peritheos uses bounded nonlinear least squares and reports covariance,
correlation, residual, and information-criterion diagnostics. Measurement
uncertainties can be supplied for every observed state variable.

## Isothermal fitting

```python
from peritheos.eos.rt import BM3
from peritheos.fitting import fit_rt_eos

result = fit_rt_eos(
    BM3,
    volume=volumes,
    pressure=pressures,
    initial={"V0": 10.0, "K0": 120.0, "K0_prime": 4.0},
    bounds={
        "V0": (9.0, 11.0),
        "K0": (50.0, 250.0),
        "K0_prime": (1.0, 10.0),
    },
    pressure_sigma=pressure_uncertainties,
    volume_sigma=volume_uncertainties,
    absolute_sigma=True,
)

print(result.parameters)
print(result.standard_errors)
print(result.correlation)
fitted_pressures = result.model.pressure(volumes)
```

Parameters in `fixed` are passed to the constructor but not optimized:

```python
result = fit_rt_eos(
    BM3,
    volumes,
    pressures,
    initial={"K0": 120.0, "K0_prime": 4.0},
    fixed={"V0": 10.0},
)
```

## Thermal fitting

The reference EOS is deliberately held fixed so the thermal and reference
isotherm regressions remain inspectable:

```python
from peritheos.eos.thermal import MieGruneisenEinstein
from peritheos.fitting import fit_thermal_eos

result = fit_thermal_eos(
    MieGruneisenEinstein,
    rt_eos=reference_eos,
    volume=volumes,
    temperature=temperatures,
    pressure=pressures,
    initial={"gamma0": 1.5, "q": 1.0},
    fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
    bounds={"gamma0": (0.0, 3.0), "q": (-2.0, 4.0)},
    pressure_sigma=pressure_uncertainties,
    volume_sigma=volume_uncertainties,
    temperature_sigma=temperature_uncertainties,
    absolute_sigma=True,
)
```

## Errors in pressure, volume, and temperature

Pressure is the dependent observation. When only `pressure_sigma` is given,
the objective contains the familiar normalized pressure residuals.

When `volume_sigma` or `temperature_sigma` is supplied, Peritheos performs an
errors-in-variables fit. For each measurement it also fits an adjusted latent
volume and/or temperature. For a P-V-T observation, the minimized components
are:

```text
(P_model(V*, T*) - P_observed) / pressure_sigma
(V* - V_observed) / volume_sigma
(T* - T_observed) / temperature_sigma
```

This is a nonlinear orthogonal-distance formulation for independent Gaussian
measurement errors. It does not merely evaluate the EOS at the measured
volume and temperature and pretend those coordinates are exact. Parameter
bounds continue to apply.

The fitted observation coordinates and corrections are available as:

```python
result.adjusted_volume
result.volume_corrections
result.adjusted_temperature
result.temperature_corrections
```

`weighted_residuals` contains the pressure components first, followed by the
volume and temperature correction components when present. `residuals`
contains the unscaled pressure residuals evaluated at the adjusted state.

## Uncertainty semantics

- With no uncertainty arguments, pressure residuals are unweighted and
  covariance is scaled by the fitted reduced chi-square.
- With uncertainty arguments and `absolute_sigma=True`, supplied standard
  deviations define the absolute covariance scale.
- With uncertainty arguments and `absolute_sigma=False`, they provide relative
  weights and covariance is rescaled by reduced chi-square.
- `sigma` remains a compatibility alias for `pressure_sigma`; do not pass both.
- Fixed parameters have a reported standard error of zero and are excluded
  from the covariance matrix.

All supplied uncertainties must be finite, positive, and scalar or
broadcastable to the pressure array. The current formulation assumes errors in
P, V, and T are mutually independent; full per-observation covariance matrices
are not yet supported.

Always inspect `correlation`, residuals versus pressure, and results under
reasonable changes in fitting range. A small residual does not make a strongly
correlated higher-order parameter physically meaningful.
