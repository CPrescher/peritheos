# API reference

## Isothermal equations of state

```python
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
```

Common methods:

- `pressure(V)`
- `bulk_modulus(V)`
- `volume(P)` and `calculate_volume(P)`

## Thermal equations of state

```python
from peritheos.eos.thermal import (
    HollandPowell2011,
    MieGruneisenDebye,
    MieGruneisenEinstein,
    Sokolova2016,
    ThermalModifiedTait,
)
```

Common methods:

- `thermal_pressure(V, T)`
- `pressure(V, T)`
- `volume(P, T)`
- `bulk_modulus(V, T)`
- `isothermal_compressibility(V, T)`
- `thermal_expansivity(V, T)`
- `molar_heat_capacity_v(V, T)` when a caloric model exists
- `molar_heat_capacity_p(V, T)` when a caloric model exists
- `adiabatic_bulk_modulus(V, T)` when a caloric model exists
- `gruneisen_parameter(V, T)` when a caloric model exists

Mie-Gruneisen models additionally expose `gruneisen_parameter()`,
`characteristic_temperature()`, and the vibrational thermodynamic methods
documented under [Thermoelastic properties](thermoelastic-properties.md).

## Fitting

```python
from peritheos.fitting import FitResult, fit_joint_eos, fit_rt_eos, fit_thermal_eos
```

`FitResult` contains the fitted `model`, parameter and uncertainty mappings,
covariance and correlation matrices, raw and weighted residuals, chi-square,
degrees of freedom, AIC, BIC, convergence status, and solver message. It also
reports `adjusted_volume`, `adjusted_temperature`, `volume_corrections`, and
`temperature_corrections` for errors-in-variables fits. See
[Fitting P-V and P-V-T data](fitting.md) for the pressure, volume, and
temperature uncertainty API.

`FitResult.summary()` returns a compact text report. `FitResult.to_dict()` and
`FitResult.to_json()` export a versioned, JSON-safe record containing the model
identity and reconstructable parameters, fit arrays, diagnostics, and solver
configuration. Passing a path to `to_json()` also writes the record to disk.

`fit_joint_eos()` estimates reference-isotherm and thermal parameters in one
regression. Reference parameters use dotted names such as `rt_eos.V0`; its
covariance includes reference/thermal cross-correlations and is directly
compatible with `FitResult.eos_uncertainty()`.

## Units

```python
from peritheos.units import (
    convert_density,
    convert_molar_volume,
    density_from_molar_volume,
    molar_volume_from_density,
)
```

## Uncertainty propagation

```python
from peritheos import EOSUncertainty, ParameterUncertainty, PredictionUncertainty
```

`EOSUncertainty` wraps a deterministic EOS and provides `pressure()`, `volume()`,
`bulk_modulus()`, and the generic `evaluate()` method. Each returns a
`PredictionUncertainty` containing the nominal value, standard error, confidence
limits, statistical assumptions, and an optional output covariance matrix.

`FitResult.eos_uncertainty()` constructs the wrapper directly from fitted
parameter covariance. See [Uncertainty in EOS calculations](uncertainty.md).
Partial parameter-error sets are supported. Parameters omitted from the error
mapping or covariance ordering are treated as exact, rather than as having an
unknown error that Peritheos will estimate.
