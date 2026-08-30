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

Constructor signatures and special requirements are:

| Class | Signature after class name | Special requirement |
|---|---|---|
| `BM2` | `(V0, K0)` | $K_0'=4$ is implied |
| `BM3` | `(V0, K0, K0_prime)` | none |
| `BM4` | `(V0, K0, K0_prime, K0_double_prime)` | `K0_double_prime` has inverse-pressure units |
| `Murnaghan` | `(V0, K0, K0_prime)` | supports the `K0_prime=0` limit |
| `NaturalStrain2` | `(V0, K0)` | $K_0'=2$ is implied |
| `NaturalStrain3` | `(V0, K0, K0_prime)` | none |
| `NaturalStrain4` | `(V0, K0, K0_prime, K0_double_prime)` | `K0_double_prime` has inverse-pressure units |
| `ModifiedTait` | `(V0, K0, K0_prime, K0_double_prime)` | rejects singular coefficient sets and volumes outside its real domain |
| `Vinet` | `(V0, K0, K0_prime)` | none |
| `Holzapfel` | `(V0, K0, K0_prime, n, Z)` | molar volume in `J bar^-1 mol^-1` |

`Holzapfel` additionally provides `bulk_modulus_derivative(V, eps=1e-6)`.
See the [equation reference](equation-reference.md#isothermal-equations) for
the mathematical definitions and coefficient domains.

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

Thermal constructor signatures are:

| Class | Parameters after `rt_eos` |
|---|---|
| `MieGruneisenDebye` | `Tr, theta0, gamma0, q, n` |
| `MieGruneisenEinstein` | `Tr, theta0, gamma0, q, n` |
| `ThermalModifiedTait` | `Tr, theta, alpha0, n` |
| `Sokolova2016` | `Tr, QE1o, mE1, QE2o, mE2, delta, t, a_0, m, g, e_0`, followed by optional `beta, QBo, d, mb, QB1o, d1, mb1` |

The Mie-Gruneisen classes accept any `EosBase` reference; thermal modified Tait
requires `ModifiedTait`, and Sokolova requires `Holzapfel`. All thermal classes
require molar volume in `J bar^-1 mol^-1`. `HollandPowell2011` is an alias for
`ThermalModifiedTait`. Exact equations and parameter roles are documented
under [Thermal equations](equation-reference.md#thermal-equations).

Common methods:

- `thermal_pressure(V, T)`
- `dac_thermal_pressure(V, T, f_dac)`
- `pressure(V, T)`
- `volume(P, T)`
- `temperature(P, V)` and `calculate_temperature(P, V)`
- `temperature_from_volumes(V_ambient, V_heated, f_dac=...)`
- `bulk_modulus(V, T)`
- `isothermal_compressibility(V, T)`
- `thermal_expansivity(V, T)`
- `molar_heat_capacity_v(V, T)` when a caloric model exists
- `molar_heat_capacity_p(V, T)` when a caloric model exists
- `adiabatic_bulk_modulus(V, T)` when a caloric model exists
- `gruneisen_parameter(V, T)` when a caloric model exists

`dac_thermal_pressure()` returns only the additional confinement term
`f_dac * thermal_pressure(V, T)`. `temperature_from_volumes()` applies the
empirical confinement model described in
[Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md);
it requires `0 <= f_dac < 1`. In this API, `f_dac` means
`(P_hot - P_ambient) / thermal_pressure(V_heated, T)`; it is not a fraction of
the cold pressure.

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
