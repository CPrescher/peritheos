# Getting started

## Room-temperature pressure and volume

```python
from peritheos.eos.rt import Vinet

eos = Vinet(V0=10.0, K0=160.0, K0_prime=4.2)
pressure = eos.pressure(8.5)
bulk_modulus = eos.bulk_modulus(8.5)
recovered_volume = eos.volume(pressure)
```

`V0` and `V` may use any consistent volume unit for an isothermal EOS. `K0`,
pressure, and bulk modulus must use the same pressure unit; Peritheos documents
and tests public pressure values as GPa.

All methods accept NumPy arrays:

```python
import numpy as np

volumes = np.linspace(8.0, 10.0, 21)
pressures = eos.pressure(volumes)
```

## Thermal pressure and properties

Thermal models require molar volume in J bar^-1 mol^-1. This equals molar
cm^3/mol divided by ten.

```python
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye

reference = BM3(V0=1.0, K0=160.0, K0_prime=4.0)
eos = MieGruneisenDebye(
    rt_eos=reference,
    Tr=300.0,
    theta0=800.0,
    gamma0=1.5,
    q=1.0,
    n=2,
)

volume = 0.9
temperature = 1800.0
pressure = eos.pressure(volume, temperature)
recovered_temperature = eos.temperature(pressure, volume)
kt = eos.bulk_modulus(volume, temperature)
alpha = eos.thermal_expansivity(volume, temperature)
cv = eos.molar_heat_capacity_v(volume, temperature)
cp = eos.molar_heat_capacity_p(volume, temperature)
ks = eos.adiabatic_bulk_modulus(volume, temperature)
```

See [Thermoelastic properties](thermoelastic-properties.md) before using
caloric or free-energy methods.

Use `temperature(P, V)` when the total pressure and volume at the heated state
are both known. If only the reference-temperature volume and heated volume are
measured, use the dedicated two-volume method below; the ambient pressure is not
the total pressure at high temperature.

## Diamond-anvil-cell thermal-pressure contribution

The ambient pressure is calculated from the volume measured at the reference
temperature. A second volume is measured during heating. The two-volume
temperature inversion solves

```text
P_ambient = P_cold(V_ambient),
P_EOS(V_heated, T)
    = P_ambient + f_DAC * Delta P_thermal(V_heated, T).
```

Because `P_EOS = P_cold + Delta P_thermal`, this reduces exactly to

```text
Delta P_thermal(V_heated, T)
    = [P_cold(V_ambient) - P_cold(V_heated)] / (1 - f_DAC).
```

The implementation solves this reduced equation directly. This avoids
re-evaluating both sides of the original expression and makes the domain clear:
for the usual positive, monotonic thermal-pressure models, a heated state needs
the heated volume to have the larger cold-compression pressure.

Here `f_dac` is defined specifically as the fraction of the EOS thermal pressure
that appears as an increase above the reference-temperature pressure:

```text
Delta P_DAC = P_hot - P_ambient
f_DAC = Delta P_DAC / Delta P_thermal(V_heated, T).
```

It is not the ratio `Delta P_DAC / P_ambient` and therefore is not a percentage
of the cold pressure.

Use the dedicated two-volume method:

```python
f_dac = 0.25
ambient_volume = 0.80000
heated_volume = 0.80001
corrected_temperature = eos.temperature_from_volumes(
    ambient_volume,
    heated_volume,
    f_dac=f_dac,
)
ambient_pressure = eos.rt_eos.pressure(ambient_volume)
thermal_pressure = eos.thermal_pressure(heated_volume, corrected_temperature)
high_temperature_pressure = ambient_pressure + f_dac * thermal_pressure
```

### Physical interpretation and limits

The two-volume inversion requires ``0 <= f_dac < 1``:

| `f_dac` | Boundary-condition interpretation | Consequence |
|---:|---|---|
| `0` | Isobaric heating | `P_hot = P_ambient` |
| between `0` and `1` | Partial confinement | Some EOS thermal pressure raises the experimental pressure |
| approaching `1` | Isochoric limit | `V_heated` approaches `V_ambient`; volume alone ceases to constrain temperature |

This interpolation is physically consistent with the isobaric-to-isochoric
range discussed by [Yen, Williams, and Kunz
(2020)](https://doi.org/10.1029/2020JB020006). It remains an experimental
boundary-condition model, not an additional term in the material EOS.

The meaning of percentages in the DAC literature must be checked before using
them as `f_dac`. [Heinz (1990)](https://doi.org/10.1029/GL017i008p01161)
estimated roughly 40--60% of the constant-volume thermodynamic thermal pressure
for a particular elastic hotspot model; this ratio is conceptually compatible
with `f_dac`. In contrast, the approximately 30% result attributed to the
Dewaele finite-element model by Yen et al. is a fraction of the *cold pressure*.
That is a different denominator and cannot be entered as `f_dac=0.3` without an
additional conversion.

### Experimental cautions

`f_dac` depends on the sample and pressure-medium strength, hotspot and sample
geometry, temperature distribution, probe position, heating duration, and
stress relaxation. It may vary spatially, with temperature, and between heating
cycles. A literature value from a different assembly is therefore best treated
as a sensitivity scenario rather than a calibrated correction.

The two measured volumes cannot determine both temperature and `f_dac`.
Assuming `f_dac` selects one temperature from a family of possible solutions.
If an independent temperature is available, an experiment-specific fraction
can instead be estimated from

```text
f_DAC = [P_EOS(V_heated, T_independent) - P_ambient]
        / Delta P_thermal(V_heated, T_independent).
```

If total high-temperature pressure is measured independently, use
`temperature(P_hot, V_heated)` directly and do not apply the two-volume model.

The reference volume also needs experimental scrutiny. Heating can cause
irreversible pressure or stress relaxation, so a post-cooling volume need not
equal the pre-heating reference. Prefer reference measurements bracketing the
heating sequence; interpolate a drifting baseline when justified. Otherwise,
define the stable interval, report how its average was selected, and propagate
its scatter rather than treating the average as exact.

For a fixed volume pair, increasing `f_dac` raises the required thermal pressure
and therefore normally raises the inferred temperature. The result is singular
as `f_dac` approaches one and can be highly sensitive to small errors in both
volumes, particularly when their difference is small. Sweep the plausible
`f_dac` range, propagate both volume errors, and check that the resulting
temperature remains inside the calibrated range of the selected EOS. Report the
uncorrected `f_dac=0` result alongside the sensitivity cases. See
[Uncertainty](uncertainty.md) and the [DAC references](references.md).

## Broadcasting P and T

Thermal state variables use NumPy broadcasting. A volume column and temperature
row produce a full state grid:

```python
volumes = np.array([0.8, 0.9, 1.0])[:, None]
temperatures = np.array([300.0, 1000.0, 2000.0])[None, :]
pressure_grid = eos.pressure(volumes, temperatures)
```
