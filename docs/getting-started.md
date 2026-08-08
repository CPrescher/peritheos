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
kt = eos.bulk_modulus(volume, temperature)
alpha = eos.thermal_expansivity(volume, temperature)
cv = eos.molar_heat_capacity_v(volume, temperature)
cp = eos.molar_heat_capacity_p(volume, temperature)
ks = eos.adiabatic_bulk_modulus(volume, temperature)
```

See [Thermoelastic properties](thermoelastic-properties.md) before using
caloric or free-energy methods.

## Broadcasting P and T

Thermal state variables use NumPy broadcasting. A volume column and temperature
row produce a full state grid:

```python
volumes = np.array([0.8, 0.9, 1.0])[:, None]
temperatures = np.array([300.0, 1000.0, 2000.0])[None, :]
pressure_grid = eos.pressure(volumes, temperatures)
```
