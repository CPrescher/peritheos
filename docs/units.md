# Units and reference states

## Public conventions

| Quantity | Convention |
|---|---|
| Pressure and bulk modulus | GPa |
| Temperature | K |
| Isothermal volume | Any internally consistent unit |
| Thermal molar volume | J bar^-1 mol^-1 |
| Heat capacity and entropy | J mol^-1 K^-1 |
| Thermodynamic energy contributions | J mol^-1 |

One J bar^-1 mol^-1 equals 10 cm^3 mol^-1 and (10^{-5}) m^3 mol^-1.

## Conversion helpers

```python
from peritheos.units import (
    convert_molar_volume,
    density_from_molar_volume,
    molar_volume_from_density,
)

thermal_volume = convert_molar_volume(11.25, "cm3/mol", "J/bar/mol")
density = density_from_molar_volume(40.304, 11.25)  # g/cm3
volume = molar_volume_from_density(40.304, density)  # cm3/mol
```

Supported molar-volume units are m3/mol, cm3/mol, L/mol, and J/bar/mol.
Supported density units are kg/m3 and g/cm3. Exponent spellings such as
`cm^3/mol` and `J bar^-1 mol^-1` are accepted.

## Reference states

Isothermal classes define `V0` at zero model pressure. Thermal pressure is a
difference relative to `Tr`, so

```python
eos.thermal_pressure(V, eos.Tr) == 0
```

for every valid volume. Consequently, `rt_eos` must represent the same
reference temperature supplied as `Tr`. Peritheos does not silently translate
an isotherm between reference temperatures.
