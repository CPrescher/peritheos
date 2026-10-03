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

The thermal-volume conversions are

\[
1\ \mathrm{J\,bar^{-1}\,mol^{-1}}
=10\ \mathrm{cm^3\,mol^{-1}}
=10^{-5}\ \mathrm{m^3\,mol^{-1}}.
\]

Two conversion factors recur in the thermal equations:

\[
1\ \mathrm{bar}=10^{-4}\ \mathrm{GPa},
\qquad
1\ \mathrm{GPa}\times1\ \mathrm{J\,bar^{-1}\,mol^{-1}}
=10^4\ \mathrm{J\,mol^{-1}}.
\]

Consequently, $\gamma E/V$ is divided by $10^4$ to return GPa, whereas
$K_TV$ is multiplied by $10^4$ when used in an energy or heat-capacity
identity. These factors are explicit in the
[equation reference](equation-reference.md#thermal-equations).

Published atom-normalized Helmholtz parameter sets require one further
conversion. For a monatomic model,

\[
1\ \mathrm{\AA^3/atom}
=N_A10^{-25}\ \mathrm{J\,bar^{-1}\,mol^{-1}}
=0.0602214076\ \mathrm{J\,bar^{-1}\,mol^{-1}},
\]

and $1\ \mathrm{eV/atom}=96485.3321\ \mathrm{J/mol}$. An inverse-volume
coefficient printed in $\mathrm{\AA^{-3}}$ is divided by `0.0602214076` for
use with a molar volume. For a formula unit containing `n` atoms, use its
molar formula-unit volume and set `n` accordingly.

## Atom count and molar-volume basis

In Debye and Einstein thermal models, `n` counts atoms in the entity that
defines one mole of energy and volume. The crystallographic $Z$ counts those
formula units in the measured unit cell. They are distinct:

\[
V_{m,\mathrm{formula}}[\mathrm{cm^3/mol}]
=\frac{N_A V_{\mathrm{cell}}[\mathrm{\AA^3}]\,10^{-24}}{Z},
\qquad n_{\mathrm{cell}}=Zn_{\mathrm{formula}}.
\]

| Material and basis | `n` | $Z$ for conversion from the conventional cell |
|---|---:|---:|
| hcp Fe, per mole of Fe | 1 | 2 |
| hcp Fe, per mole of two-atom cells | 2 | 1 |
| B1 MgO, per mole of MgO | 2 | 4 |
| B1 MgO, per mole of eight-atom cells | 8 | 1 |

The $Z=1$ rows treat the entire crystallographic cell as the molar entity;
they do not redefine the material's ordinary chemical formula or crystallography.

For hcp Fe with $V_{0,\mathrm{cell}}=22.81$ Å$^3$:

| Basis | $V_0$ in cm$^3$/mol (EosFit input) | $V_0$ in J/bar/mol (Peritheos thermal input) | `n` |
|---|---:|---:|---:|
| Mole of Fe atoms | 6.86825 | 0.686825 | 1 |
| Mole of two-atom cells | 13.73650 | 1.373650 | 2 |

Both describe the same physical full-MGD pressure surface when all molar
volumes, reference volumes and molar energies use the matching basis.
Doubling only `n` at unchanged $V$ doubles the Debye thermal-pressure
amplitude at fixed $\gamma$ and $\Theta$. In a full-MGD refit, halving
$\gamma_0$ does not exactly undo this change because $\gamma$ also determines
$\Theta(V)$. A source-reproduction setting therefore needs its own explicit
normalization provenance; agreement with author coefficients alone does
not verify physical normalization. See the
[Miozzi normalization replay](literature-reproductions/miozzi-2020-iron.md#author-confirmed-full-mgd-with-n2-2026-10-03).

If `Cvmax` replaces the usual $3nR$ heat-capacity limit, scale it with the
chosen molar entity as well. For example, a consistent per-cell Fe model has
twice the molar `Cvmax` of the per-Fe model.

## Conversion helpers

```python
from peritheos.units import (
    cell_volume_to_molar_volume,
    convert_molar_volume,
    convert_pressure,
    convert_temperature,
    density_from_molar_volume,
    molar_volume_to_cell_volume,
    molar_volume_from_density,
)

thermal_volume = convert_molar_volume(11.25, "cm3/mol", "J/bar/mol")
density = density_from_molar_volume(40.304, 11.25)  # g/cm3
volume = molar_volume_from_density(40.304, density)  # cm3/mol
pressure = convert_pressure(100.0, "kbar", "GPa")
temperature = convert_temperature(25.0, "C", "K")

# MgO B1 has four formula units in its conventional cubic cell.
mgo_molar_volume = cell_volume_to_molar_volume(74.698, 4)
mgo_cell_volume = molar_volume_to_cell_volume(mgo_molar_volume, 4)
```

Supported molar-volume units are m3/mol, cm3/mol, L/mol, and J/bar/mol.
Supported density units are kg/m3 and g/cm3. Exponent spellings such as
`cm^3/mol` and `J bar^-1 mol^-1` are accepted.

`cell_volume_to_molar_volume()` and its inverse use conventional-cell volume
in angstrom cubed and return or consume molar volume per mole of formula units.
The required `formula_units_per_cell` argument is the crystallographic `Z`,
not the number of atoms in the cell. Scalars and NumPy arrays are accepted by
all conversion helpers.

Pressure and temperature conversions formerly imported from
`peritheos.utils` now live in `peritheos.units`. The old imports remain as
deprecated compatibility wrappers.

## Reference states

Isothermal classes define `V0` at zero model pressure. Thermal pressure is a
difference relative to `Tr`, so

\[
\Delta P_{\mathrm{th}}(V,T_r)=0
\]

for every valid volume. Consequently, `rt_eos` must represent the same
reference temperature supplied as `Tr`. Peritheos does not silently translate
an isotherm between reference temperatures.

With `Tr=None`, the double-Debye Helmholtz classes consume a 0 K
motionless-ion cold curve (`Vinet`, `BM2`, `BM3`, or `BM4`) and add absolute
ionic (including zero-point) and anharmonic free energies. Their cold-curve
`V0` therefore must not be interpreted as an ambient-temperature zero-pressure
volume. With a numeric `Tr`, they instead treat the supplied curve as the
complete isotherm at that temperature and subtract the non-cold contribution
there.

## Shock Hugoniot units and initial states

`LinearUsUpHugoniot` uses density in g/cm³, particle and shock velocity in
km/s, and pressure in GPa. In that combination the momentum jump condition
`P - P0 = rho0 * Us * up` needs no numerical conversion factor. Specific
internal-energy changes are returned in MJ/kg. Hugoniot `V0` and `V` may use
any consistent volume unit. EOSMAT records expose represented-phase
conventional-cell Å³ and require a `volume_basis` stating how many formula
units both `V` and `V0` represent and the molar mass of one formula. Peritheos
checks these values against `rho0` and `V0`. This keeps transformed branches
mass-normalized when precursor and product unit cells differ.

The initial temperature belongs to `initial_state` metadata. It is not an
independent argument to a Hugoniot evaluation.
