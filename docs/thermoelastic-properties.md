# Thermoelastic properties

Every `ThermalEOS` provides pressure-volume inversion and numerical
thermoelastic derivatives. Models with a caloric description additionally
provide heat capacities and adiabatic properties.

## Mechanical derivatives

The isothermal bulk modulus and compressibility are

\[
K_T=-V\left(\frac{\partial P}{\partial V}\right)_T,
\qquad
\kappa_T=\frac{1}{K_T}.
\]

```python
kt = eos.bulk_modulus(V, T)                 # GPa
kappa = eos.isothermal_compressibility(V, T)  # GPa^-1
```

Thermal expansivity follows the exact identity

\[
\alpha=\frac{1}{K_T}\left(\frac{\partial P}{\partial T}\right)_V.
\]

```python
alpha = eos.thermal_expansivity(V, T)  # K^-1
```

These derivatives use centered finite differences with steps relative to the
state variables. They therefore include both reference-isotherm and thermal
pressure contributions.

## Heat capacities and adiabatic modulus

Caloric thermal models implement `molar_heat_capacity_v()`. Peritheos then uses

\[
C_P-C_V=\alpha^2 K_TVT,
\qquad
K_S=K_T\frac{C_P}{C_V}.
\]

The conversion from GPa and J bar^-1 mol^-1 to J mol^-1 is handled internally.

```python
cv = eos.molar_heat_capacity_v(V, T)  # J mol^-1 K^-1
cp = eos.molar_heat_capacity_p(V, T)  # J mol^-1 K^-1
ks = eos.adiabatic_bulk_modulus(V, T) # GPa
gamma = eos.gruneisen_parameter(V, T) # dimensionless
```

`Sokolova2016` does not yet expose a complete caloric potential and therefore
raises `NotImplementedError` for heat capacity and adiabatic modulus. Its
mechanical derivatives remain available.

## Vibrational thermodynamic contributions

The Debye and Einstein Mie-Gruneisen models expose:

```python
e = eos.thermal_internal_energy(V, T)
s = eos.thermal_entropy(V, T)
f = eos.thermal_helmholtz_free_energy(V, T)
h = eos.thermal_enthalpy(V, T)
g = eos.thermal_gibbs_free_energy(V, T)
```

Energy values are in J mol^-1 and entropy in J mol^-1 K^-1. Enthalpy and Gibbs
energy use the unreferenced `vibrational_pressure()`, whereas the public
`thermal_pressure()` is referenced to `Tr`. These methods omit static formation
energies and zero-point offsets. They are vibrational
contributions suitable for differences within the same reference convention;
they are not absolute chemical potentials.
