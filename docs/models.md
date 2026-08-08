# EOS models

## Isothermal models

| Import | Parameters after `V0`, `K0` | Notes |
|---|---|---|
| `BM2` | none | Second-order Birch-Murnaghan; implied `K0_prime = 4` |
| `BM3` | `K0_prime` | Standard third-order Eulerian finite strain |
| `BM4` | `K0_prime`, `K0_double_prime` | Fourth-order Eulerian finite strain |
| `Murnaghan` | `K0_prime` | Assumes bulk modulus linear in pressure |
| `NaturalStrain2` | none | Hencky strain; implied `K0_prime = 2` |
| `NaturalStrain3` | `K0_prime` | Third-order Poirier-Tarantola logarithmic EOS |
| `NaturalStrain4` | `K0_prime`, `K0_double_prime` | Fourth-order natural strain |
| `ModifiedTait` | `K0_prime`, `K0_double_prime` | Analytically invertible Tait family |
| `Vinet` | `K0_prime` | Rydberg-Vinet interatomic-potential form |
| `Holzapfel` | `K0_prime`, `n`, `Z` | High-compression form with strict molar units |

All isothermal models implement:

```python
pressure(V)
bulk_modulus(V)
volume(P)
```

### Natural strain

With positive compressive strain

\[
f = \frac{1}{3}\ln\left(\frac{V_0}{V}\right),
\]

the implemented family is

\[
P = 3K_0\frac{V_0}{V}(f + Af^2 + Bf^3).
\]

For third order, `A = 3(K0_prime - 2)/2` and `B = 0`. Fourth order
sets `B` from `K0_double_prime`. At the reference volume the supplied values
of `K0`, `K0_prime`, and `K0_double_prime` are recovered by differentiation.

## Thermal models

| Import | Reference EOS | Thermal parameters | Caloric model |
|---|---|---|---|
| `MieGruneisenDebye` | any `EosBase` | `Tr`, `theta0`, `gamma0`, `q`, `n` | Debye |
| `MieGruneisenEinstein` | any `EosBase` | `Tr`, `theta0`, `gamma0`, `q`, `n` | Einstein |
| `ThermalModifiedTait` | `ModifiedTait` | `Tr`, `theta`, `alpha0`, `n` | Holland-Powell Einstein pressure |
| `Sokolova2016` | `Holzapfel` | mode, Gruneisen, anharmonic, electronic parameters | Multi-mode |

`HollandPowell2011` is an alias for `ThermalModifiedTait`.

### Mie-Gruneisen volume dependence

\[
\gamma(V)=\gamma_0\left(\frac{V}{V_0}\right)^q,
\qquad
\gamma=-\frac{\partial\ln\Theta}{\partial\ln V}.
\]

The second identity is enforced analytically, including the continuous `q=0`
limit. Thermal pressure is referenced to `Tr`:

\[
\Delta P_{th}=\frac{\gamma(V)}{V}
\left[E(V,T)-E(V,T_r)\right].
\]

### Selecting a model

- Use BM3 or Vinet for a conventional three-parameter compression fit.
- Use modified Tait or fourth-order natural strain when reliable second
  pressure-derivative information exists.
- Use Mie-Gruneisen-Debye when low-temperature heat-capacity behavior matters.
- Use Mie-Gruneisen-Einstein for an inexpensive single-frequency approximation.
- Use thermal modified Tait with Holland-Powell-style datasets.
- Use Sokolova2016 for its published diamond-oriented multimode formulation.

Do not choose a higher-order model solely because it has more parameters. Check
parameter correlations and extrapolation behavior with the fitting diagnostics.
