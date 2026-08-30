# EOS models

## Isothermal models

| Import | Parameters after `V0`, `K0` | Notes |
|---|---|---|
| [`BM2`](equation-reference.md#birch-murnaghan-family) | none | Second-order Birch-Murnaghan; implied `K0_prime = 4` |
| [`BM3`](equation-reference.md#birch-murnaghan-family) | `K0_prime` | Standard third-order Eulerian finite strain |
| [`BM4`](equation-reference.md#birch-murnaghan-family) | `K0_prime`, `K0_double_prime` | Fourth-order Eulerian finite strain |
| [`Murnaghan`](equation-reference.md#murnaghan) | `K0_prime` | Assumes bulk modulus linear in pressure |
| [`NaturalStrain2`](equation-reference.md#natural-strain-family) | none | Hencky strain; implied `K0_prime = 2` |
| [`NaturalStrain3`](equation-reference.md#natural-strain-family) | `K0_prime` | Third-order Poirier-Tarantola logarithmic EOS |
| [`NaturalStrain4`](equation-reference.md#natural-strain-family) | `K0_prime`, `K0_double_prime` | Fourth-order natural strain |
| [`ModifiedTait`](equation-reference.md#modified-tait) | `K0_prime`, `K0_double_prime` | Analytically invertible Tait family |
| [`Vinet`](equation-reference.md#vinet) | `K0_prime` | Rydberg-Vinet interatomic-potential form |
| [`Holzapfel`](equation-reference.md#holzapfel) | `K0_prime`, `n`, `Z` | High-compression form with strict molar units |

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
The complete coefficient definitions are in the
[equation reference](equation-reference.md#natural-strain-family).

## Thermal models

| Import | Reference EOS | Thermal parameters | Caloric model |
|---|---|---|---|
| [`MieGruneisenDebye`](equation-reference.md#mie-gruneisen-debye-and-einstein) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `q`, `n` | Debye |
| [`MieGruneisenEinstein`](equation-reference.md#mie-gruneisen-debye-and-einstein) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `q`, `n` | Einstein |
| [`ThermalModifiedTait`](equation-reference.md#thermal-modified-tait) | `ModifiedTait` | `Tr`, `theta`, `alpha0`, `n` | Holland-Powell Einstein pressure |
| [`Sokolova2016`](equation-reference.md#sokolova2016) | `Holzapfel` | mode, Gruneisen, anharmonic, electronic parameters | Multi-mode |

`HollandPowell2011` is an alias for `ThermalModifiedTait`.

`Sokolova2016` accepts the complete optional `beta` and generalized Bose-mode
terms. Its defaults (`beta=mb=mb1=0`) disable those additions and preserve the
reduced configuration used by existing Peritheos and simplified Dioptas-style
parameter sets. Temperature inversion, including `temperature_from_volumes`,
uses the same active terms and works with both configurations.

`Sokolova2016` is intentionally compatible with the calculation path and
numerical results of the accompanying Excel workbook. It is not a literal
transcription of the equations typeset in the journal article, which differ
from the workbook in several consequential details. Do not use the printed
paper equations alone to reproduce Peritheos values; see
[Paper versus spreadsheet](equation-reference.md#paper-versus-spreadsheet).

### Mie-Gruneisen volume dependence

\[
\gamma(V)=\gamma_0\left(\frac{V}{V_0}\right)^q,
\qquad
\gamma=-\frac{\partial\ln\Theta}{\partial\ln V}.
\]

The second identity is enforced analytically, including the continuous `q=0`
limit. Thermal pressure is referenced to `Tr`:

\[
\Delta P_{th}=10^{-4}\frac{\gamma(V)}{V}
\left[E(V,T)-E(V,T_r)\right].
\]

The factor $10^{-4}$ converts bar to GPa for the public thermal-volume and
energy units. See the [equation reference](equation-reference.md#thermal-equations)
for the characteristic-temperature and oscillator-energy definitions.

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
