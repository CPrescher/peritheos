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

## Thermal models

A thermal EOS wraps an isothermal `rt_eos`, which defines pressure on the
reference isotherm at `Tr`. The thermal model adds the pressure change away
from that temperature, so the combined model evaluates
$P(V,T)=P_{ref}(V)+\Delta P_{th}(V,T)$ with
$\Delta P_{th}(V,T_r)=0$. Pressure, volume, and temperature inversion all use
this combined relation. The reference and thermal parameters can also be fitted
together with [`fit_joint_eos`](fitting.md#joint-reference-and-thermal-fitting).

The allowed reference EOS depends on the thermal formulation. The
Mie-Gruneisen models accept any isothermal `EosBase` model;
`ThermalModifiedTait` requires `ModifiedTait`; and `Sokolova2016` requires
`Holzapfel`. Because thermal pressure is calculated from molar energy divided
by molar volume, the wrapped reference EOS must use the thermal-model volume
convention described under [Units and reference states](units.md).

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

### Selecting a model

- Use BM3 or Vinet for a conventional three-parameter compression fit.
- Use modified Tait or a fourth-order natural-strain EOS when reliable second
  pressure-derivative information exists.
- Use Mie-Gruneisen-Debye when low-temperature heat-capacity behavior matters.
- Use Mie-Gruneisen-Einstein for an inexpensive single-frequency approximation.
- Use thermal modified Tait with Holland-Powell-style datasets.
- Use Sokolova2016 for the workbook-compatible diamond-oriented multimode
  formulation accompanying Sokolova et al. (2016).

Do not choose a higher-order model solely because it has more parameters. Check
parameter correlations and extrapolation behavior with the fitting diagnostics.
