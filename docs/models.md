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

Most thermal EOS classes wrap an isothermal `rt_eos`, which defines pressure on the
reference isotherm at `Tr`. The thermal model adds the pressure change away
from that temperature, so the combined model evaluates
$P(V,T)=P_{ref}(V)+\Delta P_{th}(V,T)$ with
$\Delta P_{th}(V,T_r)=0$. Pressure, volume, and temperature inversion all use
this combined relation. The reference and thermal parameters can also be fitted
together with [`fit_joint_eos`](fitting.md#joint-reference-and-thermal-fitting).
`DoubleDebyeHelmholtz` is the full-free-energy exception described below.

The allowed reference EOS depends on the thermal formulation. The
Mie-Gruneisen, multi-oscillator, and constant linear thermal-pressure models
accept any isothermal `EosBase` model;
`LogVolumeThermalPressure` additionally requires `V0`, and
`ThermalModifiedTait` requires `ModifiedTait`.
Because thermal pressure is calculated from molar energy divided
by molar volume, energy-based models require the volume convention described
under [Units and reference states](units.md). The volume-independent linear
correction instead inherits the reference EOS volume convention.

| Import | Reference EOS | Thermal parameters | Caloric model |
|---|---|---|---|
| [`DoubleDebyeHelmholtz`](equation-reference.md#double-debye-helmholtz) | `Vinet` **0 K cold curve** | `Vp`; three sets of `theta_*0`, `a_*`, `b_*`; optional `n`, `alpha0`, `Ve`, `kappa`, `phi0` | double Debye + $T^2$ |
| [`MieGruneisenDebye`](equation-reference.md#mie-gruneisen-debye-and-einstein) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `q`, `n`; optional `debye_temperature_law` | Debye |
| [`MieGruneisenEinstein`](equation-reference.md#mie-gruneisen-debye-and-einstein) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `q`, `n` | Einstein |
| [`Tange2009Debye`](equation-reference.md#tange-2009-mgo-thermal-model) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `a`, `b`, `n` | Debye |
| [`LinearThermalPressure`](equation-reference.md#linear-thermal-pressure) | any `EosBase` | `Tr`, `alpha_KT` | none |
| [`LogVolumeThermalPressure`](equation-reference.md#logarithmic-volume-linear-thermal-pressure) | reference EOS exposing `V0` | `Tr`, `alpha_KT_ref`, `dK_dT_V` | none |
| [`ThermalReferenceStateEOS`](equation-reference.md#temperature-dependent-reference-state) | reference EOS exposing `V0`, `K0` | `Tr`, `alpha0`, `dK_dT`; optional `alpha1`, `thermal_expansion_law`, `reference_volume_law` | none |
| [`ThermalModifiedTait`](equation-reference.md#thermal-modified-tait) | `ModifiedTait` | `Tr`, `theta`, `alpha0`, `n` | Holland-Powell Einstein pressure |
| [`MultiOscillatorGruneisenThermalEOS`](equation-reference.md#multi-oscillator-gruneisen-thermal-pressure) | any `EosBase` | mode, Gruneisen, anharmonic, electronic parameters plus `n` | Multi-mode |

`HollandPowell2011` is an alias for `ThermalModifiedTait`.

`MultiOscillatorGruneisenThermalEOS` accepts the complete optional `beta` and generalized Bose-mode
terms. Its defaults (`beta=mb=mb1=0`) disable those additions and preserve the
reduced configuration used by existing Peritheos and simplified Dioptas-style
parameter sets. Temperature inversion, including `temperature_from_volumes`,
uses the same active terms and works with both configurations.

The Sokolova catalog composition of `MultiOscillatorGruneisenThermalEOS` is intentionally compatible with the calculation path and
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
  Its `debye_temperature_law` defaults to `integrated_gruneisen`; select
  `variable_exponent` only for scales that explicitly publish
  $\Theta_D=\Theta_0(V/V_0)^{-\gamma(V)}$.
- Use Mie-Gruneisen-Einstein for an inexpensive single-frequency approximation.
- Use `Tange2009Debye` for the publication-specific Tange MgO Gruneisen law.
- Use linear thermal pressure only when the primary calibration reports a
  constant fitted $\alpha K_T$ product over the intended range.
- Use thermal modified Tait with Holland-Powell-style datasets.
- Use `MultiOscillatorGruneisenThermalEOS` for the multimode formulation; pass
  `n` explicitly and select the reference isotherm independently. Catalog
  entries pin the exact source-validated combinations.
- Use `DoubleDebyeHelmholtz` when the source supplies one thermodynamically
  complete Vinet cold curve plus volume-dependent double-Debye and $T^2$
  free-energy terms. Its `rt_eos` argument is a motionless-ion 0 K cold curve,
  not a room-temperature reference isotherm.

Do not choose a higher-order model solely because it has more parameters. Check
parameter correlations and extrapolation behavior with the fitting diagnostics.

`Sokolova2016` remains a compatibility alias. New code and serialized model
identifiers use neutral mechanism-based names; paper names belong to catalog
record identifiers and provenance.
