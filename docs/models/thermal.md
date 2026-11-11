# Thermal models

Select a model name in the inventory below for its complete equations and limits.

## Reference EOS and volume conventions

Most thermal EOS classes wrap an isothermal `rt_eos`, which defines pressure on the
reference isotherm at `Tr`. The thermal model usually adds the pressure change away
from that temperature, so the combined model evaluates
$P(V,T)=P_{ref}(V)+\Delta P_{th}(V,T)$ with
$\Delta P_{th}(V,T_r)=0$. Pressure, volume, and temperature inversion all use
this combined relation. The reference and thermal parameters can also be fitted
together with [`fit_joint_eos`](../fitting.md#joint-reference-and-thermal-fitting).
`SecondOrderTaylorThermalPressure` instead adds an absolute pressure polynomial
to a cold curve; its `Tr` is an expansion coordinate, not the temperature of
the reference EOS.
`AsymptoticDebyeTabulatedPressure` likewise uses an absolute 0 K cold curve,
with a separate ambient volume for its phonon law and bounded electronic and
residual pressure tables. See the [Yokoo PVT model](../literature-reproductions/yokoo-2009-pvt-library.md)
for its reference conventions and derived reconstruction limits.
`DoubleDebyeHelmholtz` and `DoubleDebyeLogMomentHelmholtz` can either represent
the absolute simulated free energy or, when their optional `Tr` is supplied,
add the simulated non-cold contribution relative to `Tr` to an experimental
reference isotherm.

Both double-Debye models support `Vinet`, `BM2`, `BM3`, and `BM4` with a
matching reference energy, so total pressure remains the negative volume
derivative of total Helmholtz energy. The same choices apply in Python, Rust,
and `.eosmat` files. Changing a published reference curve creates a new
composition whose accuracy must be checked separately.

The allowed reference EOS depends on the thermal formulation. The
Mie-Gruneisen, multi-oscillator, and constant linear thermal-pressure models
accept any isothermal `EosBase` model;
`LogVolumeThermalPressure` additionally requires `V0`,
`HollandPowellThermalPressure` requires `K0`, and `ThermalModifiedTait`
requires `ModifiedTait`.
Because thermal pressure is calculated from molar energy divided
by molar volume, energy-based models require the volume convention described
under [Units and reference states](../units.md). The volume-independent linear
correction instead inherits the reference EOS volume convention.

For Debye and Einstein models, `n` is the number of atoms per formula unit,
with the same formula unit defining the molar volume. It is not $Z$, the
number of formula units in a crystallographic cell. For hcp Fe, the standard
choice is `n=1` with volume per mole of Fe; $Z=2$ is used to convert its
two-atom cell volume. See the
[MGD equations, derivations, and worked calculation](../equations/thermal/mie-gruneisen.md)
and [normalization examples](../units.md#atom-count-and-molar-volume-basis).

## Models and parameters

| Import | Reference EOS | Thermal parameters | Caloric model |
|---|---|---|---|
| [`AsymptoticDebyeTabulatedPressure`](../literature-reproductions/yokoo-2009-pvt-library.md) | cold `EosBase` exposing `V0` | `Tr`, `theta0`, `gamma0`, `a`, `b`, `n`, `cold_volume_ratio`; pressure tables and numerical bounds | pressure only |
| [`DoubleDebyeHelmholtz`](../equations/thermal/double-debye.md) | `Vinet`, `BM2`, `BM3`, or `BM4`: 0 K cold curve, or reference isotherm when `Tr` is set | `Vp`; three sets of `theta_*0`, `a_*`, `b_*`; optional `n`, `alpha0`, `Ve`, `kappa`, `phi0`, `Tr` | double Debye + $T^2$ |
| [`DoubleDebyeLogMomentHelmholtz`](../equations/thermal/double-debye.md#logarithmic-moment-double-debye-variant) | `Vinet`, `BM2`, `BM3`, or `BM4`: 0 K cold curve, or reference isotherm when `Tr` is set | `Vp`; cutoff A, cutoff B, and `theta_0` parameter triples; optional `n`, `anharmonic_a`, `phi0`, `Tr` | logarithmic-moment double Debye + $T^2$ |
| [`MieGruneisenDebye`](../equations/thermal/mie-gruneisen.md) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `q`, `n`; optional `Cvmax`, `debye_temperature_law`, `thermal_pressure_reference` | Debye |
| [`MieGruneisenEinstein`](../equations/thermal/mie-gruneisen.md) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `q`, `n` | Einstein |
| [`Tange2009Debye`](../equations/thermal/tange-2009.md) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `a`, `b`, `n` | Debye |
| [`AsymptoticPowerLawMieGruneisenDebyeExcess`](../equations/thermal/asymptotic-debye-excess.md) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `a`, `b`, `n`, `beta0`, `m` | Debye plus volume-dependent $T^2$ excess |
| [`SoundVelocityDebyeHelmholtz`](../equations/thermal/sound-velocity-debye.md) | any `EosBase` with `bulk_modulus` | `Tr`, molar mass, `n`, and longitudinal/shear velocity intercepts and slopes | sound-derived Debye |
| [`Dewaele2006`](../equations/thermal/dewaele-2006.md) | any `EosBase` | `Tr`, `theta0`, `gamma0`, `gamma_inf`, `beta`, anharmonic and electronic terms, `n` | Debye + $T^2$ |
| [`DorogokupetsOganov2007`](../equations/thermal/dorogokupets-oganov-2007.md) | any `EosBase` | Four oscillator modes, `gamma0`, `gamma_inf`, `beta`, anharmonic, electronic, and defect parameters | generalized Bose + Einstein |
| [`LinearThermalPressure`](../equations/thermal/linear-thermal-pressure.md) | any `EosBase` | `Tr`, `alpha_KT` | none |
| [`SecondOrderTaylorThermalPressure`](../equations/thermal/second-order-taylor.md) | reference EOS exposing `V0` | `Tr`, `eta0`, `c0`--`c5` | none |
| [`LogVolumeThermalPressure`](../equations/thermal/log-volume-thermal-pressure.md) | reference EOS exposing `V0` | `Tr`, `alpha_KT_ref`, `dK_dT_V` | none |
| [`ThermalReferenceStateEOS`](../equations/thermal/reference-state.md) | reference EOS exposing `V0`, `K0` | `Tr`, `alpha0`, `dK_dT`; optional `alpha1`, `thermal_expansion_law`, `reference_volume_law` | none |
| [`HollandPowellThermalPressure`](../equations/thermal/holland-powell.md) | reference EOS exposing `K0` | `Tr`, `theta`, `alpha0`, `n` | Holland-Powell Einstein pressure |
| [`ThermalModifiedTait`](../equations/thermal/holland-powell.md) | `ModifiedTait` | `Tr`, `theta`, `alpha0`, `n` | Holland-Powell Einstein pressure |
| [`MultiOscillatorGruneisenThermalEOS`](../equations/thermal/multi-oscillator.md) | any `EosBase` | mode, Gruneisen, anharmonic, electronic parameters plus `n` | Multi-mode |

## Compatibility and source-specific conventions

`ThermalModifiedTait` is the modified-Tait compatibility form;
`HollandPowell2011` remains its alias.

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
[Paper versus spreadsheet](../equations/thermal/multi-oscillator.md#paper-versus-spreadsheet).

[Equation index](../equation-reference.md) · [Choosing a model](../models.md#selecting-a-model)
