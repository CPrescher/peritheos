# Diamond-anvil-cell thermal-pressure contribution

These advanced workflows predict a confined hot state from a measured cold
pressure or estimate temperature from separate cold and hot volumes. They
introduce an experiment-specific boundary-condition parameter and should only
be used when the total high-temperature pressure is not measured independently.
Review the [thermal model](models.md#thermal-models), [unit](units.md), and
[uncertainty](uncertainty.md) conventions before applying them.

## Forward prediction from a cold pressure

For phase-line prediction, first establish the pressure on the
reference-temperature isotherm and then choose the hot temperature and retained
fraction. Peritheos solves

\[
P_{\mathrm{EOS}}(V_{\mathrm{hot}},T)
=P_{\mathrm{cold}}
+f_{\mathrm{DAC}}\,\Delta P_{\mathrm{thermal}}(V_{\mathrm{hot}},T).
\]

The thermal increment depends on the unknown hot volume, so evaluating it once
at either the cold or isobaric volume is not self-consistent. The forward API
solves for that volume directly:

```python
cold_pressure = 40.0  # GPa at the reference temperature
hot_temperature = 2000.0  # K
f_dac = 0.25

heated_volume = eos.volume_with_dac_confinement(
    cold_pressure,
    hot_temperature,
    f_dac=f_dac,
)
thermal_pressure = eos.thermal_pressure_increment(heated_volume, hot_temperature)
confinement_pressure = eos.dac_thermal_pressure(heated_volume, hot_temperature, f_dac)
total_hot_pressure = cold_pressure + confinement_pressure
```

`thermal_pressure` is the full EOS heating increment, while
`confinement_pressure` is only the fraction retained by the DAC. The latter is
what raises the experimental pressure. An `EOSRecord` exposes the same three
methods and accepts conventional-cell volumes, which lets diffraction
applications predict line positions without performing molar-volume conversion.

At `f_dac=0`, the result is the ordinary isobaric `volume(cold_pressure, T)`.
Increasing the fraction reduces the predicted thermal expansion and moves the
volume towards its cold value.

## Two-volume temperature inversion

The ambient pressure is calculated from the volume measured at the reference
temperature. A second volume is measured during heating. The two-volume
temperature inversion solves

\[
\begin{aligned}
P_{\mathrm{ambient}} &= P_{\mathrm{ref}}(V_{\mathrm{ambient}}), \\
P_{\mathrm{EOS}}(V_{\mathrm{heated}}, T)
&= P_{\mathrm{ambient}}
   + f_{\mathrm{DAC}}\,\Delta P_{\mathrm{thermal}}(V_{\mathrm{heated}}, T).
\end{aligned}
\]

Here

\[
P_{\mathrm{ref}}(V)=P_{\mathrm{EOS}}(V,T_r),\qquad
\Delta P_{\mathrm{thermal}}(V,T)=P_{\mathrm{EOS}}(V,T)-P_{\mathrm{ref}}(V).
\]

The volume-pair equation therefore reduces exactly to

\[
\Delta P_{\mathrm{thermal}}(V_{\mathrm{heated}}, T)
= \frac{P_{\mathrm{ref}}(V_{\mathrm{ambient}})
       - P_{\mathrm{ref}}(V_{\mathrm{heated}})}
      {1-f_{\mathrm{DAC}}}.
\]

For the usual thermal models, `P_ref` is the stored cold/reference isotherm and
their existing `thermal_pressure` already vanishes at $T_r$. For
the double-Debye Helmholtz classes, `P_ref` is their complete 300 K isotherm and subtraction
removes the zero-point and 300 K thermal baseline only from the heating
increment—not from either total pressure. This makes the volume-pair method
compatible with their absolute Helmholtz formulation.

Here `f_dac` is defined specifically as the fraction of the EOS thermal pressure
that appears as an increase above the reference-temperature pressure:

\[
\begin{aligned}
\Delta P_{\mathrm{DAC}} &= P_{\mathrm{hot}} - P_{\mathrm{ambient}}, \\
f_{\mathrm{DAC}} &=
\frac{\Delta P_{\mathrm{DAC}}}
     {\Delta P_{\mathrm{thermal}}(V_{\mathrm{heated}}, T)}.
\end{aligned}
\]

It is not the ratio `Delta P_DAC / P_ambient` and therefore is not a percentage
of the cold pressure.

## Using the two-volume method

```python
f_dac = 0.25
ambient_volume = 0.80000
heated_volume = 0.80001
corrected_temperature = eos.temperature_from_volumes(
    ambient_volume,
    heated_volume,
    f_dac=f_dac,
)
ambient_pressure = eos.pressure(ambient_volume, eos.Tr)
confinement_pressure = eos.dac_thermal_pressure(
    heated_volume, corrected_temperature, f_dac
)
high_temperature_pressure = ambient_pressure + confinement_pressure
```

If total high-temperature pressure is measured independently, use
`temperature(P_hot, V_heated)` directly and do not apply the two-volume model.
See the [thermal API reference](api.md#thermal-equations-of-state) for both
inversion methods.

## Physical interpretation and limits

The two-volume inversion requires `0 <= f_dac < 1`:

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

## Experimental cautions

`f_dac` depends on the sample and pressure-medium strength, hotspot and sample
geometry, temperature distribution, probe position, heating duration, and
stress relaxation. It may vary spatially, with temperature, and between heating
cycles. A literature value from a different assembly is therefore best treated
as a sensitivity scenario rather than a calibrated correction.

The two measured volumes cannot determine both temperature and `f_dac`.
Assuming `f_dac` selects one temperature from a family of possible solutions.
If an independent temperature is available, an experiment-specific fraction
can instead be estimated from

\[
f_{\mathrm{DAC}} =
\frac{P_{\mathrm{EOS}}(V_{\mathrm{heated}}, T_{\mathrm{independent}})
      - P_{\mathrm{ambient}}}
     {\Delta P_{\mathrm{thermal}}
      (V_{\mathrm{heated}}, T_{\mathrm{independent}})}.
\]

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
