# Published Pt pressure scales checked against Fei

## Scope and equation checks

This comparison uses published coefficients without refitting or changing any
catalog default or validation status. The models are Fei (2004), Fei (2007),
Matsui (2009), and Zhu's released-v3 optimizer/property parameterization.

The supplied Matsui primary PDF was checked directly, including Equations
(5)-(10), Table II, Table III and the Fei comparison on printed page 013505-6.
Matsui's thermal model is the integrated constant-q Debye law plus the
Tsuchiya-Kawamura tabulated electronic pressure, referenced to 300 K. Linear
interpolation between electronic table nodes remains an explicit Peritheos
choice. The independent Table I/III audit passes; its largest calculated-table
pressure discrepancy is 0.008583 GPa. This is equation reconstruction, not a
replay of the original optimization.

Zhu's check uses the previously recovered released-v3 optimizer/property
constants, not the different preprint-v1 table or standalone calculator.
A fresh download of the original v3 release was unavailable during this check.
An additional SI implementation derives the excess pressure from
`F_ex = -beta0/2 (V/V0)^m T^2` and integrates the stated gamma law for theta.
It uses theta0=240 K, gamma0=2.75, a=0.39, b=5.1,
beta0=0.002145 J/mol/K^2 and m=0.65. The 85-row released thermal refit was
also rerun: gamma0=2.75224510 and b=5.10989878, reproducing the stored source
precision. This check does not resolve inconsistencies between source versions.

Independent pressure evaluation agrees with the production records over a
20-state grid (V/V0=0.7-1.0, T=300-3000 K):

| Model | Largest implementation difference (GPa) |
|---|---:|
| Fei 2004 | 4.33e-10 |
| Fei 2007 | 4.55e-13 |
| Matsui 2009 | 4.55e-13 |
| Zhu v3 | 3.74e-10 |

Independent integration also checks `ln(theta/theta0) = -integral gamma dlnV`
for Matsui and Zhu to below 1e-12. Both gamma functions decrease with
compression. Fei 2007 retains its printed variable-exponent theta expression.
For gamma=gamma0(V/V0)^q, a positive q below 1 can still make gamma/V and
the high-temperature phonon pressure slope increase with compression; that
does not imply that gamma itself increases.

## Fei measurement comparison

The RT comparison uses all 36 Dewaele (2004) measurements with the revised
ruby pressures. The hot comparison uses all 35 hot Fei (2004) paired Au/Pt
measurements at 1473, 1673 and 1873 K, spanning approximately 12.70-29.29 GPa
after Fei (2007) Au recalibration. Au pressures are recalculated independently
from measured Au volumes; Pt predictions never set the target pressure.
Original Fei (2004) pressure columns are compared separately.

| Model | RT RMS residual (GPa) | Hot RMS, 2007 Au recalibration (GPa) | Hot RMS, original 2004 pressures (GPa) |
|---|---:|---:|---:|
| Fei 2004 | 1.202 | 0.795 | 0.468 |
| Fei 2007 | 0.262 | 0.588 | 0.474 |
| Matsui 2009 | 0.321 | 0.666 | 0.696 |
| Zhu v3 | 0.640 | 0.721 | 0.540 |

These are descriptive, equally weighted pressure residuals. They are not
uncertainty-normalized goodness-of-fit statistics or proof of absolute pressure
accuracy. The calibration changes the comparison; neither Matsui nor Zhu
improves on Fei 2007 for the recalibrated hot subset. Matsui has a near-zero
mean hot residual (+0.029 GPa), but larger scatter. A diagnostic replacing
each cold curve with Fei's fixed 2007 RT curve gives hot RMS values 0.643 GPa
for Matsui and 0.750 GPa for Zhu, versus Fei's 0.588 GPa. Those diagnostic
hybrids are not published EOSs.

## Fei Figure 2 and thermal differences

The preserved raster samples of Fei Figure 2 are calculated curves, not
additional observations. At volumes around 48-55 A^3, RMS discrepancies for
the 1873-minus-1473 K curve spacing are:

| Model | Thermal-spacing RMS discrepancy (GPa) |
|---|---:|
| Fei 2004 | 0.051 |
| Fei 2007 | 0.032 |
| Matsui 2009 | 0.093 |
| Zhu v3 | 0.353 |

Raster spacing has readout sensitivity; small differences do not establish
exact numerical parity. Zhu predicts a lower thermal increment here than
Fei. Cold-curve differences partially compensate in the total pressure.
At V=50 A^3 and T=1873 K:

| Model | Total pressure (GPa) | P(1873)-P(300) (GPa) |
|---|---:|---:|
| Fei 2007 | 96.140 | 12.750 |
| Matsui 2009 | 95.107 | 12.015 |
| Zhu v3 | 95.659 | 11.026 |

This state is outside the pressure range of the recovered hot measurements;
it compares model predictions, not measured validation at 96 GPa.

Matsui's paper states that its pressures differ from Fei 2007 by less than
2 GPa up to 250 GPa and 3000 K. A dense volume sweep at its Figure 8
temperatures 300, 1000, 2000 and 3000 K reproduces that statement: maxima
are 1.552, 0.730, 1.278 and 1.909 GPa, respectively, where Matsui pressure
is between 0 and 250 GPa. This sampled check is not a continuous-domain proof.

The result supports the existing equation implementations but does not
reproduce Fei's original fitting procedure, establish either alternative as
the universally best Pt EOS, or change Fei's `parity_not_achieved` disposition.

## Reproduction and sources

### Isotherm offsets at 100-150 GPa

An additional calculation inverts Fei 2007 separately at each pressure and
temperature, then evaluates Matsui and Zhu at the same volume and temperature.
The following values are alternative minus Fei pressure, in GPa; they include
both cold and thermal differences, with no refitting:

| T (K) | Matsui at 100 GPa | Matsui at 150 GPa | Zhu v3 at 100 GPa | Zhu v3 at 150 GPa |
|---|---:|---:|---:|---:|
| 300 | -0.219 | +0.188 | +1.665 | +3.158 |
| 1000 | -0.729 | -0.557 | +0.665 | +2.000 |
| 1873 | -1.061 | -1.214 | -0.423 | +0.670 |
| 2000 | -1.076 | -1.278 | -0.564 | +0.490 |
| 3000 | -1.103 | -1.705 | -1.512 | -0.795 |

Zhu's higher cold pressure partly cancels its smaller thermal pressure as
temperature increases. Consequently, differences in total pressure can be
small even when thermal increments differ. At equal pressure and temperature,
the checkpoint volume differences are at most 0.21% for Matsui and 0.35% for
Zhu. These hot predictions extend beyond the recovered Fei measurement range.

Run `PYTHONPATH=. .venv/bin/python scripts/plot_pt_isotherm_differences_100_150gpa.py`.
All inverted states are checked to reproduce their target pressure within
1e-9 GPa. The [numerical report](../data/pt-isotherm-differences-100-150gpa.json),
[checkpoint CSV](../data/pt-isotherm-differences-100-150gpa.csv), and
[plot](../data/pt-isotherm-differences-100-150gpa.png) preserve the comparison
convention and source-version qualification. Ruff passes and the plot was
visually inspected.

Run `PYTHONPATH=. .venv/bin/python scripts/compare_pt_thermal_eos_with_fei.py`.
The report preserves input hashes, per-model residual summaries, Figure 2
checks and example pressures. No source data or EOS coefficients are modified.

- [Numerical report](../data/pt-thermal-eos-fei-crosscheck.json)
- [All 78 comparison rows](../data/pt-thermal-eos-fei-crosscheck.csv)
- [Comparison plot](../data/pt-thermal-eos-fei-crosscheck.png)
- [Matsui 2009 primary paper](https://doi.org/10.1063/1.3054331)
- [Fei 2007 primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890468/)
- [Zhu released-v3 dataset](https://doi.org/10.17632/6kxnhc2g73.3)
- [Matsui source audit](matsui-2009-platinum.md)
- [Zhu source/version audit](zhu-2025-pressure-standards.md)
- [Fei fitting audit](fei-2007-platinum.md)

Validation: Matsui primary-table audit passes; Zhu Pt released thermal refit
passes; 18 existing Matsui/Zhu/Fei tests pass; Ruff passes for the comparison
script. The plotted PNG was visually inspected.
