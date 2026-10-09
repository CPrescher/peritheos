# Au PVT validation: published equations and derived reconstructions

The pressure-only audit is reproduced with
`.venv/bin/python -m scripts.audit_gold_pvt --check`. Its numerical results,
input hashes and comparison points are archived in
[gold-pvt-validation.json](../data/gold-pvt-validation.json).
Entropy, heat capacity and electronic energy are not acceptance criteria
for this PVT audit.

## What is validated

| Model | PVT outcome | Evidence limit |
|---|---|---|
| `gold_fei_2007_vinet_2` | Published complete PVT equation verified | Original optimization and covariance are not recovered |
| Yokoo (2009) published BM3 and Vinet 300 K records | Published reference-isotherm equations verified | These fitted branches are distinct from the full thermal equation |
| `gold_yokoo_2009_pvt_reconstruction` | Numerical Table III reconstruction checked; scientific status `not_reproduced` | Explicit-selection diagnostic with altered phonon coefficients; published analytical PVT remains unreproduced |
| Yokoo (2009) full equation with printed coefficients | Not reproduced | Cold Vc is unreported and printed phonon normalization leaves a pressure mismatch |

These outcomes describe different evidence. Equation implementation,
graphical source agreement, calibration-dependent measured-data comparison,
and recovery of an author's original optimization are assessed separately.

## Fei (2007): complete published Au PVT

The primary PDF of [Fei et al. (2007)](https://doi.org/10.1073/pnas.0609013104)
was inspected at Table 1, Equations 2–3 and Figure 1 on page 9183, and the
Au thermal-refit paragraph continuing onto page 9184.

| Parameter | Value | Printed error width |
|---|---:|---:|
| V0 (Å³/four-atom cell) | 67.850 | 0.004 |
| K0 (GPa) | 167 | not reported; adopted |
| K0′ (Vinet) | 6.00 | 0.02 |
| θ0 (K) | 170 | not reported |
| γ0 | 2.97 | 0.03 |
| q | 0.6 | 0.3 |

The thermal error widths are now retained in the catalog. Confidence level
and covariance remain unspecified. The pressure equation is

`P = P_Vinet + gamma(V) [E_D(V,T) − E_D(V,300 K)] / Vm`,

with `gamma(V)=2.97 (V/V0)^0.6` and the **printed direct exponent**
`theta(V)=170 (V/V0)^(-gamma(V))`. This source choice differs from the
integrated constant-q Debye-temperature law. The independent audit uses
adaptive Debye quadrature, SI R, and `Vm=Vcell N_A 10^-30/4` in m³/mol Au
atoms. The electronic pressure table is not added to this published Fei
equation, which has no separate electronic term.

The independent pressure replay checks 119 states at 17 volumes and seven
temperatures. Pressure-to-volume inversion, temperature inversion and
`.eosmat` reconstruction are tested. This numerical test domain does not
assert a measured validity rectangle.

Figure 1 is also checked independently of implementation parity. The
source PDF contains vector paths for the three calculated isotherms.
[fei-2007-gold-source-curves.json](../data/fei-2007-gold-source-curves.json)
retains the PDF checksum, axis tick positions, colors, line width and all
37 vertices of each path. Linear calibration uses the printed tick
coordinates, without fitting EOS coefficients. The pressure RMS/max
differences are:

| Isotherm (K) | RMS (GPa) | Maximum absolute difference (GPa) |
|---:|---:|---:|
| 300 | 0.0983 | 0.2186 |
| 1473 | 0.0314 | 0.0865 |
| 2173 | 0.0310 | 0.0863 |

Every point is within one full printed line width in the pressure direction
(0.278 GPa). This is agreement at graphical precision, not exact numerical
parity. Line width is a graphical comparison scale, not an experimental
uncertainty or confidence interval. These are calculated curves, not extra
observations.

Measured-data comparisons use all 26 Au–MgO rows in Fei (2004) Table 1
and all 37 Au rows in Dewaele (2004) Table I. Fei Table 1 was visually
checked on page 519, including its MgO-calibration footnote. The original
source values and uncertainty columns remain unchanged:

| Comparison | Source pressure calibration | Pressure RMS/max difference (GPa) |
|---|---|---|
| 26 hot Fei (2004) rows | Speziale (2001) MgO | 0.4359 / 0.8625 |
| 37 cold Dewaele (2004) rows | Revised ruby scale | 0.3185 / 0.9059 |

The cold source series was measured at 298 K; this comparison explicitly
evaluates the paper's 300 K reference isotherm. The hot rows cover
8.64–25.56 GPa and 1273–2173 K; the cold rows reach 93.6 GPa. These are
marginal envelopes, not complete hot coverage to 93.6 GPa. They are source
fit/comparison data, not independently withheld observations. The paper's
six digitized new cold observations remain explicitly plot-only. Original
fit selections, weights, unrounded inputs and covariance are not recovered,
so no claim of exact author-fit reproduction is made.

## Yokoo (2009): qualified pressure reconstruction

The [Yokoo Au audit](yokoo-2009-gold.md) and
[thermal reconstruction](../data/yokoo-2009-gold-thermal-reconstruction.json)
retain the published coefficients separately from the diagnostic refit.
The [PVT library audit](yokoo-2009-pvt-library.md) documents the registered
Au and Pt thermal reconstruction records, bounded inversions and interchange.
Table III on page 104114-4 was visually checked, including the six
parenthesized first-liquid states and six omitted cells.

The derived model fits 156 unmarked calculated states with RMS **0.008923
GPa** and maximum **0.021133 GPa**. The 74 states withheld on alternate
isochores give RMS **0.009076 GPa**, maximum **0.018672 GPa**. All six
first-liquid markers are held out of the primary fit and assessed separately.
An adaptive SI quadrature check agrees with the evaluator to below 10^-10
GPa. Multiple-start convergence, domain bounds and corrected electronic
nodes are covered by the reconstruction tests.

This verifies the executable **derived pressure reconstruction** under its
documented assumptions. It does not validate the reconstructed coefficients
as independently measured material parameters. In particular, gamma0 changes
from 2.96 to 2.923265 (−1.24%), beyond printed decimal rounding, and Vc is
inferred from model outputs. The residual maximum exceeds the table's 0.005
GPa half-step. The cause of the discrepancy with printed coefficients remains
unestablished, so the full published Yokoo PVT model is not marked reproduced.
The [pressure-convention audit](yokoo-2009-pressure-conventions.md) isolates
the mismatch using thermal increments, where cold-curve and constant
reference offsets cancel, and tests explicit last-digit rounding intervals.

The supplied Tsuchiya–Kawamura electronic **pressure** table is sufficient
for the bounded pressure reconstruction, with explicitly chosen linear
interpolation and volume independence. Missing electronic energy concerns
the original shock-temperature reduction and caloric calculations; it is
not a requirement for evaluating P(V,T). The represented 0.6–1 volume ratios
and 0–3000 K rectangle do not establish solid-phase validity beyond the
source phase annotations.
