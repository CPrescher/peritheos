# Fei et al. (2007): Pt thermal pressure scale and joint-fit audit

## Investigation disposition

On 2026-10-09 the user paused further Fei Pt investigation pending a reply
from Fei to the drafted fitting-method inquiry. Published Pt EOS coefficients,
defaults and existing validation distinctions remain unchanged. Do not continue
Pt refitting or attempt to resolve the discrepancy until the author response
is available and the user resumes the work. Other Fei 2007 materials may use
the published Pt calibration unchanged where required by their source.

## Accepted published equation

The canonical record `platinum_fei_2007_vinet_mgd` implements the Pt Vinet
and Mie–Grüneisen–Debye scale from [Fei et al. (2007), *Toward an internally
consistent pressure scale*](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890468/),
Table 1 and Equations (2)–(3). The primary full text was inspected on
2026-10-08. This is the published parameterization; the diagnostic fits below
do not replace its coefficients. The existing
`platinum_fei_2007_vinet_300k` remains a separate reference-isotherm record.
The compatibility identifier `pt_fcc_fei_2007` now resolves to the canonical
thermal record and preserves the legacy numerical pressure calculation.

Volumes are conventional fcc-cell Å³ with four Pt atoms per cell. The thermal
energy is per mole of atoms, so cell volume is converted to molar volume using
`V_m = V_cell N_A / 4`. The reference temperature is 300 K.

| Parameter | Published value | Printed error width |
|---|---:|---:|
| V0 (Å³/cell) | 60.38 | 0.01 |
| K0 (GPa) | 277 | not reported |
| K0′ | 5.08 | 0.02 |
| θ0 (K) | 230 | not reported |
| γ0 | 2.72 | 0.03 |
| q | 0.5 | 0.5 |

The Pt value `q=0.5(5)` is identical to Fei (2004), Section 3.2 and Table 3,
where it appears among optimized Pt parameters. Fei (2007) describes fitted
model parameters but does not explicitly enumerate q as free in its updated
optimization. A newly optimized versus retained q therefore remains uncertain;
its one-decimal rounding and reported error do not establish that it was fixed.
Our free-q diagnostic is not a verified replay of the author's parameter choices.

The pressure is `P(V,T) = P_Vinet(V) + γ(V)[E_D(V,T) − E_D(V,300 K)]/V_m`.
The source prints `γ(V) = γ0 (V/V0)^q` and
`θ_D(V) = θ0 (V/V0)^(−γ(V))`. This last expression is the direct
variable-exponent law already supported by Peritheos; it must not be replaced
by the integrated constant-q law. This is an implementation of the printed
pressure scale, without a claim that every caloric derivative follows a
thermodynamically consistent free-energy model. Missing confidence levels and
parameter covariance remain null. The adopted 277 GPa bulk modulus is retained
as fixed in the record; the complete original optimization constraints are
not supplied by the article.

The recovered inputs span 0–93.6 GPa and 300–1873 K. These are marginal
data envelopes: the high-temperature measurements cover substantially lower
pressures than the cold-compression series. They do not establish a fully
measured rectangle extending to 93.6 GPa at 1873 K.

## Recovered measurements and Au recalculation

The Pt paragraph states that the fit combines the room-temperature data of
Dewaele et al. (2004) with the P–V–T data of Fei et al. (2004), recalculated
using this study's Au scale. The audit therefore combines:

- All 36 rows in `platinum-dewaele-2004-table1-compression.csv`, using the
  revised ruby pressures and four times the atomic volume. The original
  298 K source provenance is preserved; this diagnostic assigns the series
  to the paper's 300 K reference branch.
- All 42 paired Au–Pt rows in `platinum-fei-2004-table2.csv`, including
  300, 1473, 1673 and 1873 K. Pressures are independently recalculated from
  the measured Au volumes using Fei (2007)'s Au Vinet-MGD coefficients:
  V0=67.85 Å³, K0=167 GPa, K0′=6.0, θ0=170 K, γ0=2.97 and q=0.6.

Pt's own predicted pressures never define the fitting target. The original
2004 pressure columns and all source measurements remain unchanged. The
42 derived Au pressures, their original pressure counterparts, and SHA-256
hashes of both source CSVs are retained in
[`fei-2007-platinum-reproduction.json`](../data/fei-2007-platinum-reproduction.json).
These are derived calibration values, not newly measured pressures. Source
error columns are preserved, but the diagnostic does not claim a propagation
of their uncertainties through the revised Au scale.

## Equation parity and qualified fit discrepancy

Run `.venv/bin/python -m scripts.reproduce_fei_2007_platinum --check`.
The independent implementation uses explicit Vinet and Debye-integral
expressions with SI R and N_A, rather than the Peritheos EOS evaluator.
The Pt equation agrees over a nine-state grid spanning V/V0=0.8–1.0 and
T=300–1873 K to below 10⁻⁹ GPa. Independently recalculated pressures for
all 42 Au observations also agree with the canonical Au evaluator below
10⁻⁹ GPa. Pressure-to-volume and pressure/volume-to-temperature inversion,
catalog lookup, alias resolution and `.eosmat` reconstruction are tested.

An additional comparison uses the author's calculated 300, 1473 and 1873 K
curves in Figure 2, independently of our implementation-to-implementation
checks, without refitting any coefficients. At 23 sampled raster positions
(cell volumes approximately 48–55 Å³), the pressure residual RMS values are
0.369, 0.178 and 0.180 GPa, respectively. Maximum absolute differences are
0.491, 0.235 and 0.233 GPa. Calculated pressures are slightly lower than the
sampled curve centers; the 300 K offset is systematic and exact numerical
parity is not established. An earlier illustrative three-pixel sensitivity
check alone is not sufficient evidence for exact parity.

The spacing between the 1873 and 1473 K curves at equal volume removes the
cold pressure contribution. With the published thermal coefficients unchanged,
the spacing residual RMS is 0.032 GPa and maximum absolute difference is
0.073 GPa, smaller than one horizontal source pixel (0.089 GPa). This supports
the published thermal coefficients as an interpretation of the drawn curves.
It does not distinguish the printed and integrated Debye-temperature laws,
whose pressure differences here are smaller than the raster readout scale.
The source image, axis calibration, sampled pixels and computed pressures
are preserved in
[`fei-2007-platinum-source-curve-check.json`](../data/fei-2007-platinum-source-curve-check.json).
The overlay and residual panels are generated with
`python -m scripts.plot_fei_2007_figure2_validation` and saved in
[`fei-2007-platinum-figure2-validation.png`](../data/fei-2007-platinum-figure2-validation.png).
The error bars combine observed horizontal half-stroke width with an
illustrative one pixel per calibrated coordinate; they are deterministic
readout sensitivities, not statistical uncertainties. These calculated curves
are not additional experimental observations or fitting targets.

The diagnostic coefficients below substantially depart from the drawn thermal
branches at compression beyond the recovered thermal inputs. On the same
Figure 2 samples, their RMS pressure differences are 3.436 GPa at 1473 K and
4.701 GPa at 1873 K. At V=50 Å³ and 1873 K, the published and diagnostic
pressures are 96.139 and 101.741 GPa, respectively. The 1873-minus-1473 K
spacing is 3.273 GPa with the published coefficients and 4.719 GPa with the
diagnostic coefficients. These comparisons use unchanged diagnostic
coefficients, not coefficients fitted to the figure.

For context, a separate Fei (2004) diagnostic using its original 42 Au–Pt
rows and pressure targets gives q=-0.977 with K0, K0′, γ0 and q free and
V0/θ0 fixed. The 2004 published Pt coefficients reproduce its Table 2 Pt
calculated-pressure column to within 0.194 GPa; that output column is not
the fitting target. In three manually reviewed raster rows of Figure 4 per
isotherm, published-coefficient pressure differences are at most 0.224 GPa.
Both the published and diagnostic curves remain visually close over Figure
4's limited pressure range, despite different fitted coefficients. This
comparison does not recover the author's fitting protocol or establish
exact rounded-coefficient/table-output parity. The two source figures,
published calculations and diagnostic curves are compared in
[`fei-2004-2007-platinum-refit-curve-comparison.png`](../data/fei-2004-2007-platinum-refit-curve-comparison.png);
the numerical results and source provenance are retained in the
[comparison report](../data/fei-2004-2007-platinum-refit-curve-comparison.json).

### Fei 2004 fitting-weight sensitivity

A subsequent diagnostic uses only the original 42 Table 2 observations and
their original pressure calibration, without Dewaele data or Fei 2007 Au
recalibration. The seven RT rows alone give K0=294.352 GPa and K0'=2.105,
consistent with the reported RT-only K0=290(10) GPa and K0'=2.7(9) within
the printed error widths, without exact coefficient parity. The full unweighted
joint fit gives K0=289.739 GPa, K0'=2.636, gamma0=2.55940 and q=-0.97720,
with RMS 0.37436 GPa versus 0.44257 GPa for the published full-model coefficients.

The sign of q is sensitive to weighting and constraints. Pressure-equivalent
weights derived from the printed Pt volume errors, evaluated at the published
EOS coefficients and then held fixed, give a joint fit with q=0.15436.
Its other coefficients do not reproduce the full published parameter set.
With the published K0=273 GPa and K0'=4.8 fixed, a fit of gamma0 and q to
only the 35 hot rows, weighted by combined Au/Pt volume-error contributions,
gives gamma0=2.68342 and q=0.23771. Both thermal values lie within the
published 2.69(3) and 0.5(5) error widths. The fixed cold coefficients are
inputs, not independently reproduced results. Original source confidence
levels and covariance remain unspecified.

These weights approximate pressure-equivalent errors by symmetric differences
at plus/minus the reported volume errors and assume independent Au/Pt errors.
They omit temperature uncertainties, calibration-parameter uncertainties and
correlated systematic errors. They are explicit sensitivity choices, not
recovered author weights. Three q starts converge to the same objective in
each variant. This does not recover the original optimization but shows that
negative q is not an unavoidable result of fitting the 2004 data. An email
should describe the negative q from the unweighted reconstruction, rather
than imply that every reasonable fit gives negative q.

Run `PYTHONPATH=. .venv/bin/python scripts/audit_fei_2004_platinum_fit_weights.py`.
All seven constraint/weight variants and the input hash are saved in
[`fei-2004-platinum-fit-weight-sensitivity.json`](../data/fei-2004-platinum-fit-weight-sensitivity.json).

An equal-pressure-weight joint least-squares diagnostic fixes V0, K0 and θ0,
and fits K0′, γ0 and q. Those constraints and weights are explicit audit
choices, not a recovered author objective. The article does not specify the
complete fitted/fixed parameter list, weights, covariance or unrounded input
values. The recovered-table diagnostic gives:

| Parameter | Published | Diagnostic |
|---|---:|---:|
| K0′ | 5.08 | 5.086960 |
| γ0 | 2.72 | 2.612403 |
| q | 0.5 | −1.660642 |

The RMS pressure residual decreases from 0.442161 GPa at the published
coefficients to 0.347249 GPa. Freeing V0 gives 60.375849 Å³ and q=−1.665203;
excluding the seven 1673 K rows, which are not shown in Figure 2, gives
q=−1.570512 on 71 observations. These sensitivities do not recover the printed
positive q. Their conditional regression standard errors are separate from
the source's unspecified confidence convention; no combined-2σ parity claim
is made.

### Nonnegative-q diagnostic

`python -m scripts.fit_fei_platinum_q_bounds` repeats the same pressure
objective with q constrained to be nonnegative. The 2007 solution reaches
the lower boundary, q≈0, with K0′=5.089434, γ0=2.736640 and RMS residual
0.389014 GPa, versus 0.347249 GPa unconstrained and 0.442161 GPa at all
published coefficients. Fixing q=0.5 while refitting K0′/γ0 gives
γ0=2.774315 and RMS residual 0.415483 GPa. Requiring q≥10⁻⁶ instead gives
q=10⁻⁶ and essentially the same residual as the nonnegative fit; it does
not select an interior positive optimum. At q=0, gamma is volume-independent.

The 2004 diagnostic similarly reaches q≈0, with K0=293.504575 GPa,
K0′=2.563001, γ0=2.606256 and RMS residual 0.390010 GPa, versus
0.374360 GPa unconstrained and 0.442570 GPa at its published coefficients.
The different cold/thermal constraints of the two diagnostics remain
explicit; the boundary does not establish the source's optimization protocol.
Starts at q=0.05, 0.5 and 2 converge to the same constrained residual, and
fixing q=10⁻⁴ slightly increases the residual in each diagnostic. Symmetric
regression errors are not assigned to these boundary solutions. The 2004
independent theta calculation uses an algebraically equivalent `expm1`
expression to avoid cancellation near q=0, with pressure parity against the
native published evaluator still better than 10⁻⁸ GPa.

Numerical results and input hashes are retained in
[`fei-platinum-q-bounded-diagnostics.json`](../data/fei-platinum-q-bounded-diagnostics.json).
Published catalog coefficients remain unchanged.

### Thermal-only fit with the RT EOS frozen

`python -m scripts.fit_fei_2007_platinum_thermal_subset` first fits K0′ to
all 43 RT rows (36 Dewaele plus seven Fei), with V0=60.38 Å³ and K0=277 GPa
fixed. It then freezes K0′=5.086992 along with V0/K0/θ0 and fits only γ0/q
to the 35 Fei (2004) hot rows: six at 1473 K, seven at 1673 K and 22 at
1873 K. These thermal measurements span 12.70–29.29 GPa after recalculation
using the published Fei (2007) Au scale. No RT residual enters the thermal
optimization and no cold parameter changes during it.

| Thermal constraint | γ0 | q | Hot-only RMS (GPa) |
|---|---:|---:|---:|
| Unconstrained q | 2.612404 | −1.660616 | 0.424133 |
| q≥0 | 2.736730 | 0 (boundary) | 0.498497 |
| q fixed at 0.5 | 2.774434 | 0.5 | 0.544075 |

These hot-only RMS values must not be compared directly with the preceding
78-row combined RMS values. The thermal inputs in the earlier joint fit were
already exclusively Fei (2004); Dewaele supplied only the RT branch. Freezing
the RT coefficients therefore leaves the negative-q discrepancy essentially
unchanged. Sensitivities freezing the Dewaele-only RT fit or the exact
published cold coefficients also give negative unconstrained q (−1.661484
and −1.666343, respectively) and boundary q for the nonnegative fit.
Different q starts converge to the same residuals, and all 43 RT pressures
are verified unchanged by every thermal fit. Results, hot run IDs, derived
Au pressures and input hashes are retained in
[`fei-2007-platinum-fixed-rt-thermal-subset.json`](../data/fei-2007-platinum-fixed-rt-thermal-subset.json).

### RT curve incorrectly treated as a 0 K baseline

`python -m scripts.audit_fei_2007_platinum_reference_temperature` tests
omitting the reference-energy subtraction while keeping the RT cold
coefficients frozen: `P=P_RT+gamma*E_D(T)/V_m`, instead of
`P=P_RT+gamma*(E_D(T)-E_D(300))/V_m`. Debye energy here excludes zero-point
energy. On the same 35 hot rows, the incorrect Pt-only reference gives
γ0=2.294911, q=−1.351253 and RMS residual 0.426924 GPa. If the same omission
is also made in the hot Au pressure recalculation, the Pt fit gives
γ0=2.644439, q=−1.360114 and RMS residual 0.425024 GPa. Both nonnegative-q
variants again reach zero. Multiple q starts give the same residuals.

This specific reference error therefore does not recover the published
positive q. With the published coefficients, it adds 1.660662 GPa at
V=V0 and 300 K. Its Figure 2 sampled-curve RMS differences become 1.194,
1.367 and 1.362 GPa at 300, 1473 and 1873 K, versus 0.369, 0.178 and
0.180 GPa with the 300 K subtraction retained. This is evidence against
the omission as an explanation of the published calculation, without a
claim about the author's unavailable code. The combined Au/Pt scenario
changes hot calibration targets while retaining the earlier frozen RT
coefficients; it is not a full recalibration of all RT observations.

A consistent conversion of the baseline to 0 K first subtracts the 300 K
vibrational pressure from the RT curve. Adding the full vibrational pressure
then reproduces the original referenced calculation, verified to below
10⁻¹⁰ GPa. The reference-temperature choice alone is therefore not an error;
inconsistent handling of the baseline is. Counterfactual results, baseline
shifts and source-input hashes are saved in
[`fei-2007-platinum-reference-temperature-sensitivity.json`](../data/fei-2007-platinum-reference-temperature-sensitivity.json).

The record's fit-ledger status is therefore `parity_not_achieved`, while its
equation implementation remains `primary_source_validated`. Equation parity
and reproduction of the author's fit are distinct results. The retained
disposition is that the recovered Pt measurements do not validate the
published thermal fit in our reconstruction; the source-curve checks do not
change that measurement-fit status. The paper clearly
identifies the physical model, input publications and pressure-calibration
sequence; this diagnostic does not establish that its description is
insufficient to reproduce the fit. The discrepancy remains unresolved in our
reconstruction. Further investigation should follow the stated calibration
and fitting sequence and check the recovered input and uncertainty
conventions. Original weights, constraints and unrounded data would help
resolve remaining ambiguities, but their necessity has not been demonstrated.
