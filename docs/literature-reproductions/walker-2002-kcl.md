# Walker (2002): Figure 1 reference and recovered Birch calibration

Audit date: 2026-10-08. The selected B1 record now uses **Figure 1's
V0=37.50 cm3/mol**, converted to 249.080860076 A3 per conventional Z=4 cell.
Tables 1/3 give the conflicting 249.53 A3 (37.57 cm3/mol); both raw CSVs,
including the derived Table 1 anchor, remain unchanged. Conversion digits
are not additional measurement precision. K0=17.7 GPa, K0'=5,
beta=0.00195 +/- 0.00005 GPa/K and Tr=296.15 K remain the published inputs.
Defaults and the Walker B2 reference are unchanged.

Run `PYTHONPATH=. python scripts/audit_walker_2002.py --check` and
`PYTHONPATH=. python scripts/reproduce_kcl_variants.py --check`.
The [numerical report](../data/walker-2002-reproduction.json) records the
independent fits, all 69 bundled NaCl input rows, their actual temperatures,
source hashes, Birch benchmarks and conditional B1/B2 pressure replay.

## Primary evidence

The final [Walker article](https://doi.org/10.2138/am-2002-0701), pages 806-810,
was inspected from the locally held original PDF. SHA256:
`01585fbd045709f3c08516aeedb30dcb553aec1d3073542c5450a72a6c3ff4b6`.
Page 806 attributes pressure for both KCl phases to the Birch (1986) NaCl-B1
BE2 thermal EOS. Page 808 states an unweighted sum of squared pressure
residuals minimized in Excel Solver; B1 fixes V0 and K0'=5. The Table 3
cold-first staging footnote belongs to B2. Figure 1 prints B1 V0=37.50 cm3/mol.
Table 1's B2-KCl volume heading is erroneous: the volumes describe B1, Z=4.

The supplied original [Birch (1986) article](https://doi.org/10.1029/JB091iB05p04949)
was inspected, including rendered pages 4951-4952. PDF SHA256:
`3b5b75419a3ceb70fc46e7c676b77cca2ad43dddb0449a4a5e3f2f9bc28c8e3a`.
Table 6's final extended 25 Celsius row supplies K0=238.8 kbar,
a=1.796 and b=-5.00. Equation 8 supplies the thermal term:

```text
r = V/V00, with V00 the zero-pressure 25 Celsius reference
f = (r^(-2/3)-1)/2
P(GPa) = 3*23.88*f*(1+2*f)^(5/2)*(1+1.796*f-5.00*f^2)
         + 0.00286*(T-298.15)
```

BE2 is a quadratic strain polynomial, not BM2. The printed polynomial
coefficients are used directly: reconstructing a from rounded K0'=5.20
would give 1.80 and degrade agreement. All eight Table 5 reference-isotherm
benchmarks agree within half a printed pressure digit (0.0005 GPa).
Temperature-specific Table 6 polynomials approximate the additive Equation 8;
they are not substituted for it. The construction covers 25-500 Celsius and
0-30 GPa. The old missing-coefficients assessment is superseded.

## Conditional B2 NaCl replay

Walker Table 2 supplies separate ambient NaCl anchors: Run 1 r34439,
a=5.6414(14) A at 23 Celsius; Run 2 r35101, a=5.6468(4) A at 24 Celsius.
The first seven bundled observations belong to Run 1; the remaining 32 to
Run 2. The audit assumes each anchor represents zero pressure at its actual
temperature. It first solves Birch's equation for that ambient ratio r0,
then evaluates the observed ratio `(a/a_anchor)^3*r0` at the measured
Celsius +273.15 temperature. No parameter is fit to Walker pressures.

All 39 B2 pressures agree with the reported values within **0.01093 GPa**.
This supports the run-normalization interpretation, but the author spreadsheet,
unrounded inputs, exact anchor corrections and uncertainty propagation remain
unavailable. Ambient 23/24 Celsius and heated 600 Celsius states are flagged
as outside Birch's construction range. The conditional replay extrapolates
Equation 8 to them; it does not recover the author's extrapolation protocol.
No ready observation-pressure reduction or executable catalog calibration
edge is registered; the catalog's
`reference_eos_not_bundled` status means that edge is absent, not that Birch's
numerical coefficients are still missing.

## Conditional B1 NaCl replay

Table 1 contains separate ambient NaCl measurements at **36 Celsius**:
r57689, a=5.6479(4) A in the mixed KCl-NaCl pellet, and r57693,
a=5.6473(4) A in the pure-NaCl spot check. Page 806 describes these two
measurement families. The audit assigns each sample or spot-check row to
its corresponding reference and retains every measured row temperature.

An explicit **23 Celsius normalization hypothesis** treats these two lattices
as zero-pressure anchors at 296.15 K despite their measured 309.15 K.
It uses the same Birch ratio normalization as B2, without fitting coefficients
to Walker pressures. For nonzero pressures the agreement is:

| NaCl family | Rows | RMSE (GPa) | Maximum difference (GPa) |
|---|---:|---:|---:|
| Mixed-pellet sample | 22 | 0.00069980 | 0.00145779 |
| Pure-NaCl spot check | 5 | 0.00063829 | 0.00103240 |

This agreement does not establish the author's temperature convention.
The two bracketed ambient pressures are preserved as imposed zeros:
evaluating their actual 36 Celsius temperatures under the hypothesis gives
**0.03718 GPa**, rather than zero. They are displayed as diagnostics and
excluded from the 27 nonzero-pressure comparison. The derived 23 Celsius
reference has no reported pressure and is not replayed.

The report also evaluates both anchors at their printed 36 Celsius
temperatures. Nonzero-pressure RMSE then rises to 0.03906 GPa for samples
and 0.03764 GPa for spot checks (maximum differences 0.04693 and
0.04200 GPa). This sensitivity makes the unresolved reference-temperature
inconsistency visible. The conditional replay is retained as the available
reproduction, with the reference-temperature convention recorded as an evidence
limit rather than further numerical work. An exact author reduction would
require the original spreadsheet or confirmation of that convention. The audit
does not alter raw observations, published KCl coefficients, source
error widths, Figure 1 V0, or catalog calibration availability.

## B1 fit with the selected Figure reference

The 23 measured KCl states are fit with fixed Figure V0, K0'=5 and
Tr=296.15 K. Six NaCl-only spot checks and the derived reference anchor
are excluded from fit residuals. With beta=alpha0*K0, the stated objective
is linear in K0 and beta. Independent linear and nonlinear solves agree.

| Quantity | Published | Joint Figure-reference fit |
|---|---:|---:|
| K0 (GPa) | 17.7 | 17.68342185 |
| beta (GPa/K) | 0.00195 +/- 0.00005 | 0.001927145871 |
| alpha0 (1/K) | 0.00011, rounded | 0.000108980371 |
| Pressure RMSE (GPa) | 0.05443363 | 0.05382736 |

K0 and alpha0 round to the published values; beta differs by 1.17%, within
its printed error width. This is qualified numerical similarity, not exact
solver/covariance reproduction. The conflicting table reference remains a
separate diagnostic (joint RMSE 0.06063608 GPa). The compression-positive
BM3 convention is explicit because the printed BE1 signs disagree with the
printed positive strain definition. No unreported elastic errors are assigned.

## B2 and Ma boundary

The independent B2 cold/thermal staged fit remains unchanged: eight 23/24
Celsius rows give K0=23.77539804 GPa, K0'=4.41587367; all 39 thermal rows
give beta=0.002766245675 GPa/K. The actual 297.15 K of the eighth cold row
is retained. The seven-23-Celsius-row sensitivity also remains visible.

The [Ma audit](ma-2024-kcl.md) applies the same source-backed run-normalization
principle to Matsui's NaCl EOS, then brings KCl to 300 K with Walker's
B2 thermal coefficient. All eight Ma pressures agree within 0.000377 GPa,
below printed-input precision. This is separate from the B1 Figure choice.
Ma's deposited reductions, joint-fit observations and published parameters
are unchanged. Exact algorithm, source weighting and covariance remain open.

Tests cover Birch's printed benchmarks, actual and assumed anchor temperatures,
separate B1 reference families, imposed-zero diagnostics, both fit
objectives, raw-source preservation, catalog round trips and generated-report
freshness. Article PDFs are not redistributed. Contributor CC0 applies to
transcription, normalization and metadata, not third-party article rights.
