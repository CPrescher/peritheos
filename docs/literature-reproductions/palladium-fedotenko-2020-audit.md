# Fedotenko (2020) Pd: corrected weights and poor fixed-K' data agreement

Audited 23 October 2026; classification clarified 24 October 2026.
**Both Pd parameterizations achieve coefficient parity with the paper under
the documented diagnostic assumptions.** The fixed-K'=4 curve remains a poor
fit to Table 1; that finding is recorded separately as `data_agreement_status:
poor`, not as failed reproduction. It is retained for explicit diagnostic
selection and is not recommended for quantitative use or pressure calibration.
Its published coefficients and errors are unchanged. Matching the published
fit reproduces its poor data agreement too; it does not establish exact author
workflow or imply that the authors claimed excellent BM2 data agreement.

**Do not use the fixed-K'=4 entry for quantitative predictions or pressure
calibration.** Its label explicitly says "poor fit; diagnostic only". Ordinary
material loading excludes it; deliberate source inspection and explicit record
selection retain it for comparisons and reproducibility.

## Source and preserved observations

Fedotenko et al., “Synthesis of palladium carbides and palladium hydride in
laser heated diamond anvil cells,” J. Alloys Compd. 844, 156179 (2020),
[DOI](https://doi.org/10.1016/j.jallcom.2020.156179).
The [author's thesis copy](https://epub.uni-bayreuth.de/id/eprint/5724/1/Thesis_Fedotenko.pdf),
Section 5.C, printed pp. 125-127 and 132, supplies methods, Equation 1,
Table 1, Figure 1 and Table 2. All 15 bundled Table 1 observations agree
with that source. Volumes are conventional fcc-cell A^3, with four Pd atoms
per cell; pressure is GPa. The range is 0-31.2 GPa at room temperature.
Neon pressure-calibrant volumes are unavailable, so pressure reduction is not
independently replayed. The source names EoSFit7-GUI, but does not establish all
constraints, confidence, covariance or error-scaling conventions.

Table 1, the main fit paragraph and Figure 1 give V0=58.868(2) A^3;
Table 2 alone gives 58.863(2). The retained value is 58.868. Substitution of
58.863 changes the BM2 pressure RMS only to 0.928243 GPa and does not resolve
the misfit. The separate fixed-K'=5.28 paragraph reports K0=190(2) GPa but
does not give V0 separately; no new catalog record is inferred from it.

## Corrected audit weights

The ambient row has a measured volume 58.868(2) A^3 and no printed pressure
uncertainty. The generic audit reader previously discarded the entire pressure
uncertainty vector because of that blank. Fitting then used unit pressure
weights, losing all 14 reported 0.2 GPa widths and producing a misleading
reduced chi-square of 0.933860 for BM2.

The dedicated Fedotenko adapter bypasses that incomplete-vector fallback:

- All 14 positive-pressure rows use their printed pressure widths of 0.2 GPa
  and volume widths of 0.003 A^3, with exact latent-volume EIV residuals.
- BM3 treats the ambient volume as a separate 58.868 +/-0.002 A^3 observation
  of V0 at reference P=0. All three EOS coefficients remain free.
- BM2 retains V0=58.868 as an explicitly diagnostic fixed constraint, with
  K'=4 intrinsic to BM2. Its ambient row contributes no independent residual
  or degree of freedom. Table 2 does not establish that the author fixed V0.
- No pressure uncertainty is invented for the ambient row. Printed widths
  are conditionally interpreted as independent absolute standard deviations.
  Covariance is inverse J^T J, without residual scaling, after accounting for
  latent volumes. This is not a recovered author uncertainty convention.

This material-specific treatment does not reinterpret incomplete uncertainty
columns or fit protocols for unrelated catalog records.

## Quantitative finding

| Curve or diagnostic reconstruction | Pressure RMS, all 15 rows (GPa) | Reduced chi-square |
| --- | ---: | ---: |
| Published BM3, V0=58.868, K0=157, K'=9.9 | 0.088533 | Not inferred |
| Corrected BM3 reconstruction | 0.075757 | 0.177875 |
| Published BM2, V0=58.868, K0=203, K'=4 | 0.934246 | Not inferred |
| Corrected fixed-V0 BM2 reconstruction | 0.933704 | 25.005483 |

The BM2 reconstruction gives K0=202.659287 +/-0.541006 GPa, agreeing with
the printed 203(3) GPa. Its weighted sum of squares is 325.071284 with 13
degrees of freedom. The published curve misses points by up to 1.438559 GPa,
7.19 times the printed pressure error width. It lies high at low pressures
and low at high pressures: the systematic curvature mismatch survives refitting
K0. The residual ratio and chi-square interpretation are conditional on the
assumed error convention and independence.

The corrected BM3 reconstruction gives V0=58.867988 +/-0.002000 A^3,
K0=157.332722 +/-2.587854 GPa, K'=9.906463 +/-0.424563. It agrees with the
published triplet under this explicit ambient treatment. This does not validate
extrapolation or resolve the unusually high K'. Neither fit claims exact
author-workflow reproduction.

Pressure residuals are P(model)-P(observed) at measured volumes; both RMS
values use all 15 input rows. These differ from latent-coordinate weighted
residuals used in optimization. A separate audit sensitivity allowing V0 to
vary with its measured ambient width still gives BM2 RMS=0.932010 GPa and
reduced chi-square=24.960047; fixing V0 does not explain the poor agreement.
No diagnostic reconstruction replaces the source parameterization.

## Reproduction and checks

The repository-only adapter is
`scripts/reproduce_fedotenko_2020_palladium.py`, called by
`scripts.validate_primary_eos_refits._fit_record` before the generic row reader.
Its regression tests check the printed-error weighting, independent BM
equation values, measured-coordinate RMS, degrees of freedom, preservation of
coefficients/errors and CSV, public explicit selection and round trips.

[Curves, coefficients and residuals](../reconstruction-reports/palladium.md#record-palladium_fedotenko_2020_bm2_2)
and the [numerical ledger](../data/primary-eos-refits.json) retain coefficient
parity separately from the poor-data finding and usage restriction.
