# CaSiO3 validation and qualified use

Qualification checkpoint: 2026-10-08. This is a metadata clarification; no EOS
coefficient, reported error, observation, pressure scale or source classification
is changed.

`scientific_validation.status = primary_source_validated` means that the stated
source representation has been checked against primary evidence. It does not
by itself establish reproduction of the original regression or accuracy as a
pressure standard. The existing extensible validation object now separately
records `reproduction_status`, `fit_reproduction_note` where applicable, and
`usage_recommendation`. These fields qualify executable records without
changing their identities, numerical behavior or published/refit classification.

| Record | Source verification | Numerical assessment | Qualified use |
| --- | --- | --- | --- |
| Sun (2016), cubic model 1 | Published equation, coefficients and all 144 printed Table 1 observations checked | Rounded parameters and errors reproduced under an inferred unweighted pressure-residual objective | Published thermal comparison within the reported range; original measurements and workbook discrepancies remain unresolved |
| Fu (2023), published joint Sun–Gréaux fit | Published equation conventions and Table S3 coefficients checked | **Not reproduced** from the 174 reconstructed candidate observations | Retain the unchanged literature parameterization with a non-reproduction warning; the discrepancy does not establish that every prediction is invalid |
| Peritheos unweighted joint candidate-data refit | Source observations and conventions checked | Independent refit reproducible for its explicit unit-weight objective | Opt-in comparison; separate coefficients, errors and identity, never a silent replacement for Fu |
| Chen (2018), tetragonal 300 K Vinet | Published coefficients and all seven I4/mcm rows checked | Qualified numerical check; original weights, pressure errors and optimizer not recoverable from the rounded table | Published comparison over 28.824–62.477 GPa; K0′ fixed at 4, V0 extrapolated, source 2-sigma errors retained |

Sun and Fu describe the high-temperature cubic branch. Their 300 K reference
states are extrapolations, not evidence of a stable cubic phase at 300 K. Chen's
parallel P4/mmm refinements are structural diagnostics and must not be counted
as extra independent observations in the stated I4/mcm fit. Theoretical EOS
records, fitted-model benchmarks and digitized observations retain their
existing classifications; none is promoted to a direct experimental table by
these qualifications.

## Sun original-data recovery closure

The user reports that Ningyu Sun no longer has the original Sun (2016) data.
This is recorded as **user-reported private author communication**, reported on
2026-10-08; the date of the author's communication is not inferred. No original
correspondence is archived or independently verified here. Further recovery is
closed. The published observations and EOS remain available unchanged.

The previously supplied `Sun et al. 2016 data.xlsx` remains separate evidence
where held by the user. It is not bundled in Peritheos and is not a replacement
for Table 1. The previous workbook audit (2026-09-29) left incomplete coverage
and row correspondence, volume-error discrepancies, zero-error meanings, a
candidate central-volume mismatch and an unidentified alternative Pt EOS
unresolved. Preserve those uncertainties; do not fabricate missing rows,
recalibrate pressures, or apply blanket corrections to volume errors.

Numerical agreement with the printed Sun table is compatible with these
unresolved measurement questions: it checks the published parameterization,
not the unavailable original measurements. Fu's unknown row selection,
objective and weights remain a separate unresolved reproduction problem.
Peritheos's formula/molar normalization check does not justify changing the
physical five-atom CaSiO3 normalization or any published coefficients.

Detailed evidence and diagnostics remain in the Sun reproduction script,
[Fu audit](fu-2023-eos.md), and [Chen audit](chen-2018-casio3.md). The aggregate
primary-source ledger repeats the concise reproduction status and usage
qualification; full evidence stays in each material's validation metadata.

## Studio presentation

Studio supplies compatible record notices for the existing pinned identities.
The notices convey these qualifications without rebuilding or repinning the
engine, rewriting archived source exports or regenerating observation assets.
The published Fu warning and independent Peritheos refit remain distinct.
