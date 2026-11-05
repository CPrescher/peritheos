# Morard cold-source tracing and Campbell table regeneration

Reviewed 5 November 2026. This follows the initial
[numerical redistribution review](numerical-dataset-redistribution.md).
At this review stage the optional selection was **415 of 453 material dataset entries**:
two previously held Morard entries are admitted, and four new material entries
expose one separately generated Campbell polynomial grid. The 38 remaining held
entries include the unchanged original Campbell integration collection.
The two Barker supporting-study entries remain deferred. The split is unreleased.

The subsequent [collection review](shen-historical-redistribution-review.md)
raises the current selection to **427 of 453 entries, with 26 held**. The
415-entry count above describes the preceding Morard/Campbell stage.

## Morard: all 21 cold inputs identified

The original author file name suggests Sata/Ohfuji, but it is not a complete
study attribution. Thirteen cold rows match Ohfuji et al. (2007), Table 2,
[DOI 10.1007/s00269-007-0151-0](https://doi.org/10.1007/s00269-007-0151-0),
reprinted in Sata et al. (2010), Table 1,
[DOI 10.1029/2009JB006975](https://doi.org/10.1029/2009JB006975).
The remaining eight match Ono and Kikegawa (2006), Table 2, p. 1943,
[DOI 10.2138/am.2006.2347](https://doi.org/10.2138/am.2006.2347).

| Author cold-row IDs | Original study | Published P (GPa) | Published cell V (Å³) |
| --- | --- | --- | --- |
| 001–013 | Ohfuji (2007), reprinted Sata (2010) | 39–214.5 | All 13 row mappings retained in the lineage receipt. |
| 014 | Ono/Kikegawa (2006), Table 2 | 43.2(1) | 82.32(9) |
| 015 | Same | 49.3(1) | 81.07(11) |
| 016 | Same | 56.3(1) | 78.59(9) |
| 017 | Same | 58.9(2) | 78.63(9) |
| 018 | Same | 64.7(1) | 77.55(6) |
| 019 | Same | 67.2(1) | 76.96(12) |
| 020 | Same | 110.6(2) | 70.04(5) |
| 021 | Same | 120.1(3) | 68.85(5) |

All 21 author molar volumes equal the published conventional-cell volumes
multiplied by the original **0.15055** factor, exactly at the preserved decimal
precision. This corresponds to the rounded Avogadro value 6.022e23 with Z=4.
It is not silently replaced by an exact-SI conversion. Twenty pressures agree;
`literature_cold_012` remains **101.2 GPa in the author file versus 101.1 GPa
in Ohfuji/Sata**. That difference is unexplained, not corrected or attributed
to a particular recalibration. The original non-monotone volume pair also remains.

The author errors (1% P, 0.5% V and 5 K) are assigned fit inputs, distinct from
the original tables' parenthesized last-digit errors. No confidence convention,
missing calibrant reading or covariance is invented. The 11 author-excluded
quenched observations remain excluded; this review does not revisit scientific
acceptance, coefficients, fits or the residual-envelope qualification.

The [lineage receipt](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/morard-cold-source-lineage.json)
records all row IDs, original tokens, differences and source bindings. It uses
the preserved exact Ohfuji primary-PDF identity receipt. Ono identification uses
[indexed text of the primary Table 2](https://rruff.info/uploads/AM91_1941.pdf);
fresh direct PDF requests returned 403/timeouts. **No fresh Ono PDF binary or
layout verification is claimed**, and no hash is invented for an unavailable PDF.
Pair identification does not reconstruct the undocumented author acquisition path.

The unchanged 167-row selected input and its overlapping 21-row literature subset
are admitted under the recorded author-held CC BY 4.0 grant plus a scoped
numerical-facts assessment of these identified primary-study measurements.
Ohfuji/Ono/Sata source terms remain unspecified; this is not a new source grant
or Debian approval. Original declarations, including their earlier unresolved
identity notes, remain intact. The updated installed notice and bound review
carry this subsequent provenance finding and the complete lineage receipt.

## Campbell: separate polynomial tables can be regenerated

The preserved final [Campbell (2009) paper](https://doi.org/10.1016/j.epsl.2009.07.022),
p. 560, Eqs. (5)–(6), was reread and visually checked using the PDF skill.
Its SHA-256 remains `e48e23a348370b1613f4e9439bd39e5d49b94c01b1e6c3d3d9a797f092ebd9f0`.
The [official workbook](https://ars.els-cdn.com/content/image/1-s2.0-S0012821X09004373-mmc2.xls)
was freshly retrieved and matches the original curation SHA-256
`b857b9e40f76389dd78da9b757c504b17607d69e4cdc55087fff6e2cb8566f30`.
The original paper/workbook and material files are preserved privately; no source
PDF or workbook is added to either distribution.

Rows 26–32 of Tables S4–S5 supply mathematical polynomial expressions,
coefficient names/values and fit summaries. The independent standard-library
script `scripts/generate_campbell_2009_polynomial_tables.py` evaluates these
expressions using only the numerical coefficients and an explicitly declared
grid: 0–100 GPa in 1 GPa steps and 1000–2600 K in 100 K steps.
`generate()` does not read the original output collection. `compare_original()`
separately compares all 3434 keyed results against the preserved originals.
The source parameters and precise formulas ship inside each dataset declaration.

| Approximation | Rows | RMS difference, log10 fO2 | Maximum absolute difference | Source stated RMS / max |
| --- | ---: | ---: | ---: | --- |
| IW, Table S4 polynomial | 1717 | 0.01757083636 | 0.05136112623 | 0.018 / 0.046 |
| NNO, Table S5 polynomial | 1717 | 0.01304040941 | 0.05389992863 | 0.013 / 0.054 |

The IW maximum exceeds the workbook's stated 0.046. That source-output
qualification is retained; no new fit repairs it, and coefficient rounding is
not asserted as the cause. These are **polynomial approximations**, not exact
original tables, independently measured data or validated replacement EOS.
The dataset identifier is `campbell_2009_independent_buffer_polynomial_grid`
under `fe_fcc`, `feo`, `nickel` and `nickel_oxide_b1`; all four entries share one
new CSV. They are not included in `fit_datasets` or used to promote any EOS.

The original `campbell_2009_tables_s4_s5_buffer_results` entries and their CSV
remain unchanged and held. A separate maintainer mathematical-evaluation
assessment admits the new grid; contributor-held output/arrangement rights are
dedicated under CC0, while the generator is MIT. No open grant is asserted for
the workbook or original collection. The source coefficients are scientific
parameters, not copied author software or publication layout. Downstream
redistributors may still require their own assessment.

## What has not been recovered

Full thermochemical reintegration of the original absolute tables is not
established. Campbell adopts Chase (1998) FeO and Knacke et al. (1991) NiO
standard-state Gibbs data with extrapolation above their solid-state input
ranges. IW additionally switches between fcc and Dewaele hcp iron at the
Shen et al. transition boundary. The existing volume-slope checks do not
recover these absolute reference/extrapolation and boundary choices.
Evaluating the explicit polynomial avoids inventing them but does not reproduce
that original integration workflow. Deferred Fe/FeO source records and all
previous scientific findings remain unchanged.

## Verification and preservation

The [follow-up receipt](../data/morard-campbell-redistribution-followup.json)
binds the before/after scientific snapshots. Removing exactly the four new
calculated declarations restores the previous complete scientific snapshot.
All original dataset declarations, numerical resource bytes and author artifacts
remain unchanged. Tests verify each original material after removing its one
new declaration, all 21 cold-volume matches, original assigned errors, the
pressure discrepancy, polynomial arithmetic, domain errors, every original grid
key and the source mismatch. Package tests retain exact allowlists, hashes and
checks preventing held resources from entering through shared files.

Repository maintenance checks:

```bash
uv run python scripts/generate_campbell_2009_polynomial_tables.py --check
uv run pytest -q -W error tests/test_morard_campbell_redistribution.py tests/test_distribution_split.py tests/test_domain_review_preservation.py
```
