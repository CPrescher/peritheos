# Baty (2024) experimental palladium BM3: unresolved reproduction

Reviewed 17 October 2026. This is a focused evidence review of
`palladium_baty_2024_bm3_1`, following the
[curation stopping rule](../adding-materials-and-eos.md#stopping-rule-unresolved-reproduction-discrepancies).

**Disposition: `not_reproduced`, explicitly not recommended for scientific
predictions or pressure calibration.** The published coefficients and all
observations should remain unchanged. The numerical formulation is resolved;
the disagreement between the published parameterization and its source table
is unresolved. A replacement fit is neither needed nor justified by this
review.

## Source authority and scope

The primary article is Baty et al., *Palladium at high pressure and high
temperature: A combined experimental and theoretical study*, J. Appl. Phys.
**135**, 075103 (2024), [DOI:10.1063/5.0179469](https://doi.org/10.1063/5.0179469).
The experimental section describes room-temperature fcc Pd in helium, with
pressure from a tungsten marker using Dewaele et al. (2004). The paragraph
discussing room-temperature EOS comparisons states that the authors' 0–80 GPa
measurements are described by BM3 with `V0=14.72(4) Å³/atom`, `K0=190(3) GPa`,
and `K0'=5.3(2)`.

The [official supplementary collection](https://doi.org/10.60893/figshare.jap.c.7021680)
provides Table S1 on PDF page 2 and Table S3 on page 3. Fresh downloads on the
review date match the previously recorded supplement provenance:

| Source | Stable download | SHA-256 |
|---|---|---|
| Final supplementary PDF | [Figshare file 44186192](https://ndownloader.figshare.com/files/44186192) | `cf2e4ba7cbcef2e9c4e437752f22a8092fa2af4db0af641a36c82995232ad35a` |
| Deposited final LaTeX source | [Figshare file 44186189](https://ndownloader.figshare.com/files/44186189) | `02d825458eff5dec00add56fb4b9705f4acd34246846cbf9169d237dd2382ada` |
| Bundled Table S1 CSV | `peritheos/data/datasets/palladium-baty-2024-table-s1-compression.csv` | `513835bba914f6fa1f02646901e0453edda808c6da7bc425fbb71c8cf1c99829` |

Table S1 contains 78 observations at 300 K, from 0 to **73.7963 GPa**, rather
than an exact 80 GPa endpoint. All 78 pressure/atomic-volume pairs were freshly
compared as exact decimal values with the official LaTeX table, retaining
pairing and multiplicity; every pair agrees. The CSV sorts observations by
pressure rather than preserving the supplement's three-block layout. Every
conventional-cell volume is exactly four times the source atomic volume.
The PDF table layout and caption were also visually checked.

The caption gives dataset-wide uncertainty limits: lattice parameter below
0.003 Å, atomic volume below 0.01 Å³/atom, and pressure uncertainty increasing
from 0.06 to 0.12 GPa. It supplies neither row-wise uncertainties nor a
confidence convention. These limits must not be turned into invented row
weights.

## Distinguish the two Baty records and the current default

| Record | Represented parameterization | Consequence of this review |
|---|---|---|
| `palladium_baty_2024_bm3_1` | Experimental 300 K BM3: conventional-cell `V0=58.88(16) Å³`, `K0=190(3) GPa`, `K0'=5.3(2)` | Affected: unresolved reproduction, not recommended |
| `palladium_baty_2024_bm3_dft_2` | Static 0 K DFT BM3: conventional-cell `V0=58.4668 Å³`, `K0=195.0(2.5) GPa`, `K0'=5.1(1)` | Outside this experimental discrepancy finding |

The public API's `get_material("palladium").default_record().identifier`
initially returned `palladium_baty_2024_bm3_dft_2` during this review.
The user initially selected Frost's published 300 K Vinet on 17 October 2026.
The subsequent [Frost audit](palladium-frost-2023-audit.md) did not pass strict
fit reproduction, so that selection has been withdrawn and **automatic Pd
equilibrium default selection is on hold**, without a replacement default.
The Baty DFT record is retained separately and explicitly nondefault at
`temperature_ref=0.0`. The experimental Baty record remains `not_reproduced`
and not recommended for quantitative use.

Table S3 explicitly labels its Baty column `BM3 (T=0)`. Its pressures are EOS
evaluations, not independent experimental observations or the original DFT
fitting grid. Agreement with that column cannot validate the experimental
triplet against Table S1. The experimental BM3 is not a Table S3 column.

## Direct evaluation of the unchanged publication

Equation (1) is standard BM3:

\[
P=\frac32 K_0(\eta^{7/3}-\eta^{5/3})
\left[1+\frac34(K'_0-4)(\eta^{2/3}-1)\right],
\qquad \eta=\rho/\rho_0=V_0/V.
\]

A fresh literal evaluation and the Peritheos BM3 evaluator agree to within
`1.96e-13 GPa` over all 78 rows; a numerical agreement tolerance of `1e-10 GPa`
comfortably separates floating-point effects from the observed discrepancy.
Atomic and conventional-cell forms have identical `V0/V`.

For `residual = P_published(V_observed) - P_observed`, direct evaluation gives:

| Metric | Result |
|---|---:|
| All-row pressure RMS | 1.273118191730 GPa |
| Mean residual | −0.592384137467 GPa |
| Maximum absolute residual | 2.718538426094 GPa |
| 0 ≤ P < 20 GPa: 15 rows, RMS | 2.196661 GPa |
| 20 ≤ P < 40 GPa: 35 rows, RMS | 0.983718 GPa |
| 40 ≤ P ≤ 73.7963 GPa: 28 rows, RMS | 0.848852 GPa |

For example, the observed pair `P=12.8983 GPa`, `V=14.04041159 Å³/atom`
gives **10.179761573906 GPa** from the published coefficients. This is an
internal source-data comparison, independent of another laboratory's EOS.
The published RMS agrees with the existing ledger within `3e-14 GPa`.

The table's decimal resolution cannot explain this difference. Using each
source pressure token's printed precision and the eight-decimal atomic
volumes, a local first-order pressure-equivalent rounding allowance is below
`5.1e-5 GPa`. This checks rounding of observations, not parameter covariance.

For scale only, `0.12 GPa + |dP/dV_cell| × 0.04 Å³/cell` ranges from 0.249 to
0.578 GPa across the table. Fifty-three residuals exceed even that sum of
coordinate uncertainty limits. This is a conservative local comparison, not
a chi-square test or a parameter confidence band: the source does not define
row-wise statistical errors or the covariance of its reported coefficients.

## What the existing diagnostics establish

These are **previously stored diagnostic results**, read from
[the primary refit ledger](../data/primary-eos-refits.json) and the
[existing Pd audit](../literature-reproductions.md#palladium-baty-2024).
No optimizer or new refit was run for this review.

| Existing calculation | Rows | V0 (Å³/cell) | K0 (GPa) | K0' | Pressure RMS (GPa) |
|---|---:|---:|---:|---:|---:|
| Unweighted pressure residuals | 78 | 59.9603 | 152.0573 | 6.09936 | 0.704337 |
| Unweighted volume residuals | 78 | 59.4701 | 178.8237 | 5.00994 | 0.777234 |
| Maximum-error EIV sensitivity check | 78 | 59.6387 | 169.0629 | 5.38318 | 0.733536 |
| Published V0 fixed, pressure residuals | 78 | 58.88 fixed | 207.4765 | 4.29636 | 0.972168 |
| Published V0 fixed, P ≥ 40 GPa | 28 | 58.88 fixed | 193.6640 | 5.02431 | 0.754310 |

Direct replay of the stored unweighted diagnostic coefficients reproduces
its `0.704336692085 GPa` RMS. That verifies the comparison, not publication
parity or a replacement EOS. The generous EIV diagnostic assigns the maximum
reported coordinate errors to every row; its stored reduced chi-square is
8.1077, conditional on treating those limits as independent standard errors.
The source does not establish that interpretation. No listed objective
jointly recovers the published triplet. The closer high-pressure subset has
no source-supported exclusion mask and cannot reproduce the stated full-range
fit.

Fresh sorting of Table S1 also reproduces 12 adjacent volume increases. The
largest is `0.04196845 Å³/atom` between 65.2306 and 65.5269 GPa. This scatter
exceeds the quoted per-point volume bound; it does not identify its physical
or reduction origin. The existing Frost comparison is secondary evidence of
an interstudy offset, not the basis for rejecting this experimental fit.

## Interpretation and defensible disposition

### Published experimental description is not a recovered fitting workflow

A close reading of printed page 4 distinguishes three statements. The
authors explicitly describe fitting the DFT P(V) results, explicitly call
the `(K0, K0')=(183.8 GPa, 5.36)` pair their own fit to Mao et al.'s data,
and say that their own experimental 0–80 GPa observations can be described
by the separate `(V0, K0, K0')=(14.72 Å³/atom, 190 GPa, 5.3)` triplet.
For this last experimental parameterization, they do not describe an
independent regression procedure or establish whether only their own rows
entered it. There is also no statement that Mao or Frost observations were
pooled with their observations. This review therefore does not claim either
an established own-data-only fit or an established combined-data fit. The
unresolved finding concerns the printed parameterization's agreement with
Table S1, regardless of how its coefficients were obtained.

The checks resolve source identity, observation transcription, pairing,
volume normalization, pressure units, 300 K reference state, and the BM3
equation. The exact regression objective, weights, fixed-parameter choices,
row selections or averaging, uncertainty calculation, and coefficient
covariance remain unspecified. The row-wise tungsten observations required
to replay pressure reduction are also unavailable in the bundled source
data.

Undocumented fitting choices, a reduction inconsistency, or a fitting/reporting
error could explain the mismatch. The available evidence does not select
among them. Neither the lower residual of an exploratory fit nor the known
hydrogen sensitivity of Pd establishes what happened in this experiment.
The supported conclusion is **unresolved reproduction of this published
experimental parameterization**, not proof that the observations or the
article's separate DFT results are invalid.

Apply `not_reproduced` and explicit non-recommendation to the experimental
record across its claimed room-temperature 0–80 GPa scope, with the directly
tested range stated as 0–73.7963 GPa. Preserve the source coefficients,
uncertainties, and observations. Keep every exploratory result diagnostic
and unsuitable for scientific prediction or pressure calibration. Reopening
would require source evidence that reconciles the published coefficients
with the observations, such as the original fit inputs/protocol or an
authoritative correction; further speculative fitting is not the next
curation step.

The initial review added only this Baty-specific note. The subsequent
user-authorized default selection changes only Pd's selection flags and
adds selection regression coverage; source coefficients, observations,
scientific statuses, shared ledgers, counts, scripts, and the main
recommendation audit remain unchanged by that selection. No author contact,
replacement-model integration, commit, push, or publication was performed.
