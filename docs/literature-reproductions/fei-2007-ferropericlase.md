# Fei et al. (2007): original ferropericlase compression EOSs

Primary source: Y. Fei, L. Zhang, A. Corgne, H. Watson, A. Ricolleau,
Y. Meng and V. Prakapenka, “Spin transition and equations of state of
(Mg, Fe)O solid solutions,” *Geophysical Research Letters* **34**, L17307,
[doi:10.1029/2007GL030712](https://doi.org/10.1029/2007GL030712).

The audit on 2026-09-08 checked the complete primary paper in local Zotero
(item `33JIXYGF`, PDF attachment `FKWA95NJ`), including visual inspection of
pages 2–3, and the three official auxiliary tables supplied by the user.
The record-level source checks contain the PDF SHA-256; each dataset contains
SHA-256 digests for its normalized CSV and its byte-preserved original file.
No copyrighted article PDF is redistributed.

## Published EOS inventory and source ownership

Before this change, no record was owned by the GRL DOI. The two fp39 records
with `solomatova_2016_fei` in their identifiers belong to Solomatova's 2016
Table 7 refit, with that paper's DOI and different coefficients. They remain
unchanged. The existing fp20 material also contains a separate Speziale 2007
EOS. Reuse those composition documents; add only the missing fp58 material.

Fei publishes exactly these five sample EOS parameterizations. All pressures
are in GPa, all volumes in Å³ per conventional B1 cell containing four formula
units, and all reference isotherms are 300 K. Parentheses below are the printed
parameter errors; confidence levels and covariance are not given.

| Composition / branch | V0 | K0 | K0′ | Published approximate fit interval | Primary location |
|---|---:|---:|---:|---|---|
| fp20, Mg0.80Fe0.20O, HS B1 compression | 76.16 | 158(3) | 4 fixed | up to 35 GPa | Figure 1; paragraph 12 |
| fp20, LS B1 compression | 74.2(1) | 170(3) | 4 fixed | 40–95 GPa | Figure 1; paragraph 14 |
| fp39, Mg0.61Fe0.39O, HS B1 compression | 77.48 | 156(2) | 4 fixed | up to 55 GPa | Figure 2; paragraph 12 |
| fp39, LS B1 compression | 73.6(1) | 170(3) | 4 fixed | 60–148 GPa | Figure 2; paragraph 14 |
| mw58/fp58, Mg0.42Fe0.58O, HS B1 compression | 79.36 | 153 | 4 | up to 42 GPa | Figure 3 caption |

Paragraphs 12 and 14 explicitly identify **third-order Birch–Murnaghan**:

`P = (3 K0 / 2) (η^7 − η^5) [1 + (3/4)(K0′ − 4)(η² − 1)]`,
with `η = (V0/V)^(1/3)`. At K0′=4 this reduces algebraically to BM2. Store the
source-named BM3 once, without duplicate BM2 entries. Figure 3 gives the same
K0′=4 parameter convention for mw58; no uncertainties are printed there.
The measured HS reference volumes of fp20/fp39 have no printed error. Only
K0′ is explicitly registered as fixed in the source-owned records; holding HS
V0 during our validation is a separate declared choice.

No separate decompression EOS is published. Figure 2 identifies compression
points as the fitted observations and plots decompression as an additional
comparison. Paragraph 8 reports reversibility and no observable hysteresis;
that does not establish a second fitted coefficient set. No full mixed-spin
EOS or population mixer is published. The relation K0=160−10XFeO is a summary
for MgO-rich HS compositions, not a license to generate additional records or
to replace the explicitly printed mw58 K0=153 with 154.2. The MgO comparison
curve in Figures 1–3 is not a new sample determination by this paper.

For mw58, the cubic–rhombohedral transition is near 44 GPa and the spin
transition near 80 GPa (paragraphs 9 and 15). Neither a rhombohedral HS EOS nor
an LS EOS is supplied. Do not apply its B1 EOS to those regimes.

## Tables, paths and pressure calibration

| Auxiliary source | Material | Rows | Path coverage |
|---|---|---:|---|
| `grl23467-sup-0002-ts01.txt`, S1 | fp20 | 46 | compression; retain original sample IDs and ordering |
| `grl23467-sup-0003-ts02.txt`, S2 | fp39 | 77 | 38 compression; 39 decompression |
| `grl23467-sup-0004-ts03.txt`, S3 | mw58 | 27 | compression, including the later rhombohedral regime |

The originals reside in `peritheos/data/datasets/fei-2007/`; normalized CSVs
are named `<material>-fei-2007-table-sN-pv.csv`. Every printed sample ID,
lattice parameter, pressure, and error is retained. The supplied previously normalized fp39 CSV has 76 rows and omits the final
`mw39c_073` observation (5.96 GPa, a=4.2188 Å). All 76 shared rows agree
exactly; this transcription restores the missing row from the official S2.
The fp39 `mw39c` block is
assigned decompression from paragraph 8, Figure 2 and the descending-pressure
block; paths are curator annotations rather than a new printed table column.
Do not reconstruct collection order from pressure sorting or treat fp20's
sample-ID prefixes as independently fitted EOS branches.

B1 sample cell volumes are calculated as `a³`, with first-order propagated
errors `3a² error(a)`. These are transformations of published lattice
observations, not independently published volume/error columns. Printed zero
errors remain zero strings in the transcription; they do not mean exact
measurements. S3 supplies no sample-lattice errors, so those cells remain blank.
At and above 44 GPa, the S3 scalar `a_mw58` does not supply the full
rhombohedral cell metric. Keep the original scalar, leave the conventional-cell
volume blank, and never fit these rows as B1 sample volumes.

Paragraph 6 names **Brown (1999), The NaCl pressure standard, J. Appl. Phys.
86, 5801–5808**, for B1 NaCl and **Fei et al. (2007), Toward an internally
consistent pressure scale, PNAS 104, 9182–9186,
[doi:10.1073/pnas.0609013104](https://doi.org/10.1073/pnas.0609013104)** for B2.
The latter DOI is not the GRL sample-EOS DOI. Their nominal boundary is 26 GPa,
but S1 includes B1 NaCl at 29.30 GPa and B2 at 26.86 GPa. The normalized phase
flag therefore uses the lattice branch (`a>4 Å`: B1; `a<4 Å`: B2), not a pressure
cutoff. B1 NaCl has Z=4, B2 has Z=1; do not mix their cell volumes.

The original `P_NaCl` values, including S2's −0.10 and −0.09 GPa rows, are
retained without recalibration. NaCl lattice observations are available, but
the exact two calibrations are not both in the shared executable EOS catalog;
recalculation is marked `reference_eos_not_bundled`, not missing observations.
Laser annealing near 1600 K does not turn the 300 K observations into P–V–T data.

## Reproduction and independent validation

`PYTHONPATH=. python scripts/reproduce_fei_2007_ferropericlase.py --write`
regenerates normalized CSVs from the bundled originals and reports residuals
using the unchanged published coefficients. The primary refit campaign also
performs explicitly labelled **independent validation refits**, with unweighted
pressure residuals, K0′=4, HS V0 held at the printed reference value, and LS V0
and K0 free. These are validation results, not Fei's published parameters.
Exact source weighting and row masks are not supplied.

Only compression rows enter validation. The stated approximate fitting windows
supply deterministic selections: HS fp20 ≤35 GPa (21 rows); LS fp20 40–95.48
GPa (19); HS fp39 −0.10–55 GPa (13); LS fp39 60–148 GPa (24); HS mw58 ≤42.75 GPa
(7). The 95.48 GPa endpoint represents the article's rounded 95 GPa maximum.
For mw58, we explicitly assume (as agreed with the user) that the caption's
approximate 42 GPa limit includes S3's 42.75 GPa B1 point. The earlier strict
42 GPa cutoff left only six points ending at 19.98 GPa; the current selection
includes all seven points through 42.75 GPa. This is a documented
rounded-endpoint assumption, not proof of the authors' exact row selection.
The rows at and above the reported 44 GPa structural transition remain excluded. Table rows near spin crossover
are not relabelled as independently measured pure-spin populations.

| Branch | Published-curve pressure RMSE, GPa | Validation V0, Å³ | Validation K0, GPa |
|---|---:|---:|---:|
| fp20 HS | 1.1724 | 76.16 held | 157.2364 |
| fp20 LS | 1.7678 | 76.3733 | 144.1865 |
| fp39 HS | 0.9371 | 77.48 held | 155.8327 |
| fp39 LS | 4.0475 | 73.5935 | 176.6534 |
| mw58 HS | 0.4879 | 79.36 held | 153.1774 |

The [fp58 HS endpoint decision](#fp58-hs-endpoint-decision) below records why
the 42.75 GPa point is included and the resulting agreement with the source.

The fp20 LS validation fails the campaign's coefficient-similarity criterion;
fp39 LS meets its combined uncertainty criterion but has a substantial
published-curve residual. Those facts do not justify replacing either
published coefficient set. Source ownership and transcription are verified;
exact replication of the original regression remains qualified. See the
[record-level refit ledger](../primary-eos-refits.md) for uncertainties,
residuals, selections and classifications. No validation curve, analytical
checkpoint or spin-transition summary has been added as a primary observation.

Focused tests cover source ownership, coefficient/error conventions, BM3
execution and volume basis, byte checksums, every original table field,
zero/missing errors, path separation, phase limits, and nonmutation of the
published records by validation. The baseline suite passed 1,860 tests before
these changes; no pre-existing test failures were found.

Published-record verification before the additional refit: 155 focused tests
and all 1,866 full-suite tests passed;
Ruff and `git diff --check` passed. Existing EOS and refit records were compared
against the starting commit and remain unchanged. An unrelated floating-point
change in regenerated Baty palladium diagnostics was excluded from this change.

## fp58 HS endpoint decision

**Decision recorded 2026-09-08:** interpret Figure 3's “up to 42 GPa” as a
rounded endpoint that includes Table S3's **42.75 GPa B1 observation,
`mw58_033`**. This is the adopted audit assumption, not confirmation of the
original authors' exact fitting mask. Exclude observations at and above
44 GPa, where the paper reports the cubic–rhombohedral transition.

The endpoint matters because the preceding observation is at only 19.98 GPa:

| Selection | Points | Actual pressure range, GPa | Refitted K0, GPa (formal 1σ) |
|---|---:|---|---:|
| Earlier strict ≤42 GPa cutoff | 6 | 9.37–19.98 | 156.38 ± 1.58 |
| **Adopted ≤42.75 GPa cutoff** | **7** | **9.37–42.75** | **153.18 ± 1.43** |
| Fei's published value | — | caption: up to 42 | 153; no uncertainty printed |

Both validation fits hold **V0=79.36 Å³ per conventional B1 cell** and
**K0′=4** fixed and minimize unweighted pressure residuals. For the adopted
seven-point selection, the published curve has pressure RMSE **0.48788 GPa**;
the validation refit has **0.48726 GPa**. Thus the published K0 is recovered
closely without changing its stored value.

The validation classification remains **`similar`**, since Fei prints no
K0 uncertainty for fp58. Formal refit errors do not include calibration
systematics. The catalog retains the source-owned
`mg042fe058o_fei_2007_hs_b1_bm3` coefficients unchanged; no separate fp58 refit
EOS entry was added. Its validity notes, `validation_selection` and
`scientific_validation.reproduction_summary` explicitly record this decision.
The numerical result is also stored in the
[machine-readable refit ledger](../data/primary-eos-refits.json).

Verification after adopting this endpoint: **145 focused checks and all 1,867
full-suite tests passed**. The regression test asserts inclusion of
`mw58_033`, seven selected points, the 42.75 GPa maximum, and exclusion of the
rhombohedral observations.

## Additional Peritheos fp20 LS refit

At the user's request, the catalog now also exposes
`mg080fe020o_fei_2007_ls_b1_bm3_refit` as an explicit `record_kind: refit`,
derived from `mg080fe020o_fei_2007_ls_b1_bm3`. There are still exactly **five
source-owned published parameterizations**; this sixth entry is a Peritheos
fit to Fei's observations. The Fei DOI identifies the data source, and the
record label, reference details, parameter provenance and audit explicitly
identify Peritheos as responsible for the refitted coefficients.

The fit uses all 19 compression points from 40.95 to 95.48 GPa, corresponding
to paragraph 14's approximate 40–95 GPa interval. It does not choose a higher
cutoff to improve agreement with Fei's coefficients. V0 and K0 are free,
K0′=4 is fixed, and unweighted pressure residuals are minimized. The executable
entry contains the exact source row numbers, sample IDs, dataset checksum,
solver and package versions, full covariance and correlation matrices, and
residual statistics. No additional dataset or synthetic observations are added.

| Parameter | Peritheos refit | Formal 1σ error |
|---|---:|---:|
| V0, Å³ per conventional B1 cell (Z=4) | 76.37327748 | 0.37404085 |
| K0, GPa | 144.18647281 | 4.17820021 |
| K0′ | 4 fixed | not fitted |

The V0–K0 correlation is −0.9927063. Covariance is calculated from the local
Jacobian and scaled by SSE/(19−2); uncertainties exclude pressure-calibration
systematics and uncertainty in phase/row selection. V0 is an extrapolated
zero-pressure reference volume, not a newly observed ambient LS volume.

On these same observations, the pressure RMSE decreases from 1.76784 GPa
for the published curve to 0.86562 GPa for the refit. This is an **in-sample
improvement**, not independent evidence of a physically superior LS EOS.
Near-transition rows may retain mixed-spin character; use within the observed
pressure interval. The published record and the existing material default
remain unchanged. The higher-cutoff sensitivity variants remain diagnostics,
not extra catalog records.

Reproduce the fit, errors, covariance and selected rows with:

```sh
PYTHONPATH=. python scripts/reproduce_fei_2007_ferropericlase.py --refit-fp20-ls
```

Select the new executable entry explicitly:

```python
from peritheos import get_eos_record

refit = get_eos_record("mg080fe020o_fei_2007_ls_b1_bm3_refit")
pressure_gpa = refit.pressure(60.0)  # conventional-cell volume in Å³, at 300 K
```

The global refit ledger labels its numerical check `stored_refit_reproduction`:
parity here is against the stored Peritheos coefficients, not Fei's published
ones. The published fp20 LS entry retains its failed coefficient-reproduction
result. A regression test independently checks covariance using the analytic
BM2 Jacobian, as well as row selection, catalog execution and unchanged source
parameters.

Verification including the additional refit: 156 focused tests and all 1,867
full-suite tests passed. Ruff and `git diff --check` also passed.
