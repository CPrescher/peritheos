# EOS recommendation dispositions and audit repairs

On 17 October 2026 the user authorized seven not-recommended designations,
a hold on the wadsleyite default, and correction of fifteen misleading audit
comparisons identified by the preceding catalog review. These changes preserve
all published coefficients, parameter errors/covariance and observations.
They introduce no replacement EOS or new material default.

## Recommendation changes

The following published records now have `scientific_validation.status:
not_reproduced`, `usage_recommendation: not_recommended_for_quantitative_use`,
and `catalog_access: explicit_selection`. They are not recommended for
scientific predictions or pressure calibration. They remain available for
explicit diagnostic inspection, with their original equation and coefficients.
Strict primary-validated material loading omits them unless the caller
explicitly accepts their scientific status.

| Record | Unresolved finding |
|---|---|
| `ca_perovskite_fu_2023_bm3_mgd_refit` | Published thermal model retains the high-temperature discrepancy after the source/normalization audit; 5.917 GPa pressure RMS on the 140 selected Sun observations |
| `palladium_baty_2024_bm3_1` | The experimental BM3 fit is not recovered from the 78 verified source rows under the completed audit checks |
| `fe09ni01_hcp_sakai_2014_type1_thermal` | About 7.1 GPa RMS on the 32 P7 hot observations, versus about 1.1 GPa reported in Table 7 |
| `fe09ni01_hcp_sakai_2014_type2_thermal` | Same source-data discrepancy |
| `fe09ni01_hcp_sakai_2014_type4_1_thermal` | Same source-data discrepancy |
| `fe09ni01_hcp_sakai_2014_type4_2_thermal` | Same source-data discrepancy |
| `fe09ni01_hcp_sakai_2014_type4_3_thermal` | Same source-data discrepancy |

The Fu identifier's historical `_refit` suffix denotes Fu's **published**
parameterization, not the separately attributed Peritheos candidate-data refit.
The Pd finding affects only the experimental BM3, not Baty's separate DFT
record or the current Pd default. Sakai's separate isotherms retain their
source-specific qualifications. Passing its Table 8 extrapolation checkpoints
does not reconcile the hot-data residuals.

See the [Fu audit](../literature-reproductions/fu-2023-eos.md),
[Pd audit](../literature-reproductions.md#palladium-baty-2024), and
[Sakai audit](../literature-reproductions/iron-fei-input-papers.md).
The Baty-specific evidence review continues in its separate user-requested chat;
it does not block the authorized recommendation designation here.

`wadsleyite_katsura_2009_bm3_1` is now `deferred` and explicitly on hold,
pending pressure-calibration and MGD reference-state/Debye mapping review.
Its unchanged coefficients give 1.89443 GPa RMS and 3.46380 GPa maximum
absolute residual on all 85 source rows. As the only wadsleyite EOS, it is
excluded from normal executable discovery; the document, data and explicit
unvalidated loading route remain. No alternative is silently selected. This
hold does not establish that the publication is wrong.

## Fifteen corrected audit comparisons

| Published records | Count | Correction |
|---|---:|---|
| Jacobsen (2002), Fe27/Fe56/Fe75 | 3 | Use the sample cubic lattice parameter and propagate its uncertainty to volume; do not use quartz calibrant volume |
| Jacobsen (2005), Mg0.73Fe0.27O and Fe0.93O B1 | 2 | Select the correct cubic sample-volume column and its source fit-inclusion flag, instead of Fe0.93O rhombohedral volume |
| Finkelstein (2017), helium/neon cubic | 2 | Use cubic sample volumes and separate medium/phase masks; do not interpret a hexagonal c/a ratio as volume |
| Dobrosavljevic (2019), B1/rhombohedral | 2 | Apply each phase's explicit source fit mask, rather than pooling both branches |
| Thompson (2017), HC epsilon-FeOOH | 1 | Use the twelve HC rows at 30–140 GPa, excluding lower-pressure HOC/transition rows |
| Gong (2004), reference isentrope | 1 | Classify direct coefficient refitting as unavailable from the raw Hugoniot states; do not compare thermodynamic paths without thermal reduction |
| Dewaele (2019), Fe Mao/Dor scale variants | 2 | Use the existing dedicated 53-row hcp source audit with explicit ruby-scale conversion |
| Solomatova (2016), Fp48 HS/LS reference branches | 2 | Use the existing coupled spin-crossover audit, rather than separately fitting each end-member to mixed-spin observations |

These are corrections to validation, not revisions to any published EOS.
The ten corrected sample/phase diagnostics compare coefficient agreement under
their stated objective. Finkelstein's diagnostic does not replay its ambient
V0/neon derivative priors, and the Dobrosavljevic diagnostic does not replay
the MINUTI volume/prior objective. Those limitations remain explicit even when
the conditional coefficients pass the ledger's agreement criteria. The Dewaele
and Solomatova routes retain their dedicated `similar` classifications; Gong
has no spurious isentrope-versus-Hugoniot pressure RMS or regression result.

The previous catalog screening and numerical comparisons are preserved under
`output/eos-recommendation-audit-2026-10-17/`. The implementation backup also
records all material equation/thermal fields, parameter uncertainties and
provenance, all dataset declarations, and checksums of every bundled dataset
file so preservation can be checked independently of recommendation metadata.

## Verification

The preservation check compares all 663 original records' EOS/thermal fields,
parameter errors/confidence/covariance, fixed parameters and parameter/fit
provenance, plus every dataset declaration. All match the saved pre-change
state. All 486 bundled dataset files also retain their exact SHA-256 checksums.

The regenerated record ledger has no extraction/solver failures. For the
fifteen repaired comparisons, ten are conditional `parity` matches, four are
`similar` under their dedicated source audits, and one is `not_refittable`
because its required thermodynamic reduction is unavailable. These labels
are not new author-workflow or physical-validation claims.

## Commit scope

The screening and preservation evidence above describes the shared checkout
at the time of this review, including concurrent work. The committed snapshot
excludes subsequent Au source-recovery changes and the separate Baty, NaCl
and ringwoodite investigations. The common ledgers are regenerated for that
snapshot; its Sakai Au record remains directly unrefittable without the
concurrent source-recovery inputs. Every existing record's EOS, thermal model,
parameter uncertainties and provenance, and every dataset declaration was
compared with the parent commit and remains unchanged. The brucite addition
contains three new isotherms and their primary data; it supplies no recommended
continuous thermal replacement.

The complete pre-change checkout backups remain local; they include unrelated
concurrent source recovery. The small committed `commit-verification.json`
records the separate comparison with the parent Git commit and the tests of
this isolated snapshot.
