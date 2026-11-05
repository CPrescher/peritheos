# Barker numerical datasets: source review and optional distribution

Reviewed 5 November 2026 against `58f2090b`. Both existing Barker resources are
admitted to the unreleased `peritheos[datasets]` selection under a scoped
numerical-content assessment. Their original CSV bytes remain unchanged.

| Resource | Preserved contents | Scientific qualification |
| --- | --- | --- |
| `argon-barker-1987-table1-liquid.csv` | Seven Table I rows; 100 K, 27.04 cm³/mol; pV/NkT and U/NkT | One experimental comparator and six theoretical results, all including quantum corrections; no reported errors or upstream origin for the experimental row |
| `argon-barker-1987-figure5-selected.csv` | 21 BFW+AT curve samples and 11 selected Ross comparison marks; 298 K | Incomplete graphical subset; solid phase unspecified; picking bounds are not experimental errors |

The package now contains **451 material dataset entries plus two supporting-study
entries**, using **392 distinct CSV resources**. The material catalog remains
456 entries, with the same five held entries referring to the original Sakai S8
and Campbell compilations. No argon phase or executable Barker EOS is invented.

## Source evidence and preservation

The current Zotero Barker PDF is byte-identical to the original user-supplied
source: SHA-256 `d8a43b6186bdfbeaa00b8113e5275efe396bf44bab2fc4664fa6a68c2d6a6b98`.
Methods, references, Table I on printed p.1510 and Figure 5 on p.1511 were
inspected; both relevant pages were rendered and visually checked. All seven
printed numerical pairs match the CSV. The Figure 5 caption confirms the solid
BFW+AT curve and Ross experimental-circle identities. The existing coordinates
and axis calibration are retained; a new independent point-picking pass is not
claimed. No additional high-pressure or ambiguous observations are fabricated.

The Ross primary PDF, DOI 10.1063/1.451346, was also inspected for experimental
origin and methods. Its current SHA-256 is
`8f5f50ab785adb5f03c3d5769d1c3b91a0ec7c3e5eb3a1374faa273112aaf9c3`.
The admitted marks remain digitized from Barker's replot; they are not replaced
with exact Ross table rows or given uncertainties inferred from another table.
Barker's liquid experimental comparator is not silently attributed to Ross's
solid experiment. No pressure recalibration is performed.

Original study metadata, the two CSVs, the diagnostic report and source binaries
are preserved privately. The canonical study receives additive typed columns,
attribution and rights qualifications only. Removing those additions restores
the original study JSON exactly. The material/EOS scientific snapshot and the
original diagnostic report remain unchanged. The study remains theoretical,
`not_reproduced`, with no self-contained potential/EOS implementation recovered.

## Rights assessment

The primary article and [IBM's author-institution record](https://research.ibm.com/publications/high-pressure-equation-of-state-for-solid-argon-from-interatomic-potentials)
identify AIP copyright. **No separate open-data grant was found** for Barker or
Ross. Source availability, a user-supplied PDF and an attribution are not grants;
no copyright-expiry or worldwide public-domain claim is made.

The admission concerns these numerical excerpts. The liquid table contains
physical numerical comparisons and ordinary potential labels, with one source
experimental comparator explicitly separated from six theoretical values.
The graphical resource contains selected numerical center/curve coordinates,
functional series identities and contributor conversions. It contains no
source image pixels, vector drawing paths, copied caption, article prose,
original table design or figure artwork.

The [CC factual-content/database guidance](https://creativecommons.org/faq/#which-components-of-databases-are-protected-by-copyright)
distinguishes individual facts from protected expression and independently
protected collections. No separately protected substantial collection investment
was established for these particular primary-study numerical excerpts. This
is a reasoned maintainer applicability assessment, not proof that all rights
or acquisition contracts are absent. A small row count or experimental origin
alone is not an exemption. No whole-publication or whole-database reuse is
admitted by these decisions, and applicable third-party rights remain binding.

New contributor-held typed metadata and review annotations are dedicated CC0;
original source works and rights are not. The source rights remain explicitly
unspecified for numerical reuse, with the published copyright notices retained.
The [bound receipt](../data/barker-1987-redistribution.json), also installed as
`barker_review`, records each original/current declaration, exact resource and
source-copy identity. The installed scoped notice explains the same boundaries.

## Installed access without inventing a material

Supporting studies now have their own typed dataset interface:

```python
from peritheos import get_study_dataset, get_study_document, list_studies

assert "argon_barker_1987" in list_studies()
study = get_study_document("argon_barker_1987")
assert study["reproduction"]["status"] == "not_reproduced"

liquid = get_study_dataset("argon_barker_1987", "argon_barker_1987_table1_liquid")
assert len(liquid) == 7
print(liquid["pv_over_nkt"])

figure = get_study_dataset("argon_barker_1987", "argon_barker_1987_figure5_selected")
experimental = figure["determination_method"] == "experimental"
assert experimental.sum() == 11
```

Core `peritheos` retains public study context and declarations, including the
scientific limitations; numerical values require the version-matched optional
package. Missing, incompatible or corrupted data use the existing actionable
`DatasetError` behavior. Repository/editable loading retains the original CSV
resources. No network access, automatic download, calibration, fit, phase
assignment or material/EOS registration is introduced. Private audit history
and source PDFs are absent from the exported study metadata.

The source study JSON is a supporting-study record, not a new EOSMAT format.
There is no change to material discovery or the material observation-view
registry. Package manifests distinguish `dataset_count` (material entries),
`study_dataset_count` and `total_dataset_count`; supporting approvals and
resources are verified separately against the same installed review ledger.
