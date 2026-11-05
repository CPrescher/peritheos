# Package contents and license inventory

This page describes **Peritheos 0.18.0**, released on 5 November 2026. Start here
for packaging decisions. The [older 0.17.3 audit](license-inventory-0.17.3.md)
is preserved separately; it inventories a much larger package before the split.

**For Dioptas, use `peritheos` without the `datasets` extra.** That installs the
calculation library, structures and EOS parameters. The separate
`peritheos-datasets` distribution supplies selected observation values and other
numerical datasets. It is optional and depends on the same core version.

## What is actually shipped?

| Content in 0.18.0 | Core `peritheos` | Optional `peritheos-datasets` |
| --- | --- | --- |
| Python/Rust calculation and fitting library | Included | Depends on core; no second implementation |
| Material documents and EOS parameter records | 234 documents, 685 records | None separately |
| Source CIFs | 10 COD files | None |
| Calculated reference peak lists | 10 lists, 1,154 lines, in the material documents | None separately |
| Material dataset descriptions | 456 descriptions, references and checksums; no observation rows | Values/declarations for 451 of those entries |
| Supporting-study datasets | Descriptions only | 2 Barker entries |
| Numerical CSV resources | 0 | 392 |
| Tables embedded in dataset JSON | 0 | 14 |
| Original source workbooks, publisher figures/PDFs, author research archives | Excluded | Excluded |
| Detailed reconstruction reports and their raw inputs | Excluded | Excluded |
| Applicable notices, license texts and review records | Included | Included |

A dataset entry is a named table or view. Several entries can use the same CSV,
and some tables are embedded in JSON. **453 entries does not mean 453 files or
453 independent experiments.** The optional data include measured observations,
numerical figure readouts, calculated tables and a source-lineage graph; their
scientific labels remain distinct.

“Core without datasets” does not mean “core without scientific numbers”. It still
contains structures, EOS coefficients, uncertainties and calculated diffraction
lines. It also retains dataset metadata so users can discover the source and
learn how to obtain missing observations.

An NC/ND label in a core dataset description still describes the absent source
table. It does not establish that the table is present in core or apply that
license to all the code and EOS parameters. Likewise, a notice may discuss
source artifacts that are not installed. The contents and manifest identify
the shipped files.

### Where to look in an installed package

| Path | What it contains |
| --- | --- |
| `peritheos/data/materials/*.eosmat` | Core structures, EOS parameters, source references and dataset descriptions. No observation rows. |
| `peritheos/data/datasets/*.cif` | The ten core COD structures. The directory name is retained for compatibility. |
| `peritheos/data/datasets/*.LICENSE.md` | Core structure, pattern and Dioptas attribution notices. |
| `peritheos/data/third_party/rust/` | Rust licenses, source attributions and standard-library notice. |
| `peritheos/data/LICENSING.md` | The installed core licensing guide. |
| `peritheos_datasets/data/datasets/` | Selected numerical resources from the optional package. |
| `peritheos_datasets/declarations/*.json` | Selected material dataset declarations, including the 14 inline tables. |
| `peritheos_datasets/study-declarations/*.json` | The two supporting-study dataset declarations. |
| `peritheos_datasets/notices/` | Source-specific dataset notices. |
| `peritheos_datasets/manifest.json` | Exact selected entries and resource checksums. |
| `peritheos_datasets/redistribution-review.json` | Why each reviewed numerical dataset was included or left out. |

The core package also supplies `eos-structure-rights-review.json`,
`calculated-pattern-inputs.json` and `datasets/cif-license-review.json` for
source-level detail. The optional package carries a copy of the EOS/structure
review for context; it does not add another structure library.

## Licenses: what the reader should learn

**The code is MIT.** Keep the Rust dependency notices and the separate Dioptas
attribution where applicable. The ten named COD CIFs have recorded CC0 source
terms. Other structures, EOS parameters and calculated patterns have their own
citations and notices. A paper citation is not a license statement, and a paper's
license does not automatically answer the rights question for every number.

The optional package has four kinds of inclusion decision:

| Why the entries are included | Entries | What this means |
| --- | ---: | --- |
| Recorded source/author terms, including one newly verified deposit license | 45 | Use the source license and attribution recorded for the particular entry. |
| Author permission combined with separately assessed literature values | 3 | Permission covers the author's contribution; the cited literature inputs are identified separately. |
| Independently generated mathematical outputs or source-lineage information | 7 | Newly generated material, with its input sources and contributor rights recorded. These are not replacements for the original author tables. |
| Assessment of numerical content, including figure readouts and collection extracts | 398 | We included the numerical transcription for the reasons in its review. This is not a new license for the source publication or collection. |
| **Total** | **453** | **451 material entries and 2 study entries.** |

These decision groups are not license names. For example, the last group includes
21 entries whose sources carry NC/ND notices. Those notices remain recorded;
our decision concerns the particular numerical extraction, not redistribution
of the original article or artwork. The
[numerical review](numerical-dataset-redistribution.md) and its source-specific
follow-ups retain the reasoning and limits.

### What the machine-readable labels mean

| Label in a record | Plain meaning |
| --- | --- |
| `open_license_declared`, `open_source_declared`, `author_permission_declared` | A source or author license is recorded; read its stated scope. |
| `mixed_author_permission` | An author contribution and literature inputs have different bases for reuse. |
| `source_terms_not_recorded`, `limited_or_unspecified_source` | No applicable source-data grant is recorded, or a dedication covers only part of the work. This alone does not establish that numerical reuse is prohibited. |
| `restricted_source` | A restrictive source notice, such as NC/ND, is recorded. Check the assessment of the exact content distributed. |
| `admit` | Selected for the optional package. It does not mean “new open license granted”. |
| `hold` | Not selected for that package. It does not necessarily mean “illegal to use”. |

## What is excluded?

Five material entries are not supplied by the optional package. They refer to
two shared resources:

| Resource | Entries not supplied | Reason and available alternative |
| --- | --- | --- |
| `sakai-2025-table-s8-calibration-lineage.csv` | `gold_sakai_2025_table_s8_lineage` | The original author-selected literature compilation is held. A separately named partial primary-source graph is included; it is not a recovered Table S8. |
| `campbell-2009-tables-s4-s5-buffer-results.csv` | `campbell_2009_tables_s4_s5_buffer_results` for `fe_fcc`, `feo`, `nickel` and `nickel_oxide_b1` | The original integration collection is held. Separately named polynomial evaluations and a relative NNO pressure integral are included, with their approximation and missing-baseline limits. |

Neither package includes the removed Fei source images, Ma publisher DOCX,
Tateno workbook, Chidester ZIP, former He thesis CIF or inherited peak lists.
Their preserved evidence and hashes are not permission to redistribute them.
The openly sourced replacement structures and independently calculated peak
lists are in core. No original papers or research workbooks are added by the
optional package.

## What should a Debian maintainer check?

For the intended **core-only package**, use the official `peritheos` source
archive, retain its code/third-party/source notices, and test it with Dioptas.
There is no need to assess the 453 optional dataset entries for that route.
The observation-table exclusions have already been made in the core archive.

The detailed EOS/structure review records questions about source compilations;
it is not a finding that seven Holland coefficient sets, or other numerical
parameters, must be removed. If Debian raises a particular source question,
the [component index](../data/license-inventory/eos-structure-review.csv)
identifies its exact material, EOS/structure component and source. We have not
established a mandatory core exclusion list or obtained a Debian acceptance.

A separate decision would be needed to package **`peritheos-datasets`**: does
Debian accept the documented basis for the numerical transcriptions, including
the entries without a source-data grant and those retaining NC/ND source notices?
If not, exclusions should name those exact entries/resources. This decision is
independent of whether Dioptas needs the optional data.

## Which source archive?

Use the uploaded release asset or PyPI sdist named below. GitHub's automatically
generated “Source code” downloads contain the full repository, including research
evidence omitted from both Python distributions.

| Published 0.18.0 source archive | Files | CSVs | CIFs |
| --- | ---: | ---: | ---: |
| `peritheos-0.18.0.tar.gz` | 477 | 0 | 10 |
| `peritheos_datasets-0.18.0.tar.gz` | 645 | 392 | 0 |

The two official sdists were inspected and checked against the published PyPI
SHA-256 digests for this guide:

- Core: `22d573ed687d7c7d07c29bcc45328f32b4d8110963c1b13e807a0bee33941277`.
- Optional datasets: `177691bbda5ac232beb0e5eae82130d13d9d28d23a6e53c6742c22f502cc97f4`.

The corresponding [GitHub release](https://github.com/CPrescher/peritheos/releases/tag/v0.18.0)
contains these files. Package metadata is also available on
[PyPI for core](https://pypi.org/project/peritheos/0.18.0/) and
[PyPI for datasets](https://pypi.org/project/peritheos-datasets/0.18.0/).

## Preserved detailed evidence

The [historical audit](license-inventory-0.17.3.md) retains the original findings,
source comparisons, file lists and later source-tree notes. In its linked CSVs,
`release` means **0.17.3**, and `local` means the **research checkout reviewed
later**, not either 0.18.0 package. The old source-rights index is not the optional
package's shipping list; use the installed manifest for that.

<a id="specific-findings-and-actions"></a>
For the historical findings referenced by older receipts, see
[specific findings and actions](license-inventory-0.17.3.md#specific-findings-and-actions).
Original archive hashes, permissions, source resources and review decisions have
not been changed by this documentation rewrite.
