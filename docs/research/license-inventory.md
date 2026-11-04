# Peritheos 0.17.3 license inventory

Reviewed 4 November 2026. This inventories the published files, their recorded
source terms, and the subsequent local licensing changes. It is a complete
accounting of this release, **not completed rights clearance or Debian
acceptance**. It does not decide whether a particular numerical fact is
protected by copyright or database rights.

All **36 official PyPI archives** were downloaded and checked against their
published SHA-256 digests: 35 wheels and one source distribution. The source
archive contains **5,226 files**. Across platforms there are **5,282 paths** and
**5,404 distinct path/content combinations**, including native binaries,
metadata and license copies. Every archive member has an inventory row.

The catalog is also enumerated independently of physical files: **449 dataset
entries**, **685 EOS records**, and **234 structure/diffraction components**.
Of the dataset entries, 435 reference 389 distinct files; 14 contain inline
data. These counts do not represent independent experiments. Coefficients,
structural facts and imported peak patterns are not omitted merely because
they have no separate CSV.

## Inventory files

| File | Contents |
| --- | --- |
| [Dataset source index](../data/license-inventory/dataset-sources.md) | Human-readable source-by-source table of all catalog datasets, including unresolved cases. |
| [Files](../data/license-inventory/files.csv) | Exact paths, hashes, archive membership, scoped terms, evidence, actions and the local overlay. Byte variants remain distinct. |
| [Catalog components](../data/license-inventory/catalog-components.csv) | Every dataset, EOS and structure/diffraction component, including citations, JSON pointers, content kinds and rights statements. Related dataset declarations are linked to EOS records without extending their grants. |
| [Sources](../data/license-inventory/sources.csv) | Source groups across all components, with affected materials, resources and evidence pointers. |
| [Rights statements](../data/license-inventory/rights-statements.csv) | Verbatim structured assertions and their locations, separately for release and local states. |
| [Rust components](../data/license-inventory/cargo-components.csv) | 47 normal/build dependency packages, declared licenses, source URLs and hashes of available license texts. |
| [Archive receipt and totals](../data/license-inventory/summary.json) | All archive URLs, digests, coverage counts, status totals, wheel comparisons and limitations. |
| [Reviewed scopes](../data/license-inventory-scopes.json) | Explicit file-level interpretations and evidence; these do not grant permissions or replace original notices. |
| [CIF license review](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/cif-license-review.json) | Exact COD source-file comparisons, current CC0 declaration and unresolved He thesis-source terms; installed with the package. |
| [License-only source bindings](../data/license-metadata-source-bindings.json) | Exact before/after licensing fields and whole-document digest chains; previous release/merge/cohort snapshots remain unchanged. |

The `release` state describes what shipped in 0.17.3. The `local` overlay
records the author-file declarations and Ma DOCX exclusion introduced in
0.17.4, followed by the Fei image removal, numerical plot replacements and
Tateno workbook/Chidester ZIP exclusions under `Unreleased`. It is a
source-tree overlay, not an archive inventory of 0.17.4 or later builds. Published 0.17.4 archives have not been rewritten.
Later declarations are not represented as notices already included in 0.17.3.

## Dataset findings

| Local status | Entries | Meaning |
| --- | ---: | --- |
| `open_license_declared` | 22 | Explicit open license in existing dataset metadata; not independent authentication of all underlying rights. |
| `open_source_declared` | 21 | A named source/deposit grant mapped to the resource, with its scope. |
| `author_permission_declared` | 1 | Fe7C3 author input: dated author/maintainer CC BY 4.0 declaration. |
| `mixed_author_permission` | 3 | Coesite/Morard author-held rights are CC BY 4.0; literature inputs remain separate. |
| `limited_or_unspecified_source` | 71 | Contributor-only dedications, source qualifications, or explicitly unspecified terms. |
| `restricted_source` | 20 | Source NC/ND terms require an applicability assessment or separate grant. Not a finding that every fact/calculated output inherits them. |
| `source_terms_not_recorded` | 311 | No applicable grant established by the component metadata and reviewed scopes. |

The 685 EOS and 234 structure/diffraction entries have no separate
component-level license assertions. Their citations and available provenance
are indexed. This does **not** mean every coefficient requires author
permission or that numerical facts are proprietary. Conversely, MIT is not
evidence of Peritheos ownership of third-party rights. Related dataset grants
are leads for checking scope, not blanket grants for an entire material.

## Specific findings and actions

**Open deposits.** Ma's named Mendeley workbook/extractions, Hu's PANGAEA
workbook/observations, Zhu's five Mendeley resources, Farla's Figshare tables,
and Morard's named public Figshare resources have scoped CC BY 4.0 declarations.
Katsura's original MATLAB files and workbook are CC BY 4.0 rather than MIT;
the separate 2009 table recovery is outside that deposit grant. The named MgO
COD/LASL inputs retain their documented CC0/public-domain basis. These grants
do not automatically cover adjacent source files or associated publications.

**New author declarations.** The named Fe7C3, Coesite and selected Morard author
files have the dated CC BY 4.0 declarations requested by Clemens Prescher.
Their notices distinguish author-held rights from embedded literature inputs.
They do not extend to publisher articles, unrelated workbooks, the separate
Fe7C3 Table S4, the unselected 178-row Morard input, or final Morard Table S1.

**Ma publisher DOCX.** `ma_2024_sources/supporting-information.docx` was present
in the released source archive and every wheel, despite the old notice saying
it was not bundled. The local exclusion corrects subsequent builds and retains
the original audit copy. The CC BY 4.0 Mendeley workbook remains separate.

**Ismailova (2016) and Sakai (2025).** The former records source CC BY-NC 4.0
terms. The latter records CC BY-NC-ND 4.0 source terms for the affected recovery.
The inventory distinguishes observations, lineage tables and calculated RS
grids: reading coefficients and computing values does not automatically create
a copyrighted adaptation. Dewaele (2018)'s separately licensed CC BY 4.0 input
is explicitly separated from the Sakai source flag. Assess which protected
elements were reproduced; where permission is needed, obtain a separate grant
or exclude the affected content from a cleared source archive.

**Tateno KCl workbook.** `kcl_variant_sources/6779TableS1 revised.xlsx` came
from the official MSA supplement `AM-19-56779.zip`, associated with
DOI `10.2138/am-2019-6779`. No CC0 declaration was established. After 0.17.4,
the original workbook is removed from the current source tree and excluded
from both distribution formats; it is preserved privately. Original manifests,
hashes and numerical CSVs remain unchanged. Numerical-data rights still need
a source-specific grant or a documented rights assessment. Fresh OpenXML
inspection uses an explicitly supplied original verified against its SHA-256;
default replay does not repeat that inspection or replace its historical report.

**Chidester supplement ZIP (excluded after 0.17.4).**
`chidester_2018_sources/AM-18-56212.zip` is a redundant source container with two
CSV tables and filesystem metadata. It is removed from the current source tree
and excluded from both distribution formats, with original bytes preserved
privately and the retrieval URL/hash retained in its unchanged manifest. Both
extracted CSVs remain byte-identical to the archive members. This removes the
redundant container without claiming to clear or restrict the numerical data.
Published 0.17.3/0.17.4 archives retain both source containers.

**Removed Fei source figures (after 0.17.4).** The four original images under
`docs/data/fei-2007-thermal-source-figures/` have been removed from the current
source tree and excluded from future builds because no redistribution grant
is recorded. Three comparison plots also embedded these rasters: the Pt
Figure 2 validation, the 2004/2007 Pt comparison and the Ne Figure 5 validation.
Their current replacements render numerical readouts and calculated curves
without source artwork. The historical inventory explicitly flags those old
overlays as source-bearing content; a generated filename is not evidence that
an image contains only contributor-owned material.

Original images and overlays are preserved in the private audit archive.
Their hashes, source citations, axis calibrations, digitized observations and
curve-readout coordinates remain traceable. The unchanged historical reports
record the earlier image-based audits. Current numerical replay uses the
[preserved pixel readouts](../data/fei-source-curve-readouts.json); it does not
perform a fresh source-image extraction or verify unavailable image bytes.
This removal does not settle the separate question of rights in numerical
transcriptions. Both 0.17.3 and 0.17.4 source archives contain the original
rasters and overlays; the correction applies to subsequent builds only.

**COD/AMCSD CIFs: current source declaration verified.** The vanadium
COD 9012770 copy is byte-identical to the Granada COD mirror; the
superhydrous-phase-B COD 9007379 copy is byte-identical to the primary COD
server. Primary/mirror variants differ only in SVN date/URL comments, with all
other bytes matching. COD's [current CC0 data declaration](https://www.crystallography.net/cod/#license)
is now recorded for these two named resources and the existing COD 1000054
brucite resource. This establishes a recorded source declaration for the CIF
data, not a publication-wide grant or independent authentication of every
underlying right. All CIF bytes and old attribution headers remain unchanged.
The `release` inventory retains the earlier unresolved legacy-header review;
the `local` overlay records this subsequent clarification. The installed
[notice](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/cod-cif-sources.LICENSE.md) and
[review receipt](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/cif-license-review.json) retain
exact URLs, original/served hashes and the comparison method.

**He (2023) brucite CIF: source terms still unresolved.** This is a Peritheos
transcription of thesis Table 4.3, not a COD input. The indexed institutional
record has no populated file license; live institutional record/API/file-detail
requests returned HTTP 503 (browser record request HTTP 429). No fresh source
reuse grant was established. Failed access and public availability are not
permission evidence, nor proof that no grant exists. Keep the contributor CC0
scope separate from original thesis/table rights and seek a dataset-specific
grant or a documented downstream rights assessment. The existing scientific
inspection, original thesis hash and all CIF values remain preserved.

**Contributor-only CC0 and missing terms.** The Brown, brucite, Ca-perovskite,
several cementite, diopside, Dorfman, ferrous-ringwoodite, Katsura-2009, phase-E
and Ye notices dedicate only contributor-held transcription/normalization/
metadata/arrangement rights. The complete index also includes older metal,
calibrant and mineral resources without standalone notices. No source grant
is inferred from public access, a citation, scientific validation, a filename,
or a file's directory.

For each source group, record an applicable file/deposit grant or a documented
rights assessment accepted by the downstream reviewer. If new permission is
needed, ask the relevant rights holders to identify the exact files/data and
authorize CC BY 4.0 or another suitable license, including redistribution and
modification by anyone, commercially as well as noncommercially. Permission
only for Clemens to use data is narrower. An entire article need not be
relicensed if the underlying author dataset can be licensed separately.

## Code and dependency notices

Project implementation and contributor-written research scripts remain MIT.
The retained Katsura MATLAB code is separately attributed CC BY 4.0.
Rust licenses include MIT/Apache alternatives, BSD-2-Clause, Unicode-3.0 and
the LLVM exception. The dependency inventory includes build tools and
all-target conditional dependencies and excludes dev-only packages. The [subsequent Rust notice audit](rust-third-party-notices.md) maps all five
wheel-target graphs, distinguishes target libraries from host tools and ships
a conservative all-47-package notice set plus the complete audited Rust
standard-library copyright document. The installed bundle records exact
archive and notice hashes. Target selection is not a claim that every byte
survives linking. The publication compiler is pinned; dependency/compiler
changes or extra bundled native libraries require a new audit. The original
release's dependency inventory remains historical evidence, not a statement
that these notices were already shipped in 0.17.3/0.17.4.

NumPy and SciPy are external `Requires-Dist` dependencies, not bundled by
virtue of that dependency. Cargo registry sources are not vendored in the
source archive. Debian copyright inventories must include the applicable
full license texts, including separately declared terms such as CC BY-SA 4.0
where those components are retained.

## Verification and reproduction

All 35 wheels contain identical numerical/source data bytes. Five
documentation/notice paths have Windows CRLF-only variants: the installed
licensing guide, material README, PAW research README, Fei supplement
`.gitattributes`, and Ma nested license notice. Exact hashes and archive
membership are retained; file hashes were not normalized. Local overlay
hashes describe this checkout and can differ on a checkout that changes text
line endings.

Tests cover scoped CC0, source NC terms, code/data boundaries, changed archives,
literal citations, complete archive coverage and inline catalog components.
Independent scientific-source preservation checks also pass; this does not
establish redistribution rights. No source observations, coefficients,
original artifacts or historical findings are rewritten by the generator.

Historical audit fingerprints include licensing metadata. The six new exact
bindings retain their historical and previous-reviewed digests and bind the
licensed declarations separately. Tests reverse only the four named license
fields in four datasets and recover the three previous whole documents.
A coefficient mutation still fails even if its new digest is repinned.
Catalog coverage changes only its material-file fingerprint, not any count,
scientific classification or ledger finding.

Repository-only reproduction, using the committed archive receipt:

```bash
mkdir -p .git/license-inventory
cargo metadata --locked --format-version 1 > .git/license-inventory/cargo.json
python scripts/inventory_distribution_licenses.py \
  --archives .git/license-inventory/archives \
  --cargo-metadata .git/license-inventory/cargo.json \
  --download --check
```

Use the maintained environment from [Contributing](https://github.com/CPrescher/peritheos/blob/main/CONTRIBUTING.md).
`--download` retrieves the exact 36 receipt URLs and verifies their digests;
it does not substitute a newer release. The dependency check requires the
release-matched dependency entries in `Cargo.lock`; workspace-only version
bumps are allowed. To record a newly reviewed local overlay, omit
`--check`; only inventory outputs are written. Original evidence is preserved.

Debian maintainers and archive reviewers decide acceptance. This report
supplies the evidence and unresolved actions for that review. See
[Debian Policy](https://www.debian.org/doc/debian-policy/ch-archive.html) and
the [copyright-file specification](https://www.debian.org/doc/packaging-manuals/copyright-format/1.0/).
