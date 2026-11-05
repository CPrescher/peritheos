# Licensing of Peritheos code and scientific data

The source checkout retains complete research evidence. Published builds now
export a runtime package without experimental observation values and a separate
optional `peritheos-datasets` package selected by `peritheos[datasets]`. Each has its
own source archive. The runtime archive installs a shorter, scope-specific
version of this guide. The source inventory below describes research resources,
not a claim that every resource is installed by either package. This split is
unreleased; previously published 0.17.4 archives are unchanged.

Peritheos code uses the **MIT (Expat) license**. The scientific data bundled
with that code have **source-specific terms**. An article citation establishes
scientific provenance; it is not, by itself, a permission to redistribute.

## What each statement covers

| Material | Applicable statement |
| --- | --- |
| Peritheos Python/Rust implementation, contributor-written research scripts in `scripts/`, and contributor-written documentation | MIT, in the root `LICENSE`, unless an explicit file-specific notice says otherwise. Citing or implementing a published equation does not apply the article's license to the script. Explicit additional grants, such as the Sakai recovery's contributor CC0 notice, remain available. |
| A source dataset explicitly released under CC BY 4.0 | The named source dataset and its attributed redistribution retain CC BY 4.0. Preserve the named creators, source, license notice and description of changes. This does not automatically cover the associated journal article. |
| A source explicitly dedicated under CC0 or identified as public domain | The scope established by that source's notice; it does not automatically extend to the associated publication or other sources. |
| A Peritheos contributor's scoped CC0 dedication | Only rights held by that contributor in the stated transcription, normalization, metadata or arrangement. It does not relicense the original observations, article, figures, publisher files or author artifacts. |
| A resource whose source terms are unspecified | No source license grant is asserted. Neither the code's MIT license nor a contributor's CC0 dedication supplies missing third-party permissions. |
| Compiled Rust dependencies and standard library | [Rust third-party notices](third_party/rust/NOTICE.md) and their full license/copyright texts; these do not change the MIT license of Peritheos code or grant scientific-data rights. |
| Retained original author code or workbooks | Their own source terms. For example, Katsura's original MATLAB files and workbook retain CC BY 4.0, even though Peritheos's independent Python diagnostics use MIT. |

Numerical coefficients, measurements and their transcription are distinguished
from article prose, figures and source-file layouts throughout the catalog.
This guide does not assert that all scientific data are free of copyright or
database rights in every jurisdiction. A source's restrictive article license
also does not, by itself, establish the rights in every numerical fact.

## Finding the terms for a resource

Start with its `.eosmat` dataset metadata and companion `*-source.json`,
manifest or README. Follow any `license`, `license_scope`, `source_rights`,
`contributor_rights`, `rights` or source-specific notice. Scientific validation
and recommendation status describe the science, not redistribution permission.
Migration from Dioptas or a beamline collection does not establish ownership
of the source data.

The following is an index of the standalone dataset notices shipped with the
package, not a list of all resources with cleared reuse terms. Their precise
file lists and qualifications control; a family name is not a blanket grant.

| Dataset or investigation | Notice |
| --- | --- |
| Brown NaCl historical recovery | [Brown rights and provenance](datasets/brown-nacl-1999.LICENSE.md) |
| Brucite numerical data and CC0 COD structure | [Brucite source terms](datasets/brucite.LICENSE.md) |
| Ca-perovskite Fu/Sun/Gréaux observations | [Contributor-only CC0 scope](datasets/ca-perovskite-fu-2023-observations.LICENSE.md) |
| Ca-perovskite Noguchi Table 1 | [Contributor-only CC0 scope](datasets/ca-perovskite-noguchi-2013-table1.LICENSE.md) |
| Cementite Hu/PANGAEA shock dataset | [Source CC BY 4.0](datasets/cementite-hu-2019.LICENSE.md) |
| Cementite Ono–Mibe | [Unspecified source terms; contributor CC0](datasets/cementite-ono-mibe-2010.LICENSE.md) |
| Cementite Sata/Li | [Unspecified source terms; contributor CC0](datasets/cementite-sata-2010.LICENSE.md) |
| Cementite Takahashi, Litasov and Prescher | [Separate source scopes](datasets/cementite-sources.LICENSE.md) |
| Coesite Bykova author inputs | [Maintainer-declared author CC BY 4.0; third-party inputs separate](datasets/coesite-bykova-2018-author-inputs.LICENSE.md) |
| COD CIFs: vanadium, superhydrous phase B and archival brucite | [Scoped COD CC0 source declaration](datasets/cod-cif-sources.LICENSE.md) |
| Diopside | [Unspecified source terms; contributor CC0](datasets/diopside.LICENSE.md) |
| Dorfman co-compression tables | [Contributor-only CC0 scope](datasets/dorfman-2012-tables-s1-s6-cocompression.LICENSE.md) |
| Ferrous ringwoodite | [Source terms and contributor CC0](datasets/ferrous-ringwoodite.LICENSE.md) |
| FeS Morard author inputs | [Maintainer-reported author CC BY 4.0; literature inputs separate](datasets/fes-morard-2026-author-inputs.LICENSE.md) |
| Forsterite Katsura workbook, MATLAB and table recovery | [Zenodo CC BY 4.0 and separate 2009 terms](datasets/forsterite-katsura-2022.LICENSE.md) |
| KCl Ma Mendeley workbook and extractions | [Source CC BY 4.0](datasets/ma_2024_sources/LICENSE.md) |
| Fe7C3 Prescher author inputs and fit export | [Author-declared CC BY 4.0; Table S4 separate](datasets/iron-carbide-fe7c3-prescher-2015-author-inputs.LICENSE.md) |
| MgO Tange input recovery | [COD CC0 and LASL public-domain inputs](datasets/mgo-tange-2009-redistributable-inputs.LICENSE.md) |
| Phase E | [Unspecified source terms; contributor CC0](datasets/phase-e-crichton-ross-2000.LICENSE.md) |
| Rubidium halide Farla supplement | [Source CC BY 4.0 and contributor scope](datasets/rubidium-halides-farla-2025.LICENSE.md) |
| Sakai gold and reused Dewaele tables | [Restricted article terms, source distinctions and contributor CC0](datasets/sakai-2025-au.LICENSE.md) |
| Ye co-compression tables | [Contributor-only CC0 scope](datasets/ye-2017-300k-cocompression.LICENSE.md) |
| Zhu Mendeley dataset | [Source CC BY 4.0 and normalization changes](datasets/zhu-2025-mendeley-data.LICENSE.md) |

For resources without a standalone notice or explicit metadata grant, source
reuse terms remain **unspecified**. Do not infer MIT or CC0 from their directory,
file format or the absence of another notice. This also applies to inline EOS
and structural data and to author inputs with limited provenance-only use.

The installed [CIF license review](datasets/cif-license-review.json) records
source URLs, hashes and the current COD data declaration for the three named
COD CIFs. It preserves legacy headers and structural values. The He (2023)
thesis-derived brucite CIF was subsequently replaced by COD 1000054 and its
original evidence preserved privately. Its source grant remains unresolved;
the replacement does not relicense that former thesis transcription.

## Installed packages and downstream redistribution

The [file and catalog license inventory](https://github.com/CPrescher/peritheos/blob/v0.17.4/docs/research/license-inventory.md)
accounts for the published 0.17.3 archives and records the author-permission
and packaging updates introduced in 0.17.4. It lists scoped grants and
unresolved source rights; it is not a blanket redistribution clearance.

This guide and the dataset notices are installed under `peritheos/data/`.
The distribution's `License-File` metadata also includes these notices, the
MIT license, complete CC BY 4.0 and CC0 texts, and the complete Rust notice
bundle. The Rust manifest records locked crate sources, target roles and
notice hashes; its standard-library document matches the pinned official
release compiler. Redistributors rebuilding with another compiler must
refresh/review the matching standard-library notices. Python package metadata uses
`LicenseRef-Peritheos-Distribution` to refer to this component-specific scope;
it does not change the MIT license of the library code or grant unknown rights.

Read the installed guide without a repository checkout:

```python
from importlib.resources import files

print(files("peritheos").joinpath("data/LICENSING.md").read_text(encoding="utf-8"))
```

Wheels contain the installed library and its data. Source distributions also
contain research scripts, reports and selected original artifacts under
`docs/data/`, including the unchanged Katsura CC BY 4.0 author files. Locally
preserved publisher full texts and private evidence excluded by the build
configuration are not covered by contributor dedications or distributed by
virtue of being retained locally. Inspect the actual archive being repackaged;
repository, source-distribution and wheel contents differ.

Ma's publisher supporting DOCX is repository-only audit evidence and is
excluded from both distribution formats. Its separate source terms are not
covered by the CC BY 4.0 grant for Ma's Mendeley workbook, which remains
distributed. The DOCX was inadvertently included in 0.17.3; the exclusion
applies to subsequent builds.

The four Fei source raster figures and three historical raster-bearing
comparison plots were included in the 0.17.4 source distribution. Subsequent
source builds omit the original rasters and replace those plots with numerical
readout/calculation plots that contain no source artwork. Original image
hashes, numerical readouts and scientific findings are retained; the image
bytes and old overlays are preserved only in the private evidence archive.
Published 0.17.4 archives remain unchanged. Rights in numerical transcriptions
remain a separate inventory question.

After 0.17.4, the original Tateno KCl workbook
`datasets/kcl_variant_sources/6779TableS1 revised.xlsx` and Chidester supplement
`datasets/chidester_2018_sources/AM-18-56212.zip` are also excluded from both
source distributions and wheels. Original bytes are preserved privately;
source manifests and hashes remain. The numerical Tateno CSVs and both
byte-identical extracted Chidester CSVs remain distributed. Excluding these
containers does not establish redistribution rights for their numerical data.

For Debian-style copyright inventories, record code and data separately and
retain the individual source notices. An unspecified source license is an
unresolved packaging question, not a finding that the resource meets Debian's
Free Software Guidelines. This guide documents the current evidence; it does
not supply a blanket redistribution clearance for the entire data directory.
