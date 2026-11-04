# Licensing of Peritheos code and scientific data

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
| Brucite, including COD and thesis structures | [Brucite source terms](datasets/brucite.LICENSE.md) |
| Ca-perovskite Fu/Sun/Gréaux observations | [Contributor-only CC0 scope](datasets/ca-perovskite-fu-2023-observations.LICENSE.md) |
| Ca-perovskite Noguchi Table 1 | [Contributor-only CC0 scope](datasets/ca-perovskite-noguchi-2013-table1.LICENSE.md) |
| Cementite Hu/PANGAEA shock dataset | [Source CC BY 4.0](datasets/cementite-hu-2019.LICENSE.md) |
| Cementite Ono–Mibe | [Unspecified source terms; contributor CC0](datasets/cementite-ono-mibe-2010.LICENSE.md) |
| Cementite Sata/Li | [Unspecified source terms; contributor CC0](datasets/cementite-sata-2010.LICENSE.md) |
| Cementite Takahashi, Litasov and Prescher | [Separate source scopes](datasets/cementite-sources.LICENSE.md) |
| Coesite Bykova author inputs | [Maintainer-declared author CC BY 4.0; third-party inputs separate](datasets/coesite-bykova-2018-author-inputs.LICENSE.md) |
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

## Installed packages and downstream redistribution

This guide and the dataset notices are installed under `peritheos/data/`.
The distribution's `License-File` metadata also includes these notices, the
MIT license and complete CC BY 4.0 and CC0 texts. Python package metadata uses
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

For Debian-style copyright inventories, record code and data separately and
retain the individual source notices. An unspecified source license is an
unresolved packaging question, not a finding that the resource meets Debian's
Free Software Guidelines. This guide documents the current evidence; it does
not supply a blanket redistribution clearance for the entire data directory.
