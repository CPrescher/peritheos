# Code, data and redistribution

The unreleased distribution split makes experimental observations optional:
`pip install peritheos` supplies the library, structures and EOS parameters;
`pip install "peritheos[datasets]"` additionally installs `peritheos-datasets`.
The two packages have separate source archives. The optional package contains
451 explicitly selected material dataset entries and two Barker supporting-study
entries, including separately assessed
numerical transcriptions. Source grants and factual-content assessments remain
distinct; unspecified source terms are not relabeled open. Five material entries
remain held: Sakai Table S8 and four views of the original Campbell integration collection. The [numerical redistribution review](research/numerical-dataset-redistribution.md)
records every decision and the collection, source and Debian qualifications.
The [Morard/Campbell follow-up](research/morard-campbell-redistribution-followup.md)
identifies all selected cold-row sources and separately calculated buffer
polynomial approximations; original Campbell integration tables remain held. See the
[installation and data-access guide](datasets.md). The
[Shen/Smith and historical collection review](research/shen-historical-redistribution-review.md)
adds three unchanged numerical resources under source-specific assessments;
original source terms and all scientific qualifications remain explicit. The
[Sakai/Ismailova/iridium review](research/restricted-source-numerical-review.md)
admits 21 numerical extractions on an applicability assessment, preserving all
NC/ND source flags. This concerns general redistribution of the optional data
package; core `peritheos` is independently installable without these values.
The [independent-generation follow-up](research/independent-lineage-buffers.md)
adds a separate partial lineage graph and a relative NNO pressure integral. It
does not establish full Campbell absolute-table regeneration or admit the held
originals. This reduces distributed
source content; it is not a blanket Debian clearance of remaining structural
facts or EOS parameter compilations. Complete research evidence remains in the
repository. The published 0.17.4 archives are unchanged.

The [EOS/structure compilation review](research/eos-structure-compilation-rights.md)
now accounts for all 685 parameter records and 234 structural/diffraction
components, with exact source/runtime bindings. It preserves the pinned Dioptas
MIT notice and separates numerical-content assessments from remaining
collection/acquisition questions: seven Holland/THERMOCALC extracts, ten legacy
peak lists and two ICSD-linked topologies. These 19 components remain present
in the runtime package and need resolution or explicit downstream exclusion
before claiming an entirely cleared archive. No source grant, scientific
reclassification or dataset selection change is made by this review.


The [0.17.3 license inventory](research/license-inventory.md) accounts for all
published archive files and catalog components, including inline scientific
data, original author artifacts and third-party code. Its
[source-by-source dataset index](data/license-inventory/dataset-sources.md)
lists both recorded grants and unresolved rights. The inventory is evidence
for downstream review, not a blanket redistribution clearance.

After 0.17.4, the original Fei source images have been removed from the source
tree and excluded from future distributions. Three plots that embedded those
images have been replaced with plots of preserved numerical readouts and
calculations. Source hashes and scientific evidence remain traceable; the
original artwork is private audit evidence. This correction does not alter
the already published 0.17.4 archives or resolve numerical-data rights.
The original Tateno KCl workbook and redundant Chidester supplement ZIP have
also been moved to private evidence storage and excluded from both package
formats. Their source manifests, checksums and unchanged numerical CSVs remain.

The [CIF license review](research/license-inventory.md#specific-findings-and-actions)
records matched COD files and the current CC0 source declaration for the named
vanadium, superhydrous-phase-B and archival brucite CIFs. He (2023)'s former
thesis structure has been replaced with CC0 COD 1000054. Original He evidence
is preserved privately, without asserting a new source grant.

The [Rust third-party notice audit](research/rust-third-party-notices.md)
documents the complete installed notice bundle, source/archive verification,
platform roles and audited release compiler. License texts, additional source
attributions and the full standard-library copyright document ship in both
package formats and their license metadata. Extra native libraries trigger a
new audit through the installed smoke check.

The [Barker review](research/barker-1987-redistribution.md) admits two unchanged
numerical resources through the supporting-study interface. No open source-data
grant was found; inclusion relies on the scoped numerical-content assessment.
Only newly added contributor-held metadata/annotations receive a CC0 dedication.
The source artwork/PDF, absent crystal-phase identity and unreproduced Monte
Carlo implementation are not supplied or relabeled by the admission.

A [Morard author-data permission template](research/morard-data-permission-template.md)
is available for file-specific confirmation. It is an unsigned draft and does
not extend the existing author-file grant until confirmation is recorded.

Peritheos implementation code and contributor-written research scripts use
the MIT (Expat) license. Bundled scientific data have source-specific terms.
Article citations in scripts identify scientific sources; they do not change
the script's code license. Original author code retained as evidence keeps its
own license.

The [installed licensing guide](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/LICENSING.md)
explains the scopes, indexes every standalone dataset notice, and describes
the different contents of wheels and source distributions. It is also available
from any installed package:

```python
from importlib.resources import files

print(files("peritheos").joinpath("data/LICENSING.md").read_text(encoding="utf-8"))
```

A contributor CC0 notice covers only the rights explicitly held and dedicated
by that contributor. It does not relicense third-party source data or an
associated paper. Some source datasets explicitly use CC BY 4.0 or CC0; others
have unspecified terms. Keep the original source attribution, notices and
descriptions of changes when redistributing data. Missing terms remain
unspecified, including for resources without a standalone license file.

The Python distribution's `LicenseRef-Peritheos-Distribution` metadata refers
to this component-specific scope. Its `License-File` entries include the code
license, installed guide, dataset notices and complete Creative Commons texts.
This documents the distributed components without asserting that all bundled
data have an established open source-data license.

The [CC BY 4.0 legal code](https://creativecommons.org/licenses/by/4.0/legalcode.en)
defines attribution and modification notices and limits the grant to rights
the licensor can grant. The
[Python packaging specification](https://packaging.python.org/en/latest/specifications/core-metadata/#license-expression)
defines license expressions for the containing distribution, rather than only
its implementation code. Downstream distributors should assess the individual
resources against their own policies, including the
[Debian Free Software Guidelines](https://www.debian.org/social_contract#guidelines).
