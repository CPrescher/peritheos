# Code, data and redistribution

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
