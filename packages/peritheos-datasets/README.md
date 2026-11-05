# Optional Peritheos datasets

This package adds selected observation datasets to the Peritheos library.
Install `peritheos[data]` to request the matching data package. Installing
`peritheos` alone supplies structures, EOS parameters and scientific provenance
without experimental observations. Neither installation downloads data at runtime.

Only the dataset declarations explicitly listed in `datasets.json` in the
source repository are admitted. They have recorded open or author-declared
source terms. Mixed author/literature, unresolved and NC/ND sources are not
included. This is not a clearance of every dataset in the research repository.
The installed `peritheos_datasets/manifest.json` records the exact selection,
declaration fingerprints and original resource checksums. All numerical source
values, error conventions, attribution and data licenses are preserved.

The data package version must match the library version. Dataset access through
Peritheos verifies resource checksums. Other source evidence is available through
the cited original deposits and links from the Peritheos material records.

Code is MIT; data retain their individual terms. Full licenses and scoped notices
are included in this distribution. A contributor dedication does not relicense
third-party sources. Please cite the original authors and report transformations.
