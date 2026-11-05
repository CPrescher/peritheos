# Optional Peritheos datasets

This package adds selected observation datasets to the Peritheos library.
Install `peritheos[datasets]` to request the matching data package. Installing
`peritheos` alone supplies structures, EOS parameters and scientific provenance
without experimental observations. Neither installation downloads data at runtime.

The explicit policy selects 451 material dataset entries plus two Barker
supporting-study entries (453 total). The material selection comprises: 44 with previously
recorded source/author terms, 363 admitted after a numerical factual-content
assessment, one newly verified CC BY 4.0 deposit extraction, three author datasets
with separately assessed numerical literature inputs, and four separately named
Campbell polynomial-evaluation entries. This follow-up also admits ten views of
the shared Shen/Smith numerical collection and two fixed historical excerpts.
Another 21 existing Sakai/Ismailova/iridium numerical entries are admitted under
source-specific applicability assessments, with all NC/ND source notices retained.
Five material entries remain held: Sakai Table S8 literature lineage and four
views of the original Campbell integration collection. These decisions concern
general redistribution, independently of whether a consumer needs datasets.
Three separate new declarations expose a 37-edge primary-resource lineage graph
and two views of a 1717-row calculated relative NNO pressure integral. These are
partial reconstructions: author selections, absolute Gibbs baselines and the
full Campbell integration tables are not recovered. Original identifiers are
not repurposed; the installed notice explains the source comparison and limits.
The three collection-follow-up resources retain their original unspecified or
all-rights-reserved terms; their installed scoped notice explains the numerical
and collection assessments, source-access limits and downstream qualifications.
The new Campbell entries are calculated approximations, not the original
integration tables; their source mismatch and scientific limitations remain explicit.
Barker Table I (seven liquid comparison rows) and selected Figure 5 readouts
(21 theoretical curve samples, 11 Ross experimental marks) are admitted after
a separate numerical-content assessment. Use `get_study_dataset()`; no crystal
phase or executable Barker EOS is inferred. Its source publication copyright
and unknown data-grant status remain explicit in the installed notice and review.

Factual-content admission is a maintainer assessment, not a new source license,
legal opinion or Debian approval. Unspecified source terms remain unspecified.
Contributor CC0 dedications do not override third-party rights. The installed
`peritheos_datasets/redistribution-review.json` identifies the exact decisions,
source references, limitations and held entries. The accompanying numerical
transcription notice explains copyright, collection and contractual caveats.
`manifest.json` records the policy binding, declaration fingerprints and original
resource checksums. Numerical values, error conventions, attribution, original
source-license fields and scientific qualifications are preserved.

The data package version must match the library version. Dataset access through
Peritheos verifies resource checksums. Other source evidence is available through
the cited original deposits and links from the Peritheos material records.

Code is MIT; data retain their individual terms. Full licenses and scoped notices
are included in this distribution. A contributor dedication does not relicense
third-party sources. Please cite the original authors and report transformations.

The separate `eos-structure-rights-review.json` and compilation notice document
the runtime package's EOS/structure assessment and unresolved collection cases.
Their presence does not add observations or extend this dataset selection.
The pinned Dioptas upstream MIT attribution is also preserved; it cannot clear
third-party publications or inherited beamline-card acquisition rights.
