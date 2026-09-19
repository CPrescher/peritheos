# Incomplete datasets and author requests

This page tracks whether obtaining additional source data would resolve a
documented reproduction gap, exactly what to request, and the outreach status.
The [paper investigation ledger](paper-investigation-ledger.md#papers-with-unavailable-direct-refits)
provides the complete inventory of catalog papers with unavailable direct
refits and their record-level blockers. Its
[withheld/deferred papers](paper-investigation-ledger.md#withheld-or-deferred-papers)
also include candidates awaiting primary evidence.

The list below contains reviewed requests. Papers absent from it have **not
yet been assessed for author outreach**; absence does not mean their data are
complete or contacting the authors would be unhelpful. Some ledger gaps concern
computed grids, adopted coefficients, or model definitions rather than missing
experimental observations.

## Reviewed requests

| Paper / material | Recommendation | What additional data would enable | Outreach status | Last reviewed |
|---|---|---|---|---|
| [Crichton et al. (2016), bcc vanadium](#crichton-et-al-2016-bcc-vanadium) | Contact authors; high priority for original-fit verification | Refit the complete 62-state experiment; resolve the two incomplete fixed-K′ alternatives | Not contacted; historical contact route available | 2026-09-19 |

### Crichton et al. (2016): bcc vanadium

Paper: *High-temperature equation of state of vanadium*,
[DOI 10.1080/08957959.2015.1123256](https://doi.org/10.1080/08957959.2015.1123256).
Record: `vanadium_bcc_crichton_2016_bm3_thermal`.
Evidence: [source audit and reproduction](literature-reproductions/crichton-2016-vanadium.md).

**Why contact is worthwhile.** The published equation reproduces the plotted
curves, but only 29 distinguishable experimental symbols could be digitized,
and their exact temperatures are missing. The nominal-temperature proxy does
not recover all published thermal coefficients. Obtaining the complete input
table would make a source-faithful refit possible; it would not guarantee
coefficient parity. The existing published-coefficient record remains usable
within its documented range.

**Minimum request for the main fit:**

- The original 62 P–V–T rows in their original units and volume normalization,
  with exact row-wise temperatures, run identifiers, and the selected/excluded
  points, including treatment of the initial stressed-foil cycle.
- The EOS-FIT5.2 input/output files or equivalent fit settings: residual
  definition, weighting, fixed/free parameters, and reference-state convention.
- Measurement uncertainties and parameter covariance, if retained. These would
  support statistical checks; their absence should not prevent sharing the
  central observations.

**Useful additional requests:**

- Paired NaCl/Au lattice measurements, the precise calibrant EOS references and
  coefficients, and PTX-Cal settings, to replay the original P–T reduction.
- Complete fitted parameter sets for the K′ = 3.5 and K′ = 4 alternatives,
  especially V0 and the thermal coefficients, to assess the withheld records.
- Clarification of the Conclusion's thermal-expansion percentages, which differ
  from calculations using the printed coefficients.
- A public repository deposit or permission to redistribute the supplied data
  with attribution, so the reproduction can be bundled and independently run.

**Access checked.** The UCL accepted manuscript and publisher abstract were
inspected. Publisher full text/supplement routes, institutional repositories,
correction searches and the local Zotero EOS collection did not yield the
complete numerical dataset. The final typeset full text was inaccessible.
See the audit for the exact routes and limitations.

**Contact route and history.** The accepted manuscript lists W. A. Crichton as
the corresponding author with `crichton@esrf.fr`. This is the paper's historical
address; current deliverability has not been verified. No request has been
sent and no response or permission has been received.

## Maintaining the list

Add one entry per paper when an audit identifies a concrete information gap.
Link the affected records and audit, distinguish essential inputs from optional
improvements, record the access routes already checked, and explain what a
response would unlock. Recommend outreach when the requested material could
resolve that gap. A coefficient mismatch alone is not evidence of missing data;
first inspect the model, units, fit selection and objective.

Use explicit outreach states: **not contacted**, **contacted** (with date),
**response received**, **data received—validation pending**, **resolved**, or
**closed without data** (with reason). Record follow-up dates and sharing terms
when known. Receipt of a file does not establish parity: preserve its provenance,
rerun the reproduction, and update the scientific ledgers before marking the
gap resolved. Keep unresolved secondary requests visible if only part of the
requested data arrives.
