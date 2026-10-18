# EOSMAT structure and usability review

EOSMAT can be made substantially easier to read and maintain without reducing its scientific quality. The main obstacle is that the current material document serves as an executable EOS definition, crystal structure, dataset catalogue, source audit, and investigation notebook simultaneously. A simpler core with concise scientific qualifications and access to complete evidence would serve people better.

The best starting point is the existing JSON format and one material document per distinct material or phase. Replacing JSON with YAML or splitting every EOS into another file would leave most of the underlying complexity intact and introduce additional maintenance work.

This review covers the local working files on 17 October 2026, including existing changes. It proposes a design; the production schema, readers and catalogue were not changed by this review.

**Measured complexity**

| Measure | Current local catalogue |
| --- | ---: |
| Material documents | 232 |
| EOS records | 663 |
| Median material file | 390 lines |
| Gold | 7,550 lines, 23 records |
| MgO | 6,827 lines, 37 records |
| Iron | 8,549 lines, 38 records |
| Dataset declarations | 423 |
| Datasets stored as external resources | 414 |
| Distinct field names at EOS record level | 83 |
| Fields explicitly described by the EOS record schema | 35 |
| Distinct field names inside scientific_validation | 68 |

About 91% of records use BM2, BM3 or Vinet as their primary equation. The common authoring task is much simpler than the full schema suggests. Thermal models and the three current Hugoniot records need additional information, but those requirements need not dominate the introduction or a basic material template.

The primary and thermal equation objects together account for about 4% of compactly serialized EOS-record field values. Scientific validation alone accounts for approximately 36%. This is a measure of text volume, not a claim that the remaining information is dispensable: uncertainties, reference states, calibration and scientific qualifications are essential. It does show that changing punctuation or shortening parameter names cannot solve the main readability problem.

The measurements and exact definitions are in [audit.json](audit.json).

**What is necessary**

A material definition needs an unambiguous identity and phase, explicit units and volume basis, a reference state, model identifiers and coefficients, parameter uncertainties and fixed coefficients, source attribution, and an honest description of its supported domain and scientific status. Distinct source parameterizations need stable identifiers. Published records and independent refits must remain distinguishable.

Conventional-cell volume, atomic volume and molar volume cannot be interchanged by visual shorthand. Thermal models need the appropriate cell contents and explicit conventions for thermal pressure, reference-temperature subtraction and zero-point terms. Shock records need their initial state, mass basis, loading path and branch domain. These are conditional requirements with numerical consequences.

Dataset declarations need their units, uncertainty meanings, selections, resource hashes and links to the records they support. Original observations and transformed views must remain distinguishable. A rectangular pressure–temperature envelope does not establish a measured joint domain or permission for extrapolation.

A reader must also retain enough information to prevent a deferred, diagnostic or explicitly restricted record from becoming a default merely because its detailed audit is elsewhere. A source parameter check, numerical curve reproduction and reproduction of the authors’ fitting workflow are different claims.

**What is avoidably complicated**

| Current pattern | Why it is difficult | Recommended treatment |
| --- | --- | --- |
| Both eos.type and eos.model identify the same equation | Authors must supply a redundant pair and readers check their consistency | Choose one stable model identifier in a future version; generate the legacy spelling for v3 export |
| Experimental ranges and validity ranges overlap | The reader must understand their relationship and precedence | Distinguish observation coverage from permitted evaluation domain; store each fact once |
| Thermal choices can occur in configuration and at established component locations | The same choice can acquire multiple representations | One authoritative configuration location; legacy export handles duplication |
| Long audit findings accompany every parameterization | The numerical definition is difficult to locate and compare | Concise status and critical qualifications in the main card; detailed investigation in linked evidence |
| Many unrelated extensions occur directly on records | Readers cannot tell which fields affect execution | Document a small core and a deliberate evidence/extension section |
| Explicit null covariance/confidence and empty optional collections are common | Boilerplate obscures the available information | Define omission as unknown or absent for these optional fields; preserve meaningful nulls |
| Model dictionaries, parameter arrays and short ranges expand across many lines | Related numerical information is visually scattered | Format short numerical groups together |
| Material identity and format markers sometimes follow large record arrays | A reader encounters detailed content before knowing the document identity | Always put identity, units and selection information first |

There are 337 records with identical pressure ranges in both validity.pressure_gpa and experimental_pressure_range_gpa. These are candidates for normalization, not grounds for deleting every validity object. The two concepts can differ, and the Python loader currently gives the validity object precedence when present.

Unknown-field preservation is valuable for interchange, but it does not make an unlimited collection of case-specific fields a good human interface. The 83 record field names versus 35 schema properties show how much meaning currently lives outside the documented core. Such extensions can remain available without requiring every reader to study them.

The primary-source audit ledger repeats information from the material cards. The inspected audit script can regenerate ledger entries from the cards, so this duplication can be a generated index with a clear owner. Merely deleting the ledger would not resolve ownership elsewhere in the curation workflow. Each scientific fact should have one authoritative authoring location, with summaries and indexes derived from it.

**A simpler file organization**

Keep one material file per distinct material or phase. For a simple material, a single file should be sufficient. For a heavily investigated material, use one optional evidence companion for the material, rather than a separate file for every coefficient, source and validation step.

A possible arrangement is:

```text
materials/
  gold.eosmat
  gold.evidence.json
  coesite.eosmat

datasets/
  gold-fei-2007.csv
  ...
```

The main card would contain identity and units, a compact optional structure, EOS records, source citations, uncertainties, calibration identity, supported domain, default selection and concise scientific qualifications. The companion would carry detailed parameter-source locations, digitization procedures, fit replay reports, historical corrections and investigation findings. Existing large observation tables would remain CSV resources.

This is a proposed future capability. EOSMAT v3 does not currently resolve a general evidence companion. The companion should be optional for numerical evaluation, but retained and resolvable for scientific inspection. Essential restrictions must remain in the main card, and missing evidence must remain visible to an audit-aware application. Full exports must embed the evidence or include the referenced files and checksums in a documented portable bundle. Existing package-relative dataset paths already mean that a material file is not always self-contained.

Avoid making two manually maintained representations of the same record. A future simplified source can generate legacy v3 exports and catalogue indexes. A generated human material sheet should also come from that source. Editing a published record in a spreadsheet and separately editing its EOSMAT JSON would recreate the synchronization problem.

**What an everyday EOS record should expose**

A human should be able to locate the following information together:

```text
SrSiO3 6H perovskite — Xiao et al. (2013)
Source: 10.2138/am.2013.4470
Origin: theoretical; published parameterization
Equation: BM3 with K0_prime fixed to 4
Reference temperature: 0 K
V0 = 311.23 A³ per six-formula-unit cell; reported error 0.10 A³
K0 = 183.8 GPa; reported error 0.6 GPa
K0_prime = 4.0, fixed; no reported numerical error
Reported calculation pressure span: 0–40 GPa
Calibration: not applicable; first-principles calculation
Source check: validated
Underlying numerical states: plot only
```

These are values from an existing record, not a new scientific fit. The current identifier contains “bm2” while its executable model is BM3 with the derivative fixed to 4; the summary deliberately exposes the stored executable model. Simplification must preserve identifiers and aliases rather than silently rename historical references.

An EOSMAT representation should put these facts near one another, with detailed evidence available by following the record’s evidence entry. It should not require reading a digitization history to find K0. A parameter comparison table generated from the cards would be particularly useful for materials with dozens of EOS choices.

For hand editing, retain familiar parameter names and short arrays. Keep a single model discriminator and consistent parameter-error conventions in a future version. A new system of elaborate objects around every coefficient would add nesting; it needs a demonstrated advantage before adoption.

**Changes that can be made without a new format**

Improve key ordering and formatting first. Both demonstration files beside this review retain every original JSON value. Gold falls from 7,550 to 4,191 lines; the small SrSiO3 material falls from 129 to 79 lines. Across all 232 current documents, the same formatting rules reduce 167,537 lines to 95,322, approximately 43%. Fewer lines alone do not prove better usability, but grouping short arrays and related parameters makes comparisons easier without a migration.

The present writer uses uniform JSON indentation, which cannot express this grouping. A human-oriented writer can, although writing through an unchanged application may restore its own formatting. This should eventually be a common writer behavior rather than a one-off catalogue reformat.

Also provide a genuinely small authoring example and generated material summaries. The introduction should teach an ordinary isothermal record first, then thermal components, datasets, shock records and evidence. Remove optional empty boilerplate only after defining the relevant omission semantics. An absent reported error must remain unknown, and an absent covariance must never imply independence. Some null values, including certain thermal-reference conventions, have model-specific meaning and cannot be stripped generically.

**Changes that require a coordinated format revision**

Removing type while retaining model is rejected by both the current Python validator and JSON Schema. Moving active metadata, changing selection rules or introducing evidence-file resolution also requires coordinated reader work. Treat these as an explicit versioned migration, not silently different v3 content.

The Python and Rust readers both participate in validation and executable model construction. Dioptas reads the shared material documents and saves through Peritheos. Studio imports and exports JSON and evaluates through a pinned Rust engine. Its export digest uses canonicalized JSON, so formatting alone is different from removing or relocating data. Other raw-file checksums and generated observation metadata must still be checked before a bulk rewrite.

A new version should have a small core with conditional sections for advanced models. Keep separate answers to whether parameters were checked, whether numerical results or fits were reproduced, and whether ordinary use is recommended. Older records need conservative explicit mappings, especially where an earlier broad validation label does not establish fit reproduction. Avoid inventing a new all-purpose “validated” flag.

A single internal model registry can supply type/model mappings and reduce repeated implementation lists where practical. Shared accepted/rejected fixtures should check schema, Python and Rust behavior. JSON Schema remains responsible for shape; executable readers still need checks for references, unit conversions, mass consistency and scientific selection rules. Replacing the current validators wholesale with shape validation would lose those checks.

**Recommended sequence**

1. Establish a human-oriented v3 writer, small templates and generated material summaries. Review the gold and small-material prototypes with real readers before bulk formatting.
2. Catalogue core execution and selection fields separately from evidence fields. State authoritative locations, range semantics and omission behavior. Keep CSV observations and their provenance intact.
3. Trial the proposed evidence companion and single-model representation on a simple isotherm, a complex thermal standard, a restricted or deferred record and a Hugoniot. Use one authored representation and a v3 exporter.
4. Require numerical parity, unit-conversion parity, unchanged default/hold behavior, complete evidence resolution and unknown-field preservation before migrating the catalogue. Update Peritheos, Dioptas and Studio together where their responsibilities require it.

This sequence delivers visible improvement before committing to a new exchange format. It also tests whether moving detailed evidence actually improves editing and review enough to justify companion files. A YAML authoring layer should remain a later option if comments and manual editing prove essential; introducing a second parser now would not address the main source of complexity.

**Verification and source pointers**

All 232 material documents passed both current Python structural validation and the distributed JSON Schema during this review. Every formatting experiment decoded to exactly the original JSON value. The existing EOSMAT test module passed 143 tests. A separate compact-v3 experiment omitted optional crystallographic metadata and null covariance/confidence for the small theoretical SrSiO3 record; both validators accepted it, and its pressures matched the original at V/V0 = 1.00, 0.95 and 0.90 to an absolute tolerance of 1e-12 GPa. That projection was a feasibility check, not a lossless archival export.

The future companion-file proposal has not been implemented or tested for migration parity. No claim of a measured authoring-time improvement is made.

Owning sources inspected: peritheos/data/eosmat-v3.schema.json; peritheos/eosmat.py; peritheos/materials.py; peritheos/datasets.py; crates/peritheos/src/eosmat.rs; docs/eosmat-schema.md; docs/datasets.md; tests/test_eosmat.py; scripts/apply_primary_source_audit.py; the local Dioptas EOS database/builder adapters; and Studio’s JSON canonicalization and refit export code.

Review artifacts: [gold with identical information](gold-same-information.eosmat), [SrSiO3 with identical information](srsio3_6h_perovskite-same-information.eosmat), and [measurements](audit.json).
