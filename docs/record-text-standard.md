# Writing standard for record text

This standard applies to the free text of bundled records in
`peritheos/data/materials/<id>.eosmat`. Standalone user files may follow it but
need not. Field layout: [EOSMAT v5 record reference](eosmat-reference.md);
label words and generic caveats: [EOS labels](eos-labels.md).

## 1. Fidelity rule (overrides every other rule here)

Nothing scientific is deleted. A rewrite of `label`, `notes`, `validity.notes`
or `thermal.notes` is applied with `scripts/apply_record_text.py` (a repository
maintenance tool), which keeps the complete original of every changed field
verbatim in the record's `superseded_record_text` block of the audit sidecar
(`peritheos/data/audits/<id>.audit.json`). Content that belongs under a
specific sidecar key (table in section 2) should also be written there.

- Never change a number, a source location a user needs, or a qualification
  that affects how the EOS may be used. A rewrite rounds nothing and adds no
  number that is not in the old text or the structured fields; the apply tool
  refuses one that does.
- If you are not sure a sentence is dispensable, keep it and flag it for the
  curator.
- If a structured field contradicts a sentence (for example
  `parameter_error_confidence: null` while the note says "one standard
  deviation"), keep the sentence and flag the field; do not delete either.

## 2. Fields

| Field | Purpose | Limit |
|---|---|---|
| `label` | Identify the record in a list row | ≤ 60 characters where possible, never over 80 |
| `notes` (record) | What this EOS is and what a user must know that no field says | ≤ 400 characters, ≤ 3 sentences (aim ≤ 250) |
| `validity.notes` | Record-specific facts about the stated range | ≤ 2 entries, ≤ 160 characters each |
| `thermal.notes` | Conventions of the thermal model not expressed by `thermal` fields | ≤ 250 characters |
| material `notes` | What the material/phase is; material-wide caveats | ≤ 300 characters |
| `check.source_note` | What was checked against the source for this record (copy of the sidecar `record_note`, else of the verbatim `finding`) | ≤ 300 characters |
| `check.warning_reason` | Why the record's warning applies (copy of the sidecar `warning_reason`) | one clause, ≤ 140 characters |

`scripts/check_record_text.py` reports every record over these limits (with
the 80-character label limit) and every label that is not unique within its
material; the test suite fails on a case that is not in its baseline file. A
`record_note` within the limit resolves a record's `source_note_length` case.

**`label`.** Pattern `Author (year), model[, distinguishing condition]`.
Author from the stored `reference.authors`: one surname; "A and B" for two;
"A et al." for three or more. Refits start `Peritheos refit of Author (year),
...`. Add a condition only if it distinguishes the record from another record
of the same material or states what the record is (0 K DFT, Hugoniot,
low-spin, pressure scale). Never: status words ("DO NOT USE", "diagnostic
only", "unvalidated", "provisional"); identifiers or snake_case; journal names;
the material name or composition; table or row codes (see section 4). A status
belongs in `check.warning` (generated from coded restrictions). A replaced
label, such as a legacy Dioptas/JCPDS form, stays searchable in
`superseded_record_text.label`.

**`notes` (record).** Belongs: the data and method in one phrase; how the
record differs from sibling records when the label cannot say it; the printed
source value when the stored value is converted (e.g. "Published V0 =
14.72(4) Å³/atom; stored per conventional cell (×4)"); reference-state and
normalization conventions; the reason for a `warning` in one sentence when
`warning` and `why` do not give it; uncertainty semantics that no field can
express. Never: coefficients, errors, fixed flags or ranges that are structured
fields; "confidence/covariance not reported" (see section 5); refit numbers,
RMS values, objectives, weights, solver details; row bookkeeping and counts of
bundled rows; repository paths, dataset or record identifiers; audit dates;
curation history ("retained", "accepted", "former", "was removed");
disclaimers of things nobody claimed ("coefficients are preserved"); any
sentence already in `source_note`, `warning`, `warning_reason` or `why`. Where such content
belongs in the sidecar (keys under `records.<record id>`):

| Content | Sidecar key |
|---|---|
| Refit coefficients, RMS, residuals, objective, weights | `evidence.reproduction.assessment.fit` / `.qualifications` |
| Why a reproduction is incomplete; unpublished inputs | `evidence.scientific_validation.unresolved` |
| Source pages, tables, figures, curves | `evidence.scientific_validation.primary_source_check.locations` |
| What data exist and are bundled | `evidence.scientific_validation.primary_data_check.finding` |
| Where a coefficient comes from, choice among printed values | `evidence.parameter_provenance.<parameter>` |
| Audit dates | `evidence.scientific_validation.audit_date` |
| Curation history, related Peritheos models, the replaced wording itself | `superseded_record_text` (kept verbatim by the apply tool) |

**`validity.notes`.** Only in a record that has a `validity` object; a
`validity` object replaces the experimental ranges as the enforced envelope
([Which range the reader uses](eosmat-reference.md#which-range-the-reader-uses)),
so a note never creates one: in a record without it, range facts go in
`notes`. Only statements about the range: why a bound differs
from the paper's headline number; where data are sparse inside the rectangle
(e.g. "Thermal data cover 92.7-148.9 GPa only"); a usage restriction on the
domain ("Not for the mixed-spin region"). Never: the range itself; the generic
envelope statement (section 5); anything that is not about the range (move it
to `notes`).

**`thermal.notes`.** Which variant of the paper's thermal model, and
integration or reference conventions. Never thermal parameters or their errors.

**Material `notes`.** Never: naming the default (that is `default_records`),
listing the records, or a record-level finding.

**`check.source_note`.** Written per record in the sidecar field
`evidence.scientific_validation.primary_source_check.record_note` (at most 300
characters) with `scripts/apply_check_text.py`, which regenerates the `check`
block; without a `record_note` it is the verbatim `finding`, which is never
edited for this purpose and keeps the complete audit statement. It says what
was checked against the paper for this record and any conversion or
correction that matters to a user (e.g. "Table 2 V0 and K0 match; V0
converted from cm³/mol"). Never: a restatement of coefficients, text about
other records ("The ten stored records are..."), or how to use the record. A
note shared by a set of records is refused unless the curator overrides it
deliberately.

**`check.warning_reason`.** One clause stating the cause of the restriction,
taken from the audit (e.g. the sidecar's `unresolved` or
`recommendation_reason`), written to
`evidence.scientific_validation.warning_reason` with
`scripts/apply_check_text.py` (at most 140 characters). Only for a record with
a `warning`; the `warning` wording itself is generated and never edited.
`record.summary()` prints `Warning: <warning> <reason>`.

**`check.warning`.** Generated, never written. For a bundled record it comes
from a coded `evidence.scientific_validation.usage_recommendation` (for
example `qualified_use` or `not_recommended_for_quantitative_use`; the codes
and their fixed wording are listed in
[EOSMAT reference](eosmat-reference.md#giving-a-bundled-record-a-warning)),
set with `scripts/apply_check_text.py` together with its `warning_reason`.
Free text in `usage_recommendation` never produces a warning; when a code
replaces it, the text is kept verbatim as `recommendation_reason`.

## 3. Notation

| Item | Write | Do not write |
|---|---|---|
| Pressure derivative | `K0'` (prose: `K0' = 4`; label: `K0'=4`) | K0_prime, K0-prime, K', Kprime, K0 prime |
| Other parameters | the key in `parameters`: V0, K0, gamma0, theta0, q0, a_0 | γ0, θ0, Greek in prose |
| Printed error | the source's parenthetical form: `14.75(4) Å³`, `195.0(2.5) GPa` | `+/-`, `±`, "plus or minus" |
| Error without a value | `error 0.4` | `+/-0.4` |
| Volume | `Å³`, `Å³/atom`, `cm³/mol`; density `g/cm³` | A^3, A3, angstrom^3, cm^3/mol |
| Ranges | `0-80 GPa`, `300-2500 K` (unit once) | `0 GPa to 80GPa` |
| Temperature | `300 K`; `23 °C` only when quoting | degC |
| Models | BM2, BM3, BM4, Vinet, MGD; in prose "Mie-Grüneisen-Debye" | Gruneisen, VINET |
| Spin state | prose "high-spin", "low-spin"; labels may use HS/LS | hs, ls |
| Source parts | Table 2, Table S1, Figure 3, Equation (B2), Section 5.4 | Tab., Fig., Eq. |
| Software | with a noun: "MINUTI program", "EosFit7 program" | bare code names |

## 4. Table, fit and internal codes

A code from the source (Fit3, Table S1 row 7, type4_1, P1-P7) or from
Peritheos (snake_case, identifiers, "source-owned", "fp20") never appears
unexplained.

- In a `label`, replace the code by its plain meaning, taken from structured
  fields or the record's own notes: the pressure scale, the fixed value, the
  branch (e.g. "LS BM3, K0'=4.5, V0(LS/HS)=0.99"). If no plain meaning is
  available, keep the source's code in the source's spelling and flag the
  label for the curator.
- In `notes`, the source's own name may appear once, as a locator: "(Fit3 in
  the paper)".
- Peritheos identifiers and paths never appear in record text.

## 5. Generic statements

A caveat that holds for every record is stated once and not repeated in
records:

- A stated range is the published fit or observation envelope, not a
  phase-stability field or an extrapolation guarantee
  ([Which range the reader uses](eosmat-reference.md#which-range-the-reader-uses)).
- An absent or `null` `parameter_error_confidence` or `parameter_covariance`
  means the source did not report it; errors without a stated confidence are
  used as 1 sigma ([Parameters](eosmat-reference.md#parameters)). Retire a record's sentence
  saying so only where both fields are absent or `null`; if either is set, the
  sentence and the field must be reconciled by the curator first.
- A published record keeps its published values and qualifications
  ([record origin](eos-labels.md#record-origin)).

Sentences these definitions retire, verbatim:

| Sentence | Defined once in |
|---|---|
| Observation envelope; extensions to inner-core conditions are extrapolations. | range statement |
| Published fit/data envelope, not a phase-stability or extrapolation guarantee. | range statement |
| Published Table I data envelope, not a phase-stability or extrapolation guarantee. | range statement |
| Bounds describe fitted observations, not a phase-stability field. | range statement |
| Printed ± confidence and covariance are not specified. | error statement |
| Parameter covariance and the confidence convention for printed ± values are not reported. | error statement |
| The confidence convention and covariance of published errors are not specified. | error statement |
| Weights, covariance, and confidence of parameter errors are not specified. | error statement (weights belong in the sidecar fit description) |
| The published coefficients and errors are preserved. | published record |
| Measured pressure span at this temperature; no temperature interpolation by an isothermal record. | `equation_kind: isothermal` |

A sentence that only resembles one of these may carry record-specific content;
review it sentence by sentence.

## 6. Applying a rewrite

1. Write the reviewed text as a JSON file of `before`/`after` values per record
   (format in the docstring of `scripts/apply_record_text.py`), and labels as
   rows of a CSV flagged `mechanical`.
2. Run `python scripts/apply_record_text.py FILE.json --labels LABELS.csv
   --check`. Every refusal is listed: a `before` text that differs from the
   file, a validity note that would create or remove the record's `validity`
   object, or a failed guard (a new number, a length limit, a duplicate label,
   a repository path). Fix the text; use `--allow RECORD:GUARD` only for a stated
   false positive.
3. Run it without `--check`, then `--verify`, and regenerate the generated
   pages.

`source_note`, `warning_reason` and a coded `usage_recommendation` are
written the same way with `scripts/apply_check_text.py` from a JSON file
`{"<material>": {"<record>": {"source_note": ..., "warning_reason": ...,
"usage_recommendation": ...}}}` (`null` removes a text; for `source_note` the
verbatim `finding` is shown again; for `usage_recommendation` the stored
original returns). A `usage_recommendation` must be one of the codes with a
fixed warning, and a change that gives a record its first warning needs a
`warning_reason` in the same file; the dry run prints each record's resulting
`check.warning`. Its guards are the two length limits, `numbers` (every number, sign
included, must occur in the record's finding, other sidecar evidence, record
fields or stored originals, or be a correct rounding of one at the precision
written with at least two significant digits, such as 23.89 for 23.89043),
`forbidden_text` and `shared_note`; override one guard for one record with
`--allow RECORD:GUARD` only for a stated false positive. A `warning_reason`
for a record without a warning is always refused. Use `--check` first, then
apply and run `--verify`; either tool writes nothing unless the whole run is
valid.

## 7. Checklist

1. Label matches the pattern, ≤ 60 characters where possible, unique within the
   material, no status, code or identifier.
2. Notes ≤ 400 characters; each sentence says something no field,
   `source_note`, `warning`, `warning_reason` or `why` says.
3. Every removed sentence is kept in `superseded_record_text`, and content with
   a sidecar key in section 2 is also written there.
4. Numbers in the new text all occur in the old text or the structured fields;
   nothing rounded.
5. Validity notes are about the range only; generic statements (section 5) are
   not repeated.
6. Notation table followed; no unexplained codes.
7. Doubtful removals kept and flagged; field/text contradictions flagged, not
   resolved silently.
