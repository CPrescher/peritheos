---
name: peritheos-eos-curation
description: "Curate or audit Peritheos material/EOS records from primary publications: recover source observations, map equations, reproduce fits, and integrate validated records. Use for adding materials, thermal EOS coverage, pressure-calibration provenance, and published-fit parity investigations in Peritheos. General literature questions without catalog work belong to scientific evidence review."
---

# Peritheos EOS curation

Produce an auditable scientific contribution at the scope requested: a source investigation, staged candidate, independent refit, or integrated catalog addition. A numerical match alone does not establish reproduction of the authors' workflow.

## Establish the current context

Work from the current Peritheos checkout, including its local instructions and changes. Read `docs/adding-materials-and-eos.md`, `docs/eos-labels.md` and the relevant parts of `docs/eosmat-v5.md` before authoring a record; these maintained files own acceptance criteria, labels and fields. Read `docs/material-eos-candidates.md`, `docs/dataset-requests.md`, and the material's existing audit only when relevant. Use existing model families whenever they exactly represent the source formulation.

Distinguish a material (composition and phase), a source parameterization, an independent refit, and a new equation family. Extend the matching material; keep scientifically distinct phases/compositions separate. For audit-only requests, report findings without changing catalog data. For requested additions, carry the work through the authorized integration or explain the specific evidence preventing acceptance.

## Recover and preserve source evidence

- Use the final primary publication, official supplements, and corrections as authority. Record exact table, equation, page, figure, worksheet, or cell locations and source versions/checksums in the audit sidecar. Secondary sources are leads.
- Recover all available observations relevant to the represented fit. Preserve original tokens/units, missing cells, row/run identities, reported coordinate errors and their meaning, and paired calibrant readings. Keep normalized or recalibrated quantities traceable to the originals.
- Distinguish measured, calculated, digitized, fixed, adopted, and fitted values. A pressure column calculated from the published EOS is not an independent pressure measurement or fit benchmark. Label plot digitization and its precision.
- Missing source licensing alone is not a reason to omit recoverable factual observations. Preserve attribution and explicit source terms; any CC0 dedication covers only contributors' transcription, normalization, metadata, and arrangement rights, never third-party article text, figures, or publisher files. Keep unspecified source terms unspecified.
- If data cannot be recovered, name the files/repositories checked and distinguish inaccessible from absent. Do not synthesize observations. A fully defined published parameterization can still be usable under the guide's independent-benchmark requirements.

## Resolve identity and formulation before fitting

Verify composition, phase, conventional/primitive/formula/atomic/molar volume basis, formula-unit count Z, units, reference pressure/temperature, and determination method. For new materials, supply diffraction-ready structure or an explicitly supported fallback, with structural provenance separate from EOS provenance. Check cell contents and volume consistency.

Map the actual equations, reference-state conventions, fixed/adopted parameters, thermal laws, atom normalization, and pressure-scale version. Record actual experimental coverage and phase limitations; marginal P/T ranges do not establish stability over the entire rectangle. Missing uncertainty is missing, not zero. When the source states the error convention, set `parameter_error_confidence` (`0.6826894921370859` for 1 sigma, `0.95` for 95%) and keep the printed errors unchanged.

## Validate at the level the evidence supports

Keep three claims separate in the audit: **source-equation validation** (the implemented formulation reproduces an independently specified source result within a justified tolerance; reference-state identities alone are insufficient), **independent refit** (available observations fit with a documented objective, selection, weights, constraints and order), and **author-workflow reproduction** (the original inputs, selections, calibration and protocol replayed). Conditional coefficient agreement does not establish the last.

Investigate discrepancies without changing the published record. Missing runs do not prove they entered the original fit. Alternative masks, weights, or fixed parameters are diagnostic hypotheses until supported by source evidence; keep the literal source replay separate from experiments to improve agreement.

### Stop at unresolved non-parity

Apply the [contributor guide's stopping rule](../../../docs/adding-materials-and-eos.md#stopping-rule-unresolved-reproduction-discrepancies). After checking transcription, equation mapping, units/reference state, and the reported fitting protocol, an unresolved failure to reproduce the published EOS or its reported data agreement is a completed finding: headline `disagrees`, source-validation `not_reproduced` (or `deferred` for incomplete mapping), and the evidence in the sidecar. Do not build a replacement from successive EOS families, thermal laws, extra parameters or protocols; that needs an explicit user request as separate research. Scope non-parity to the affected parameterization/domain. When a separate refit is requested and accepted, store it as `record_kind: refit` with complete `fit_provenance`; do not invent a catalog parent for a direct dataset refit.

## Write the record and the audit

Set the headline from the outcome code as defined in `docs/eos-labels.md`: `parity` → `reproduced`, `similar` → `close`, `parity_not_achieved`/`not_reproduced` → `disagrees`, otherwise `untested`. `untested` is not a failure. `check.compared` (what the headline refers to), `why`, `source_check` and `source_note` follow mechanically from the assessment and source-validation status (`docs/eosmat-v5.md`); the split command derives them. Add a `warning` only for an explicit usage restriction the reader must see, for example a reproduced fit that describes its data poorly. Write the full source-check result in `primary_source_check.finding` (sidecar). When it is long, restates coefficients or is shared by several records, add a per-record `primary_source_check.record_note` (≤ 300 characters: what was checked against the paper for this record and any conversion or correction a user needs; never coefficients or other records); a restricted record gets a one-clause `scientific_validation.warning_reason` (≤ 140 characters) taken from the audit. Both are copied into `check.source_note` / `check.warning_reason`; set them on split materials with `scripts/apply_check_text.py` (`--check`, apply, `--verify`), never by editing `check`. Do not change a default (the material's `default_records`); propose it to the maintainer.

The record holds what a reader needs to use the EOS; the per-material sidecar `peritheos/data/audits/<id>.audit.json` holds scope, basis, conditions, fixed parameters, refit numbers and source-check locations. Field layout is in `docs/eosmat-v5.md`.

- State each finding once. Do not repeat it across basis, finding, notes and qualifications.
- Do not restate a label's generic caveat; `docs/eos-labels.md` holds it.
- No solver diagnostics or full-precision refit numbers in the record.
- No hard-coded catalog counts in prose, and no sections that only say nothing is recorded.
- Prefer a short specific statement over defensive qualification.
- Labels, notes, validity notes and thermal notes follow `docs/record-text-standard.md`. Rewrite existing record text only with `scripts/apply_record_text.py`, which keeps every replaced original in the sidecar's `superseded_record_text`.

Bad (record): "Parity achieved. This does not establish experimental accuracy, phase stability or author-workflow reproduction. Parity here means … Refit V0=263.929881998 …", with the same sentence again in notes.
Good (record): `reproduction: reproduced`, `warning: null`. Good (sidecar entry): scope `parameter_agreement`; basis "within printed source widths; K0 fixed at 210 GPa"; refit coefficients with errors; source "Table 2, pp. 64-65".

## Integrate and report

Use material-specific datasets, scripts, tests, and audit notes. Respect concurrent ownership; coordinate shared manifests and ledgers only within the authorized scope. Staging is appropriate when requested or needed to keep a concurrent integration stable.

Check transcription, independent literature values, schema/public loading, evaluation/inversion, and applicable Rust/Python parity, following the contributor guide; distinguish new failures from existing ones. Do not equate schema validity or a passing test suite with scientific validation.

Report briefly: source authority, recovered data, equation mapping, headline and any warning, changes, checks, and remaining gaps. Draft author questions when requested; sending them requires explicit authorization.
