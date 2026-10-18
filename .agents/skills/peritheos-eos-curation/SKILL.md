---
name: peritheos-eos-curation
description: "Curate or audit Peritheos material/EOS records from primary publications: recover source observations, map equations, reproduce fits, and integrate validated records. Use for adding materials, thermal EOS coverage, pressure-calibration provenance, and published-fit parity investigations in Peritheos. General literature questions without catalog work belong to scientific evidence review."
---

# Peritheos EOS curation

Produce an auditable scientific contribution at the scope requested: a source investigation, staged candidate, independent refit, or integrated catalog addition. A numerical match alone does not establish reproduction of the authors' workflow.

## Establish the current context

Work from the current Peritheos checkout, including its local instructions and changes. Read `docs/adding-materials-and-eos.md` and the relevant parts of `docs/eosmat-schema.md` before authoring a record; these maintained files own acceptance criteria, schema, and contributor checks. Read `docs/material-eos-candidates.md`, `docs/dataset-requests.md`, and the material's existing audit only when relevant. Use existing model families whenever they exactly represent the source formulation.

Distinguish a material (composition and phase), a source parameterization, an independent refit, and a new equation family. Extend the matching material; keep scientifically distinct phases/compositions separate. For audit-only requests, report findings without changing catalog data. For requested additions, carry the work through the authorized integration or explain the specific evidence preventing acceptance.

## Recover and preserve source evidence

- Use the final primary publication, official supplements, and corrections as authority. Record exact table, equation, page, figure, worksheet, or cell locations and source versions/checksums where useful. Secondary sources are leads.
- Recover all available observations relevant to the represented fit. Preserve original tokens/units, missing cells, row/run identities, reported coordinate errors and their meaning, and paired calibrant readings. Keep normalized or recalibrated quantities traceable to the originals.
- Distinguish measured, calculated, digitized, fixed, adopted, and fitted values. A pressure column calculated from the published EOS is not an independent pressure measurement or fit benchmark. Label plot digitization and its precision.
- Missing source licensing alone is not a reason to omit recoverable factual observations. Preserve attribution and explicit source terms; any CC0 dedication covers only contributors' transcription, normalization, metadata, and arrangement rights, never third-party article text, figures, or publisher files. Keep unspecified source terms unspecified.
- If data cannot be recovered, name the files/repositories checked and distinguish inaccessible from absent. Document the precise missing inputs; do not synthesize observations. A fully defined published parameterization can still be usable under the current acceptance guide's independent-benchmark requirements.

## Resolve identity and formulation before fitting

Verify composition, phase, conventional/primitive/formula/atomic/molar volume basis, formula-unit count Z, units, reference pressure/temperature, and determination method. For new materials, supply diffraction-ready structure or an explicitly supported fallback, with structural provenance separate from EOS provenance. Check cell contents and volume consistency.

Map the actual equations, reference-state conventions, fixed/adopted parameters, thermal laws, atom normalization, and pressure-scale version. Record actual experimental coverage and phase limitations; marginal P/T ranges do not establish stability over the entire rectangle. Missing uncertainty is missing, not zero.

## Validate at the level the evidence supports

Keep these claims separate in the report and metadata:

1. **Source-equation validation:** the implemented formulation and coefficients reproduce an independently specified source result within a tolerance justified by printed precision, measurement uncertainty, or digitization error. Reference-state identities alone are insufficient.
2. **Independent refit:** the available observations are fit with an explicitly documented objective, row/pair selection, weights, constraints, and fitting order. Report coefficient and curve/residual differences from the source.
3. **Author-workflow reproduction:** the original input set, selections, calibration, constraints, objective, and relevant fitting details are recovered and replayed. Conditional coefficient agreement does not establish this stronger claim.

Investigate discrepancies without silently changing the published record. Missing runs do not prove they entered the original fit or explain non-parity. Alternative masks, weights, or fixed parameters are diagnostic hypotheses until supported by source evidence. Preserve the literal source replay separately from experiments intended to improve agreement.

### Stop at unresolved non-parity

Apply the [contributor guide's stopping rule](../../../docs/adding-materials-and-eos.md#stopping-rule-unresolved-reproduction-discrepancies). After checking transcription, equation mapping, units/reference state, and the reported fitting protocol, an unresolved failure to reproduce the published EOS or its reported data agreement is a completed curation finding: **document the discrepancy and do not recommend the affected EOS**. Preserve source coefficients, observations, quantitative comparisons, justified tolerances, affected domain, and unresolved inputs. Use `not_reproduced` for a demonstrated unresolved reproduction problem, or `deferred` for incomplete source/model mapping; do not promote it as validated or a recommended default. Keep catalog descriptions and scientific ledgers consistent. Any retained executable diagnostic must be explicitly diagnostic and not recommended for scientific predictions or pressure calibration.

Do not turn failed reproduction into successive alternative EOS families, thermal laws, extra parameters, or fitting protocols to construct a replacement. Existing exploratory refits remain diagnostic evidence; improved residuals or replay of their own coefficients do not establish publication parity or a validated replacement. Replacement-model development requires an explicit user request as a separate research task. A source-protocol refit used to check reproduction remains curation work.

Scope non-parity to the affected parameterization/domain; independently verified isotherms do not validate an unresolved continuous EOS. Missing exact author-workflow details alone do not establish numerical failure if the published parameterization passes independent source benchmarks. Unavailable observations follow the guide's parameterization-only validation route.

When a separate independent refit is explicitly requested and accepted for catalog integration under the guide, store it as `record_kind: refit` with complete `fit_provenance` and the required identifier convention. Preserve source coefficients and distinguish its recommendation status from the original publication's unresolved status. Do not invent a catalog parent for a direct dataset refit.

## Integrate and report

Use material-specific datasets, scripts, tests, and audit notes. Respect concurrent ownership; coordinate shared manifests, ledgers, and totals only within the authorized integration scope. Staging is appropriate when requested or needed to keep a concurrent integration stable, not a reason to leave an otherwise authorized addition unfinished.

Check transcription, independent literature values, schema/public loading, evaluation/inversion, and applicable Rust/Python parity. Follow the current contributor guide for required checks and applicable generated-ledger updates; distinguish new failures from existing ones. Do not equate schema validity or a passing test suite with scientific validation.

Deliver a concise account of source authority, recovered data, equation mapping, reproduction level, results, changes, checks, and remaining evidence gaps. Keep a non-parity finding useful even if exact reproduction is impossible. Draft author questions when requested; sending them requires explicit authorization.
