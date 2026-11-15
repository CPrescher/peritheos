# EOS labels and reconstruction status

A label describes **one EOS record and one check**. It does not rate a material,
publication or library: one material can hold a published EOS, an independent
refit and a diagnostic record with different labels.

The caveats on this page apply to every record that carries a label. They are
stated here once; individual records and reports do not repeat them.

## Headline: reproduction

Each curated EOS record has a `check` block whose `reproduction` is one word:

| Word | Meaning | What it does not establish |
|---|---|---|
| **reproduced** | A refit of the original data gave the published coefficients under the stated criterion. | Experimental accuracy, phase stability, recommended use or replay of the authors' full workflow. |
| **close** | A refit gave similar coefficients but did not meet the full criterion. | Statistical equivalence or exact coefficient recovery. |
| **disagrees** | A refit was attempted and did not match. | The cause, an error in the paper, or permission to replace its coefficients. |
| **untested** | There are no data to refit, or only the source equation was rebuilt. | Failure. `untested` is not a negative result and must never be shown as one. |

## Warning

`warning` is null for almost every record. When set, it states an explicit usage
restriction, such as "do not use for pressure calibration". A positive headline
never hides a warning. A reproduced fit can still describe its data poorly;
reproducing a poor fit is still reproduction. Fedotenko (2020) Pd with fixed
K'=4 is such a case; see its
[audit](literature-reproductions/palladium-fedotenko-2020-audit.md).

A default record is `reproduced` where the material has one, else `close`, else
`untested`. A record that `disagrees` or carries a warning is never a default.
Default selection is a maintainer decision. A default is not globally best,
valid at arbitrary conditions or experimentally proven.

The exact field layout is defined in [EOSMAT v5](eosmat-v5.md). Standalone user
files need no `check` block; loading a file never gives it a positive label.

## Detailed outcome codes

The per-material audit sidecar (`peritheos/data/audits/<id>.audit.json`), the
[reconstruction ledger](primary-eos-refits.md) and older files keep the detailed
codes. Keep them unchanged when exchanging existing records.

| Code | Headline | Meaning |
|---|---|---|
| `parity` | reproduced | Every free coefficient, or the identified source calculation, met the stated criterion. |
| `similar` | close | Every compared free coefficient passes the similarity screen below, but parity is not met. |
| `parity_not_achieved` | disagrees | Ledger code: a reconstruction failed its criterion, or a dedicated audit found an unresolved published-coefficient mismatch. |
| `not_reproduced` | disagrees | Native EOSMAT code: the named reconstruction has a documented unresolved problem. It applies to that reconstruction, not to every EOS of the material. |
| `not_refittable` | untested | Observations or reduction inputs are insufficient for a full refit. A subset or equation check may exist. |
| `source_reconstruction` | untested | A source calculation was rebuilt without refitting the full coefficient set. Calculated curves are not measured observations. |
| `not_assessed` | untested | No result for this check. |

## What "parity" compares

For the generic fit comparison, **every free coefficient** must satisfy both the
combined-error comparison and the numerical similarity screen:

\[
|p_{\mathrm{refit}}-p_{\mathrm{source}}|
\leq 2\sqrt{s_{\mathrm{source}}^2+s_{\mathrm{refit}}^2}.
\]

The error values are the comparison quantities recorded by the audit. Source
confidence convention, error model, covariance, weighting, fixed constraints and
row selection still matter; the rule is not a universal 95% confidence claim.
Fixed/adopted coefficients do not become independently recovered coefficients
because a comparison keeps them fixed.

| Coefficient | Numerical similarity limit |
|---|---|
| `V0` | 5% relative difference |
| `K0` | 15% relative difference |
| `K0_prime` | 1.0 absolute **or** 20% relative difference |
| `gamma0`, `q`, `theta0` | 25% relative difference |
| `alpha...`, `dK_dT` | 30% relative difference |
| Other coefficients | 20% relative difference |

**Similar** means each free coefficient passes either the combined-error
comparison or its similarity limit; **parity** requires both. A dedicated source
audit may use a stricter or explicitly documented alternative criterion, such as
agreement within printed source error widths; its recorded basis takes
precedence. Agreement within printed widths without a reconstruction covariance,
or at printed precision, is not a formal combined-error result.

## Qualified results and comparison scope

The audit records what was compared:

| Scope | What was checked |
|---|---|
| `parameter_agreement` | Coefficients agree under the recorded constraints and criterion. |
| `author_workflow` | Replay of the author's inputs, selections, calibration and fitting procedure. A parameter match alone is insufficient. |
| `independent_refit` | A separately defined Peritheos analysis reproduces under its own inputs and objective. This is not source-author parity. |
| `source_equation` | The source equation or calculated checkpoints were reproduced; this is not an independent data refit. |

Qualifiers such as **conditional parity** or **within reported parameter errors**
narrow a result; they do not remove its limits. Matching central values does not
reproduce uncertainties, calibration or weighting. A cold-branch reconstruction
does not establish a thermal fit. Reproducing a stored refit does not reproduce
the original paper.

[Reconstruction reports](reconstruction-reports/index.md) show each comparison
and coefficient graph; a report is not source-author or physical validation.
Older v4 assessment fields are described in
[EOSMAT v4](eosmat-v4.md#a-consistent-reconstruction-assessment).

## Source-validation status

Source review is a separate axis from reconstruction:

| Status | Meaning |
|---|---|
| `primary_source_validated` | Equation, coefficients, conventions and qualifications were traced to the primary source; for a refit, its inputs and analysis are traceable. Not fit parity or physical accuracy. |
| `not_reproduced` | A documented unresolved reconstruction problem is retained with its affected scope and use restrictions. Executable diagnostic access is not a recommendation. |
| `deferred` | Source or model mapping is unresolved; the record is kept as evidence and excluded from ordinary executable loading. |
| `pending_primary_source_check` | An unreviewed user or staged record. Structural validity and numerical construction do not change this. |

See the [validation guide](validation.md) and
[contribution requirements](adding-materials-and-eos.md) for eligibility. An
explicit **do not use** restriction remains binding whatever the other labels say.

## Record origin

| `record_kind` | Meaning |
|---|---|
| `published` | A literature parameterization with its published values and qualifications. |
| `refit` | Independently fitted coefficients with their own observations, objective, constraints and provenance. |
| `derived` | A declared construction from other equations or source relationships, with its own provenance and checks. |
| `diagnostic` | An inspection, teaching or exploratory record; not a recommendation for predictions or pressure calibration. |

A published record is not automatically validated, a refit is not automatically
more accurate, and a diagnostic result never replaces the published record.
