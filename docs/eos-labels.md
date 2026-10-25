# EOS labels and reconstruction status

These labels describe an **individual EOS record and a particular check**.
They are not a rating of an entire material, publication or library. One material
can have a published EOS, an independent refit and a diagnostic reconstruction
with different findings. Always read the record's scope, basis and restrictions.

## Reconstruction outcomes

The numerical reconstruction ledger uses the following labels. The native
EOSMAT assessment uses `not_reproduced` for a demonstrated unresolved failure;
the ledger's `parity_not_achieved` is its failed numerical-comparison outcome.
Keep the original codes and their scopes when exchanging existing records.

| Label | What it means | What it does not establish |
|---|---|---|
| **Parity** (`parity`) | The stated numerical comparison criterion was met for the varied coefficients or explicitly identified source calculation. | Experimental accuracy, phase stability, recommended use or recovery of the author's full workflow. |
| **Similar** (`similar`) | All compared free coefficients meet the ledger's broader agreement screen, but the complete parity criterion is not met. | Statistical equivalence, exact coefficient recovery or a reproduced author fit. |
| **Parity not achieved** (`parity_not_achieved`) | A numerical reconstruction was compared and failed its parity/similarity criterion, or a dedicated audit found an unresolved published-coefficient mismatch. | A diagnosis of the cause, proof of an error in the paper, or permission to replace its coefficients. |
| **Not reproduced** (`not_reproduced`) | The explicitly named reconstruction has a documented unresolved reproduction problem. Read the report to identify the affected coefficients, equation or workflow. | A blanket claim about every phase, isotherm, source observation or other EOS for that material. |
| **Not directly refittable** (`not_refittable`) | Available observations or source-reduction inputs are insufficient for a full independent refit. A qualified subset or equation check may exist. | A failed fit or proof that the published parameterization is wrong. |
| **Source reconstruction** (`source_reconstruction`) | A source calculation or composition of source equations was reconstructed without independently fitting the complete EOS coefficient set. | Independent author-fit reproduction. Calculated curves are not measured fit observations. |
| **Not assessed** (`not_assessed`) | No result is supplied for this check. New user files must not acquire a positive status merely because they load. | Either success or failure. |

These are separate outcomes, not a numerical star score. In particular,
`not_refittable` and `not_assessed` must not be displayed as failed parity.

## What "parity" compares

For the ledger's generic fit comparison, **every free coefficient** must satisfy
both the combined-error comparison and the numerical similarity screen:

\[
|p_{\mathrm{refit}}-p_{\mathrm{source}}|
\leq 2\sqrt{s_{\mathrm{source}}^2+s_{\mathrm{refit}}^2}.
\]

Here the error values are the comparison quantities recorded by the audit.
The source confidence convention, the refit's error model, covariance, weighting,
fixed constraints and row selection still matter. The ledger rule is not a
universal 95% confidence claim. Fixed/adopted coefficients do not become
independently recovered coefficients merely because a comparison keeps them fixed.

A reproduced parameterization can fit its observations poorly. Record coefficient
parity and poor data agreement separately: reproducing the paper's poor fit is
still reproduction. Usage recommendations may remain negative. Fedotenko (2020)
Pd with fixed K'=4 is such a case; see its
[audit](literature-reproductions/palladium-fedotenko-2020-audit.md).

The generic similarity limits are:

| Coefficient | Numerical similarity limit |
|---|---|
| `V0` | 5% relative difference |
| `K0` | 15% relative difference |
| `K0_prime` | 1.0 absolute **or** 20% relative difference |
| `gamma0`, `q`, `theta0` | 25% relative difference |
| `alpha...`, `dK_dT` | 30% relative difference |
| Other coefficients | 20% relative difference |

A generic **similar** result means each free coefficient passes either the
combined-error comparison or its similarity limit; **parity** requires both
for every free coefficient. A dedicated source audit can impose stricter limits,
cap an outcome, or use an explicitly documented alternative criterion. Its
recorded basis and qualification take precedence over this generic description.
See the [record-by-record ledger](primary-eos-refits.md) for actual comparisons.

The [25 October 2026 uncertainty follow-up](literature-reproductions/uncertainty-reconstruction.md)
uses an explicit printed-source-error-width criterion for 61 conditional
parameter comparisons, and the combined-error/numerical-similarity criterion
for Chen (2018) and Litasov (2007) HTBM/T03. Chen's source two-sigma widths are
converted to standard errors for that comparison; Litasov's unspecified-confidence
widths remain diagnostic comparison quantities. Shared-data cross-fit covariance
is unknown, so these are conditional numerical comparisons rather than formal
confidence tests. It retains all calibration/workflow limits.
Farla (2025) RbCl achieves conditional parameter parity while retaining separate
unresolved calibration replay, source benchmark and author-workflow findings.
The Datchi diamond and Fei preliminary Pt comparisons are explicitly conditional
on a source-unconfirmed selection/weighting hypothesis and a diagnostic fixed-V0
mask, respectively. Litasov siderite's parameter agreement retains its unresolved
source-RMS discrepancy. These conditions and findings do not cap the limited
`parameter_agreement` outcome; no hypothesis, author protocol or recommendation
is validated by promoting that numerical outcome. All 63 reviewed records meet
their stated conditional parameter criterion.

The [remaining-similar scope review](literature-reproductions/similar-parameter-review.md)
checks the other 96 similar findings. It promotes 21 conditional parameter
comparisons, two registered independent-refit replays and one source-equation
checkpoint check, while 72 retain similarity. Source-width comparisons can
establish central-value agreement without a reconstruction covariance; that
missing covariance stays unknown and no combined-error claim is inferred.
Verified half-last-digit tolerances establish printed-precision agreement,
not statistical uncertainty parity. Priors and fixed anchors remain explicit
conditions rather than independently recovered coefficients. Data-agreement,
calibration, original-publication and recommendation findings remain separate.

## Qualified results and comparison scope

A parity label must identify what was reproduced. Native EOSMAT's optional
`evidence.reproduction.assessment` requires an explicit `scope`, and an assessed
result also requires `basis`:

| Scope | What was checked |
|---|---|
| `parameter_agreement` | Coefficients agree under the recorded constraints and comparison criterion. |
| `author_workflow` | The assessment concerns replay of the author's actual inputs, selections, calibration and fitting procedure. A parameter match alone is insufficient. |
| `independent_refit` | A separately defined Peritheos analysis is reproducible under its own declared inputs and objective. This is not source-author parity. |
| `source_equation` | The source equation or calculated checkpoints were reproduced; this is not an independent data refit. |

Terms such as **conditional parity**, **partial reproduction**,
**within reported parameter errors**, **central-value parity** and
**independent refit reproduced** qualify a result; they do not erase its limits.
For example:

- Agreement inside printed parameter error widths can be a documented numerical
  criterion without establishing a formal combined-two-sigma result.
- Matching one fitted central value at printed precision does not reproduce its
  uncertainty, adopted coefficients, calibration network or original weighting.
- A cold-branch reconstruction does not establish full thermal-fit reproduction.
- Reproducing a stored independent refit does not reproduce the original paper.
- Dewaele (2008) Table VI source-equation parity checks identified Figure 4
  curves within declared graphical reading and coefficient-printing bounds.
  It does not claim recovery of the raw PAW fit or parameter-refit parity;
  the four unresolved Fe LDA sets remain separately qualified.

The [Dewaele PAW research collection](research/dewaele-2008-paw.md) preserves
its 54 assessments and full checks in its own `reproduction.json`. Those files
are outside the main catalog and its reconstruction totals; explicit file loading
does not register them as database records.

Historical records may express these findings through `parity_basis`,
`reproduction_status`, `summary_status`, scientific-validation notes or dedicated
reports. Those fields are preserved. The standard assessment is optional and
explicit; the reader does not manufacture one from a plot or source-check flag.

Bundled files also provide the standardized assessment, copied from existing
ledger findings without changing their classifications. Read `status`, `scope`,
`basis` and `qualifications` together. `source_status` preserves the ledger's
original code, including `parity_not_achieved`. The `parameters` comparison table
retains actual reconstructed coefficients and reported errors even when parity
failed. The main EOS coefficients remain unchanged. Fixed constraints and the
available protocol are in `fit`; they do not become independently recovered
coefficients. See [reconstructed coefficients in EOSMAT](eosmat-v4.md#actual-reconstructed-coefficients).

Use the [EOS reconstruction reports](reconstruction-reports/index.md) to open each
record's explanation and available coefficient graph. A report's existence, or a
coefficient graph, does not establish source-author or physical validation. The
stored `report_url` opens Read the Docs latest; `report_path` also supports local
and version-specific builds.

## Source-validation status

Source review is a separate axis from reconstruction:

| Status | Meaning |
|---|---|
| `primary_source_validated` | The represented equation, coefficients, conventions and qualifications were traced to the primary source. For a refit, its source inputs and declared analysis are traceable. This alone does not establish fit parity or physical accuracy. |
| `not_reproduced` | A documented unresolved reconstruction problem is retained with its affected scope and use restrictions. Executable diagnostic access is not a recommendation. |
| `deferred` | Source or model mapping remains unresolved; the record is retained as evidence and is excluded from ordinary executable loading. |
| `pending_primary_source_check` | An unreviewed user or staged record. Structural validity and numerical construction do not change this status. |

The [validation guide](validation.md), [contribution requirements](adding-materials-and-eos.md)
and the individual record/report determine eligibility and restrictions.
A recommended or default EOS is a separate selection decision. A default does
not mean globally best, valid at arbitrary conditions, or proven experimentally.
An explicit **do not use** restriction or recommendation hold remains binding;
a positive source or numerical label must not hide it.

## Record origin

| `record_kind` | Meaning |
|---|---|
| `published` | A literature parameterization, with the published values and source qualifications retained. |
| `refit` | Independently fitted coefficients with their own observations, objective, constraints and provenance. |
| `derived` | A declared construction from other equations or source relationships; it needs its own composition provenance and checks. |
| `diagnostic` | An inspection, teaching or exploratory record. Its presence is not a recommendation for scientific predictions or pressure calibration. |

A published record is not automatically validated. A refit is not automatically
more accurate. A diagnostic result must not replace the published record.
