# Frost (2023) palladium: fit reproduction and recommendation hold

Audited 17 October 2026. **Neither published Pd fit passes the strict
coefficient-reproduction check when the reported coordinate errors are treated
as absolute standard deviations.** This is a conditional diagnostic conclusion:
the paper does not specify that covariance convention, and residual-scaled
refit errors make both parameter comparisons pass the same two-width screen.
An unconditional failure of the authors' fitting protocol has not been established.
Both records are `not_reproduced`, explicitly not recommended for quantitative
use or pressure calibration, and available for diagnostic selection. Pd has
no automatic equilibrium default pending resolution or a separately audited
replacement. The source coefficients, errors, datasets and separate Baty 0 K
DFT record are preserved. No diagnostic refit is installed as an EOS.

This is an unresolved reproduction finding, not proof that the paper's EOS is
physically wrong. The numerical coefficient classification remains `similar`;
that screening category is not scientific acceptance or author-fit reproduction.

## Primary evidence recovered

M. Frost, D. Smith, E. E. McBride, J. S. Smith and S. H. Glenzer,
“The equations of state of statically compressed palladium and rhodium,”
J. Appl. Phys. **134**, 035901 (2023),
[DOI](https://doi.org/10.1063/5.0161038).

The final typeset article was recovered from the local library during the
user-requested protocol recheck. All five printed article pages were read,
including Methods on p. 035901-2, Table I, Figure 1 and caption, Results,
data availability and the references. Its fitting description and coefficient
table agree with the previously read public
[accepted manuscript](https://pubs.aip.org/aip/jap/article-pdf/doi/10.1063/5.0161038/18050563/035901_1_5.0161038.am.pdf).
The final PDF is stored at
`/Users/clemens/Zotero/storage/85R42HTW/Frost et al. - 2023 - The equations of state of statically compressed palladium and rhodium.pdf`.
It contains the publisher cover and five printed article pages.

The [official supplement](https://pubs.aip.org/jap/article-supplement/2903099/zip/035901_1_5.0161038.suppl_material/)
was independently downloaded again. Its archive and PDF checksums exactly
match the previously recorded source checksums:

| File | SHA-256 |
|---|---|
| Final typeset article, local library copy | `c80f3e91b4c6deb121819f3c72bcc70ebb49f9240e55899be03b0746955f8f71` |
| Accepted manuscript | `8de7587ecdf2b9ce635bd76d958047f834d93794924efffee6893006cfc24c1f` |
| Official supplement ZIP | `7f03394af010c3eaa6b4346ac6ed0e2273ff86f28ba1cf910471d564873ae015` |
| Supplementary PDF | `7e47e13580ecc7faec0a372347217e860aee65c6e32cc279a2c6e66651653f3b` |
| Bundled CSV, unchanged | `92e4ece39614b1a9f980bff32ae3986083b06bd9673ef0e615f4d27bd39fa75c` |

Every printed Pd pressure, lattice parameter, lattice error, loading identifier
and marker in Tables I-III was compared again with the CSV. **All 93 agree.**
The independent transcription retains original decimal tokens in
[source rows](../data/frost-2023-palladium-source-rows.json). Rh Table IV and
hydrogen-contaminated runs are excluded by material identity, not by fit quality.

| Table / loading | Rows | Pressure (GPa) | Marker |
|---|---:|---:|---|
| I / BX18 | 26 | 7.20–60.32 | W, Dewaele (2004) |
| II / BX32 | 20 | 2.99–62.74 | W, Dewaele (2004) |
| III / HP | 47 | 29.10–181.79 | Au, Dorfman (2012) |

All are hydrogen-excluded fcc Pd in Ne at room temperature. The record's
300 K value represents that room-temperature isotherm. A conventional fcc
cell contains four Pd atoms. Source normalization is `V=a^3`, with
`sigma(V)=3a^2 sigma(a)`; plots use `V/4` in Å³/atom and pressure horizontally.
The actual compression table starts at 2.99 GPa; the paper's 182 GPa endpoint
is rounded. Pressure-marker identities are recovered, but paired W/Au volumes
are not published here. Thus the pressure calibration cannot be independently
recomputed from this supplement.

## What is and is not specified about the source fit

Methods, final p. 035901-2, states nonlinear Marquardt-Levenberg least squares
accounting for pressure errors estimated as 0.1 GPa in all runs and volume
errors derived from Le Bail refinement. Table I supplies all three coefficients
and their errors. The source does not identify its software, exact objective,
whether residual-coordinate linearization was used, uncertainty rescaling,
parameter covariance or confidence convention. No fixed V0 constraint or row
exclusion is reported. All three coefficients therefore vary in this replay;
this tests the published parameterization rather than imposing a new ambient
constraint. The original author workflow remains unrecovered.

Supplement II.A separately reports `a=3.89088(2) Å` at ambient pressure “in
addition to Tables I, II and III.” It is displayed separately and is not silently
added to the 93 compression rows. Whether it participated in the original
regression is not stated. The fitted V0 is an extrapolation, not that independent
ambient measurement. Previous ambient-inclusion diagnostics did not recover
the published triplets; they do not establish author row selection.

### Full-protocol recheck: reported steps versus our assumptions

| Stage | What the final article specifies | Boundary of our replay |
|---|---|---|
| Sample and loading | Hydrogen-excluded Pd in Ne; three DAC loadings, up to 182 GPa | All 93 hydrogen-free supplementary Pd rows retained; contaminated-material comparisons excluded |
| Diffraction reduction | Dioptas integration followed by JANA06 Le Bail fitting | Published lattice parameters/errors used; raw images and Le Bail projects not replayed |
| Pressure reduction | W EOS, Dewaele (2004), for Runs 1/2; Au EOS, Dorfman (2012), for Run 3 | Published pressures retained; paired calibrant observations unavailable |
| EOS model | Vinet and BM3, with Table I coefficients and errors | Both equations implemented literally; no alternate model introduced |
| Regression | Nonlinear least-squares Marquardt-Levenberg accounting for both P and V uncertainty | Exact orthogonal and effective-variance objectives are our explicit interpretations, not recovered author code |
| Input errors | Estimated 0.1 GPa in all runs; volume uncertainty from Le Bail fitting | `sigma(V)=3a^2 sigma(a)` is our propagation of the supplementary lattice errors |
| Parameter constraints and ambient point | No fixed V0, prior, exclusion mask or ambient-fit inclusion is explicitly stated | All three coefficients free and ambient point omitted is a diagnostic selection, not an established author manifest |
| Parameter uncertainty | Table I error widths, without confidence/covariance prescription | `absolute_sigma=True` is our assumption; residual-variance scaling is not ruled out by the paper |

The final article contains no additional regression formula, software name for
the EOS fit, or fitting-protocol reference beyond these statements. JANA06 is
identified for the Le Bail stage; that does not establish it as the EOS fitter.
The protocol is sufficient for a documented uncertainty-aware independent
refit, but not an exact author-workflow replay. The earlier conclusion must
therefore remain conditional on our explicit choices.

## Reproducible calculations

Reproduction script: `scripts/reproduce_frost_2023_palladium.py`.
The [full numerical report](../data/frost-2023-palladium-audit.json) is bundled.
Run `python scripts/reproduce_frost_2023_palladium.py`; optional
`--source-text` checks freshly extracted supplementary text and `--plot-dir`
creates comparison plots. Factual source rows are bundled, so a new download
is not needed for numerical replay.

Vinet is evaluated literally as
`P=3K0(1-x)x^-2 exp[1.5(K0'-1)(1-x)]`, `x=(V/V0)^(1/3)`.
BM3 is evaluated as
`P=(3K0/2)(eta^7-eta^5)[1+(3/4)(K0'-4)(eta^2-1)]`,
`eta=(V0/V)^(1/3)`. Independent equations agree with the native evaluator
within 1e-10 GPa. Units are GPa and conventional-cell Å³ throughout fitting.
The existing cross-publication checkpoint, Baty supplementary Table S3,
gives 10.04 GPa for Frost's Vinet at 14 Å³/atom and is still reproduced within
its 0.005 GPa printed-pressure rounding width. This verifies evaluation of
the published curve; those calculated Table S3 values are not Frost fit
observations and do not resolve the regression discrepancy.

The explicit unweighted baseline minimizes `sum[P(Vobs)-Pobs]^2`. The
uncertainty-aware replay minimizes the exact orthogonal objective

`Q=sum{[(P(v_i;b)-Pobs_i)/0.1]^2 + [(v_i-Vobs_i)/sigma(V_i)]^2}`,

with three global coefficients `b` and 93 latent volumes `v_i`. There are
`186-96=90` residual degrees of freedom. No exclusions, robust loss, adjustable
error multipliers, alternate EOS families or additional coefficients are used.
The same objective is independently solved by SciPy ODR using literal equations;
coefficients agree to better than 4e-6 GPa for K0 and 1e-7 for K0'. A second
check now explicitly uses `least_squares(method="lm")` on the 186 joint
coordinate residuals and 96 unknowns. It reaches the same orthogonal minimum
and coefficients (K0 difference below 1.3e-5 GPa). Thus the choice of optimizer
does not resolve the disagreement under this assumed objective. A separate
source-consistent interpretation minimizes the first-order effective-variance
objective `sum[(P(Vobs)-Pobs)^2/(sigmaP^2+(dP/dV sigmaV)^2)]`. It also fails to
recover the published central values. These are documented interpretations of
the prose, not a claim of recovered author software.

| Model / objective | V0 (Å³/cell) | K0 (GPa) | K0' |
|---|---:|---:|---:|
| Published Vinet | 58.678(73) | 189.3(3.0) | 5.473(63) |
| Vinet orthogonal replay | 58.766518 | 184.190605 | 5.603712 |
| Vinet effective variance | 58.766861 | 184.172488 | 5.604565 |
| Vinet unweighted baseline | 58.592415 | 193.436574 | 5.381540 |
| Published BM3 | 58.588(76) | 197.5(3.3) | 4.996(66) |
| BM3 orthogonal replay | 58.684581 | 191.476721 | 5.147382 |
| BM3 effective variance | 58.685012 | 191.450808 | 5.148579 |
| BM3 unweighted baseline | 58.477437 | 202.947706 | 4.880678 |

Formal absolute-coordinate-error standard errors for the Vinet replay are
`(0.014476, 0.612883, 0.013813)` and for BM3
`(0.014085, 0.621865, 0.013925)`. Complete covariance/correlation matrices are
in the report. They are conditional on independent, correctly specified errors,
and omit pressure-scale, run correlations and model inadequacy. These formal
errors are not a scientific uncertainty budget for a replacement EOS.

## Correcting the earlier RMS comparison

The older audit reported a pressure RMS of 0.276 GPa for the Vinet refit and
0.293 GPa for BM3. Those numbers use `P(v_i)-Pobs_i` at **adjusted latent volumes**;
the published RMS used measured volumes. That was an unfair comparison.
It is corrected here and in the maintained ledger, retaining the adjusted-volume
component under an explicitly named diagnostic field.

| Quantity, all 93 rows | Vinet | BM3 |
|---|---:|---:|
| Published RMS `P(Vobs)-Pobs` (GPa) | 0.771992 | 0.845419 |
| Orthogonal replay RMS `P(Vobs)-Pobs` (GPa) | 0.825323 | 0.915380 |
| Unweighted baseline RMS `P(Vobs)-Pobs` (GPa) | 0.763097 | 0.829299 |
| Published Q, latent volumes optimized with coefficients fixed | 1833.641 | 2161.045 |
| Free-coefficient replay Q | 1714.418 | 1987.165 |
| Replay reduced chi square, Q/90 | 19.049 | 22.080 |

An orthogonal fit minimizes Q, not unweighted pressure RMS. Its larger measured-
volume pressure RMS is therefore compatible with its smaller Q. The source
does not report a numerical RMS or chi-square against which to claim that a
published fit statistic is contradicted. The plot shows small overall volume
differences with structured residuals; its coordinate errors alone do not
describe the observed scatter.

## Acceptance decision and its limits

The coefficient differences remain outside the ledger's strict conditional
combined-two-error-width check for K0': 2.027 combined widths for Vinet and
2.244 for BM3. V0 and K0 pass that check. All coefficients pass the looser
numerical `similar` screening thresholds. The strict Vinet miss is marginal,
and source error confidence is unspecified; it must not be exaggerated.
The combined-width rule is a catalog screening criterion, not a calibrated
hypothesis test: the two fits share observations and their cross-covariance
is unavailable.

Scaling the replay covariance by reduced chi-square increases its errors and
would make the marginal two-error-width intervals overlap. This changes the
comparison convention, not the recovered central coefficients or objective.
The source does not document its covariance convention. We report both error
conventions without claiming either as the recovered author convention. Correlated runs
and systematic errors make a naive chi-square significance test inappropriate.

| Replay covariance convention | Vinet K0' shift / combined error width | BM3 K0' shift / combined error width | All three coefficients inside two widths? |
|---|---:|---:|---|
| Absolute coordinate errors | 2.027 | 2.244 | No, for either model |
| Covariance scaled by replay residual variance | 1.499 | 1.629 | Yes, for both models |

Covariance scaling leaves the fitted central coefficients unchanged. Thus it
does not establish recovery of the original regression, but it demonstrates
that an unconditional claim that the published parameters fail uncertainty
agreement is unsupported. The source's uncertainty convention must first be
resolved. The current `not_reproduced` designation means the original fit has
not been recovered; the recommendation/default hold remains pending that
clarification. It must not be presented as proof that the authors' protocol
was replayed and failed.

Source-coordinate rounding contributes at most about 0.047 GPa per row under
a first-order bound (`P` rounded to .01 GPa, `a` to .0001 Å), much less than the
largest 2.575 GPa Vinet residual. This does not prove all regression-weight
rounding is negligible: the lattice errors themselves have limited printed
precision. The previous categorical claim that rounding was ruled out as an
explanation of the coefficient mismatch was too strong. No unsupported
uncertainty rescaling is adopted here.

The user-requested recommendation/default hold is precautionary pending
resolution of the source protocol. The absolute-error diagnostic alone does
not establish an unconditional failure of the publication or its protocol;
missing exact author-workflow details do not by themselves invalidate an EOS.
Both source records are retained unchanged scientifically, with explicit
`not_reproduced` status denoting the unrecovered original fit.
No better-fitting diagnostic is promoted as a validated replacement. This
decision is limited to these two Pd parameterizations; it does not audit Rh
or reject the publication's entire experimental work.

Pd carries `recommendation_hold.status=on_hold` and
`scope=equilibrium_default`. Python and Rust return no equilibrium default
while this hold is active, including after round trip or explicit record
selection. Other source records remain selectable. Baty's DFT stays at 0 K
and does not become an implicit 300 K substitute. Resolution requires the
actual author fit inputs, masks/ambient treatment, objective, uncertainty
transformation and covariance convention, or a justified independent
validation of another source. No author has been contacted by this audit.

## Software verification

All 3,000 Python tests passed across the complete run and corrected rerun:
the full run passed 2,988 tests and its 12 reporting-path failures passed
on rerun with warnings treated as errors. Those failures came from an
intermediate reporting-code variable-placement error, corrected while the
full run was in progress. The dedicated 205-test selection also passed with
warnings treated as errors. All 160 Rust workspace tests passed with all
features enabled. Repository lint/format checks, Rust formatting and strict
Clippy, and mypy for all 56 package source files passed.

The complete 663-record numerical-ledger consistency check passed, as did
the paper, source-index and catalog-coverage consistency checks. A final
comparison confirms that all seven Pd source equations, parameter errors,
fixed choices, reference temperatures, parameter/fit provenance and pressure
calibration metadata, plus all dataset declarations, are preserved. The
Frost CSV checksum is unchanged. Python and Rust tests verify the hold,
explicit diagnostic selection, 0 K DFT separation and round trips.

The strict documentation build remains blocked by seven existing broken
links in the separate NaCl and ringwoodite audits. The new Pd audit links
pass. These software checks establish reproducibility and catalog behavior;
they do not turn the unresolved scientific fit into a validated EOS. Changes
remain local and uncommitted.
