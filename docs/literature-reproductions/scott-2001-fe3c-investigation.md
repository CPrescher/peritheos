# Scott (2001) Fe3C: why the reconstruction misses the published coefficients

Investigated 30 October 2026. This is a diagnostic source investigation, with no
replacement EOS or change to the published catalog record.

**The equal-pressure-weight fallback does not reproduce the source's stated
weighted fit. Weighting and the graphical input precision demonstrably move the
coefficients by the observed amount. A defensible conditional comparison can
reach both printed source coefficient widths; the original author workflow and
parameter uncertainties remain unrecovered. No tested equation, unit conversion,
parameter-order or optimizer error explains the original mismatch.**

There is also direct positive evidence for the published parameterization:
isolated portions of the drawn third-order curve agree with the published BM3
within about 1.1 native pixels RMS. This finding is more informative than simply
stating that original weights are missing. It does not establish experimental
accuracy, extrapolation validity, or the original fitting protocol.

## Source evidence and constraints

Authority is the final four-page article by Scott, Williams and Knittle,
[GRL 28, 1875-1878, DOI 10.1029/2000GL012606](https://doi.org/10.1029/2000GL012606).
The institutional PDF is retained under
`output/cementite-eos-audit-2026-10-25/sources/scott-2001.pdf`, SHA256
`abca7f6489523e6c6148cc721e3a79fa942fd444bae873bd9abd1c52871ae7e9`.
Source locations and recovery/access history in
`peritheos/data/research/cementite/scott-2001-source-audit.json` distinguish
inspected evidence from unavailable online material.

| Location | Established information | Unresolved detail |
| --- | --- | --- |
| p. 1875, Experimental | Independently measured ambient V0=155.26 +/- 0.14 Å³; conventional orthorhombic cell contains 12 Fe and 4 C, Z=4 | The article does **not explicitly say** V0 was fixed, free, constrained by its error, or part of the fit input |
| p. 1876, EOS paragraph | Weighted least squares; BM3 K0=175.4 +/- 3.5 GPa and K0′=5.1 +/- 0.3; 300 K isotherm | Numerical rows, weights, pressure/volume residual direction, error confidence/scaling and covariance |
| p. 1876, additional fits | BM2 K0=186.5 +/- 2.9 GPa; constrained K0′ range 2.5-7; approximate relation K0=226.7-9.97 K0′ GPa | Exact constrained-fit grid and objective |
| p. 1877, Fig. 3/caption | Data from two runs; third-order BM and Vinet curves indistinguishable; open circle quenched from >1500 K at 16 GPa | Per-point run membership; explicit inclusion/exclusion of the open point from fitting |
| Abstract versus body | Abstract gives K0=175 +/- 4 and K0′=5.2 +/- 0.3; body gives the values above | Reason for the different printed derivative is unspecified; do not silently substitute 5.2 for 5.1 |
| Fig. 1 versus title | Highest diffraction pattern labelled 73.2 GPa; title/abstract round to 73 GPa | Digitized maximum need not equal the original numerical maximum |

The pressure scale explicitly cites **Mao et al. (1978)** ruby fluorescence.
The compression medium was 16:3:1 methanol-ethanol-water. Above its reported
15.4 GPa hydrostatic limit, the DAC was externally heated to about 500 K between
measurements. Separate laser-heated experiments used no medium. No original
ruby readings or calibration recalculation were recovered or invented.

The complete final article contains no numerical P-V table. Public searches by
DOI, title, correction/supplement keywords and author affiliation did not recover
an original dataset, supplement, correction or different primary version.
Crossref metadata has an empty relation object and no `update-to` field.
Live Wiley full-page/abstract fetches were inaccessible (abstract fetch 403),
although the official abstract was indexed in search. These access and metadata
limits do not establish that no additional evidence exists. No authors were
contacted. Citing papers were leads, not replacement Scott observations.

## Independent equation and fitting audit

The direct replay uses SciPy and the explicit conventional-cell BM3 equation:

\[
z=(V_0/V)^{2/3},\qquad
P=\frac32 K_0(z^{7/2}-z^{5/2})
\left[1+\frac34(K'_0-4)(z-1)\right].
\]

Pressure is in GPa; V and V0 share the same volume units and normalization.
Using conventional-cell volumes or dividing both volumes by four is identical.
The finite-strain form provides a second algebraic check. Direct evaluation
agrees with Peritheos to 1.8e-13 GPa on 101 compression states. Inversion is
performed independently by scalar root finding.

For fixed V0, BM3 is linear in K0 and K0(K0′-4). An ordinary linear least-squares
solution therefore checks the optimizer globally for this pressure-residual
objective. It gives **161.341424227 GPa / 6.044124284**, agreeing with both the
direct SciPy nonlinear replay and Peritheos. Optimizer starting-point tuning is
unnecessary. The discrepancy is not a local-minimum failure.

The original 21-row CSV remains byte-for-byte unchanged, SHA256
`bf198807ec30f0b978d73e6d31791e07677ba8c8d393ff706e5875d6039e5b19`.
Thirteen rows have recoverable plotted P and V bars. The generic reader returns
`None` for each sigma array containing any missing/nonpositive entry, so both
arrays are discarded. It retains all 21 observations and performs equal pressure
residual fitting. That is a well-defined fallback diagnostic, **not a replay of
the stated weighted fit**. Missing bars are not zero errors. Digitization bounds
are not experimental sigmas and are not substituted for them.

The unweighted parameter covariance is residual-scaled inverse information:
SSE/(21-2) times the inverse normal matrix. Independently calculated errors are
8.597975 GPa and 0.635573; correlation is about -0.985. Calling the current API
with `absolute_sigma=True` or `False` gives the same covariance when no sigmas
are supplied, as required by its `not absolute_sigma or not uncertainty_supplied`
condition. The historical flag alone therefore does not denote absolute measured
errors. Changing an error scaling cannot change these central coefficients or
establish parity. Source confidence and shared-data covariance are unknown.

## What changes the coefficients

Every row below is a diagnostic with its assumptions stated. Except for the
explicit free-V0 rows, V0 is fixed at 155.26 Å³; this is a conditional treatment,
not a recovered source instruction. Errors and full covariance for every case
are in the [numerical report](../data/scott-2001-parity-investigation.json).

| Input and objective | K0 (GPa) | K0′ | Interpretation |
| --- | ---: | ---: | --- |
| Original 21 markers, equal pressure residuals | 161.3414 | 6.04412 | Existing fallback; both outside printed source widths |
| Same 21, equal volume residuals | 168.0311 | 5.55589 | Residual direction matters; does not recover both widths |
| Same 21, free V0, pressure residuals | 145.0716 | 6.76160 | V0=156.1215 Å³; freeing V0 does not fix the discrepancy |
| Same 21, free V0, volume residuals | 161.0944 | 5.86646 | V0=155.5520 Å³; still outside both widths |
| Original 13 with visible bars, equal pressure residuals | 161.3987 | 6.09304 | Subset selection alone barely changes K0 |
| **Same 13**, pressure bars as relative numerical weights | **171.9824** | **5.27476** | Both within source widths; demonstrates the weighting effect, without recovering a justified author row mask |
| Same 13, P and V bars, latent-volume errors-in-variables | 166.7568 | 5.63693 | Treating both graphical widths as independent coordinate errors gives a different answer |
| Refined pixels, 21 markers, equal pressure residuals | 160.9509 | 6.22592 | Better axis calibration alone does not restore parity |
| Original 20 nonambient markers, relative-pressure residuals | 177.8557 | 4.78408 | K0′ misses its lower printed width by about 0.016 |
| **Refined 20 nonambient markers, relative-pressure residuals** | **173.6247** | **5.21387** | Conditional agreement within both printed widths; original weights unconfirmed |
| Original 20 nonambient markers, equal normalized-stress F weights | 178.0799 | 4.83599 | Another transformation-dependent diagnostic |

Relative-pressure loss means minimizing sum[(Pmodel-Pobs)/Pobs]². It is tested
because the 13 recoverable pressure-bar half-widths are 9.34-14.43% of pressure,
with median 10.70%. This supplies a source-based reason to test that loss; it does
**not** establish missing experimental errors or the original weighting rule.
The ambient marker is omitted from this objective because division by zero is
undefined. Its reference volume is held fixed. No positive-pressure marker is
selected or discarded to obtain the agreement.

The refined relative-pressure case differs from the body values by -1.7753 GPa
and +0.11387, inside 3.5 GPa and 0.3 respectively. Its residual-scaled errors are
4.180 GPa and 0.3563; these do not reproduce the printed errors. Its curve differs
from the published BM3 by at most **0.0543 Å³** over 0-73.2 GPa. Agreement is at
the printed-error-width criterion, **not identical central-value rounding** and
not a formal confidence test. No source covariance or source workflow is inferred.

Further constrained checks bound the explanation:

- Varying fixed V0 only over 155.12-155.40 Å³ gives 164.0776/5.93499 through
  158.6295/6.15536. The printed ambient width alone does not restore the pair.
- Removing the original exact (0 GPa, V/V0=1) marker leaves fixed-V0 central
  values unchanged: its residual and parameter derivatives vanish. It affects
  residual-scaled error degrees of freedom. An omitted ambient datum is not
  the mechanism here.
- Including an approximate, separately identified quenched marker gives
  162.9062/5.93370 on the original 21-row set plus that point. Its uncertain
  membership does not account for the discrepancy under equal-pressure fitting.
- Source-supported magnetic/pressure subdivisions are not specified for this
  fit; arbitrary range masks were not used to obtain agreement.

The additional published fits offer cross-checks beyond the free BM3 pair:

| Protocol | BM2 K0 (GPa) | Approximate constrained-fit K0 = intercept + slope K0′ (GPa) |
| --- | ---: | --- |
| Published | 186.5 +/- 2.9 | 226.7 - 9.97 K0′ |
| Original equal pressure | 193.8739 | 263.7819 - 16.9638 K0′ |
| Visible-bar subset, pressure widths | 190.4099 | 249.4219 - 14.4061 K0′ |
| Original relative pressure | 186.8140 | 231.7554 - 11.0873 K0′ |
| Refined relative pressure | 186.9162 | 230.7691 - 10.8265 K0′ |

Our constrained grid uses K0′=2.5:0.1:7; the source gives the range without its
sampling. Relative-pressure weighting brings BM2 inside its printed width and
reduces the discrepancy in the tradeoff relation. It does not exactly recover
the reported relation, whose linear approximation and original grid are unknown.
This additional evidence supports protocol sensitivity rather than a claim that
one guessed rule has recovered the authors' procedure.

## Digitization and what the figure can resolve

All four PDF pages contain a single one-bit scanned image and no vector drawing
paths. Figure 3 therefore has no hidden vector markers or numerical coordinates.
Its native page image is 2381 by 3387 pixels. The original digitization's
two-point calibration is preserved; a **separate** extraction refines its marker
seeds and calibrates against seven volume-ratio ticks and nine pressure ticks.
Pixel coordinates, calibration, drawing samples and qualifications are retained
in `peritheos/data/research/cementite/scott-2001-figure3-audit.json`.
This is not an exhaustive independent point census: the seed
refinement cannot discover hidden or overlapping additional measurements or
prove the run assignments. The visual inspection did not establish an omitted
separate filled marker.

The refined ambient marker is allowed to retain its small calibrated offset;
it is not forced to (0,1): the result is -0.0615 GPa and V/V0=0.999491.
The maximum calibrated pressure is 73.2846 GPa, not an additional measurement.
Maximum coordinate changes from the old extraction are 0.531 GPa and 0.000885
in V/V0. Pixel skew, shared calibration, marker/stem overlap
and subpixel glyph-centre choices remain systematic limitations. Relative loss
puts substantial leverage on the few low-pressure markers, which explains why
small graphical changes can affect that comparison appreciably.

For an independent source-equation check, three visually selected windows
isolate portions of the drawn continuous third-order line, avoiding dashed BM2
and error-bar stems. Their 23 trace samples are correlated drawing pixels,
**not 23 measured states**. The published BM3 gives x-coordinate residual RMS
1.074 pixels and maximum 2.942 pixels, consistent with a scanned line several
pixels wide. The original equal-pressure refit gives 5.186 pixels RMS and
7.737 pixels maximum on the same samples. Fitting the traced drawing gives
K0=173.2663 GPa / K0′=5.21203, within both printed source widths. The caption's
overlapping Vinet line prevents using these pixels to discriminate BM3 from
Vinet. Trace-fit covariance describes drawing residuals, not experimental
parameter uncertainty.
The traced portions span about 19.93-77.56 GPa, including the drawn curve's
short extension above the 73.2 GPa experimental maximum. Those extension pixels
check the authors' calculated drawing and do not extend the experimental domain.

For the observational digitization, independent uniform +/-1-pixel perturbations
in both coordinates give a 2000-draw central sensitivity envelope of
K0=160.11-162.66 GPa and K0′=5.945-6.141 for the original equal-pressure fit.
Using the archival +/-0.8 GPa and +/-0.0035 ratio bounds expands that envelope
to 151.60-171.58 GPa and 5.302-6.836. These are explicitly artificial sensitivity
experiments, not confidence intervals or probabilities. Shared calibration
errors are absent from these independent perturbations.

A complementary adversarial local score calculation finds coordinated movements
bounded by about 2.384 native pixels can make the source pair stationary for
equal-pressure fitting. A full nonlinear replay after those linearized
movements gives 175.695/5.079. These constructed movements are **not recovered
observations** and are not used as a dataset or fit benchmark. They show why
graphical identifiability cannot be assessed from random jitter alone. The
published-versus-original-refit maximum curve separation, 0.4419 Å³ (about
0.309%), is around 7.4 native x pixels and lies below the archival ratio
digitization bound of approximately 9 pixels.

## Assessment and remaining inputs

| Claim | Finding |
| --- | --- |
| Source transcription/formulation | Body values, normalization and BM3 mapping checked; abstract derivative retained separately |
| Numerical implementation | Independent linear solution, SciPy and Peritheos agree; no demonstrated implementation error |
| Source-equation validation | Published BM3 consistent with separately traced source curve at scanned-line precision; does not test experimental accuracy |
| Original 21-row equal-pressure coefficient parity | Fails; its existing qualified-similarity assessment remains unchanged |
| Conditional coefficient agreement | Achieved within printed widths using refined markers plus relative-pressure loss; also using the recovered pressure bars on a 13-row diagnostic subset |
| Published parameter-error reproduction | Not established; residual scaling, error confidence and covariance differ or remain unknown |
| Author-workflow reproduction | Not achieved |
| Physical/data validation and recommendation | Not established by this investigation; no recommendation or catalog classification change |

The most specific supported explanation is that the existing reconstruction
estimates a different weighted/transformed problem from the paper's fit, on a
finite-precision graphical dataset with strong K0/K0′ tradeoff. Weighting has a
demonstrated large effect; calibration and graphical precision affect the
conditional agreement. The exact contribution of the authors' actual weights
cannot be identified from the evidence recovered. The paper need not be wrong
for this particular unweighted reconstruction to miss its coefficients.

Exact replay still needs the original numerical P-V rows and run membership,
ambient/quenched selection, V0 constraint, residual or transformed objective,
weights and their definitions, and covariance/error scaling. Ruby readings
would additionally be required to audit pressure calibration itself. Inflating
refit errors or silently assigning errors to the missing bars cannot supply
those inputs.

## Reproducibility and scope

Repository-only maintenance command (Poppler `pdfimages` required):

```bash
uv run --with pillow python -m scripts.audit_scott_2001
uv run pytest -q -W error tests/test_scott_2001_audit.py tests/test_cementite_sources.py
uv run ruff check scripts/audit_scott_2001.py tests/test_scott_2001_audit.py
```

The script imports Peritheos only after its independent SciPy calculations. It
writes the two Scott-specific numerical/extraction reports and diagnostic plots
under `output/scott-2001-parity-2026-10-30/`; it does not write material records,
original digitization, shared fit code or catalog-wide ledgers. The source PDF
remains a local, unmodified third-party document. Protocol comparison plots put
pressure on x and volume on y. Focused tests cover equation/volume convention,
linear-versus-nonlinear fitting, covariance scaling, missing sigma behavior,
source-value and input-array preservation, conditional comparisons and the
ambient-point constraint. Passing them is numerical verification, not
experimental validation.

Verification: 14 focused tests passed with warnings treated as errors; Ruff
check and formatting passed for the new script and test. The complete plotting
command completed, and the protocol plot and native-source overlay were visually
inspected. A strict-warning run of the numerical `run()` function also passed.
The retained replay uses Poppler extraction and SVD of the Jacobian for
covariance. The broader repository suite was not run for these isolated audit
artifacts. No commit, push, release or deployment was made.
