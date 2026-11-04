# Sata et al. (2010): original Fe3C static EOS and pressure-basis audit

Audited and integrated **1 November 2026**. The original published static curve
is available as `cementite_nonmagnetic_sata_2010_bm3` in
`cementite_nonmagnetic`. It is a **nondefault published** record, restricted to
**300 K and 27.2–186.6 GPa**, including the source's explicitly cited Li points.
**Qualified use:** some quoted point errors are smaller than the overall
experimental scatter. This does not establish that the published coefficient
errors are wrong. The published EOS is retained; the recalculated-pressure
refit remains a separate diagnostic. Conditional parameter parity and the
numerical data-agreement assessment remain unchanged.
Existing records, defaults, structures and datasets are preserved.

The original source is Sata, Hirose, Shen, Nakajima, Ohishi and Hirao,
*Journal of Geophysical Research* **115**, B09204,
[DOI 10.1029/2009JB006975](https://doi.org/10.1029/2009JB006975).
The final publisher HTML supplies Tables 1, 4 and 6, Equations 1–2 and
paragraphs 5–6, 12, 14–15 and 24–27. It is the primary authority.

## What was recovered and verified

The available local PDF, `sata_JGR_SE.pdf`, is a **2009 manuscript**, not the
final publication: metadata title `sata_JGR_5th_last`, created 11 September
2009, 51 PDF pages. Its Fe3C Table 1, Eq. 1 and relevant fit prose were rendered
and read. All **nine Fe3C rows / 72 numerical tokens**, including every
parenthesized error and paired MgO reading, match the final publisher HTML.
This correspondence does not turn the manuscript into the version of record.
Its coefficient-summary Table 4 became final Table 6; unrelated FeS entries
differ and were not imported. The publisher PDF and official Table 1/Table 6
text-download routes returned 403; the final HTML tables were accessible.

| Data | Retained content | Role |
|---|---|---|
| Sata final Table 1 | Nine Fe3C states in two runs (2 + 7), all a/b/c/V/P error tokens, paired MgO a/V/errors and VC/VFe | Original measurements; new-state coverage 49.87–186.60 GPa |
| Li (2002), Table 1 T1-12/T1-13/T1-14 | Unchanged 27.2, 30.3 and 30.5 GPa states, volumes and pressure/volume widths | The three measured points above 25 GPa required by Sata paragraph 15 |
| Reconstructed source selection | All nine Sata states plus those three Li states, 12 rows total | Declared joint-input replay, with source identities and original pressure bases retained |
| Sata final Table 4 | All 12 hkl, measured d/error, source-calculated d and printed delta-d entries | Diffraction provenance only, never EOS fit inputs |

The three Li observations already exist in `cementite_li_2002_table1`; the
selected table does not make them new independent measurements. No additional
ambient, quenched or calculated EOS point is inserted. Sata's annealing below
1300 K precedes room-temperature XRD; those heating temperatures are not hot
P–V–T observations.
Li supplies no run/measurement labels; those fields remain empty/null in the
joint table, while its original T1 row IDs identify the observations.

The packaged source manifest,
`peritheos/data/datasets/cementite-sata-2010-source.json`,
retains version/normalization decisions, source hashes, retrieval limits,
overlap with Li and source-rights qualifications. Original tokens are in
`peritheos/data/datasets/cementite-sata-2010-table1-tokens.csv`; normalized
Table 1 and the selected input are reusable typed datasets. Source numerical
licensing is unspecified; the contributor CC0 dedication covers only rights
in transcription, normalization and arrangement, not article text or figures.

## Material, magnetic branch and volume convention

Sata identifies high-pressure cementite and assigns its elastic branch to
the nonmagnetic phase using earlier magnetic-transition evidence near 25 GPa.
The XRD data alone are not a magnetic-state measurement or phase diagram.
Minor Fe7C3 and Re peaks were also present at the highest pressure; the
indexed Fe3C peaks were refined separately. No mixed composition is invented.

Table 1 explicitly gives **Z=4**, hence Fe12C4 and **16 atoms per conventional
P cell**. Table 6 prints V0=9.341 Å³/atom:

```text
V0_cell = 9.341 * (4 atoms/formula) * (4 formula units/cell)
        = 149.456 Å³/cell.
```

The final prose's 149.46 Å³/cell is rounded. The catalog retains the atomic
source precision via 149.456, with no invented V0 error. This is an extrapolated
nonmagnetic branch parameter, not a measured stable ambient cell.

The installed card retains its separately sourced Prescher (2012) high-pressure
Pnma atomic model. Sata's axis order is `(standard c, standard a, standard b)`;
the volume and Z are invariant. Its hkl and lattice columns are retained in
their own setting, not silently assigned to the Prescher coordinate system.
Using the final Table 4 cell reproduces its calculated d-spacings within
4.8×10⁻⁶ Å, consistent with their printing precision. Observed and calculated
d remain separate. The Table 4 title says **193 GPa**, while paragraph 14 and
Figure 3 say **187 GPa**, and Table 1 measurement **193** has **186.60(33) GPa**.
Both labels are preserved; no pressure is repaired from an inferred cause.

## Published equation forms: one curve, not several independent EOSs

For the modified third-order Birch–Murnaghan equation, define
`f=((Vr/V)^(2/3)-1)/2`. Sata Eq. 1 is

```text
P = (1+2f)^(5/2) * [Pr + (3Kr-5Pr)f
                    + (4.5Kr(Kr'-4)+17.5Pr)f²].
```

Paragraph 15 reports the selected reference Vr=9.341 Å³/atom,
**Pr=0.0±1.6 GPa, Kr=290±13 GPa and Kr′=3.76±0.18**. Vr is selected;
Pr is a fitted output, not a source-fixed zero. With its **printed central
Pr=0**, this is exactly ordinary BM3, as stated in paragraph 12 and Eq. 2.
Final Table 6 explicitly publishes that zero-pressure representation.

The sole installed source record uses Table 6: V0=149.456 Å³/cell,
K0=290 GPa and K0′=3.76. The reference fit, including its Pr error, is retained
in `published_reference_fit`. No source covariance or V0 error is supplied;
Pr's error is not silently converted into a V0 error. The source record has
no inferred fixed-V0 fitting flag. The fixed-V0 comparison below is an
explicit diagnostic constraint.

**Litasov (2013) is separate.** Its Table 2 later recalibrates earlier pressure
standards and lists a Sata-attributed **Vinet** comparison at 298 K with
V0=148.5 Å³/cell, K0=316(6) GPa and K0′=3.5(1). Those values are neither
Sata's original Table 6 nor another original Sata EOS variant. They remain
in the completed [Litasov investigation](litasov-2013-fe3c.md#every-table-2-comparison-row).
No Litasov coefficient, pressure or 298 K convention enters this integration.

## Conditional fit replay and residuals

All calculations use the literal printed volumes, not lattice products
substituted for the volume column. The independent reference-form replay is
linear least squares in `(Pr, Kr, Kr*Kr′)` at fixed Vr. A second two-coefficient
linear solution holds the published zero-pressure V0 fixed and fits
`(K0, K0*K0′)`; its covariance is transformed to physical coefficients.
Neither path calls Peritheos's fitting implementation.

| Comparison | Pr (GPa) | K or Kr (GPa) | K′ or Kr′ | Unweighted pressure RMS (GPa) |
|---|---:|---:|---:|---:|
| Published source curve on all 12 selected rows | 0 (printed) | 290±13 | 3.76±0.18 | 1.106733 |
| Literal mBM3, joint 12 rows, equal pressure weights | 0.01±1.5 | 290±12 | 3.76±0.18 | 1.080402 |
| Ordinary BM3, published V0 fixed, same 12 rows | 0 by model | 290±4 | 3.76±0.08 | 1.080405 |

Full-precision diagnostic values are retained in the
[numerical report](../data/sata-2010-fe3c-audit.json): the mBM3 solution is
Pr=0.0117001 GPa, Kr=290.407438 GPa and Kr′=3.76370512; the fixed-V0 BM3
solution is K0=290.494528 GPa and K0′=3.76249517. All three reference-form
central values and both fixed-V0 free coefficients lie within their printed
source widths. The common ledger therefore records **conditional parameter
parity** for the explicitly named fixed-V0 comparison. This is not a recovered
author workflow, uncertainty-confidence match or statistical equivalence test.

The nine Sata rows alone instead give Kr≈296.97 GPa and Kr′≈3.68357;
omitting Li is not a replay of the reported joint selection. Pressure-width
and effective-P/V-width sensitivities on the same 12 rows are separately
retained, with no preferred weighting inferred from numerical resemblance.
Effective-width diagnostics assume uncorrelated P/V errors and omit shared
calibration covariance. Diagnostic covariance is RSS/dof-scaled and conditions
on fixed inputs; the source error confidence and covariance remain unknown.

Residuals mean **P(model) − printed P**, in GPa at measured volume. The unchanged
published curve has maximum absolute residual **2.591677 GPa** on the joint
selection, and 1.255813 GPa RMS on the nine new Sata points. The observations
do not all lie within their small printed coordinate widths. The source gives
no strict residual envelope or exact regression objective; coefficient
agreement does not turn the scatter into experimental validation.

The coordinate-error comparison uses
`s_eff=sqrt(sP²+(K(V)*sV/V)²)`, with the published curve's bulk modulus
and all printed pressure/volume error widths. This is a first-order diagnostic
assuming uncorrelated coordinates. At the printed central coefficients:

| Sata pressure (GPa) | P(model) − P(table) (GPa) | Printed P error width (GPa) | Combined P/V width (GPa) | Absolute residual / combined width |
|---|---:|---:|---:|---:|
| 69.21 | +2.20049 | 0.39 | 0.41069 | 5.36 |
| 77.32 | −0.82361 | 0.20 | 0.24070 | 3.42 |
| 100.78 | −2.59168 | 0.23 | 0.28723 | 9.02 |
| 186.60 | −1.08254 | 0.33 | 0.39073 | 2.77 |

**Six of nine Sata observations exceed their combined coordinate-error width**;
seven exceed their pressure-error width alone. All three included Li points
are within their larger combined widths. Per-row calculations and the complete
12-row comparison are retained under `data_agreement` in the numerical report.
These ratios are descriptive, **not formal sigma significances**: confidence
conventions, coordinate/calibration covariance and fitted-parameter covariance
are unavailable. No predictive confidence band or reduced chi-square is inferred.
The published parameter errors are retained; their unknown covariance is not
replaced by an independent-coefficient assumption to construct a wider band.

The numerical data-agreement finding is **poor**, and the usage recommendation
is **qualified use** of the original published EOS within its stated limits.
The record remains a named nondefault source option with this explicit caveat:
**some quoted point errors are smaller than the overall experimental scatter.** The primary-source status
means the coefficients, equation and conventions are source-traceable;
conditional `parameter_agreement: parity` remains unchanged. Neither status
establishes physical validation or adequate agreement with the measurements.
This follows the [separate coefficient/data-agreement assessments](../eos-labels.md#what-parity-compares)
and the [acceptance requirements](../adding-materials-and-eos.md#acceptance-requirements):
the source supplies no quantitative residual criterion whose replay has been
shown to fail, so poor observed agreement is not relabelled as failed
coefficient reproduction.

The generated [coefficient, curve and residual comparison](../reconstruction-reports/cementite_nonmagnetic.md#record-cementite_nonmagnetic_sata_2010_bm3)
uses all 12 rows on their original pressure bases. Native pressure evaluation,
inversion and reference normalization are tested independently of fitting.

## Pressure scales and unresolved reduction precision

All nine new Fe3C states are paired with **MgO**. The final methods cite
Speziale et al. (2001), although the Table 1 footnote labels it 2000. The
printed reduction is BM3 with V0=74.71 Å³/cell, K0=160.2 GPa and K0′=3.99;
the bundled `mgo_speziale_2001_bm3_2` has that 300 K isotherm. Its thermal
parameters do not contribute at the 300 K reference.

Replay from the **measured MgO volumes** differs from all nine printed pressures,
with maximum absolute difference **0.352639 GPa** (model lower). Recomputing
marker volume as a³ does not remove it. A one-coefficient reducer diagnostic
at the printed V0/K0 gives K0′≈3.99681, but this does not recover the authors'
actual coefficient precision, software or calculation. No revised MgO record
is created, and the original observations and printed 3.99 remain unchanged.
Separately labelled derived pressure coordinates are displayed below. This is a
separate pressure-reduction limit,
not a failure to recover the Fe3C source's conditional coefficients.

### Recalculated Fe3C pressure coordinates

For transparency, the nine Sata Fe3C volumes are also shown at pressures
**recalculated from their paired measured MgO volumes using the printed
Speziale BM3 coefficients**. These are derived pressure coordinates, not newly
measured points or established corrections to Sata's original reduction.
The two versions of each point share the same Fe3C and MgO observations and
are not independent data. Li's three ruby/neon points are excluded from this
comparison and retain their original pressures.

![Nine original Sata Fe3C points and the same volumes at recalculated MgO pressures, with pressure shifts and unchanged-source-curve residuals](assets/sata-2010-fe3c-pressure-recalculation.svg#only-light)
![Nine original Sata Fe3C points and the same volumes at recalculated MgO pressures, with pressure shifts and unchanged-source-curve residuals](assets/sata-2010-fe3c-pressure-recalculation-dark.svg#only-dark)

| Sata state | Measured Fe3C V (Å³/cell) | Published P (GPa) | Recalculated P (GPa) | Change in P (GPa) |
|---|---:|---:|---:|---:|
| Run 1, 29 | 130.804 | 49.87 | 49.836 | −0.034 |
| Run 1, 45 | 125.238 | 69.21 | 69.151 | −0.059 |
| Run 2, 56 | 124.076 | 77.32 | 77.247 | −0.073 |
| Run 2, 83 | 119.583 | 100.78 | 100.661 | −0.119 |
| Run 2, 111 | 114.794 | 124.81 | 124.626 | −0.184 |
| Run 2, 130 | 112.311 | 141.37 | 141.156 | −0.214 |
| Run 2, 151 | 109.874 | 159.13 | 158.869 | −0.261 |
| Run 2, 172 | 107.347 | 177.77 | 177.452 | −0.318 |
| Run 2, 193 | 106.449 | 186.60 | 186.247 | −0.353 |

The [downloadable derived-coordinate CSV](../data/sata-2010-fe3c-recalculated-points.csv)
retains source row/run/measurement identities, both measured volumes and their
printed widths, original pressure/error, recalculated pressure, shifts and
residuals. Extra calculated digits support numerical replay and do not imply
extra experimental precision. Original error bars are shown for the original
coordinates only. No independent total uncertainty or coefficient/calibration
covariance is inferred for the recalculated coordinates.

This coordinate-only comparison keeps the **published Fe3C curve unchanged** and performs
**no Fe3C refit**. Its nine-point pressure RMS changes from **1.255813 to
1.225289 GPa**, and its maximum absolute residual from **2.591677 to
2.472864 GPa**. Recalculation therefore leaves the larger observation scatter
in place. The published record retains qualified use with explicit uncertainty
and data-scatter caveats. The original-data conditional parameter-parity assessment remains
unchanged. No recalculated points or replacement fit are installed in the catalog.

Repository-only reproduction:

```sh
uv run python -m scripts.audit_sata_2010_mgo_pressure --write
uv run python -m scripts.audit_sata_2010_mgo_pressure --check
uv run python -m scripts.plot_sata_2010_pressure_recalculation
```

The displayed-coordinate follow-up passes **27** focused source, recalculation,
uncertainty and documentation checks. The light/dark SVGs were rendered and
visually inspected. Snapshot hashes preserve every catalog/material and source
observation file, the original Fe3C numerical audit and all preceding MgO
diagnostics; the new calculated rows, CSV and figures remain separate evidence.

### Separate diagnostic refit using recalculated Sata pressures

At the user's explicit request, a **separate diagnostic** fits all nine Sata
Fe3C volumes at the recalculated MgO pressures together with the **three
unchanged Li points**. The original-data fits and the published record remain
unchanged. No source observation is replaced, no new EOS is installed, and
Litasov's later recalibration is excluded.

The primary comparison uses the **same method on both pressure versions**:
ordinary BM3, V0 fixed diagnostically at the printed 149.456 Å³/cell, and equal
pressure weights at the measured Fe3C volumes. The objective is the sum of
squared pressure residuals. This fixed anchor is not a newly measured ambient
volume or an inferred source constraint. Each fit uses all 12 rows and has ten
residual degrees of freedom. The recalculated-input envelope is
**27.2–186.247361 GPa at 300 K**; this does not change the published record's
original limits.

| Curve / diagnostic | K0 (GPa) | K0′ | RMS on original pressures (GPa) | RMS on recalculated Sata + unchanged Li (GPa) |
|---|---:|---:|---:|---:|
| Published BM3, unchanged | 290 | 3.76 | 1.106733 | 1.080767 |
| Equal-weight diagnostic on original pressures | 290.494528 | 3.762495 | 1.080405 | 1.093770 |
| Equal-weight diagnostic on recalculated Sata + unchanged Li | 290.524468 | 3.752552 | 1.095609 | 1.078541 |

Extra calculated digits expose the paired coefficient shift; they do not imply
experimental accuracy. Changing the Sata pressure reduction shifts **K0 by
+0.029940 GPa** and **K0′ by −0.009944** relative to the matched original-pressure
diagnostic. The recalculated fit's conditional standard errors are about
**4.20 GPa and 0.0839**, versus 4.21 GPa and 0.0842 for the original-pressure
diagnostic. Both use RSS/dof-scaled covariance and condition on the fixed V0.
Coordinate errors, calibrant-parameter errors and calibration covariance are
not propagated. The pressure widths retained in the input CSV belong to the
original observations and are **not used as regression weights** or assigned
to the recalculated pressures. No recovered author-confidence convention or
predictive confidence band is claimed. The two fits share their observations;
their coefficient shift is not an independent-fit significance test.

The original-reference **mBM3** method provides a matched three-coefficient
check at selected Vr=149.456 Å³/cell, with Pr fitted:

| Input pressures | Pr (GPa) | Kr (GPa) | Kr′ | Pressure RMS (GPa) |
|---|---:|---:|---:|---:|
| Original Sata + original Li | 0.011700 | 290.407436 | 3.763705 | 1.080402 |
| Recalculated Sata + unchanged Li | 0.033908 | 290.272067 | 3.756051 | 1.078511 |

These are reference-volume coefficients, not a new zero-pressure parameter
list when Pr differs from zero. The mBM3 fits have nine residual degrees of
freedom; complete conditional covariance and error estimates are archived.

![Published curve and matched original-pressure and recalculated-pressure BM3 diagnostic refits, evaluated on both pressure coordinate sets](assets/sata-2010-fe3c-recalculated-pressure-refit.svg#only-light)
![Published curve and matched original-pressure and recalculated-pressure BM3 diagnostic refits, evaluated on both pressure coordinate sets](assets/sata-2010-fe3c-recalculated-pressure-refit-dark.svg#only-dark)

**The refit leaves substantial scatter.** On its own recalculated inputs, the
fixed-V0 BM3 has RMS **1.078541 GPa** and maximum absolute residual
**2.385979 GPa**, compared with **1.080405 and 2.393862 GPa** for the matched
original-pressure fit on its own inputs. The original-data conditional
parameter-parity finding and the **poor-data-agreement / uncertainty
qualifications** of the published record are preserved.
This new diagnostic does not establish the historical pressure reduction or a
validated replacement EOS.

The separate [diagnostic report](../data/sata-2010-fe3c-recalculated-pressure-refit.json)
stores the complete fit protocol, covariance, input selection, source hashes,
paired coefficient changes and all three models' evaluations on both pressure
sets. The [12-row fit-input CSV](../data/sata-2010-fe3c-recalculated-refit-inputs.csv)
retains original and diagnostic pressures and their calibration identities;
the [residual CSV](../data/sata-2010-fe3c-recalculated-refit-residuals.csv)
retains all curve/data combinations.

Repository-only reproduction:

```sh
uv run python -m scripts.refit_sata_2010_recalculated_pressure --write --plot
uv run python -m scripts.refit_sata_2010_recalculated_pressure --check
uv run pytest -q -W error tests/test_sata_2010_recalculated_pressure_refit.py
```

Verification: **30** focused refit/source/calibration/documentation checks pass.
A separately written nonlinear BM3 pressure objective agrees with both linear
fixed-V0 solutions. Report/CSV freshness, scoped Ruff/format checks and source
audit/index checks pass; the comparison figure was rendered and inspected.
Snapshot hashes preserve every preceding material, observation, numerical-audit
and ledger file, including the original source record, conditional parity,
usage qualifications, Li data and separate Litasov evidence. The new
diagnostic is not committed, published or installed in the catalog.

### Conditional K–K′ ellipses and P–V uncertainty bands

These plots use the **stored diagnostic refit covariances**, not an invented
covariance for Sata's published coefficients. The main comparison uses the
recalculated Sata pressures and unchanged Li points. An
[original-pressure version](assets/sata-2010-fe3c-confidence-original.svg)
is supplied separately; its shapes and band widths are nearly identical.
Neither calculation changes any fit or catalog record.

![Conditional joint K–Kprime confidence ellipses and pointwise P–V mean-curve bands for fixed-V0 BM3 and free-Pr mBM3 diagnostic fits](assets/sata-2010-fe3c-confidence-recalculated.svg#only-light)
![Conditional joint K–Kprime confidence ellipses and pointwise P–V mean-curve bands for fixed-V0 BM3 and free-Pr mBM3 diagnostic fits](assets/sata-2010-fe3c-confidence-recalculated-dark.svg#only-dark)

Both **68.3% and 95% nominal confidence levels** are shown. They are conditional
on the equal-weight regression model: independent Gaussian pressure residuals
with a common unknown variance, measured volumes treated as fixed, and V0/Vr
fixed at 149.456 Å³/cell. The residual variance is estimated as RSS/(N−p).
These are explicit diagnostic assumptions; the source's one-/two-sigma
convention and full experimental uncertainty remain unresolved.

For each joint K–K′ ellipse, the displacement d from the fitted center obeys

```text
dᵀ C_pair⁻¹ d ≤ 2 F(2, ν; coverage).
```

Here ν=10 for fixed-V0 BM3 and ν=9 for mBM3. The mBM3 pair is **Kr, Kr′ at Vr**,
with Pr fitted; the ordinary BM3 pair is **K0, K0′ at zero pressure**. The pair
covariance is extracted from the full mBM3 covariance, allowing Pr to vary.
It is not the smaller covariance conditional on a fixed Pr. Since K′ is a ratio
of fitted linear coefficients, the displayed ellipses are local covariance
approximations in physical parameter coordinates. The finite-sample F scaling
accounts for the estimated residual variance. Joint confidence regions do not
have the same interpretation as separate one-/two-standard-error bars; see
[Boggs and Donaldson, NISTIR 89-4102, section 3](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir89-4102.pdf#page=9)
for linearized-region conventions.

The shaded **pointwise confidence intervals for the fitted mean P(V)** use

```text
P̂(V) ± t(ν; (1+coverage)/2) sqrt(g(V)ᵀ C_full g(V)).
```

Pressure is horizontal and volume vertical; the bands therefore show pressure
intervals at each fixed volume. The **full 3×3 covariance**, including Pr and
all its cross terms, is propagated for mBM3. The Student-t scaling follows the
[NIST fitted-mean confidence-interval convention](https://itl.nist.gov/div898/handbook/pmd/section5/pmd511.htm).
Pointwise curve coverage and pair-joint ellipse coverage are different
probability statements. These bands are **not simultaneous whole-curve bands
or prediction intervals for another observation**; no extra residual variance
is added as a future-measurement term. Point-coordinate error bars are omitted
from this confidence comparison, rather than assigned zero uncertainty.

| Recalculated-pressure diagnostic | Conditional K standard error (GPa) | Conditional K′ standard error | K–K′ correlation | 95% pointwise P half-width over measured volumes (GPa) |
|---|---:|---:|---:|---:|
| Fixed-V0 BM3 | 4.20 | 0.0839 | −0.9779 | 0.658–1.658 |
| mBM3, Pr fitted | 12.05 | 0.1790 | −0.9860 | 1.081–1.911 |

The strong anticorrelation allows sizeable changes in K and K′ to compensate
along a narrow family of compression curves. This explains why the mBM3
coefficient ellipse is much larger while the mean P–V band remains relatively
narrow over the sampled volumes. The calculation spans only
**106.449–138.120 Å³/cell**; no ambient extrapolation is included. The bands
omit fixed-anchor, measured-volume, pressure-calibrant and systematic uncertainty
and do not establish physical validation or remove the usage qualifications.

The [confidence-method report](../data/sata-2010-fe3c-confidence.json) records
the covariance matrices, correlations, degrees of freedom, nominal coverage
factors and source-refit checksum. Download the complete
[ellipse coordinates](../data/sata-2010-fe3c-confidence-ellipses.csv) and
[pointwise mean bands](../data/sata-2010-fe3c-confidence-bands.csv), including
both original- and recalculated-pressure fits.

Repository-only reproduction:

```sh
uv run python -m scripts.plot_sata_2010_confidence --write --plot
uv run python -m scripts.plot_sata_2010_confidence --check
uv run pytest -q -W error tests/test_sata_2010_confidence.py
```

Verification: **22** focused covariance/refit/source/documentation checks pass.
Ellipse points satisfy the specified F-scaled Mahalanobis radius. Analytic
physical-parameter derivatives agree with finite differences, and propagated
mean-curve variance agrees with an independent calculation in the linear
coefficient basis, including all mBM3 Pr terms. Report/CSV freshness and scoped
Ruff/format checks pass. The light/dark figures were rendered and inspected;
snapshot hashes preserve every preceding material, observation, refit and audit
data file. No catalog record, classification, recommendation or coefficient was
changed, and nothing was committed or published.

### Original Speziale coefficients, rounding and formulation follow-up

The original [Speziale et al. (2001) author-hosted final PDF](https://duffy.princeton.edu/sites/g/files/toruqf616/files/speziale_et_al-2001-jgrse.pdf)
was checked directly: section 4.2 (p.518), Table 2 (p.520), section 5.3
(pp.520–521), and equations 1–2 (p.516). The PDF SHA-256 is
`c0431bdddf8f4bc8269a7e0d1b79d30d92a38b4a778d9edddad9cb4c45b1f94b`,
matching the completed Speziale audit. Its preferred coefficients are exactly
V0=74.71(1) Å³/cell, fixed isothermal K0=160.2 GPa (adopted uncertainty
0.2 GPa), and K0′=3.99(1). No additional coefficient digits are printed there.
The Table 2 value is explicitly the constrained BM3 fit. The alternative
unconstrained BM3 fit in section 4.2 is V0=74.53, K0=170 and K0′=3.59.

The reproducible follow-up is in
[`sata-2010-mgo-pressure-check.json`](../data/sata-2010-mgo-pressure-check.json),
generated by the repository-only command
`uv run python -m scripts.audit_sata_2010_mgo_pressure --check`.
All comparisons below use the unchanged nine Sata MgO volume/pressure pairs.

| Calculation | Pressure RMS difference (GPa) | Maximum absolute difference (GPa) |
|---|---:|---:|
| Preferred published BM3 | 0.21018 | 0.35264 |
| Same BM3, volume calculated from printed a³ | 0.20646 | 0.34522 |
| BM2 / K0′=4 diagnostic | 0.09877 | 0.15184 |
| Original unconstrained BM3 alternative | 6.80224 | 12.01327 |
| Original Vinet alternative, K0=160.2 and K0′=4.12 | 2.66940 | 4.98543 |
| Original logarithmic alternative, K0=160.2 and K0′=4.21 | 5.26935 | 9.54820 |
| K0′-only BM3 diagnostic, approximately 3.99681 | 0.00586 | 0.00929 |

Section 5.3 does not restate V0 for the Vinet/logarithmic fits; those controls
explicitly adopt the printed 74.71 Å³/cell. They do not recover an unstated
reference volume. Substituting BM3's 3.99 directly into Vinet makes the maximum
discrepancy still larger, about 9.96 GPa. The published alternatives therefore
give no explanation for Sata's specified BM3 pressure column at its high
compressions. The constant-q, final variable-q, and previously audited
local-power gamma controls all reduce to the same BM3 at **300 K**, because
the thermal-energy difference is zero. The historical thermal-integration
question cannot generate this room-temperature offset.

**Ordinary display rounding is too small for a complete nine-row replay.**
At fixed printed coefficients, MgO-volume ±0.0005 Å³ and pressure ±0.005 GPa
rounding allow at most 0.01326 GPa, far below 0.35264 GPa. Coefficient rounding
has a larger high-compression effect: at the last row, increasing V0, K0 and
K0′ individually by half their printed last digits changes pressure by
approximately 0.0521, 0.0581 and 0.2522 GPa. An independently chosen coefficient
box for each row can cover every target, so magnitude alone does not exclude
coefficient precision as a cause.

A bounded minimax diagnostic instead requires **one common coefficient set**
for all nine rows. It permits V0=74.705–74.715, K0=160.15–160.25 and
K0′=3.985–3.995, plus each row's P/V half-digit intervals. This is deliberately
liberal, including K0 although the original fit fixed it. Three starting points
converge to the same remaining maximum interval gap, **0.01156 GPa**. Thus
even this joint nearest-rounding box does not fully reproduce the printed
table. These are deterministic display bounds, not parameter-confidence or
measurement-error intervals; the small remaining gap is not a physical
pressure-scale validation test.

The maintained replay uses complex-step derivatives for the analytic SLSQP
objective and constraints to avoid finite-difference cancellation on the minimum
SciPy stack. The objective, bounds and three starting points are unchanged;
every start must report convergence and satisfy the original feasibility check.
The archived report is retained, and the replay reproduces its gap within
`1e-7 GPa`. This numerical portability fix does not recover a source reducer.

The K0′-only diagnostic **3.99680594** lies inside the published ±0.01 error
width. Its maximum central-value discrepancy is 0.00929 GPa, and every row
overlaps after allowing P/V display rounding. But this value would normally
round to **4.00**, not 3.99. Truncation or another slightly different coefficient
set is plausible, but neither is documented by the sources. Consequently the
result supports a small BM3 coefficient-precision/reduction difference, without
identifying the historical reducer or authorizing replacement of the printed
scale. The MgO card, Sata observations, original Fe3C record and separate
Litasov evidence are unchanged by this follow-up.

Verification: all **22** focused Sata/Speziale tests, the report's `--check`
and scoped Ruff checks pass. Snapshot hashes confirm that all material cards,
the catalog manifest and the preceding Sata/Speziale numerical reports remain
unchanged. This follow-up adds diagnostic evidence only.

The three Li observations retain their **original ruby/neon** pressure basis.
The existing [Li audit](li-2002-fe3c.md#observations-and-calibration) preserves
their neon readings and the limitations of the unrecovered ruby/reduction
workflow. They are never relabelled as MgO-calibrated or silently reduced
onto Litasov's later scale. The overall calibration remains `partially_resolved`
with no complete exact-author pressure-recalculation route asserted.

Sata's separate PM#1/PM#2/PM#3 experiments compare Au, Ar, B2-NaCl, Pt and MgO.
They contextualize relative scale consistency; those are not additional Fe3C
sample measurements or alternative Fe3C EOS parameterizations. No new NaCl,
Au, Pt or Ar catalog record is created by this audit.

## Disposition and reproducibility

The published source BM3 is retained as a **qualified, nondefault static** record on the
existing source-assigned nonmagnetic card. Its validity is restricted to the
combined measured envelope **27.2–186.6 GPa at 300 K**; there is no thermal
extrapolation, ambient nonmagnetic stability claim or default change. Original
record values, datasets and the completed Litasov finding are preserved.
Its uncertainty and scatter qualifications remain explicit; the earlier blanket
negative recommendation is superseded by **qualified use**. The separate
recalculated-pressure refit is a diagnostic and does not replace this source EOS.
The catalog caveat is: **some quoted point errors are smaller than the overall
experimental scatter.** This changes the usage
disposition, not any numerical finding. All 22 focused source/refit/covariance/
documentation checks and the dedicated Rust test pass after this update;
snapshot checks preserve the coefficients, observations, all numerical results,
conditional parity, existing default and all other audit entries. Nothing was
committed or published.

This installed-package example selects the source explicitly:

```python
from peritheos import get_eos_record

# Select the published EOS with its stated validity and uncertainty caveats.
record = get_eos_record("cementite_nonmagnetic_sata_2010_bm3")
volume = record.volume(100.0, 300.0, check_validity=True)
assert abs(record.pressure(volume, 300.0) - 100.0) < 1e-8
```

Repository-only maintenance commands:

```sh
uv run python -m scripts.audit_sata_2010_fe3c --check
uv run python -m scripts.audit_sata_2010_fe3c --check \
  --manuscript output/sata-2010-fe3c-2026-11-01/sources/sata-2010-manuscript.pdf
uv run pytest -q -W error tests/test_sata_2010_fe3c.py
```

`--write` regenerates only the numerical report and two normalized data tables
from retained original tokens/selected Li rows. The append-only integration
script `scripts/integrate_sata_2010_fe3c.py` rejects an already-present record;
it is not an installed-package API or permission to overwrite existing data.

## Integration verification

The initial targeted Python integration run passed **455 tests**, covering the new
source audit, existing cementite/Takahashi records, catalog inventories, native
schema and round trips, reconstruction metadata and documentation. After the
final correction to represent unavailable Li run/measurement identifiers as
missing values, all **13** source, coverage and documentation checks passed.
The dedicated Rust Sata test passed. A freshly built wheel was installed outside
the checkout and verified for explicit record selection, the 300 K and
27.2–186.6 GPa limits, inversion, retained observations and the unchanged default.
Source-token and numerical-report checks, generated inventory/index checks,
scoped report/figure freshness and reconstruction metadata checks passed.
The generated coefficient and pressure–volume/residual figures were rendered
and visually inspected.

Snapshot comparisons preserve every preceding material card, both existing
nonmagnetic EOS records and their datasets, structure/default metadata, and all
682 preceding entries and their order in each audit ledger. The separate
Litasov evidence and MgO card are unchanged.

Broader repository checks are **not all green**. The full Python run stopped
after 466 passes at the existing missing determination-method rationale for
`coesite_iv_v_bykova_2018_experimental_bm2_reference`; full coverage was not
completed. The Rust workspace stopped at the existing brucite retained-isotherm
expectation. Repository-wide lint/format checks also retain unrelated failures;
the changed Python code passes scoped Ruff checks and Python 3.9 mypy passes.
The full refit-ledger freshness check still detects the preceding stored-order
versus generator-order discrepancy. Strict documentation building stops on the
existing Litasov link to a repository script; it reports no new Sata warning.
These checks were recorded without changing unrelated scientific findings or
reordering existing ledger entries. No commit, publication or release was made.

### Earlier recommendation follow-up verification, 1 November 2026 (superseded)

The earlier record, source-integration definition, numerical report, reconstruction
ledger, primary-source audit, material/reconstruction summaries and catalog
notes stated **poor data agreement** and
`not_recommended_for_quantitative_use`. Conditional parameter parity remains
unchanged. The record remains published and nondefault, with explicit
source-curve reproduction/comparison as its only supported use; catalog
discoverability was unchanged. The blanket negative usage recommendation is
superseded by the current qualified-use disposition above; the numerical
findings and all published values remain unchanged.

All **51** focused source, uncertainty-propagation, recommendation and existing
Takahashi checks passed. The dedicated Rust test and a fresh installed-wheel
smoke passed, including the recommendation metadata, preserved parity, source
coefficients, observations, limits, inversion, round trip and existing default.
Scoped Ruff/format checks, Python 3.9 mypy, numerical-report checks, reconstruction
metadata and generated audit/coverage/index checks passed. Snapshot comparisons
preserve all preceding numerical fit/replay results, every source observation,
other material cards and all other 682 entries and their order in both ledgers.
The separate Finger, Litasov and MgO evidence is unchanged.

A broader targeted run had **480 passes and four failures**: three plotting
tests stopped on the installed Matplotlib/pyparsing deprecation warning, and
one global report-freshness test detected the unchanged `cementite.md` report.
The three plotting tests pass with only that dependency-warning category
filtered. The Sata report is current; unrelated reports and dependencies were
left untouched. These results do not claim a clean full repository suite.
The strict documentation build completed with the single existing Litasov
repository-script link warning and no new Sata warning.
