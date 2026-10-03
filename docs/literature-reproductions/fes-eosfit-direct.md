# FeS: direct refits with the official EosFit console

Run date: 2026-10-03. **The data can be refitted in the actual EosFit7c
program. All six runs below converged.** This is separate from an independent
implementation or an extracted CrysFML calculation harness.

The official Mac executable from
[Ross Angel's EosFit distribution](https://www.rossangel.com/text_eosfit.htm)
identifies itself as EosFit7c 7.60, version date 13-May-2021. The
[run manifest](../data/fes-eosfit7c/manifest.json) pins its SHA256 and every
input/output file. The EOS executable was not modified. Its X11 dependencies
were extracted from the official XQuartz package to temporary folders; an
authenticated temporary virtual display was used for its startup. Auxiliary
X11 dependency paths were relocated; no system-wide installation was performed.
The executable emits a font-config warning at startup; fitting and saving
completed. It restores its saved working directory, so the batch runner uses
absolute file names instead of changing the user's directory preference.

## Inputs and selected model

Morard Section III.D describes a staged procedure: first fit the cold BM3
curve to Sata/Ohfuji literature data, then refine only gamma0 for the thermal
component, fixing theta0 and q. It does not describe simultaneous refinement
of all four cold/thermal coefficients against a combined table. Our
thermal-only gamma fits follow that staged parameter selection; the joint
four-parameter runs below are additional diagnostics. Whether the authors
kept the cold rows loaded during their gamma refinement is unspecified.
At the 300 K reference, their thermal pressure is zero, so those cold rows
cannot constrain gamma0 when the cold coefficients are fixed.

Three input selections were fitted with equal weights and, separately, all
supplied errors:

- Cold: all 13 FeS-VI rows from Sata (2010) Table 1 at 300 K; volumes in
  Å³/cell. FeS-VII rows are excluded by phase.
- Morard: all 146 unchanged supplied Table S1 observations; cell volumes and
  errors converted to cm³/mol FeS using Z=4.
- Combined diagnostic: those 146 thermal and 13 cold rows, with all volumes
  in cm³/mol FeS. No inter-study pressure rescaling or new selection is applied.

All runs use third-order Birch–Murnaghan. Thermal runs explicitly select full
MGD, answering **N** to q-compromise, with theta0=417 K, n=2, q=1 and Tr=300 K
fixed. Thermal-only runs refine gamma0 with the printed Morard cold curve
fixed. Combined runs refine V0, K0, K′0 and gamma0. Cold runs refine V0, K0 and
K′0 with no thermal model. Scale factors remain fixed to one.

The installed manual's input prompt defaults to N for q-compromise. Its full
MGD law is
\(\gamma=\gamma_0(V/V_0)^q\) and
\(\theta=\theta_0\exp[(\gamma_0-\gamma)/q]\), with the q=0 limiting law.
The q-compromise approximation instead keeps theta and gamma/V constant and
has no refinable q. Manual pages and equation-image hashes are in the
[manual evidence](../data/eosfit-mgd-manual-audit.json).
These manual defaults do not establish which option the authors selected.

The installed GUI manual adds an important workflow distinction: its MGD and
Einstein parameter **estimation** always uses q-compromise for stability, and
the user can subsequently switch to the general model. Thus the console
prompt's N default does not establish what model a GUI workflow would retain.
This is not evidence of the authors' actual final selection. Setting q=1 in
the full MGD model is also not equivalent to q-compromise: the full model still
has a volume-dependent Debye temperature. The manual evidence additionally
pins the relevant GUI, refinement, input and dataset-scaling documentation.

## Results

The [cold audit](sata-2010-fes-vi.md) gives the source-reference comparison,
including the transformation of correlated errors. The EosFit equal-weight
cold result is V0=98.93066 ± 2.80286 Å³/cell, K0=148.56732 ± 34.40257 GPa and
K′0=4.52374 ± 0.45851; RMS=1.980916 GPa. With P/V errors it becomes
K0=138.01448 ± 16.68440 GPa and K′0=4.67266 ± 0.28835, with weighted
χ²/dof=12.2004. Convergence does not imply that the quoted data errors explain
the cold scatter.

| Thermal selection and objective | V0 (cm³/mol FeS) | K0 (GPa) | K′0 | gamma0 | Pressure RMS |
|---|---:|---:|---:|---:|---:|
| 146 hot rows, equal weights | 15.40 fixed | 115.5 fixed | 4.99 fixed | 2.28211 ± 0.03231 | 1.803098 GPa |
| 146 hot rows, P/V/T errors | 15.40 fixed | 115.5 fixed | 4.99 fixed | 2.29464 ± 0.03830 | 1.804022 GPa |
| 159 joint rows, equal weights | 16.20995 ± 0.14340 | 74.65259 ± 5.89003 | 5.87565 ± 0.20583 | 2.51506 ± 0.07845 | 1.317203 GPa |
| 159 joint rows, P/V/T errors | 15.31635 ± 0.07709 | 116.04768 ± 5.38969 | 5.04478 ± 0.11438 | 2.58429 ± 0.05638 | 1.727727 GPa |

Errors are the program's reported esds. Full saved covariance and output
precision are retained; the covariance is conditional on the selected model,
rows and error treatment. EosFit's iterative weighting is distinguished from
minimizing a parameter-dependent effective-variance objective in SciPy.
The thermal equal-weight gamma result is close to the independent
2.28198 ± 0.03231 replay. Differences at the final displayed digits are
retained rather than forced to bitwise equality.

The initial unchanged Morard coefficient replay in EosFit has RMS=1.910838 GPa,
residual range −3.77819 to +6.14860 GPa, and **15/146 rows exceeding ±3 GPa**.
It agrees with the independent audit's approximately 1.911067 GPa RMS and
6.149320 GPa maximum residual. The residual issue therefore persists in the
actual EosFit executable.

After the gamma-only refit, 18 thermal rows still exceed ±3 GPa under either
weighting. The equal-weight joint fit leaves four thermal rows outside that
bound, plus four cold rows; the weighted joint fit leaves 14 thermal and two
cold rows. None of these diagnostics reproduces the paper's full stated
residual bound.

## Qualifications and warnings

### Same-process staged refit: recovered cold data first

A subsequent check follows the user's requested sequence exactly within one
EosFit process: fit all 13 recovered cold VI points to BM3 in molar units,
retain that fitted curve in memory, load all 159 cold/hot rows, add full MGD,
fix V0/K0/K′0 and refine only gamma0. Theta0=417 K, n=2, q=1 and Tr=300 K
remain fixed. No rounded file reload or manually substituted cold coefficient
intervenes between stages. Both stages converge under both objectives.

| Objective in both stages | Cold V0 (cm³/mol FeS) | Cold K0 (GPa) | Cold K′0 | Stage-two gamma0 | Hot-row RMS | Hot rows beyond ±3 GPa |
|---|---:|---:|---:|---:|---:|---:|
| Equal pressure weights | 14.89429 ± 0.42194 | 148.57274 ± 34.40418 | 4.52367 ± 0.45847 | 2.44764 ± 0.04153 | 2.442846 GPa | 33/146 |
| Supplied errors | 15.01981 ± 0.19732 | 138.00906 ± 16.68608 | 4.67274 ± 0.28837 | 2.52046 ± 0.05190 | 2.159749 GPa | 22/146 |

The largest absolute hot-row residuals are 6.53106 and 6.12060 GPa,
respectively. The cold rows' residuals are exactly unchanged, at printed
precision, from stage one to stage two. Keeping the cold rows loaded therefore
does not add gamma constraints at Tr, though their residuals and row count
affect the reported error scaling/statistics. Stage-two gamma errors condition
on the cold solution and do not propagate its substantial covariance.
The weighted all-data stage again warns about missing cold temperature esds.

The equal-weight gamma is close to Morard's 2.42 ± 0.03, but this recovered
cold table produces different central cold coefficients from their
V0=15.40, K0=115.5 and K′0=4.99. Thus the staged check is a closer procedural
comparison, without reproducing the complete published solution or residual
bound. Complete original cold input/selection and weighting remain unresolved.

The [staged manifest](../data/fes-eosfit7c-staged/manifest.json) retains both
stage results, their actual covariance, inputs, macros, saved models and
output hashes. Rerun with `python -m scripts.run_fes_eosfit_staged`, using
the same `--executable` and new `--output` arguments described below.

### Interpretation of program warnings

The weighted joint run reports χ²/dof=1.5619 and warns that some temperature
esds are zero. These are the cold rows: their source gives a 300 K isotherm
without temperature errors. The input retains zero only as the explicit
isotherm assumption used by this diagnostic; no measured zero uncertainty is
claimed. The raw source manifest keeps the missing uncertainty semantics.
This warning qualifies the joint weighted statistics and errors, even though
the numerical refinement converged.

The thermal setup emits a unit warning before the scales are entered. The
macros then set GPa and cm³/mol before any calculation; every saved model
records the correct scales. This setup warning is retained in the logs.
Saved EOS files can contain Inf in an unused density field because no mass
was supplied. The essential EOS coefficients and refined covariance are
finite; the files are retained exactly as EosFit wrote them.

The joint result is an explicitly assembled diagnostic across studies with
unverified pressure-scale alignment. It is not a recovered author fit and
does not justify replacing the source database parameters. The Morard source
record remains deferred. These initial runs preceded recovery of author input
files; see the [author-input follow-up](#author-input-follow-up-2026-10-03)
for the subsequent supplied files and their distinct uncertainty assignments.

## Q-compromise rerun after the GUI clarification

The user reports that the authors used the GUI. Because its estimation stage
uses q-compromise, all direct and staged procedures were rerun with that
specific model selected. GUI use alone does not establish the final fit model:
the GUI manual allows a subsequent switch to the general model. These are
controlled hypothesis tests, not recovered author settings.

The same official EosFit7c 7.60 executable and exact input observations are
retained. In every thermal setup the q-compromise prompt is explicitly Y.
Saved EOS files record the q-compromise flag in parameter 14; q is undefined
and is neither input nor refined. Theta0=417 K, n=2 and Tr=300 K remain fixed.
All six direct runs and both staged runs converge, with the cold results
unchanged. No original full-MGD outputs are overwritten.

| Thermal calculation | Full-MGD thermal RMS | Q-compromise thermal RMS | Q-compromise largest absolute residual | Q-compromise points outside ±3 GPa |
|---|---:|---:|---:|---:|
| Published coefficients, unchanged | 1.910838 | 2.175382 | 6.74120 | 21/146 |
| Published cold curve; gamma refit, equal weights | 1.803098 | 1.895020 | 5.63200 | 19/146 |
| Published cold curve; gamma refit, supplied errors | 1.804022 | 1.895509 | 5.67390 | 21/146 |
| Recovered cold first; gamma refit, equal weights | 2.442846 | 2.555084 | 6.83235 | 36/146 |
| Recovered cold first; gamma refit, supplied errors | 2.159749 | 2.271071 | 6.22990 | 26/146 |
| Joint four-parameter fit, equal weights | 1.183076 | 1.196703 | 4.07440 | 4/146 |
| Joint four-parameter fit, supplied errors | 1.689198 | 1.758777 | 5.52070 | 18/146 |

All pressure residual values in the table are GPa and use the identical 146
thermal observations. Q-compromise does not reproduce the ±3 GPa bound and
has larger thermal RMS than full MGD in every comparison shown. This result
does not establish which option the authors actually used.

| Parameter | Published | Joint q-compromise, equal weights | Joint q-compromise, supplied errors |
|---|---:|---:|---:|
| V0 (cm³/mol FeS) | 15.40 ± 0.45 | 16.33092 ± 0.16316 | 15.34238 ± 0.08150 |
| K0 (GPa) | 115.5 ± 27.91 | 69.97496 ± 6.39064 | 114.38818 ± 5.60624 |
| K′0 | 4.99 ± 0.51 | 6.01866 ± 0.24339 | 5.07349 ± 0.12062 |
| gamma0 | 2.42 ± 0.03 | 2.43663 ± 0.07514 | 2.48631 ± 0.05209 |

All-row joint RMS values are 1.337304 GPa (equal weights) and 1.791735 GPa
(supplied errors). The staged gamma fits give 2.36739 ± 0.04026 and
2.42977 ± 0.04899, respectively. The latter is close to the reported gamma,
but proximity of one coefficient does not reproduce the cold curve and
thermal residual claim. Weighted runs retain the warning about unreported
cold temperature errors. Parameter errors remain conditional on the selected
model and objective; staged gamma errors do not propagate the cold covariance.

Actual inputs, outputs, macros and saved covariance are retained in the
[direct q-compromise manifest](../data/fes-eosfit7c-q-compromise/manifest.json)
and [staged q-compromise manifest](../data/fes-eosfit7c-staged-q-compromise/manifest.json).
Rerun the existing console/staged runners with `--q-compromise` and a new
output directory.

An additional Peritheos cross-check composes its native BM3 and native Debye
energy evaluated at the reference volume, giving exactly the q-compromise
pressure expression. This is an audit-local composition, not an added native
MGD mode or catalog model. The maximum pressure difference from EosFit's
printed residuals is 0.001079 GPa, with identical rows outside ±3 GPa in every
case. An independent equal-weight joint fit through Peritheos's public fitting
API gives V0=16.330915 ± 0.163166, K0=69.974975 ± 6.390661,
K′0=6.018669 ± 0.243400 and gamma0=2.436493 ± 0.075135; all-row RMS is
approximately 1.337309 GPa. Its residual-scaled covariance is finite and full
rank. As for full MGD, small constant/precision differences do not change the
scientific result. The
[comparison JSON](../data/fes-q-compromise-peritheos-comparison.json) retains
the independent fit and pointwise comparisons. Reproduce it with
`python -m scripts.compare_fes_q_compromise` (or verify with `--check`).

## Reproduction

### Peritheos evaluation and independent joint fit

The public Peritheos 0.11.0 API was also run on the exact retained EosFit
input files. Forward evaluation uses the native Rust BM3 and full integrated
MGD model, with the same n=2, Tr=300 K, theta0=417 K and q=1. EosFit's
cm³/mol volumes are divided by ten to obtain Peritheos's J/bar/mol units.
The conventional-cell catalog adapter agrees with this direct API evaluation
to 4.5×10⁻¹³ GPa. Catalog validation is bypassed explicitly for this audit of
the deferred source record; its disposition and coefficients are unchanged.

All eight retained direct/staged cases were evaluated at their initial and
final coefficients. The largest pointwise difference from the actual EosFit
printed residuals is **0.001060 GPa**, and every case gives exactly the same
rows outside ±3 GPa. Peritheos's native and thermal Python fallback evaluations
agree within 4.3×10⁻¹⁴ GPa.

| Model, evaluated at EosFit coefficients | EosFit RMS, all 159 rows | Peritheos RMS, all 159 rows | EosFit thermal RMS | Peritheos thermal RMS |
|---|---:|---:|---:|---:|
| Published Morard coefficients | 1.947373 GPa | 1.947577 GPa | 1.910838 GPa | 1.911067 GPa |
| Joint equal-weight solution | 1.317203 GPa | 1.317205 GPa | 1.183076 GPa | 1.183078 GPa |
| Joint supplied-error solution | 1.727727 GPa | 1.727758 GPa | 1.689198 GPa | 1.689230 GPa |

Default Peritheos uses R=8.31446261815324 J/mol/K. A separate diagnostic
setting Cvmax=6×8.314 reduces the largest published thermal-replay difference
from approximately 0.0010 to 0.00011 GPa, consistent with the rounded gas
constant found in the audited CrysFML source. This diagnostic does not replace
the default Peritheos result or establish the exact constants in the binary.
Printed-coefficient/output precision and implementation details prevent
claiming bitwise agreement.

An independent equal-weight four-parameter refit through
`peritheos.fitting.fit_joint_eos` also converges from the published starting
coefficients. Both programs use the same unweighted pressure-residual
objective, with the cold and thermal parameters released together:

| Parameter | Actual EosFit refit | Independent Peritheos refit |
|---|---:|---:|
| V0 (cm³/mol FeS) | 16.209950 ± 0.143400 | 16.209968 ± 0.143416 |
| K0 (GPa) | 74.652590 ± 5.890030 | 74.651843 ± 5.890344 |
| K′0 | 5.875650 ± 0.205830 | 5.875675 ± 0.205857 |
| gamma0 | 2.515060 ± 0.078450 | 2.514903 ± 0.078451 |

The independently fitted Peritheos model has all-row RMS=1.317205 GPa and
thermal RMS=1.183077 GPa. Its four-parameter covariance is finite and full
rank, with 155 degrees of freedom; errors are conditional and scaled from
residual scatter. This verifies numerical evaluation and the equal-weight
refit, not independent physical validation or recovery of the author's inputs.
The weighted EosFit solutions were forward-evaluated only: Peritheos's
errors-in-variables fitter adjusts latent V/T coordinates, whereas EosFit's
iterative weighting uses a different objective. Their weighted parameter
refits are not claimed equivalent.

The [comparison JSON](../data/fes-peritheos-eosfit-comparison.json) retains
every pressure comparison, input hashes, native extension hash, fit covariance
and convergence. Reproduce it with
`python -m scripts.compare_fes_peritheos_eosfit`, or use `--check` to verify
the retained output with the same build and numerical environment.

Each case directory contains `input.dat`, `run.mcr`, `initial.eos`,
`fitted.eos`, `run.log`, `stdout.txt` and `stderr.txt`. The manifest retains
all residuals, program convergence, final-cycle coefficients, esds, covariance,
warnings and byte hashes. Logs use Pobs−Pcal; the manifest explicitly reverses
this to Pcal−Pobs to match the independent audits.

With a runnable official console and its display/library environment set:

```sh
python -m scripts.run_fes_eosfit_console \
  --executable /absolute/path/to/eosfit7c \
  --output /absolute/path/to/a/new/run-directory
```

The output directory must be new. The runner generates fresh absolute paths,
runs the real console and stops on a failed save, unexpected parameter set
or non-convergence. Archived macros contain the original absolute run paths;
regenerate them with the runner when moving the files.

Focused tests verify observation preservation, phase selection, derived
unit conversions, exact input/output hashes, saved model settings, covariance
extraction and actual convergence. The native library and independent
finite-pressure benchmark are tested separately.

## Author-input follow-up (2026-10-03)

The user subsequently supplied `Supplementary_Materials_Rev_140526.docx`,
`FeS6.eos`, `FeS6EOSfittot-SataOhfuji.dat` and
`FeS6EOSfittot-SataOhfuji-300K.dat`. This materially changes the available
evidence. The [author-source manifest](../data/fes-author-eosfit/manifest.json)
retains byte-identical EOS/data files, checksums, all actual console outputs,
refit macros, final coefficients, covariance and row residuals. Original files
were not edited. These are supplied author inputs, but their final publication
status, fit mask and fitting history are not encoded in the files.

### Saved model resolves the q-compromise hypothesis

`FeS6.eos` identifies EosFit7-GUI version 20210609 and a save date of
12 March 2024. It contains BM3 plus **full MGD**, parameter 14=0,
q=1, theta0=417 K and n=2. Thus GUI use did not imply that this saved model
used q-compromise. The earlier q-compromise runs remain hypothesis tests.

| Setting | Printed article / earlier replay | Supplied author EOS |
|---|---:|---:|
| V0 (cm³/mol) | 15.40 | 15.3700 |
| K0 (GPa) | 115.5 | 115.490 |
| K0′ | 4.99 | 4.99000 |
| gamma0 | 2.42 | 2.41671 |
| Reference temperature (K) | 300 in initial replay; 300 K isotherm shown in Figure 6 | Tref=298 explicitly |

The paper does not explicitly print an EosFit Tref setting. Thus 300 K is
the reference adopted in the initial audit from the cold-data isotherm, not
a directly recovered software setting. K0 and gamma0 in the author file agree
with the article at its printed rounding; V0=15.3700 does not round to 15.40.

Only the saved gamma variance is nonzero: parameter (18,18)=0.0016139,
corresponding to an esd of 0.04017. This is consistent with a conditional
gamma-only refinement with fixed cold parameters. It does not supply their
uncertainties or prove the exact data selection. The saved model predates the
May 2026 supplement revision, so it is not automatically the final paper fit.

### Both input files contain combined hot and cold data

| Author input | Thermal rows | Literature 300 K rows | Additional 300 K rows | Total |
|---|---:|---:|---:|---:|
| FeS6EOSfittot-SataOhfuji.dat | 146 | 21 | 11 | 178 |
| FeS6EOSfittot-SataOhfuji-300K.dat | 146 | 21 | 0 | 167 |

On 2026-10-03 the user clarified that Guillaume Morard confirmed in personal
communication that the eleven additional cold observations were excluded from
the paper's final fit because their quenched values were unsuitable. This is
user-reported author communication, not a statement extracted from the article.
The 167-row variant is therefore the author-confirmed selection for the
exclusion comparison, while the 178-row variant is an inclusion sensitivity
test. This resolves the exclusion of these eleven; it does not establish the
final weighting or every saved model setting.

The second file removes precisely the eleven additional 300 K rows; it is not a
cold-only input. Their pressures span 62.269–146.495 GPa, sigmaP=0.0056 GPa
and sigmaT=5 K. These errors are retained as supplied, without endorsing their
physical uncertainty. EosFit's combined weights also propagate volume and
temperature errors; the pressure error alone is not the effective fit weight.

Both files retain the same 146 thermal pressure/temperature observations as
Table S1, in the same order, including its duplicate seventh/eighth rows.
Volumes use the rounded factor 0.15055 cm³/mol per Å³ cell, rather than the
exact NA/(4×10²⁴). The largest difference is 0.0003248 cm³/mol and does not
resolve the residual discrepancy.

Thirteen of the 21 literature rows match the FeS-VI cell volumes originally
reported in Ohfuji (2007) Table 2 and re-reported in Sata (2010). Direct access
to the original publisher PDF on 2026-10-03 confirms all thirteen original
run/experiment identities. The other eight cold observations were absent from our earlier
reconstruction; their individual source assignment is not encoded in the
combined input. The author file uses 101.2 GPa for one matched row whose Sata
table pressure is 101.1 GPa. All 21 rows are assigned sigmaP=1% of P and
sigmaV=0.5% of V, with sigmaT=5 K. These are different from the original Sata
table error tokens used in the earlier weighted diagnostics.

### Identity of the eleven additional cold observations

All eleven match 300 K FeS-VI lattice refinements in the public `FeS-GSAS.csv`
deposited with Morard's study. This identifies where the values appear; it
does not by itself establish who originally measured them. The paper explicitly
states that the cold fit uses literature data from Sata (2008)/Ohfuji (2007)
(Section III.D), and the Figure 6 caption identifies its 300 K observations
as literature data. The additional eleven author-input points are not
individually attributed there. Their original-study attribution is therefore
kept unconfirmed; the earlier label "own cold points" was too strong.

A direct check of [Ohfuji (2007) Table 2](https://link.springer.com/article/10.1007/s00269-007-0151-0/tables/2)
on 2026-10-03 rules out a direct copy of its published cold FeS-VI observations:
the table contains thirteen 300 K FeS-VI rows, and all thirteen match the separate
literature block after the author's volume conversion. None of the additional
eleven matches even a published cold volume within its 0.01 Å³ printing precision;
the smallest difference is 0.060996 Å³, greater than the 0.005 Å³ rounding bound.
Changing only the pressure calibration therefore cannot turn these eleven into
the published Ohfuji cold rows. The additional volumes also match none of the
eleven hot FeS-VI rows in Table 2, excluding simple relabelling of their temperatures.
This is an observation-identity check, not a test
of statistical consistency or proof of independent experiments. Re-refinement
of unpublished diffraction data cannot be excluded. The original methods and
Table 2 footnote specify MgO for runs 1/4 and Ar for runs 2/3, whereas the eleven
additional rows match paired KCl markers in the Morard deposit.
The [numerical comparison](../data/fes-ohfuji-2007-identity-check.json) retains
original error tokens, run IDs, source-PDF checksum, and each comparison. The
Sata/Ohfuji (2008) full paper remained unavailable through the current De Gruyter
institutional licence; its authors' SPring-8 report and the Sata (2010) re-reporting
were consulted. No broader exclusion of unpublished or differently refined
Ohfuji data is asserted.

Matching lattice
products reproduces their author molar volumes within 5×10⁻⁹ cm³/mol using
the author's factor 0.15055. Paired KCl markers reproduce their pressures
within 3×10⁻⁸ GPa using the same numerically reconstructed calibration law
as the hot Table S1 rows. The additional-cold dataset's `provenance.gsas_crosscheck`
retains scan identifiers, source lines, lattice values and pressure residuals.

| GSAS scan | Author-input P (GPa) | V (cm³/mol) |
|---|---:|---:|
| #115 | 62.26919926 | 11.82844615 |
| #116 | 63.92562005 | 11.77763301 |
| #117 | 64.11165667 | 11.76845913 |
| #118 | 66.77626827 | 11.70991616 |
| #20 | 117.626888 | 10.52028399 |
| #21 | 118.2310956 | 10.51060996 |
| #22 | 122.9684092 | 10.42257027 |
| #23 | 123.16448 | 10.42009851 |
| #26 | 126.2360703 | 10.38420076 |
| #48 | 146.4952805 | 10.06219656 |
| #49 | 146.4682671 | 10.0611115 |

All temperatures are 300 K, with assigned input errors sigmaT=5 K and
sigmaP=0.0056 GPa; volume errors are retained separately per point. In the
GSAS export these groups follow hot scans (e.g. #113 at 2540 K precedes
#115 at 300 K). A post-heating interpretation is plausible from that ordering,
but it is an inference, and neither authorship nor the full measurement
protocol is explicitly established for these input points.
These points are absent from the thermal-only Table S1 and the 167-row
author variant.

Archived EosFit run-directory keys such as `without-own-cold` retain their
original audit names to preserve exact program outputs. They denote the
167-row input variant and do not assert authorship of the removed rows.

### Actual EosFit replay and refits

Unmodified official EosFit7c 7.60 loaded the author EOS and original input
files. Native Peritheos agrees within **0.001066 GPa** on either complete
input. Calculated-minus-observed residuals for the 146 hot rows have
RMS=**1.6317545 GPa**, maximum absolute residual=**5.0579 GPa**, and
**14 rows outside ±3 GPa**. The whole-input RMS is 1.904239 GPa for 178 rows
and 1.669359 GPa for 167 rows. Even the hot subset at P≥40 GPa has seven
out-of-bound rows, so the stated lower pressure limit alone does not explain
the discrepancy. The earlier printed-coefficient replay had hot RMS=1.910838
GPa and maximum=6.1486 GPa on the exact-SI Table S1 conversion.

Cold-only BM3 fits exclude the thermal term and refine V0, K0 and K0′.
Thermal fits load the supplied full-MGD EOS with theta0, q, n and Tref fixed.
The gamma-only diagnostics hold the saved cold coefficients fixed; the joint
diagnostics additionally refine V0, K0 and K0′. All 14 refits converged.

| Cold observations / weighting | V0 (cm³/mol) | K0 (GPa) | K0′ |
|---|---:|---:|---:|
| 21 literature / equal pressure weights | 15.40351 ±0.45392 | 113.26134 ±27.82332 | 5.03534 ±0.50650 |
| 21 literature / supplied errors | 15.37078 ±0.52385 | 115.43999 ±34.78356 | 4.99168 ±0.67891 |
| 32 including additional 300 K / supplied errors | 15.32105 ±0.37797 | 135.04871 ±27.91409 | 4.43684 ±0.42528 |

The 21-row error-weighted cold refit returns coefficients close to the supplied
EOS; the equal-weight errors are close to the printed cold-parameter errors.
This supports the importance of the previously missing cold rows and weight
assignments. Neither observation identifies a unique original fitting procedure.

| Joint fit | V0 (cm³/mol) | K0 (GPa) | K0′ | gamma0 | Hot RMS (GPa) | Hot max abs (GPa) | Hot rows outside ±3 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 167 rows, equal pressure weights | 16.20277 ±0.14188 | 74.61976 ±5.85635 | 5.88503 ±0.20293 | 2.53466 ±0.06710 | 1.181603 | 4.0074 | 4 |
| 167 rows, supplied errors | 16.20742 ±0.26487 | 73.81297 ±11.31856 | 5.93481 ±0.42103 | 2.56652 ±0.10897 | 1.183081 | 4.0241 | 4 |
| 178 rows, equal pressure weights | 16.25331 ±0.17028 | 76.56376 ±7.03096 | 5.71226 ±0.22453 | 2.26228 ±0.07035 | 1.233817 | 4.2082 | 4 |
| 178 rows, supplied errors | 16.69034 ±0.27664 | 66.15215 ±9.42460 | 5.92847 ±0.33029 | 1.67085 ±0.04621 | 1.622611 | 4.6327 | 10 |

The gamma-only supplied-error fit gives 2.46317 ±0.03949 for 167 rows and
2.46943 ±0.21942 for 178 rows. Neither exactly reproduces the saved gamma0.
Volume-only and pressure-only weighting diagnostics are also retained in the
manifest. All errors are EosFit's conditional esds, with its saved covariance;
they do not propagate study calibration uncertainty or establish physical
validation. The fixed-width saved covariance has limited precision; highly
correlated cold matrices can acquire small negative eigenvalues from that
rounding. They are retained as written, without positive-definite repair.
A numerical fit succeeds, but these selected procedures do not
reproduce every printed coefficient or the strict residual bound.

### Supplement inspection

All five rendered pages were inspected. The DOCX contains Figures S1–S5 and
their captions, no numerical tables, and no Table S2 despite its reference in
the S2 caption. Figure S5 explicitly plots calculated minus observed pressure.
Its raster visibly includes residuals near +3.6 and −3.2 GPa, already beyond
a strict ±3 GPa claim. Its cold-point pressures overlap the additional 300 K group
in the 178-row input. That overlap does not prove the complete plotting mask,
model coefficients or which input version generated the figure. The exact
author EOS/input replay still does not reproduce this narrower displayed range.
The [inspection record](../data/fes-author-eosfit/supplement-inspection.json)
retains source metadata, checksum and the visual qualification.

The source catalog entry stays deferred and its published parameters are
preserved. Three additional [typed reference datasets](../datasets.md#fes-reference-data-and-author-input-variants)
expose the selected 167-row author input, its 21 literature cold observations and
original Sata VI transcription. The user reports Guillaume Morard's personal
communication confirming exclusion of the 11 unsuitable quenched observations.
Neither those 11 rows nor the 178-row variant containing them is a database dataset.
Original selected author bytes are packaged in
`peritheos/data/datasets/morard_2026_author_sources`; rejected-input originals and
historical inclusion/exclusion comparisons remain only in the audit archive.
Shared observation IDs and uncertainty metadata identify overlapping subsets.
The paper coefficients stay on the deferred source record; saved author
coefficients are additional provenance. Regenerate/check packaging with
`python -m scripts.bundle_fes_reference_datasets` / `--check`.

The following command regenerates the **historical** two-variant diagnostic,
including the rejected observations. For the current selected-data audit use
`python -m scripts.run_fes_author_staged` below.

Reproduce these additional diagnostics with:

```sh
python -m scripts.audit_fes_author_inputs \
  --executable /absolute/path/to/eosfit7c \
  --eos /absolute/path/to/FeS6.eos \
  --inputs /absolute/path/to/FeS6EOSfittot-SataOhfuji.dat \
           /absolute/path/to/FeS6EOSfittot-SataOhfuji-300K.dat \
  --output /absolute/path/to/a/new/author-audit-directory
```


### Selected-data cold-then-thermal fits (2026-10-03)

The user reports Guillaume Morard's personal communication confirming exclusion
of the 11 unsuitable quenched observations and a cold-first, fixed-cold thermal
procedure. This audit uses only `FeS6EOSfittot-SataOhfuji-300K.dat`: 21 literature
cold observations at 300 K and the unchanged 146 thermal observations. The
cold rows and hot rows are separated without changing numeric tokens or order.
The thermal stage fits only the 146 hot rows. No quenched observations enter
either stage or the database datasets.

Actual official **EosFit7c 7.60**, with full BM3-MGD (`Thermal=7`, `Param14=0`),
performs both stages in one process. The published-cold branch fixes the printed
V₀/K₀/K′ coefficients directly; the reproduced-cold branch first refines these
three parameters on the 21 cold rows and carries the internal solution into
the thermal stage without saving/reloading rounded parameters. Both branches
refine **only γ₀**, holding θ₀ = 417 K, q = 1 and n = 2 fixed as in Section III.D.
Tr = 300 K is the primary comparison, consistent with the cold isotherm;
Tr = 298 K is an explicitly separate saved-author-file control.

| Cold curve fixed during thermal fit | Weighting | V₀ (cm³/mol) | K₀ (GPa) | K′ | Fitted γ₀ ± EosFit esd | Thermal pressure RMS (GPa) | Max absolute ΔP (GPa) | Rows outside ±3 GPa |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Published cold coefficients | Supplied P/V/T errors | 15.40000 | 115.50000 | 4.99000 | 2.29252 ± 0.03829 | 1.80631 | 5.5598 | 18/146 |
| Reproduced cold fit | Supplied P/V/T errors | 15.37080 | 115.43890 | 4.99170 | 2.46418 ± 0.03953 | 1.63361 | 5.2858 | 12/146 |
| Published cold coefficients | Equal pressure weights | 15.40000 | 115.50000 | 4.99000 | 2.28010 ± 0.03235 | 1.80539 | 5.5066 | 19/146 |
| Reproduced cold fit | Equal pressure weights | 15.40339 | 113.26873 | 5.03520 | 2.44884 ± 0.02904 | 1.60110 | 5.1461 | 10/146 |

The paper reports **γ₀ = 2.42 ± 0.03**. The weighted staged reproduction differs
by +0.04418; the equal-weight staged reproduction differs by +0.02884. Holding
the published cold coefficients instead gives −0.12748 (weighted) or −0.13990
(equal weights). These differences are conditional comparisons: γ errors do
not propagate cold-parameter covariance or inter-study calibration uncertainty.
The RMS column is the unweighted pressure-residual RMS, even for weighted fits;
EosFit's actual objective uses the selected uncertainty contributions and its
parameter-dependent effective weights.

The 21-point weighted cold fit gives V₀ = 15.37080 ± 0.52385 cm³/mol,
K₀ = 115.43890 ± 34.78351 GPa and K′ = 4.99170 ± 0.67891. The equal-weight
cold fit gives V₀ = 15.40339 ± 0.45392 cm³/mol, K₀ = 113.26873 ± 27.82267 GPa,
K′ = 5.03520 ± 0.50650. Equal weighting more closely recovers the printed
cold parameter errors; this numerical resemblance does not establish the
publication's regression settings.

Using the published cold coefficients and γ₀ = 2.42 directly, the same 146-row
EosFit replay at Tr = 300 K gives RMS **1.91610 GPa**, maximum **6.1620 GPa**
and **15/146** rows outside ±3 GPa. The previously reported RMS 1.63175 GPa
and maximum 5.0579 GPa belong to the different **saved author model**
(V₀ = 15.37, K₀ = 115.49, K′ = 4.99, γ₀ = 2.41671, Tr = 298 K), not the
printed-coefficient model. The staged weighted reproduced-cold fit therefore
approaches the saved model's residual RMS, but does not reproduce the
paper's claimed complete residual envelope.

At Tr = 298 K, the weighted published-cold branch gives γ₀ = 2.29004 ± 0.03824,
RMS 1.80576 GPa, maximum 5.5594 GPa and 18 outliers. The weighted
reproduced-cold branch gives γ₀ = 2.46155 ± 0.03948, RMS 1.63273 GPa,
maximum 5.2840 GPa and 12 outliers. Changing Tr by 2 K does not resolve
the mismatch. None of these staged fits is promoted to a source EOS record.

Independent native Peritheos full-MGD replay agrees with all six EosFit fits
within **0.00112 GPa** for every pressure residual, using the final printed
coefficients. This verifies numerical pressure evaluation, not independent
physical validation. The source coefficients remain unchanged and deferred.

The [machine-readable audit](../data/fes-author-staged/manifest.json) preserves
source and executable SHA-256 values, all six macros and console outputs,
full saved covariance matrices, cold solutions, conditional thermal esds,
individual residuals and native replay checks. Regenerate into a new directory:

```sh
python -m scripts.run_fes_author_staged \
  --executable /absolute/path/to/eosfit7c \
  --output /absolute/path/to/a/new/staged-audit-directory
```


### GUI 20210609 versus the audited console: version check (2026-10-03)

The supplied `FeS6.eos` explicitly identifies **EoSFit7-GUI Program (20210609)**.
The console used for the official fits is **EosFit7c 7.60, dated 13-May-2021**.
The current official Mac installer also contains a GUI binary whose embedded
banner is **20211014**. These are different application/build labels within
the 2021 release family, not evidence that the audit used a recent CrysFML
numerical engine against an old author model.

| Release/change | Documented numerical relevance to this FeS audit |
|---|---|
| June 2017, v7.4 | MGD thermal-pressure EOS introduced |
| August 2019, v7.5 | GUI MGD fitting bugs and unweighted parameter-esd bug fixed; calculator MGD parameter editing/saving fixes also predate the author's build |
| Summer 2021, v7.6 | Adds q-compromise MGD, improves volume inversion and refinement termination; both the author's GUI and audited console belong to this release family |
| Late 2021 Mac GUI | Bundled manual says the re-release is identical to summer 2021 apart from a graphics-library crash fix |
| August 2023 Windows update | Documented Windows 11 calculator compatibility fix; no changed MGD pressure law reported |

The v7.6 numerical bug fixes concern fourth-order **linear** EOS moduli and
thermal-pressure refinement to **T–V-only** data. Neither describes the volumetric
BM3 fit to the supplied P–V–T observations. The manual warns that termination
changes relative to older versions can shift fitted parameters by small fractions
of an esd; this is not a documented June-to-October 2021 change.
The manual's Mac re-release label is November 2021; the current website calls
the graphics fixes December 2021. The executable's October build banner is
reported separately rather than treating these dates as identical.
Sources: [official release page](https://www.rossangel.com/text_eosfit.htm),
[bundled Version Info](../data/fes-eosfit-version-audit/VersionInfo.txt).

The public CrysFML repository snapshots immediately before June 9 and October 14,
2021 have **byte-identical** `CFML_EoS.f90`, `CFML_Math_General.f90`,
`CFML_Optimization_LSQ.f90` and `CFML_LSQ_TypeDef.f90`. No commits to these
files are recorded between those dates in the inspected repository history.
The June-era source already contains full integrated MGD and q-compromise,
with parameter 14 selecting the latter. The supplied author EOS has parameter
14 = 0 and therefore selects full MGD; GUI estimation's q-compromise default
does not override this explicit saved model.

In addition, the current distributed GUI and audited console have matching
address-normalized disassembly for `Get_DebyeT`, `Get_Grun_V`, `Pthermal`
and `EthDebye`. This supports the same compiled MGD pressure calculation.
Normalization removes instruction addresses and hexadecimal operands; this
comparison does **not** prove identical constants or whole-program behavior.
The GUI and console refinement-control routines differ, and the console has
additional group-scale fitting support. Agreement of pressure routines does
not establish identical weighting defaults, selected data, or GUI fit controls.
The [GUI methodology paper](https://doi.org/10.1107/S1600576716008050) describes
shared pressure least-squares code and selectable weighting, rather than an
alternative GUI EOS definition.

**Conclusion:** no documented change from GUI 20210609 to the later 2021
release explains the FeS residual mismatch. This is strong version-history
and pressure-engine evidence, but the exact original **20210609 GUI executable
has not been run**. A GUI-specific regression/default/selection issue therefore
remains an untested possibility. Later public CrysFML changes, including June
2022, must not be assumed to be compiled into the distributed 2021 programs.

The [version audit manifest](../data/fes-eosfit-version-audit/manifest.json)
records the binary/manual/source hashes, repository commit IDs, normalized
instruction artifacts, chronology and limitations. Original author files and
database coefficients remain unchanged; the 11 rejected quenched observations
remain excluded from current database datasets and fits.


### Author-used inputs confirmed by the user (2026-10-03)

The user confirms that the data and EOS files provided here are all inputs used
by Guillaume Morard. The audit treats them as the **complete author-used inputs**;
requesting those same files again is not an appropriate next step. The current
accepted selection still excludes the 11 unsuitable quenched observations.

The remaining discrepancy is numerical/source-output reproducibility: the
supplied saved model and printed parameters differ, the saved-model replay
exceeds the stated residual envelope, and the staged console fits do not recover
all reported outputs. Convergence and native pressure agreement establish
successful fitting and pressure evaluation; they do not establish the exact
GUI fit workflow. The empirical original-GUI versus console comparison remains
pending. The release/source/binary comparison above supports equivalent MGD
pressure evaluation, but does not resolve GUI weighting, group selection,
estimation or refinement controls.

The database retains the published source coefficients as deferred, the
accepted selected datasets, and explicit confirmation of the author-used input
provenance. No replacement fit is promoted automatically and no further
missing-data request is implied by the deferred status.


### Parameter agreement within reported error bars (2026-10-03)

The reproduced-cold staged fits are **similar within reported uncertainties**.
This is a more accurate parameter-level conclusion than an unqualified statement
that the fit was not reproduced. The weighted and equal-weight cold refits put
all three cold coefficients inside the published parameter error widths.

| Parameter | Published | Weighted reproduced-cold staged fit | Equal-weight reproduced-cold staged fit |
|---|---:|---:|---:|
| V₀ (cm³/mol) | 15.40 ± 0.45 | 15.37080 ± 0.52385 | 15.40339 ± 0.45392 |
| K₀ (GPa) | 115.5 ± 27.91 | 115.43890 ± 34.78351 | 113.26873 ± 27.82267 |
| K′ | 4.99 ± 0.51 | 4.99170 ± 0.67891 | 5.03520 ± 0.50650 |
| γ₀ | 2.42 ± 0.03 | 2.46418 ± 0.03953 | 2.44884 ± 0.02904 |

The equal-weight γ₀ estimate is inside the published interval [2.39, 2.45].
The weighted γ₀ central estimate is outside that interval, but its conditional
error interval [2.42465, 2.50371] **overlaps** the published interval. Its difference
from the published estimate is 0.04418, compared with an error quadrature of
0.04963. Because the estimates share observations and their cross-covariance
is unknown, this quadrature comparison is descriptive, not an independent
statistical significance test. The confidence convention of the source errors
is also unspecified. Thermal esds do not propagate the cold-fit covariance.

EosFit reports weighted χ²/dof = **0.6274** for the reproduced-cold thermal stage.
The residual maxima above 5 GPa do not imply inability to fit the supplied data
or numerical failure; their conflict with the stated complete ±3 GPa envelope
is a **separate residual-reporting discrepancy**. Exact GUI/source output parity
remains unestablished. Conversely, holding the printed cold coefficients fixed
and refitting only γ₀ gives 2.29252 ± 0.03829, whose error interval does not overlap
the published γ₀ interval; that branch should not be conflated with the
cold-then-thermal reproduced-cold fits.

The [parameter comparison](../data/fes-author-staged/parameter-comparison.json)
records interval overlap, central-value containment, descriptive error ratios
and the separate residual diagnostics. The fit reproduction is classified **`similar`** in both the material metadata
and primary refit ledger. The separate residual-claim status remains
`not_reproduced`. The source card retains its existing deferred availability
status; parameter agreement and the remaining source qualifications are
recorded separately. Neither qualification means EosFit cannot fit the data.
