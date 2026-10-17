# Farla 2025 B2 rubidium halides

Initial audit: 2026-10-15; discrepancy follow-up: 2026-10-16.
Three published preferred BM3-MGD records are added:
`rbcl_b2_farla_2025_bm3_mgd_cscl`,
`rbbr_b2_farla_2025_bm3_mgd_cscl`, and
`rbi_b2_farla_2025_bm3_mgd_cscl`. RbCl's Campbell-Heinz record and default
selection are preserved. RbBr and RbI are separate diffraction-ready B2
materials. No independent-fit coefficients are substituted into these records.

**Fit parity remains incomplete.** The follow-up recovers historical CsCl
calibration code and defaults, reproduces 104 of 109 distinct gauge pressures
within printed pressure/volume rounding, and improves the independent refits
to within 1.63 reported parameter standard errors. Original EoSFit7 files and
exact fitted-coefficient parity remain unavailable. The subsequent discrepancy
audit localizes RbCl's remaining benchmark inconsistency to 5 GPa/1500 K;
the other 53 entries admit one rounding-compatible coefficient vector.
Direct EosFit7c 7.60 refits on 16 October confirm the independent refitting
results; using the original fitting program does not restore publication-digit parity.
The separate constrained Case 3 fits now agree within 0.085 reported thermal
parameter standard errors. Their remaining last-digit differences can be
accounted for by precision-compatible fixed DAC anchors, as demonstrated below.

## Assessment and follow-up decision (2026-10-16)

The remaining differences are small numerical-reproduction limits. They do
not, by themselves, establish a substantive error in the published equations
of state or justify a correction request. Exact author-fit parity remains
unachieved; close agreement does not recover the original inputs or covariance.

| Comparison | Remaining difference | Interpretation |
|---|---|---|
| Preferred Case 1 | Largest individual coefficient differences are 0.40, 1.21 and 1.63 published marginal standard errors for RbCl, RbBr and RbI, respectively. | Modest differences between reconstructed and published fits. These marginal comparisons are not a joint significance test; original covariance, unrounded inputs and row flags remain unavailable. |
| Constrained Case 3 | Every fitted thermal coefficient differs by less than 0.085 published marginal standard errors. | Negligible relative to reported parameter errors. Illustrative fixed DAC anchors within printed rounding intervals make all four thermal pairs round to the publication; those anchors are not recovered author values. |
| CsCl calibration | 104/109 distinct readings match pressure/volume rounding. The five remaining central-value differences reach 0.07090 GPa, against reported pressure errors of 0.07 GPa. | Differences are comparable to the stated uncertainty scale; all readings agree within reported pressure errors plus rounding. The exact original calibration session remains unrecovered. |
| RbCl Table S5 benchmark | At 5 GPa/1500 K, a joint reconstruction consistent with the other 53 ratios gives about 0.92828 (rounding to 0.9283), versus printed 0.9282. | Approximately 0.01% in relative volume. This is a last-digit consistency diagnostic under stated constraints, not a confirmed table typo or a recovered author model. |

**Current decision: document the qualifications and defer author contact.**
The user chose on 16 October 2026 not to contact the authors for now. No
message has been sent and no outreach is planned under the current scope.
Original EoSFit files, full-precision fixed anchors, observation flags/errors,
covariance and saved CsCl inputs remain evidence gaps, rather than pending
author requests. Reconsideration would require a new decision to pursue
original-input recovery or new evidence of a material scientific discrepancy.

Published coefficients, source observations and the original S5 value remain
unchanged. Diagnostic fits, covariance and rounding examples remain separately
labelled. Equation implementation agreement, exact fit reproduction, calibration
replay and independent physical validation retain their distinct meanings.

## Primary evidence and access

R. Farla, A. Néri, M. Pöppelbaum and K. Glazyrin (2025), *Thermal equations
of state of B2-structured rubidium halides RbCl, RbBr, and RbI*, Journal of
Applied Physics 137, 115901, [doi:10.1063/5.0248905](https://doi.org/10.1063/5.0248905).

The source equation and parameter audit uses the official DESY
[final accepted manuscript](https://bib-pubdb1.desy.de/record/616759/files/Final%20version%20before%20publication.pdf),
listed as “Final version before publication,” deposited 17 March 2025.
The web reader exposed its complete text. Direct downloads returned an HTML
challenge. The publisher article/PDF returned HTTP 403; the restricted DESY
version-of-record file could not be retrieved. **A final publisher PDF was not
obtained, and this audit does not imply a page-by-page comparison of versions.**

The official, publicly downloadable
[AIP Figshare supplement v1](https://doi.org/10.60893/figshare.jap.28515917.v1)
was retrieved through the public Figshare API, which links its resource DOI
to this article. File 52698776 is a DOCX containing seven tables, figures and
Debye-temperature discussion, licensed CC BY 4.0 in the API metadata.
It was parsed from OOXML and rendered for visual checks of table headings,
errors, blank cells and benchmark/melting footnotes. The DESY supplement
requires login; the public AIP file supplies the numerical tables.

The source hash, URLs, source-version limitations, and hashes of every
normalized product are recorded in
`peritheos/data/datasets/rubidium-halides-farla-2025-source.json`.
Raw numerical cells from all seven supplement tables are retained alongside
the CSV files. A reader can regenerate every normalized row with
`scripts/recover_rubidium_halides_farla_2025.py` and the official DOCX.

The article says GSAS-II refinements, raw diffraction and model/input data
files are supplied. The public Figshare article lists **one DOCX only**, with
no embedded files; DESY's file inventory has no accessible numerical/project
archive. Original GSAS-II/EoSFit7 input files, ruby wavelengths and covariance
matrices were not recovered. Numerical tables remain usable independently.

## Equation and parameter choices

Accepted manuscript Eqs. (1), (6)-(9) define

\[
P(V,T)=P_{\mathrm{BM3}}(V;V_0,K_0,K'_0)
 +\frac{\gamma(V)}{V}[E_D(V,T)-E_D(V,300\,\mathrm K)],
\quad \gamma(V)=\gamma_0(V/V_0)^q,
\quad \theta(V)=\theta_0\exp[(\gamma_0-\gamma(V))/q].
\]

This is the **integrated Gruneisen Debye-temperature law**, with reference
energy evaluated at the same compressed volume and 300 K. It is not the
variable-exponent `theta0*(V/V0)^(-gamma(V))` law or an isentropic energy
reference. Use `MieGruneisenDebye`, `integrated_gruneisen`,
`thermal_pressure_reference=reference_temperature`. The model already exists
in Python and Rust; no new equation family is introduced.

Each formula unit contains two atoms. Each Pm-3m conventional cell contains
one formula unit (Z=1), so the Debye atom count is **n=2**. Native molar energy
uses the existing angstrom-cubed-to-J/bar/mol scale, `NA*1e-25` for Z=1.
The ordered ideal B2 basis is Rb at 1a (0,0,0) and Cl/Br/I at 1b
(1/2,1/2,1/2), each fully occupied. This is the CsCl prototype explicitly
identified in the paper, not a refinement of occupancies. The new cards'
lattice parameters are cube roots of their hypothetical reference volumes;
they are not ambient B2 measurements. Existing RbCl structure metadata is
preserved.

Preferred Table I case 1 uses CsCl LVP pressures and ruby-calibrated DAC data:

| Quantity | RbCl | RbBr | RbI |
|---|---:|---:|---:|
| Fixed V0, cm3/mol | 36.25(3) | 41.74(4) | 51.11(2) |
| Printed V0, A3/cell | 60.20(2) | 69.31(7) | 84.87(2) |
| Converted canonical V0, A3/cell | 60.1945411851 | 69.3109006638 | 84.8701517233 |
| K0, GPa | 19.89(8) | 16.28(4) | 13.69(4) |
| K0' | 5.00(2) | 5.28(2) | 4.95(1) |
| gamma0 | 1.96(4) | 2.18(14) | 2.21(7) |
| q | 1.05(9) | 1.52(24) | 1.42(10) |
| Fixed theta0, K | 177 | 145 | 118 |

V0 is first constrained by an independent 300 K DAC BM3 fit and then fixed
in the combined fit; K0, K0', gamma0 and q are free in case 1. Theta0, n and
Tr remain fixed. The molar reference volumes are the values used by EoSFit7
and named in benchmark table captions. They are converted with exact SI NA;
the separately rounded cell-volume row is retained as source metadata.
The printed molar/cell **uncertainties are inconsistent after conversion**:
for RbCl 0.03 cm3/mol gives 0.0498 A3, rather than printed 0.02 A3; for RbI
0.02 cm3/mol gives 0.0332 A3, rather than 0.02 A3. Canonical errors are
explicitly propagated from the molar row; both printed forms remain available.
The fixed V0 uncertainties are not set to zero. Parameter errors are described
as EoSFit7 estimated standard deviations; confidence level and covariance are
not supplied. No errors are invented for theta0, n or Tr.

The official supplement derives theta0 from **DFT elastic constants** of
Haddadi et al. (2009), with Anderson's sound-speed procedure. Consequently the
complete parameterization is classified **hybrid**, although the free fitted
coefficients are experimental. Main-text reference 41 instead names Kumara
Swamy et al. (1996); the supplement explicitly identifies the calculation's
DFT source. Its printed sound-speed/temperature formulas contain apparent
negative-exponent/unit typos. An independent SI calculation using the physical
positive square/cube roots and the supplement's printed K, G, density and
molar mass gives theta=177.0221, 144.9719 and 117.9674 K, respectively;
its mean sound speeds agree with the reported values within 0.4 m/s.
These checks reproduce the actual fixed temperatures after rounding and
identify the formula typos without changing any adopted coefficient.

Other source parameterizations are distinct and not silently blended into
case 1. Table I case 2 (Mo-weighted) has tuples (K0, K0', gamma0, q):
RbCl (19.16(23), 5.22(7), 2.03(6), 1.75(13));
RbBr (16.28(5), 5.28(2), 2.25(14), 2.41(23));
RbI (13.61(8), 4.98(3), 2.38(14), 2.42(20)). Case 3 fixes DAC-derived
static parameters: RbCl (19.93,4.99), RbBr (16.32,5.26), RbI (13.68,4.96);
its (gamma0,q) are (1.96(4),1.04(9)), (2.17(13),1.49(23)),
(2.22(7),1.43(10)). An alternative RbBr case 3 fixes V0=41.55(4) cm3/mol,
K0=17.55 and K0'=4.98, giving gamma0=2.11(18), q=1.50(32).
These alternatives and the Table II HT-BM3 model are recorded as audit
context, not additional executable records. Only case 1 has independent
tabulated model benchmarks S5-S7 in this supplement.

## Recovered observations and uncertainties

All **376 printed observation rows** are retained: S1 has 16 RbCl, 17 RbBr
and 16 RbI DAC observations; each compound has 109 LVP rows from S2-S4.
One RbI LVP row has no sample or Mo volume and remains present with blank
fields. There are 375 available sample-volume observations, not 376 measured
sample volumes. Every reported P, V, temperature, density, parenthetical
error and paired CsCl/Pt/Mo volume/pressure is retained. No recalculated
calibrant pressure overwrites a source pressure.

The source's S3/S4 sample-volume headings incorrectly say `VRbCl-B2`.
Compound identity follows each table title and consistent reported density;
raw headings are archived. S1 contains a malformed token `1.19(0` and
several printed pressure errors `(0)`. The missing parenthesis is not filled
in the raw table or token column. The parser expands the printed zero, and
does not claim physically exact pressure. Underlying unrounded errors are
missing. DAC error confidence convention is not specified separately.
LVP volume errors are GSAS-II estimated standard deviations; the **30 K**
temperature error is an author-assigned uncertainty for temperature gradients.
The pressure errors were calculated with EosCross PVT correlation propagation
and exclude EOS-parameter uncertainty. The source correlation matrices are
unavailable.

Three 300 K LVP rows per compound are retained even though the source figures
describe 300 K fitted observations as DAC data. Diagnostic fits report both
including and excluding those cold LVP rows. Source author row flags are
unavailable. The S1 RbI lowest row is labelled B2 at 0.54 GPa, whereas the
Figure S7 caption describes an RbI B1 pattern at 0.6 GPa. The row is retained
without an invented exclusion or phase reassignment; exact near-transition
selection needs original diffraction/project evidence.

## Calibration and phase coverage

The preferred fit uses reported **CsCl** LVP pressures. S2-S4 footnote a
specifies Dewaele (2020) isothermal plus Decker et al. (1973) thermal, but
main-text reference 2 is Decker (1971). Historical public EosCross code now
identifies the composite Vinet plus integrated Debye model and its coefficients.
Its GUI defaults reproduce 104/109 distinct gauge rows within pressure/volume
rounding; saved runtime settings and five discrepancies remain unresolved.
No existing Campbell-Heinz CsCl record is linked as this calibrant. Mo comparison pressures use
Huang et al. (2016); Pt footnotes say Matsui (2007), while the main reference
30 identifies Matsui (2009). These bibliographic discrepancies remain explicit.
Pt is not used as the preferred fitting pressure. Case 2 uses Mo/CsCl averages;
its exact averaging/correlation inputs are not reproduced here.

DAC pressure uses IPPS-Ruby2020 (Shen et al., 2020). Pressure/errors are
available, but raw ruby wavelengths are absent. Full recalibration from raw
observations is therefore unavailable. No shared calibrant cards or registry
entries are changed by this work.

The authors observe Mo and especially Pt apparent pressure deficits relative
to CsCl above about 1100 K (up to roughly 1 and 3 GPa) and prefer CsCl.
They also flag possible CsCl-model problems in BT815 near 17-19 GPa at high T.
The thermocouple is type C and its emf is not pressure-corrected; the paper
discusses possible systematic temperature shifts up to about 70 K at 1800 K.
Those issues are not absorbed into invented parameter errors or corrections.

DAC bounds are RbCl/RbI 0.54-26.05 GPa, RbBr 1.19-26.47 GPa at 300 K.
LVP CsCl pressures span 2.43-20.52 GPa, with temperatures to 1800 K.
The cards' rectangular ranges are **marginal observation bounds**, not a
claim that every P-T combination was measured or that the phase is stable
throughout a rectangle. Hot coverage reaches about 21 GPa, not 26 GPa.
B2 is a high-pressure phase with a hypothetical zero-pressure V0; the
source's approximate “above 0.5 GPa” statement is not a measured full
transition boundary for each compound. No B1/solidus phase boundary is
invented. Verify phase experimentally before pressure-marker use.

## Independent numerical reproduction and fits

Run `uv run python scripts/reproduce_rubidium_halides_farla_2025.py`.
The checked result is `docs/data/rubidium-halides-farla-2025-reproduction.json`.
The reference implementation evaluates BM3 explicitly, integrates the Debye
energy by 48-point Gauss-Legendre quadrature, and inverts by Brent's method.
It does not call Peritheos to calculate reference pressures/volumes.

All **162** model ratios in S5-S7, from 5-100 GPa and 300-2500 K, are
retained. Each table flags the 5 GPa/2500 K sample as molten. These are
equation benchmarks, not additional observations or a validity extension.
Native Peritheos ratios agree with the independent calculation within
1.1e-14. Maximum absolute differences from published ratios are 0.0006994
(RbCl), 0.0014902 (RbBr), and 0.0002401 (RbI). All states through 20 GPa /
1500 K differ by less than 0.00025 in V/V0. Universal last-digit table
reproduction is **not** achieved with the rounded Table I coefficients.

Every benchmark difference lies within a documented rounding-sensitivity
envelope: 0.00005 for the ratio table, plus the summed volume-ratio changes
under independent half-last-digit perturbations of K0, K0', gamma0, q,
theta0 and molar V0. This is a coefficient-precision consistency check,
not a covariance-based confidence interval or evidence of one common hidden
coefficient vector. It is stricter than simply allowing reported parameter
standard deviations, and is evaluated state by state in the report.
The stronger common-vector check in the follow-up succeeds for RbBr and RbI;
RbCl remains unresolved under its stated fixed reference-volume convention.

The initial independent diagnostic fit fixes source V0, theta0, n and Tr,
fits K0/K0'/gamma0/q, uses all DAC rows and nonblank LVP rows above 300 K,
and minimizes pressure residuals with effective variance
`sP^2 + (dP/dV*sV)^2 + (dP/dT*sT)^2`. It assumes diagonal observation
errors and uses a **diagnostic 0.005 GPa floor** for printed-zero DAC
pressure errors. This floor is a transparent sensitivity choice, not a
source uncertainty. It is not written into the primary tables. All-data
and unweighted fits, plus separate DAC BM3 fits, are also reported.

| Preferred diagnostic | RbCl | RbBr | RbI |
|---|---:|---:|---:|
| Rows fitted | 122 | 123 | 121 |
| Source curve P RMSE, GPa | 0.15443 | 0.15949 | 0.12781 |
| Source curve maximum absolute P residual, GPa | 0.40325 | 0.45728 | 0.39524 |
| Fitted K0, GPa | 19.72509 | 16.35776 | 13.59121 |
| Fitted K0' | 5.04813 | 5.24374 | 4.98811 |
| Fitted gamma0 | 1.98428 | 2.16435 | 2.23133 |
| Fitted q | 1.04570 | 1.47180 | 1.41350 |

These source residuals are consistent with the source's description that
most differences are within about 0.3 GPa and that outliers exist. Static
coefficient deviations reach about 3.8 reported standard deviations in
these initial conditional refits; unweighted results differ more. The
follow-up below improves the weighted comparison by freezing the weights
inside each least-squares cycle and removing the pressure-error floor.
Original row masks, unrounded errors, full fitting covariance and staging
remain unavailable. **Equation validation,
rounded-source curve consistency, independent refitting, original fitted
coefficient parity and author-workflow parity are distinct results.** The
last two are not achieved. No diagnostic coefficients enter the catalog.

## Recovery and fit parity follow-up

The reproducible follow-up is
`uv run python -m scripts.investigate_rubidium_halides_farla_2025 --output docs/data/rubidium-halides-farla-2025-follow-up.json`.
It uses the bundled numerical source facts in
`peritheos/data/datasets/rubidium-halides-farla-2025-author-inputs.json`.
Source commits and file hashes are pinned there; public software is not
silently vendored into the material records.

### Recovered CsCl inputs

The author's public EosCross v1.3 predates the publication. Commit
`9518e36112be23504b84c4f595eab0b0207a3539` (23 August 2024) supplies
[the CsCl calibrant](https://gitlab.desy.de/robert.farla/eoscross/-/blob/9518e36112be23504b84c4f595eab0b0207a3539/EosCross/calibrants/Dewaele_2020.py):

| Input | Value |
|---|---:|
| Static equation | Vinet |
| Module V0 | 42.179 cm3/mol |
| K0 | 16.74 GPa |
| K0' | 5.703 |
| theta0 | 151 K |
| gamma0 | 1.99 |
| q | 1.18 |
| n, reference T, reference P | 2, 300 K, 0 GPa |

The thermal term subtracts Debye energy at 300 K, with
`gamma=gamma0*(V/V0)^q` and the integrated theta law. The dependency file
pins BurnMan 1.2.0; its wheel and relevant source files were checksum audited.
The code attributes the thermal inputs to Decker (1971), so the supplement's
1973 attribution discrepancy persists at the bibliographic level.

The [single-calibrant calculation](https://gitlab.desy.de/robert.farla/eoscross/-/blob/9518e36112be23504b84c4f595eab0b0207a3539/EosCross/single/CsCl1_Dewaele2020.py)
overrides module V0 with the GUI lattice parameter. The
[GUI default](https://gitlab.desy.de/robert.farla/eoscross/-/blob/9518e36112be23504b84c4f595eab0b0207a3539/EosCross/MBase_CsCl%2BPt.py)
is a0=4.1221 Angstrom; the code uses NA=6.0221367e23/mol. Together these
give V0=42.1799615564 cm3/mol. This identifies historical defaults, not a
saved author session or proof that every original run used them.

The 109 CsCl observations are repeated identically across S2-S4, so there
are 109 distinct gauge comparisons rather than 327. Using the recovered
GUI defaults, 102 match pressure rounding alone and 104 match pressure plus
volume rounding. Pressure allowance is 0.005 GPa plus the maximum pressure
change under V +/-0.0005 Angstrom3, with the printed T held fixed.
All 109 agree within the reported pressure error plus this allowance.
The RMSE is 0.01104 GPa and the maximum difference is 0.07090 GPa.
Independent diagonal propagation of the printed volume errors and 30 K
temperature error matches all 109 pressure errors at their printed precision.
The historical code itself starts with diagonal V,T errors and derives a
PVT covariance; the original regression's treatment of those correlations
is still unknown.

Five readings exceed the pressure/volume rounding allowance:

| Run | T, K | CsCl V, Angstrom3 | Printed P, GPa | Replayed P, GPa | Printed pressure error, GPa |
|---|---:|---:|---:|---:|---:|
| BT792 | 1100 | 60.468 | 5.54 | 5.53346 | 0.07 |
| BT792 | 1200 | 60.901 | 5.49 | 5.51667 | 0.07 |
| BT792 | 1300 | 61.389 | 5.44 | 5.48106 | 0.07 |
| BT792 | 1400 | 61.773 | 5.44 | 5.51090 | 0.07 |
| BT795 | 300 | 52.947 | 10.03 | 10.10028 | 0.07 |

The report records implied temperature shifts as diagnostics only. Neither
those shifts nor corrections to any source observation are adopted. Using
the module V0 directly reproduces fewer readings, so the GUI override is a
material implementation detail rather than a negligible unit convention.

### Weighting and remaining coefficient differences

[Angel et al. (2014)](https://www.rossangel.com/Download/2014_Angel_etal_EosFit7.pdf)
describe updating effective-variance weights on each least-squares cycle.
The initial diagnostic instead minimized residual/sigma(parameters) directly
and introduced a 0.005 GPa floor for printed-zero DAC pressure errors.
Both choices are now tested explicitly. The follow-up holds the weights
fixed during each inner least-squares solve and recomputes them between
cycles. The denominator remains
`sP^2+(dP/dV*sV)^2+(dP/dT*sT)^2` with diagonal errors. A printed-zero sP
does not create infinite weight because the measured volume error supplies
finite effective variance. Source table tokens remain unchanged.

The comparison below uses literal printed errors, all DAC rows and nonblank
LVP rows above 300 K, fixed source molar V0, and the more precise fixed theta
values in the corrected author code (177.04, 144.96, 117.98 K). It also tests
R=8.314 J/mol/K from the locally audited 2021 CrysFML EthDebye source; this is
a convention test, not a recovered executable version used by the authors.

| Compound | K0, GPa | K0' | gamma0 | q | Maximum absolute difference / reported parameter esd |
|---|---:|---:|---:|---:|---:|
| RbCl | 19.92194 | 4.99518 | 1.96402 | 1.05922 | 0.40 |
| RbBr | 16.31594 | 5.25587 | 2.16660 | 1.47419 | 1.21 |
| RbI | 13.75511 | 4.93389 | 2.19264 | 1.41792 | 1.63 |

These are marginal parameter comparisons, not a joint significance test:
the fitted parameter covariance is not recovered. None of the three fits
rounds to all four published coefficients. Using the original rounded theta
and current gas constant barely changes this result. Restoring the 0.005 GPa
floor increases the maximum marginal discrepancies to 1.88, 2.30 and 3.03
reported esds. Including cold LVP rows or substituting reconstructed CsCl
pressures for the LVP targets also fails to restore exact parity. These
sensitivity cases and the original explicit unweighted pressure-residual
baseline remain available in the two reports.

### Common coefficients for the rounded benchmark tables

A stronger check than the initial state-specific envelope seeks one common
K0/K0'/gamma0/q vector within +/-0.005 of each printed coefficient, holding
the printed molar V0 and precise author theta fixed. RbBr and RbI each admit
such a vector for all 54 ratios within +/-0.00005. This establishes compatible
calculated-table reconstructions, not unique coefficients or author-fit
reproduction; calculated ratios never enter the observational regressions.

For RbCl the initial numerical search leaves two states beyond ratio rounding
(5 GPa/1500 K and 100 GPa/1500 K, maximum difference 0.00007576).
The optimizer did not establish infeasibility, so this remains an unresolved
diagnostic under the stated normalization and fixed-input assumptions.
The 16 October audit below supersedes this localization: excluding a single
5 GPa/1500 K entry also resolves the apparent 100 GPa/1500 K mismatch.

### Search for original files

The public search inspected the AIP Figshare article and collection versions,
related DataCite records, all available EosCross branch/tag histories, the
default-branch file inventories of 13 public author GitLab projects, public
EosCross releases/packages/wiki, and relevant Figshare and Zenodo searches.
The Figshare collection exposes only the same version-1 DOCX. No original
EoSFit7 or GSAS project archive was identified within that scope. URLs and
inventory scope are recorded in the author-input facts file. No authors were
contacted; public search results cannot prove that no other deposit exists.

The [first public rubidium evaluator](https://gitlab.desy.de/robert.farla/eoscross/-/blob/7b34d091878ae3952c3222d51256fa367ff0e91d/EosCross/calibrants/Farla_2025.py)
contains more decimal places but different coefficients, then was
[corrected six minutes later](https://gitlab.desy.de/robert.farla/eoscross/-/blob/ed3b7109042f43211812cd64cd1dc1c176ca08b6/EosCross/calibrants/Farla_2025.py)
to the published rounded coefficients. The earlier values are preserved as
historical evidence, not promoted to unrounded final author fits. The corrected
code recovers precise theta values, but not additional fitted-coefficient digits.

## Discrepancy audit on 16 October 2026

Literal parity with the five printed CsCl pressure readings and the complete
RbCl benchmark table is still not achieved. The calculations now distinguish
an implementation difference from a discrepancy in the published inputs or
calculated table. No primary observation, calculated-table cell or published
EOS coefficient has been changed.

The reproducible diagnostic is
`uv run python -m scripts.audit_rubidium_halides_farla_2025_discrepancies --output docs/data/rubidium-halides-farla-2025-discrepancies.json`.
Its independent historical execution evidence is
`docs/data/rubidium-halides-farla-2025-burnman-replay.json`.

### Historical calibration execution and printed inputs

The recovered CsCl module was executed with its pinned BurnMan 1.2.0 in an
isolated Python environment. Its 109 pressures differ from the independent
calculation by at most **6.27e-11 GPa**. The historical module is also byte
identical to the corresponding module at public main commit `f5b2249b`,
whose remote identity was rechecked on 16 October. A change between these
source versions does not explain the outliers. Replaying the original
module on the printed volumes and temperatures preserves all five discrepancies.

To repeat that external-dependency check, obtain the pinned module from
the source URL above, then run
`uv run --no-project --python 3.11 --with burnman==1.2.0 python scripts/replay_rubidium_halides_farla_2025_burnman.py --source /path/to/Dewaele_2020.py --output docs/data/rubidium-halides-farla-2025-burnman-replay.json`.
The script verifies the module's SHA256 before execution. It does not install
BurnMan into the project environment or require it for the normal tests.

The five original entries were visually rechecked in the supplement's S2
rendered pages 17-18 and agree with the extracted cells. The same inputs
are repeated across S2-S4. The historical batch calculation cubes the
lattice parameter, whereas the supplement prints volume; however, cubing an
unrounded lattice parameter consistent with the printed volume cannot exceed
the already allowed +/-0.0005 Angstrom3 volume rounding interval.

The following inverse inputs would reconcile the printed pressures with the
recovered GUI defaults. Temperature intervals allow both pressure and volume
rounding; volume shifts instead hold the printed temperature and pressure
centers fixed. They are diagnostic requirements, not recovered measurements.

| Run | Printed T, K | Required temperature shift interval, K | Required volume shift at printed T, Angstrom3 |
|---|---:|---:|---:|
| BT792 | 1100 | +0.55 to +5.17 | -0.01128 |
| BT792 | 1200 | -13.94 to -9.33 | +0.04808 |
| BT792 | 1300 | -20.18 to -15.59 | +0.07771 |
| BT792 | 1400 | -33.13 to -28.56 | +0.13969 |
| BT795 | 300 | -35.17 to -29.96 | +0.05941 |

Integer-temperature rounding of +/-0.5 K does not resolve any of the five.
Figure S1's caption, rendered page 4, lines 71-72, states that the type-C
thermocouple controlled each target within 1 K. Permitting that +/-1 K
control tolerance makes the 1100 K reading compatible with pressure/volume
rounding, but the other four remain outside it. The separate 30 K temperature
uncertainty includes broader measurement limitations; it is not evidence
that these particular shifts were used in the pressure calculations.

Changing a common reference lattice parameter cannot repair all readings
with the recovered coefficients and printed temperatures: their required
a0 intervals have an empty intersection, with lower bound 4.12215169 and
upper bound 4.11897595 Angstrom. The limiting readings are BT792 at 1100
and 1400 K. The corresponding test of a single common temperature offset
also has an empty intersection. Thus neither a simple shared GUI setting
nor a uniform temperature offset supplies exact source-pressure parity.
Original lattice parameters, original temperatures and saved per-run settings
remain necessary to distinguish an input/export inconsistency from a
calculation performed with other settings. No particular cause is established.

### RbCl benchmark localization

The new search minimizes the largest absolute ratio discrepancy, using
pressure constraints at the exact ends of each +/-0.00005 ratio interval.
It converges numerically with a maximum discrepancy of 0.00006436 under
the four coefficient half-digit bounds, fixed molar V0=36.25 cm3/mol and
theta0=177.04 K. Testing theta0=177 K, R=8.314, and half-digit variations
of V0/theta0 with either a shared or printed normalization volume does not
restore all 54 printed intervals. These local numerical optima are not
certified global infeasibility proofs.

All 54 single-entry omissions were then tested under the original four
coefficient bounds. **Only omitting 5 GPa/1500 K produces a common vector
consistent with the other 53 printed ratios.** One such diagnostic vector
is (K0, K0', gamma0, q) = (19.88765988 GPa, 5.00349408, 1.96362865,
1.04699136), with fixed V0/theta0 as above. All four coefficients round
to Table I and all 53 retained ratios differ by at most 0.00004809.

That vector predicts **0.92828693**, rounding to **0.9283**, for the omitted
entry, while S5 prints **0.9282**. Numerical minimization/maximization of
this omitted ratio subject to the other 53 rounding intervals yields
0.92827583-0.92829062; both ends round to 0.9283. Those are numerical
extrema under the stated assumptions, not certified interval bounds.
Changing just that diagnostic target to 0.9283 makes all 54 entries
rounding-compatible with one coefficient vector.

This suggests a one-digit inconsistency in the benchmark entry, but does
not establish an author correction or recover the unrounded fitted model.
The 0.9283 alternative comes from the joint table reconstruction: direct
evaluation of the rounded published coefficients instead gives 0.9280678,
rounding to 0.9281. Those two calculations answer different questions.
The original 0.9282 remains in the bundled benchmark CSV, and none of these
reconstructed coefficients enters the catalog or an observational regression.

Exact original-fit parity therefore remains false. The stronger outcome is
an executed-software cross-check for CsCl and localization of the RbCl
calculated-table inconsistency to one entry. No authors were contacted.

## Direct EosFit console refits

The official EosFit7c **7.60**, version date **13 May 2021**, was run on the
printed observation tables on 16 October 2026. This is an execution of the
official fitting program, rather than another implementation of its driver.
The [official installation package](https://www.rossangel.com/text_eosfit.htm)
was downloaded and extracted in the task workspace; no application or system
package was installed. The scientific executable was unmodified, SHA256
`6b48e46e5bb01c2dd13060bc29b4a9f7c67e70c7161866b77f8de3f238445e6d`.
The installer SHA256 matches the prior official-package audit. Its legacy
display dependency ran on a temporary authenticated virtual display, which
was stopped after the fits. Runtime provenance is retained in the manifest.

All **18 cases converged in both refinement passes**. Each fits K0, K0',
gamma0 and q with BM3 plus the full integrated MGD model, n=2 and Tr=300 K.
V0 and theta0 are fixed; the q-compromise is explicitly disabled. Observed
Z=1 cell volumes and their errors are converted to cm3/mol without rounding
the conversion. The six cases per material are:

- Unit weights and printed theta0, providing the explicit unweighted baseline.
- Reported P,T,V errors and printed theta0.
- Reported errors and precise theta0 from the corrected author evaluator.
- The preceding case with all nonblank 300 K LVP rows included.
- The precise-theta case with the diagnostic 0.005 GPa pressure-error floor.
- The precise-theta case converting observations with the legacy Avogadro constant.

The principal selections contain 122 RbCl, 123 RbBr and 121 RbI observations:
all DAC rows and nonblank LVP rows above 300 K. Including the cold LVP rows
gives 125, 126 and 124. This tests an explicit selection; it does not recover
the authors' original use/omit flags.

### Native weighted results

Using reported errors and precise fixed theta0, EosFit gives the following
parameters and its own rescaled standard errors:

| Material | K0, GPa | K0' | gamma0 | q |
|---|---:|---:|---:|---:|
| RbCl | 19.92227 +/-0.08601 | 4.99511 +/-0.02382 | 1.96401 +/-0.04225 | 1.05929 +/-0.08973 |
| RbBr | 16.31590 +/-0.04081 | 5.25587 +/-0.01514 | 2.16658 +/-0.13955 | 1.47412 +/-0.24376 |
| RbI | 13.75502 +/-0.03245 | 4.93392 +/-0.01107 | 2.19266 +/-0.06105 | 1.41792 +/-0.09112 |

The largest differences from the authors' coefficients, divided by their
reported marginal standard errors, are **0.40, 1.21 and 1.63**, respectively.
The new native covariance matrices describe these diagnostic refits; they
are not the missing covariance of the published coefficients.

With matched fixed anchors, R=8.314 and the same rows/errors, the largest
native-versus-independent-IRLS differences are **0.000326 GPa in K0** and
**0.000074 in the dimensionless parameters**. Thus the original program
confirms the independent fit to substantially finer precision than the
publication's coefficient digits. It does not remove the differences from
the published coefficients. Native weighted chi-squared per degree of
freedom is 2.4678 for RbCl, 17.5511 for RbBr and 2.1750 for RbI; EosFit's
reported parameter errors and saved covariance include its own scaling.

### Native unweighted results and sensitivities

The unweighted pressure-residual baseline gives:

| Material | K0, GPa | K0' | gamma0 | q |
|---|---:|---:|---:|---:|
| RbCl | 19.32076 | 5.16197 | 2.01802 | 1.03984 |
| RbBr | 15.82327 | 5.40644 | 2.21527 | 1.40171 |
| RbI | 13.44061 | 5.04076 | 2.27423 | 1.44071 |

Native final weights are exactly one for these cases. For the weighted cases
all effective weights are finite and positive, despite the program's warning
that some reported P/T esds are zero. These printed-zero errors are retained
without silently replacing them. Infinite entries in native deltaP/esdP
diagnostic columns refer to division by the printed-zero pressure error;
they are not infinite effective regression weights.

The 0.005 GPa floor substantially changes the native static coefficients;
its RbCl K0 is 19.76405 GPa. Including cold LVP rows or using the legacy
Avogadro conversion produces smaller changes. None of the 18 cases recovers
all four published coefficient digits. The source rows, unrounded errors,
fixed anchors and saved author settings therefore remain the relevant
reproducibility questions. This does not establish which particular source
input or setting caused the difference.

### Retained native evidence and rerunning

`docs/data/rubidium-halides-farla-2025-eosfit7c/manifest.json` contains every
case, source hashes, parameter differences, native errors and covariance,
warnings and runtime identity. Each case directory retains `input.dat`, a
source-row map, `run.mcr`, raw `stdout.txt`/`stderr.txt`, native `run.log`,
`initial.eos` and `fitted.eos`. The original bytes are protected by a scoped
Git attributes rule; the generated macros contain their original absolute
output paths.

Run `uv run python -m scripts.run_rubidium_halides_farla_2025_eosfit --executable /path/to/eosfit7c --output /path/to/new-results`
with a runnable external console and its required display/runtime environment.
The output directory must be new, so earlier evidence is not overwritten.
The script generates all inputs and macros with paths for that run and
performs two refinement passes per case. It does not install software.
For clones, compare archived macros using the absolute folder recorded in
each macro, rather than substituting the current checkout path.

The transient native thermal-model unit warning occurs before the macro
explicitly sets GPa and cm3/mol; the initial and final saved models contain
the correct units. All warnings remain in the preserved logs. No native-fit
coefficients or covariance have been substituted into the published records,
and no authors were contacted.

## Constrained Case 3 EosFit reproduction

The same unmodified official EosFit7c 7.60 was used for **24 additional
converged fits**, covering four source variants and six input/weighting
choices. These are separate from Case 1. Table I and its footnotes fix V0,
theta0, K0 and K0' in Case 3 and fit only gamma0 and q. RbBr has two explicitly
reported alternatives, with different V0 and K0 as well as K0'; they cannot
be represented merely by changing K0' in one common static model.

The transcribed inputs are
`peritheos/data/datasets/rubidium-halides-farla-2025-case3-inputs.json`:

| Source variant | Fixed V0, cm3/mol | Fixed K0, GPa | Fixed K0' | Fixed theta0, K |
|---|---:|---:|---:|---:|
| RbCl Case 3 | 36.25 | 19.93 | 4.99 | 177 |
| RbBr Case 3, DAC K0' | 41.74 | 16.32 | 5.26 | 145 |
| RbBr Case 3, selected K0'=4.98 | 41.55 | 17.55 | 4.98 | 145 |
| RbI Case 3 | 51.11 | 13.68 | 4.96 | 118 |

Figure 1's caption independently identifies the RbBr alternative
V0=41.55 cm3/mol and K0=17.55 GPa for selected K0'=4.98. The quoted errors
of the fitted DAC anchors are not propagated into the two-parameter native
covariance; these runs condition on the fixed values as the source describes.

### Fits with the literal printed anchors

The primary comparison uses reported P,V,T errors, no pressure-error floor,
all DAC rows and nonblank LVP rows above 300 K. Each fit was repeated after
updating its parameters; both passes converged. The results are:

| Source variant | Published gamma0 | EosFit gamma0 | Published q | EosFit q | Largest difference / published esd |
|---|---:|---:|---:|---:|---:|
| RbCl | 1.96(4) | 1.96029 +/-0.04110 | 1.04(9) | 1.04704 +/-0.08855 | 0.078 |
| RbBr, DAC K0' | 2.17(13) | 2.17459 +/-0.13941 | 1.49(23) | 1.50152 +/-0.24233 | 0.050 |
| RbBr, selected K0'=4.98 | 2.11(18) | 2.11447 +/-0.17655 | 1.50(32) | 1.50637 +/-0.32106 | 0.025 |
| RbI | 2.22(7) | 2.21407 +/-0.06127 | 1.43(10) | 1.42931 +/-0.09030 | 0.085 |

These differences are substantially smaller relative to the reported errors
than in Case 1. They remain marginal comparisons, not joint significance
tests or a recovery of the authors' covariance. With the printed anchors
held exact, each variant still has one thermal coefficient that rounds
differently from the source.

The other retained cases are an explicit unweighted baseline, precise
theta0 values from the corrected author evaluator, inclusion of all nonblank
cold LVP rows, and exclusion of the DAC rows. Adding or removing 300 K rows
does not change gamma0/q at native printed precision, as expected because
the reference-subtracted thermal pressure and its thermal-parameter
derivatives vanish at Tr. Those rows can still change chi-squared/error
scaling. The tests therefore keep row selection, coefficients and conditional
uncertainties distinct.

### Consistency with unprinted DAC digits

Four native examples test the precision of the fitted DAC anchors. Every
example preserves all printed fixed parameters on rounding; only the entry
shown below changes. The deliberately selected RbBr K0'=4.98 remains
**exactly 4.98**, so its example varies the fitted V0 instead.

| Source variant | Illustrative fixed-anchor value | Native gamma0 | Native q | Rounded thermal pair |
|---|---|---:|---:|---|
| RbCl | K0'=4.989 instead of printed 4.99 | 1.95930 | 1.04291 | 1.96, 1.04 |
| RbBr, DAC K0' | K0'=5.2575 instead of printed 5.26 | 2.17065 | 1.48899 | 2.17, 1.49 |
| RbBr, selected K0'=4.98 | V0=41.549 instead of printed 41.55 cm3/mol | 2.11426 | 1.50362 | 2.11, 1.50 |
| RbI | K0'=4.9606 instead of printed 4.96 | 2.21559 | 1.43327 | 2.22, 1.43 |

All four thermal pairs now round to their published Case 3 values. These
examples establish that the remaining differences **can** arise from the
unprinted precision of the fixed inputs. The anchors were chosen as
diagnostic examples after examining the discrepancies; they are not estimates
of the authors' actual unrounded values, evidence of uniqueness, or recovered
author files. No original observation, published coefficient, deliberately
selected constraint, or canonical EOS record is replaced by these examples.
This Case 3 consistency result does not establish Case 1 fit parity.

### Retained Case 3 evidence

`docs/data/rubidium-halides-farla-2025-case3-eosfit7c/manifest.json` stores
the 24 cases, native parameters/errors/covariance, warnings, file hashes,
source row maps and executable identity. Each directory preserves the
generated input and macro, native initial/final model, console output and
log. The full saved covariance verifies that all four anchor parameters
remain fixed and only gamma0/q are refined.

To rerun, use
`uv run python -m scripts.run_rubidium_halides_farla_2025_case3_eosfit --executable /path/to/eosfit7c --output /path/to/new-case3-results`
with the external executable's runtime/display environment. The output must
be new, and all earlier Case 1 evidence remains unchanged. The temporary
display and authorization used here were removed after execution. No
authors were contacted and original author-fit recovery remains false.

## Verification and integration

Material-specific tests cover raw-to-normalized transcription, hashes,
observations and blanks, JSON Schema and public loading, structure and cell
contents, preserved RbCl default, metadata round trips, reference state,
integrated theta law, source benchmarks, independent fits, pressure/volume
round trips, broadcasting and validity. A material-specific Rust integration
test exercises loading, source ratios, inversion and metadata serialization.
The combined original and follow-up Python suites pass all 95 tests; the
dedicated Rust test passes. Original RbCl fields, records and datasets remain
equal to their pre-addition values. Lint and formatting pass for the new files.
Full command results and shared-ledger work are recorded in
`docs/research/rubidium-halides-farla-2025-integration.md`.
