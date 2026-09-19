# Fischer et al. (2014): experimental Fe–FeSi equations of state

Audit date: 2026-09-19. **Seven published records and two independent B2 refits accepted on five
materials. Both source-published B2 fits remain withheld because their
equation/normalization discrepancy is unresolved.**
All eight computational EOS and all three source-owned hcp+B2 aggregate fits
are excluded from this experimental homogeneous-crystal addition.

## Sources and version control

Rebecca A. Fischer, Andrew J. Campbell, Razvan Caracas, Daniel M. Reaman,
Dion L. Heinz, Przemyslaw Dera, and Vitali B. Prakapenka (2014), “Equations of
state in the Fe–FeSi system at high pressures and temperatures,” *Journal of
Geophysical Research: Solid Earth* **119**(4), 2810–2827,
[doi:10.1002/2013JB010898](https://doi.org/10.1002/2013JB010898).
The publication of record's HTML, Tables 1–4, and methods were read; the
[author's accepted manuscript](https://geosci.uchicago.edu/~campbell/Papers/FischerJGR2014preprint.pdf)
was used as an additional check. The final online date is 15 April 2014;
15 March is the accepted-manuscript date. Focused correction searches and
Crossref relation/update metadata disclosed no correction.

The user supplied all seven official supplementary workbooks after direct
Wiley retrieval failed. Their filenames, SHA-256 digests, and exact parameter
cells are preserved in the [source audit register](../data/fischer-2014-source-audit.json).
Text S1 was also recovered from Downloads, extracted, and its Vinet equation
visually checked. These sources supersede the previous supplement-access hold.
Original workbooks and paper PDFs are not redistributed. No numerical EOS
coefficients were taken from external software catalogs.

The source-owned Fe–16Si reanalysis requires the earlier observations:
[Fischer et al. (2012)](https://doi.org/10.1016/j.epsl.2012.09.022), *EPSL*
**357–358**, 268–276. Its official
[mmc1.xls](https://ars.els-cdn.com/content/image/1-s2.0-S0012821X12005183-mmc1.xls)
was recovered and its Table S2 is bundled. The author manuscript establishes
the inherited staged fitting protocol and sample composition. The
[2013 companion phase study](https://doi.org/10.1016/j.epsl.2013.04.035)
establishes the measured Fe–9Si composition and ordering interpretation.

## Complete source-owned candidate disposition

| Candidate | Source | Disposition |
|---|---|---|
| B20 FeSi BM3 and Vinet | Table 2 / S6 | Accept `fesi_b20_fischer_2014_bm3_1` and `_vinet_2` |
| D03 Fe–9Si BM3 and Vinet | Table 1 / S4 | Accept `fe084si016_d03_fischer_2014_bm3_1` and `_vinet_2` |
| hcp Fe–9Si BM3+MGD and Vinet+MGD | Table 1 / S4 | Accept `fe084si016_hcp_fischer_2014_bm3_1` and `_vinet_2` |
| D03 Fe–16Si Vinet+MGD | Table S5 | Accept `fe073si027_d03_fischer_2014_vinet_1` |
| B2 FeSi BM3+MGD and Vinet+MGD | Table 2 / S6 | Withhold both: physical Debye normalization does not reproduce observations |
| hcp+B2 Fe–9Si BM3+MGD and Vinet+MGD | Table 1 / S4 | Exclude homogeneous-crystal representation |
| hcp+B2 Fe–16Si Vinet+MGD | Table S5 | Exclude homogeneous-crystal representation |
| hcp Fe | Table 3 / S3, column B | Exclude computational EOS |
| hcp Fe11Si | Table 3 / S3, column C | Exclude computational EOS |
| Ordered hcp Fe5Si | Table 3 / S3, column D | Exclude computational EOS |
| Disordered hcp Fe5Si, two configurations | Table 3 / S3, columns E/F | Exclude both computational EOS |
| Ordered hcp Fe3Si | Table 3 / S3, column G | Exclude computational EOS |
| Disordered hcp Fe3Si | Table 3 / S3, column H | Exclude computational EOS |
| hcp FeSi | Table 3 / S3, column I | Exclude computational EOS |

The temperature-corrected computational Fe5Si application is also excluded.
Comparison rows attributed to other authors in Tables 1/2 are not Fischer-owned
fits. Tables 4/S7 are derived core-composition applications, not additional EOS.
The source supplies no fcc+B2 or B20+B2 mixture EOS. The hcp+B2 fits average
phase volumes and enlarge their uncertainties without known modal fractions
or Si partitioning; no single crystal structure or composition can faithfully
represent those aggregates.

BM3 remains the preferred default for the three paired materials. The Vinet
alternatives are source-published fits, even though the paper prefers BM3 for
more consistent core extrapolation. The Fe–16Si Vinet fit is a **published
2014 reanalysis**, not a Peritheos-authored `record_kind: refit`. The B2-only Peritheos refits are separate records derived from observations;
they do not replace or claim to correct the published B2 coefficients.

## Equation, parameters, and volume identity

All V0 values below are cm³/mol **atoms**. Table S4–S6's abbreviated cm³/mol
labels inherit that convention; the tabulated lattice parameters and the main
text independently establish it. Equation (1) is ordinary BM3. Text S1 defines
`P=3*K0*(1-y)/y²*exp[1.5*(K0_prime-1)*(1-y)]`, `y=(V/V0)^(1/3)`.
Thermal branches add `gamma(V)/V*[E_D(V,T)-E_D(V,300)]`, with
`gamma=gamma0*(V/V0)^q` and `theta=theta0*exp[(gamma0-gamma)/q]`.
This maps to `MieGruneisenDebye`, `integrated_gruneisen`, `Tr=300 K`, with no
explicit electronic or anharmonic terms.

| Phase / equation | V0 | K0 (GPa) | K0′ | gamma0 | q |
|---|---|---|---|---|---|
| B20 BM3 | 6.803 ± .008 | 192.2 ± 1.6 | 5.03 ± .17 | — | — |
| B20 Vinet | 6.803 ± .008 | 193.1 ± 1.6 | 5.03 ± .17 | — | — |
| Fe–9Si D03 BM3 | 6.961 ± .012 | 183 ± 4 | 5.59 ± .51 | — | — |
| Fe–9Si D03 Vinet | 6.961 ± .012 | 182.0 ± 3.9 | 5.59 ± .51 | — | — |
| Fe–9Si hcp BM3 | 7.203 ± .054 | 129.1 ± 1.4 | 5.29 ± .08 | 1.14 ± .14 | 1 fixed |
| Fe–9Si hcp Vinet | 7.283 ± .025 | 111.4 ± 1.2 | 6.08 ± .08 | 1.15 ± .15 | 1 fixed |
| Fe–16Si D03 Vinet | 6.799 ± .004 | 193.9 ± 4.8 | 4.91 ± .59 | 1.90 ± .07 | 1 fixed |
| B2 BM3 (withheld) | 6.414 fixed | 230.6 ± 1.8 | 4.17 fixed | 1.30 ± .04 | 1.7 ± .2 |
| B2 Vinet (withheld) | 6.414 fixed | 236.1 ± 1.8 | 4.17 fixed | 1.17 ± .04 | 1.41 ± .21 |

All thermal theta0 values are fixed at 417 K, adopted from iron. B20 and D03
V0 values are measured and fixed even though measured errors are reported.
The hcp V0/K0/K0′/gamma0 are fitted together. Fe–16Si inherits the 2012 staged
protocol: K0/K0′ from 32 room-temperature rows, followed by gamma0 from 66
heated rows with the published static coefficients fixed. Errors on the static
coefficients therefore remain stored first-stage errors, not thermal-stage
free-parameter errors. Error confidence and covariance are unspecified and
stored as null. No missing uncertainty is replaced by zero.

Section 5.1's claim that all q values equal one conflicts with Table 2 and the
B2 fit methods; Table 2 governs. Text S1's q=1.73 comparison likely retains
more digits than Table 2's 1.7; no unreported precision is substituted. Text
S1 cites Ono “2002” in one sentence but lists the 2007 paper in its references;
Table S6 and the main text agree on 2007. The B2 adopted V0=6.414 also differs
from the nearby Sata comparator value 6.435. It is preserved without silently
repairing that adoption lineage.

## Crystallographic models

For A atoms/cell, `Vcell[Å³]=Vm[cm³/mol atoms]*A/0.602214076`.
B20 FeSi has Z=4 and A=8; the fractional alloy formula has one atom, so D03
has Z=A=16 and hcp has Z=A=2. Thermal alloy records use n=1 atom per fractional
formula. Structure cells are separate measured references, not forced to EOS V0.

- **B20:** P213 (No. 198). The independent refinement by
  [Wartchow, Gerighausen and Binnewies (1997)](https://doi.org/10.1524/ncrs.1997.212.1.320)
  gives a=4.495(2) Å at 298 K, Fe 4a x=.13650(2), Si 4a x=.84262(5),
  all coordinates (x,x,x), full occupancy. The
  [primary scanned page](https://d-nb.info/1215420978/34) was visually checked.
- **Fe–9Si D03:** the companion study measures 8.75 ± .40 wt% Si, represented
  as Fe0.84Si0.16. Fm-3m (No. 225), a=5.697387235035957 Å from S1 row 4.
  The explicit maximal-order approximation places Si occupancy .64 on 4a,
  the remainder Fe, with Fe on 4b/8c. It enforces composition but is **not a
  published occupancy refinement**. Superlattice intensities are approximate.
- **Fe–9Si hcp:** P63/mmc (No. 194), Fe .84/Si .16 on 2c. Reference axes
  use S1 row 19: a=2.3845569007078415 Å, c/a=1.6169972685976075 at
  59.3958006581 GPa/1186.915 K. This compressed cell is a diffraction reference;
  the volumetric EOS does not predict the pressure dependence of c/a.
- **Fe–16Si D03:** the measured 15.9 wt% Si sample is Fe0.73Si0.27, not
  ideal Fe3Si. The ambient cell derives from 2012 S2 row 4's atomic molar
  volume. The explicit maximal-order approximation fills 4a with Si and puts
  .08 Si/.92 Fe on 4b, with Fe on 8c. This preserves composition and D03
  superlattice peaks; order-dependent intensities are approximate, not refined.

Atomic multiplicities/occupancies reproduce every stated cell composition.
The Fe–9Si experimental alloy is distinct from theoretical Fe5Si. No aggregate
or computational crystal is installed under an experimental alias.

## Primary observations and calibration

The three CSV resources retain all 576 observation rows:

| Source | Total rows | Source phase groups |
|---|---:|---|
| 2014 S1 | 211 | 15 D03, 76 hcp, 96 hcp+B2, 13 fcc+B2, 11 fcc+hcp+B2 |
| 2014 S2 | 180 | 27 B20, 114 B2, 39 B20+B2 |
| 2012 S2 | 185 | 98 D03, 78 with B2 and hcp volumes, 9 with only B2 volume tabulated |

Full workbook numeric precision, source row numbers, axes, axial ratios,
pressure/temperature errors and KBr readings are retained. Cubic cell volumes
are a³ and hcp volumes are sqrt(3)/2*a³*(c/a). Covariance of a and c/a is not
published, so no fabricated independent-error volume sigma is introduced.
Asterisks are encoded in `source_flags`, with the associated numeric field
blank: * means insufficient peaks for an uncertainty, ** denotes phase detected
without a reliable lattice parameter. S2 row 62 places ** in both phase columns;
this source anomaly is preserved without guessing missing values. Literal zero
reference uncertainties remain zero. Derived 2012 phase labels describe volume
availability; `B2_volume_only` is not a single-phase claim.

Selected marginal coverage is 0–36.2793 GPa at 300 K for B20;
0–28.7885 GPa at 300 K for Fe–9Si D03; 46.4699–197.0644 GPa and
1128.715–2566.255 K for single-phase hcp; and 0–63.0356 GPa and
300–2395.535 K for Fe–16Si D03. B2 observations extend to 146.9149 GPa
and 3146.8 K. Overall main-text rounded coverage is not a rectangular
single-phase stability bound. Slow kinetics permit metastable room-temperature
compression; hcp strain effects near 200 GPa below 2000 K are discussed in
Section 4.1. Coexistence observations never enter the accepted single-phase fits.

The principal calibrant is the **Fischer (2012) B2 KBr thermal EOS**, not the bundled
Dewaele KBr Vinet EOS. Exact reference calibration is not installed, so records
are `partially_resolved` / `reference_eos_not_bundled`. All compressed rows
retain KBr lattice readings. The 2012 source additionally tabulates upstream,
downstream, mean surface, sample and marker temperatures with errors. The 2014
source reports only corrected sample temperatures. Those already include the
approximately 3% downward axial-gradient correction and must not be corrected
again; reconstruction of colder KBr temperature is a separate procedure.
Room-temperature experiments additionally use the Mao (1986) ruby scale;
the 2012 study also uses Dewaele (2008) neon as a secondary standard. Ruby
wavelengths and neon lattice readings are not tabulated, and per-row marker
attribution is not supplied. The available ruby calibration is linked without
claiming it replaces KBr. KBr lattice errors are standard errors across indexed
d-spacings; pressure errors propagate marker lattice and temperature errors.
The 2012 source also calls sample lattice uncertainties standard errors. The
2014 temperature errors combine a 100 K analytical term, differences between
the two sample sides and thickness correction.
No observation-level pressure recalibration is claimed. Unspecified fitting
weights and correlated pressure/temperature errors prevent a claim of exact
recovery of the authors' regression covariance.

## Independent reproduction and refits

Run `uv run python scripts/reproduce_fischer_2014_fesi.py`. Independent NumPy/
SciPy finite-strain/Vinet calculations use Debye energy per physical atom with
Boltzmann's constant; the production adapter uses molar formula-unit energy.
All accepted source curves agree with the native evaluator to better than
3e-10 GPa across all selected source rows, and invert to their input volumes.
The diagnostic fits minimize unweighted pressure residuals because no numerical
weights are published. The native refit ledger independently returns parity
for all seven records; this is the repository's combined-error/numerical
criterion, not proof of identical regression weights.

| Phase / model | Rows | Published pressure RMSE (GPa) | Diagnostic refit RMSE (GPa) |
|---|---:|---:|---:|
| fesi_b20_bm3 | 27 | 0.384246 | 0.369848 |
| fesi_b20_vinet | 27 | 0.392922 | 0.372248 |
| fe084si016_d03_bm3 | 15 | 0.371001 | 0.360421 |
| fe084si016_d03_vinet | 15 | 0.385134 | 0.360041 |
| fe084si016_hcp_bm3 | 76 | 1.364450 | 1.353861 |
| fe084si016_hcp_vinet | 76 | 1.383903 | 1.355083 |

The Fe–16Si Vinet curve gives 1.762368 GPa RMSE across all 98 D03 rows.
Its independent static diagnostic gives K0=198.3597 GPa and K0′=4.3711,
within the respective printed errors of 193.9(4.8) and 4.91(.59). With the
published static coefficients fixed, the heated-data fit gives gamma0=1.8612195,
within the printed 1.90(.07), with 2.072356 GPa RMSE. This staged fit is not
replaced by an unconstrained simultaneous three-parameter regression.

Source checkpoints away from V0 include B20 S2 row 17 at 20.5639181 GPa
(predicted BM3 21.0227007 GPa), D03 S1 row 11 at 17.7805214 GPa
(predicted 17.7155824 GPa), and hcp S1 row 57 at 124.9443953 GPa
(predicted 126.4652491 GPa). These residuals sit within the source's stated
residual bands. Source observed-minus-calculated extrema are approximately
−.9/+.7 GPa (B20), −1.1/+.5 (D03) and −2.9/+2.7 (hcp). The independent
calculation uses the opposite residual sign; small extrema differences reflect
printed coefficient precision. The tests preserve that sign distinction.

## B2 FeSi: demonstrated mismatch and an unresolved normalization hypothesis

Both B2 coefficient sets fail with physical n=2 atoms per FeSi formula and
volumes converted consistently from molar atomic volume. Across all 114 rows:

| B2 source curve | Physical normalization RMSE | Physical residual range, calculated minus observed | Doubled thermal-amplitude RMSE |
|---|---:|---|---:|
| BM3 + MGD | 7.147574 GPa | −16.529816 to −1.060646 GPa | 2.409232 GPa |
| Vinet + MGD | 6.919428 GPa | −16.378184 to −.593479 GPa | 2.360069 GPa |

Holding gamma-dependent theta unchanged but doubling only the vibrational
energy amplitude yields BM3 residuals **−5.999806/+3.924683 GPa**. Reversing
the residual sign reproduces the main paper's **−3.9/+6.0 GPa** limits almost
exactly. The demonstrated result is that our reconstruction needs an extra
thermal amplitude to recover the reported residual limits. An energy/volume
or atom-count normalization mismatch is a leading hypothesis; its location
and cause are not established. The same pattern in BM3 and Vinet points toward
the shared thermal treatment, but the two checks use the same reconstructed
Debye convention and are not independent evidence of an author-side coding
error. Transcription and phase-selection checks reduce those possibilities;
they do not eliminate an undocumented source convention or an error in our
interpretation. Matching two residual extrema does not recover the authors'
point-by-point predictions or fitting implementation.
Using an extra factor of two would correspond to n=4 per FeSi formula in the
current adapter; the paper does not define such an unphysical atom count.
Changing gamma0 instead is not exactly equivalent because it also changes theta.
A physical diagnostic refit gives gamma0 about 2.6091 (BM3) or 2.3551 (Vinet),
far outside the printed 1.30(.04) or 1.17(.04), while the other fitted
coefficients remain close. The several-GPa discrepancy is much larger than expected from the printed
coefficient rounding. Regression weights do not change predictions once the
published coefficients are fixed, although the undisclosed fitting procedure
still limits reconstruction of the published result.

The numerical factor of two is established as a successful diagnostic; a
factor-of-two error in the original publication is not established. Resolving
that attribution requires source code or author-generated cold/thermal
pressure components at specified source observations, with explicit energy,
volume and atom-count conventions.

The inferred extra amplitude is retained only as a named diagnostic. Both
published coefficient sets and all observations are preserved for review.
The published B2 parameter sets remain non-executable source evidence while
the equation/normalization discrepancy is unresolved. An authoritative source
clarification is needed before claiming a correction to the published model;
it is not a prerequisite for an independently derived fit to the observations.

The two executable B2 records are explicitly identified as Peritheos refits of
Fischer's measurements, with `record_kind: refit`, distinct `_refit` identifiers,
and complete `fit_provenance`. They preserve the measured pressures, physical
normalization, fixed assumptions and exclusions. The source observations are
cited directly; there is no fictitious executable published parent. Agreement
with these stored refit coefficients does not mean reproduction of Fischer's
published B2 coefficients.

The six Fe–9Si/B20 fits and the Fe–16Si fit are independently accepted; this
B2-specific hold does not change their acceptance.

## Accepted independent B2 refits

Run `uv run --python 3.9 python -m scripts.refit_fischer_2014_b2`.
The [complete evidence](../data/fischer-2014-b2-refits.json) includes full joint
covariance, correlations, five alternative starting guesses, pressure-weight
sensitivity, reference-parameter stress tests and three pressure-block holdouts.
All 114 exact `phase=B2` rows are used; no mixed-phase observations enter.
The objective is unweighted pressure residuals, chosen explicitly because the
source does not publish fitting weights or the P-T covariance. The pressure
errors already incorporate temperature uncertainty, so a joint independent-error
likelihood must not be inferred by adding those terms twice.

| Record | K0 (GPa) | gamma0 | q | Pressure RMSE (GPa) |
|---|---:|---:|---:|---:|
| `fesi_b2_fischer_2014_bm3_refit` | 230.948 ± 2.427 | 2.6091 ± .1812 | 1.7288 ± .5663 | 2.40747 |
| `fesi_b2_fischer_2014_vinet_refit` | 236.437 ± 2.468 | 2.3551 ± .1623 | 1.4168 ± .5627 | 2.35912 |

Errors are conditional local **1-sigma** estimates from
`cov=(J.T J)^-1*SSE/(114-3)`, treating pressure residuals as independent with
equal variance. Both Jacobians have full rank and the solutions lie inside
the broad bounds. Multistart fits recover the same minimum. The covariance
includes K0-gamma0-q cross terms; the material API uses `rt_eos.K0` for the
reference-EOS parameter so uncertainty propagation actually resolves that term.

Fixed V0=6.414 cm³/mol atoms, K0'=4.17, theta0=417 K and Tr=300 K are
inherited conditional assumptions. Physical n=2 and Z=1 are used for B2 FeSi.
The fixed-input uncertainties, calibration systematics, experimental run
correlations and model discrepancy are **not** included in the formal errors.
Pressure uncertainty propagated from this covariance describes parameter
uncertainty of the mean curve, not a predictive interval that includes the
roughly 2.4 GPa residual scatter.

The sensitivity results set the interpretation:

- K0-gamma0 and especially gamma0-q are correlated; the latter is about
  0.886 (BM3) and 0.869 (Vinet).
- Weighting by pressure errors alone moves gamma0 to about 2.208 and 2.026;
  these are sensitivity fits, not a recovered source likelihood.
- Changing K0' by an analyst-selected ±0.5 gives gamma0 about 2.00–3.14 for
  BM3 and 1.87–2.84 for Vinet at comparable residual scales. These bounds
  are stress tests, not uncertainty intervals for the source K0'.
- The source-table comparator V0=6.435 and analyst-selected theta0=300/600 K
  cases are retained separately. No preferred reference is silently changed.
- Holding out each pressure third gives prediction RMSE about 3.68–5.51 GPa
  for BM3 and 3.50–4.71 GPa for Vinet. Run identifiers are absent, so this is
  a pressure-coverage test, not validation on independent experiments.

These observations support executable **conditional empirical refits over the
measured range**, not a uniquely measured gamma0 or validated extrapolation to
core conditions. Gamma0 is extrapolated to the zero-pressure reference; over
actual source volumes gamma(V) is approximately 1.55–2.23 for BM3 and
1.54–2.07 for Vinet. Acceptance of the refits does not resolve the source
normalization hypothesis.

The new `fesi_b2` material has ideal fully occupied CsCl sites, Pm-3m (No. 221),
Fe 1a and Si 1b. Its reference a=2.666351276708074 Å is measured in S2 row 70
at 42.74768 GPa/1176.245 K, independently of extrapolated EOS V0. Both records
are visibly labelled Peritheos refits. BM3 is the default only within this new
refit-only material; there is no executable published B2 default to replace.
Native curves, inverses, covariance propagation and independent optimizer
agreement are tested. The generic ledger's B2 `parity` means stored-refit
reproduction and explicitly does not claim published-coefficient parity.

## Integration and validation

The catalog now contains 222 material documents, 221 deduplicated materials,
598 executable records, and 242 primary EOS papers (343 investigated papers
including the existing nonproduction register). The manifest, primary-source
and refit ledgers, paper ledger, source index, references, candidates and
changelog include this addition. The paper ledger explicitly notes the B2
withholding so seven accepted records cannot be mistaken for complete
acceptance of every source-owned candidate.

Tests cover workbook row counts/flags/markers, resource digests, occupied cell
contents, published coefficients and constraints, source residual bands,
independent/native agreement, inverses, both Fe–16Si stages and the B2 hold.
Full project check results are recorded in the task handoff.
