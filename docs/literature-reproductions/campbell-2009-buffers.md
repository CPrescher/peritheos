# Campbell et al. (2009): Fe–FeO and Ni–NiO thermal EOS

The publication of record is Campbell, Danielson, Righter, Seagle, Wang and
Prakapenka, *Earth and Planetary Science Letters* **286**, 556–564,
[doi:10.1016/j.epsl.2009.07.022](https://doi.org/10.1016/j.epsl.2009.07.022).
Audit date: 2026-09-20. The publisher PDF was recovered from Zotero item
`REIAHWUZ`, attachment `XUI8LMTT`; the publisher snapshot `U6TGWFQH` identifies
both official supplements. The complete source inventory and SHA-256 hashes
are in [the source manifest](../data/campbell-2009-source-manifest.json).

The latest Seagle follow-up reconstructs **79 of 81 individual pressures**
using the recovered Funamori, Boehler and Basinski sources alongside Seagle's
hcp calibration. The reconstruction remains conditional; the expanded
84-fit audit does not recover all three Campbell coefficients within their
printed errors. See the [fcc follow-up](#fcc-follow-up-with-the-recovered-primary-papers).
A subsequent [mismatch diagnosis](campbell-2009-mismatch.md) locates the main
FeO tension in Seagle’s hcp rows and tests five fcc thermal conventions;
improved pressure checkpoints still do not restore the coefficients.

## Decisions and material identity

Four experimental thermal parameterizations from Table 1 were investigated.
**The fcc-Fe and FeO records are deferred and non-executable** because the
combined-data coefficient mismatch remains unresolved. Their published
parameters and diagnostic results are retained for source inspection.

| Material | Record | Action |
|---|---|---|
| Pure fcc Fe | `fe_fcc_campbell_2009_bm3_mgd` | Deferred source record; excluded from executable catalog |
| Iron-saturated B1 FeO | `feo_campbell_2009_bm3_mgd` | Extend `feo` |
| fcc Ni | `nickel_campbell_2009_bm3_mgd` | Extend `nickel`; retain its existing default |
| NiO, source cubic-average B1 model | `nickel_oxide_b1_campbell_2009_bm3_mgd` | New card, separate from resolved rhombohedral `nickel_oxide` |

Section 4.1 explicitly declines to provide a new hcp-Fe parameterization and
adopts Dewaele et al. (2006). That existing record is not duplicated under
Campbell's name. Buffer integrations, polynomial representations and Tables
S4–S5 are derived results, not independent EOS fits. No computational EOS or
nonstoichiometric, rhombohedral or B8 FeO branch is added. No new opt-in refit
record is justified by these protocol-sensitive diagnostics.

Sections 2.1, 3 and 5 distinguish the **Fe1−xO starting powder** from the
**Fe1.00O product coexisting with Fe**. The authors select the stoichiometric
field of Stølen and Grønvold's model B: above 5 GPa at 900 K, with boundary
slope below −120 K/GPa. This is a phase-equilibrium inference, not a chemical
assay of each observation. All reported oxide volumes are B1; neither
rhombohedral FeO nor B8 was observed in the selected conditions. The existing
FeO card's misleading starting-powder label is clarified without changing its
Fischer coefficients or default.

Campbell calls all NiO observations B1. Roth's primary neutron study,
[Phys. Rev. **110**, 1333–1341 (1958)](https://doi.org/10.1103/PhysRev.110.1333),
Introduction, Fig. 1 and Table I, establishes the ideal NaCl nuclear structure,
the cubic paramagnetic state and the small rhombohedral antiferromagnetic
distortion (ambient Néel temperature 523 K). The new card explicitly represents
Campbell's cubic-average volume model, including its low-temperature MAP
observations. It does not claim those observations prove cubic symmetry,
resolve magnetic ordering or locate a high-pressure Néel boundary. This is a
single B1 structural approximation, not a metal/oxide aggregate. The occupied
ideal sites are Ni 4a and O 4b, each occupancy 1, in Fm-3m, Z=4. Its cell edge
is normalized to Campbell's reference volume; this is documented rather than
presented as an independently refined ambient cubic cell.

The fcc Fe diffraction cell is separately sourced from Nishihara et al.,
[Am. Mineral. **97**, 1417–1420 (2012)](https://doi.org/10.2138/am.2012.3958),
primary crystallographic deposition AMCSD 0019436,
[COD 9014789](https://www.crystallography.net/cod/9014789.cif): M620 at
17.59 GPa and 1273 K, a=3.521 Å, Fm-3m, Z=4, Fe at 4a with occupancy 1.
That structural cell intentionally differs from the extrapolated EOS
zero-pressure, 295 K cell. No Nishihara EOS is added in this audit.

## Exact equation, reference and parameters

Equation (5), page 560, is

\[
P(V,T)=P_{295}^{BM3}(V)+\gamma(V)[E_D(\theta(V),T)-E_D(\theta(V),295)]/V,
\]

with `gamma=gamma0*(V/V0)^q` and **`theta=theta0*(V/V0)^(-gamma(V))`**.
This is the existing `MieGruneisenDebye` **`variable_exponent`** option,
not the integrated-Grüneisen relation used by Fischer (2011). The gas constant
and volume units are converted explicitly in the independent implementation.
There is no electronic or anharmonic term, and the reference is **295 K**.

| Phase | V0, cm³/mol (fixed) | K0, GPa (fitted) | K0′ | theta0, K (fixed) | gamma0 (fitted) | q |
|---|---:|---:|---:|---:|---:|---:|
| fcc Fe | 7.076 | 133 ± 3 | 5 fixed | 470 | 1.95 ± 0.04 | 1.6 ± 0.6 fitted |
| B1 FeO | 12.256 | 146.9 ± 1.3 | 4 fixed | 380 | 1.42 ± 0.04 | 1.3 ± 0.3 fitted |
| Ni | 6.587 | 179 ± 3 | 4.3 ± 0.2 fitted | 415 | 2.50 ± 0.06 | 1 fixed |
| NiO | 10.973 | 190 ± 3 | 5.4 ± 0.2 fitted | 480 | 1.80 ± 0.04 | 1 fixed |

The italic entries in the rendered Table 1 are fixed. Its footnotes attribute
metal/NiO V0 to JCPDS (no card numbers), FeO V0 to McCammon and Liu (1984),
fcc-Fe K0′ to Funamori et al. (1996), FeO K0′ to Fei (1996), and Debye
temperatures to Chase (1998) or Knacke et al. (1991). Atom counts are 1 for
metals and 2 for monoxides per molar formula unit. Conventional-cell volumes
use `Vcell=Vmolar*4*1e24/NA`, with exact SI Avogadro constant. Missing errors
remain null; the source gives neither confidence convention nor covariance.
FeO's fixed K0′=4 makes its BM3 reference equivalent to BM2 without discarding
the source's explicit coefficient.

## Observations, pressure scale and temperature provenance

- [Official mmc1.pdf](https://ars.els-cdn.com/content/image/1-s2.0-S0012821X09004373-mmc1.pdf)
  contains Figure S1 and the supplementary table index; it has no extra EOS.
- [Official mmc2.xls](https://ars.els-cdn.com/content/image/1-s2.0-S0012821X09004373-mmc2.xls)
  contains all four worksheets. S2 has 25 Fe–FeO rows, including 15 fcc and
  10 hcp metal volumes. S3 has **101** Ni–NiO rows. The separate metal,
  oxide and NaCl molar volumes, their errors, metal/calibrant phases, paired
  volume differences, sample/NaCl temperatures and errors, pressures and
  errors are retained. Source row numbers identify each workbook observation.
- All **3434** tabulated buffer results (17 temperatures × 101 pressures × 2
  systems) are bundled as `published_model_output`. They are excluded from
  `fit_datasets` and used only as independent curve checkpoints.
- All **81** Seagle paired V–T observations already in the repository are
  linked to the two Fe-system fits. The primary article and official supplement
  were reread. They do not provide a row-wise pressure column or NaCl readings.
  Their original pressure methods are hcp Fe (Seagle 2006), fcc Fe (Funamori
  1996) and FeO (Murakami 2004); these must not silently be replaced by the
  Campbell fit being audited. Campbell does not specify its detailed
  reprocessing of those older rows.

A fresh mmc2 download differs in binary SHA-256 from the previously archived
copy used to transcribe S2; **every cell of every worksheet compares equal**.
Both hashes are recorded. CSV values retain the source numeric values, not
recalculated pressures or temperatures. Article/supplement redistribution
licensing is not asserted by the audit.

MAP pressures use Decker (1971) B1 NaCl; DAC pressures use that scale or Fei
et al. (2007) B2 NaCl. The workbook's “26 GPa” footnote is approximate: explicit
B2 phase rows below 26 GPa are preserved. The exact combined Decker/Fei
calibration route is not fully registered in the catalog, so pressure
recalculation metadata says `reference_eos_not_bundled` despite the complete
current-study calibrant observations.

The article prescribes the DAC correction
`Tsample=295+0.97*(Tmeas-295)`. Stored temperatures already include the correction. The correction is applied to the excess above
295 K, not to the whole absolute temperature. NaCl's effective temperature is
`(3*Tmeas+295)/4` in the article. Some Ni–NiO DAC workbook runs instead use
a **293 K ambient baseline** consistently in both expressions; others use
295 K. For example, S3 worksheet row 8 has Tmeas=2524 K, Tsample=2457.07 K,
TNaCl=1966.25 K and sigma_TNaCl=557.75 K, exactly the 293 K-baseline
values. Equation (4) prints the NaCl uncertainty as
`(Tmeas-295)/2`; the hot-DAC workbook values instead follow
`(Tmeas-Tambient)/4`. The EOS reference remains 295 K. Also, the
preceding prose prints a temperature difference where an average is intended.
The numerical workbook values are retained and these discrepancies are
reported without silently repairing the source. Temperature errors are
inclusive estimates; lattice errors originate from the scatter among peaks.
No pressure/temperature/volume covariance matrix is supplied.

Marginal coverage of S2 is 3.995–57.561 GPa and 873–2096.29 K; its fcc subset
ends at 19.974 GPa and 1812.08 K. Seagle extends fcc coverage to approximately
27 GPa and 2160 K, and B1 FeO to approximately 93 GPa and 2590 K. Heated/MAP
Ni observations span 2.336–65.9 GPa and 298–2457.07 K. Nine additional quenched
DAC rows lie at 293–295 K. These are observation extents, not rectangular
phase-stability domains, and the buffer grids extend beyond the fitted data.

## Independent reproduction and conditional refits

Run:

```bash
uv run python scripts/reproduce_campbell_2009_buffers.py --check
uv run python scripts/validate_primary_eos_refits.py --check
```

The dedicated script uses a separate 48-point Gauss–Legendre Debye integral,
BM3 expression in molar cm³ units and SciPy root inversion. Tests compare its
individual pressures with the executable public/native records, including
near the high-pressure, high-temperature end of the observations.

Equation (6) relates the derivative of the official buffer tables to the
metal–oxide volume difference. Five-point derivatives at 10 GPa for IW
(1500/2400 K) and 10/60 GPa for NNO recover the source's independent output.
The largest discrepancy is **0.00715 cm³/mol** (about 0.14%). A **0.02 cm³/mol**
tolerance allows the printed parameter uncertainties and precision; it does
not claim exact recovery of the authors' unrounded integration. Three- versus
five-point differences are separately recorded (at most 0.000251 cm³/mol).
These checks use the buffer tables as predictions, never as new observations
or EOS fits. Source-observation residuals supply additional individual-phase
checks. A paired-volume benchmark alone does not identify individual EOS.

No fitting sequence beyond least squares is stated, so all three free
coefficients are fitted simultaneously, with the adopted parameters held
fixed exactly. Main diagnostics minimize unweighted pressure residuals. For
Ni/NiO they exclude the nine quenched DAC rows, motivated by Section 2.2's
explicit residual-stress warning. This is an **audit selection**, not a claim
that the paper publishes that exact mask. The complete all-row refit and
pressure-error/effective-variance alternatives are also reported.

| Phase | Main rows | Free parameters (source → conditional refit) | Published/refit RMS, GPa | Ledger outcome |
|---|---:|---|---:|---|
| fcc Fe | 15 | K0 133 → 116.554; gamma0 1.95 → 2.03958; q 1.6 → 1.58052 | 0.7961 / 0.4481 | similar, subset only |
| FeO | 25 | K0 146.9 → 152.012; gamma0 1.42 → 1.35015; q 1.3 → 0.75080 | 1.3161 / 0.5150 | similar, subset only |
| Ni | 92 | K0 179 → 184.942; K0′ 4.3 → 4.36995; gamma0 2.50 → 2.23469 | 1.8571 / 1.6637 | similar |
| NiO | 92 | K0 190 → 209.235; K0′ 5.4 → 4.31613; gamma0 1.80 → 1.66634 | 1.3086 / 1.2411 | coefficient parity not achieved |

NiO is particularly sensitive to weighting: pressure-error weights give
`192.559, 5.22730, 1.78637`, while the fixed effective-variance alternative gives
`184.415, 5.81393, 1.77896`. This supports a regression-protocol explanation
for the discrepancy, but does not establish the authors' exact objective.
All-row Ni/NiO fits shift further because quenched states carry residual stress.
Published central coefficients remain unchanged.

The Seagle sensitivity reconstruction assigns all paired V–T rows to the
nearest of the seven printed Table 1 thermal-expansion lines, using their
nominal 20, 27, 50, 55, 72, 89 and 93 GPa series pressures and quoted errors.
The complete row assignments are saved. These are **reconstructed nominal
isobars**, not recovered individual pressures. Combined 29-row fcc and 106-row
FeO refits are disclosed in the JSON and show that this approximation does
not recover the original combined regression. In particular, the fcc q fit
reaches its lower bound. The original nominal-isobar sensitivity remains available alongside the
pointwise reconstruction below. No circular pressure reconstruction from
Campbell's own fitted Fe EOS is used to claim parity.

The [machine-readable reproduction](../data/campbell-2009-reproduction.json)
retains coefficients, local residual-scaled covariance, row counts, RMS values,
all weight/selection alternatives and solver results. These uncertainty
estimates belong to the audit, not to the publication. The main ledger applies
its normal numerical/combined-error comparison, but never promotes the
incomplete Fe subset to full-source parity.

## Seagle original-pressure reconstruction (2026-09-20)

The original hcp calibration was recovered from Seagle et al. (2006),
*JGR* **111**, B06209,
[doi:10.1029/2005JB004091](https://doi.org/10.1029/2005JB004091),
publication-of-record PDF, equations (1)–(4), Table 2 and paragraph 18.
The independent implementation is `scripts/reconstruct_seagle_2008_pressures.py`.
It evaluates BM3 with K0=164.8 GPa and K0′=5.33, gamma0=2.4, q=1.2,
theta0=380 K and a 300 K reference. Unlike Campbell's equation, its Debye
temperature follows the integrated Grüneisen relation
`theta=380*exp((2.4-gamma)/1.2)`. The two-atom hcp cell is converted to
molar volume with `Vmolar=Vcell*NA*1e-24/2`.

Two source inconsistencies are explicit audit interpretations:

- Table 2 prints V0=6.687 cm³/mol, whereas paragraph 18 states rho0=8.30 g/cm³.
  The latter gives **V0=55.845/8.30=6.728313253 cm³/mol**. Both versions are
  retained. This reference volume is derived from the stated density, never
  optimized against Campbell coefficients or Seagle benchmark pressures.
- Table 2 misprints the electronic coefficient's exponent and units. The
  implemented beta is **9.10e-8 kJ/g/K²**, with dimensionless electronic
  exponent 1.34. The electronic pressure in GPa is
  `gamma*(55.845/Vmolar)*beta/2*(Vmolar/V0)^1.34*(T²−300²)`.
  This follows the printed energy/pressure equations, including their
  multiplication by gamma(V). It is a documented dimensional and numerical
  interpretation, not an independently recovered Boness parameter table.

All **41 Table 3 iron pressure benchmarks** are transcribed in
[the checkpoint CSV](../data/seagle-hcp-pressure-checkpoints.csv).
Hot upstream/downstream temperatures are averaged; printed Fe volumes,
pressures and their parenthetical errors are preserved. The density-based
interpretation gives **0.2441 GPa RMS**, maximum difference **0.6058 GPa**,
and all 41 predictions within the reported pressure errors. Using the literal
Table 2 volume gives **2.4360 GPa RMS**, with only 16 of 41 inside those errors.
Five unambiguous hcp matches between Seagle (2008) Table 2's lower-eutectic
bounds and its supplement independently give **0.3463 GPa RMS** and all five
within reported errors. Two other article rows carrying hcp pressures have
fcc-sized supplement volumes (lower-eutectic rows 4 and 11); they are excluded
from this benchmark comparison without relabelling phases or halving volumes.
Melting-bound benchmarks are never added as subsolidus fit observations.

The [derived pressure CSV](../data/seagle-2008-reconstructed-pressures.csv)
retains all **81 source row IDs**. Its **65 hcp rows** contain reconstructed
individual pressures under both reference-volume conventions. The other
**14 fcc rows and two rows without an Fe volume retain blank individual
pressures**, with explicit status fields. Original observation CSVs are
unchanged. A separate column preserves the nominal Table 1 series assignment.
Measurement-only pressure errors propagate the printed Fe volume error and
global ±150 K temperature estimate by independent-input quadrature; missing
volume errors remain missing. These are not published pressure errors, omit
EOS systematic uncertainty and unknown correlations, and are not used as
regression weights.

The hcp series means agree closely with the approximate Table 1 isobars:

| Nominal series, GPa | hcp rows | Reconstructed mean, GPa | Individual range, GPa |
|---:|---:|---:|---:|
| 50 | 29 | 49.423 | 45.149–56.721 |
| 55 | 14 | 54.825 | 52.914–56.741 |
| 72 | 8 | 71.353 | 69.534–72.430 |
| 89 | 10 | 88.819 | 85.195–98.017 |
| 93 | 4 | 92.417 | 92.134–92.796 |

The larger individual excursions are retained, not clipped to nominal
isobars. Series assignment uses the previously documented nearest FeO
thermal-expansion line and is itself an audit inference.

At this hcp-only stage, the Funamori/Boehler full texts were inaccessible.
The user subsequently supplied both papers, and the reference-lattice source
was recovered independently. The [fcc follow-up below](#fcc-follow-up-with-the-recovered-primary-papers)
supersedes that access limitation and provides 14 conditional individual
pressures. The following hcp-only sensitivities are retained as historical
controls; the two rows without iron volumes still cannot be recalculated
through the Fe marker.

The conditional combined FeO results are:

| Input selection | Total rows | K0, GPa | gamma0 | q | Published/refit RMS, GPa |
|---|---:|---:|---:|---:|---:|
| S2 + 65 reconstructed hcp, density reference | 90 | 130.320 | 2.09544 | 0.01 (bound) | 3.6893 / 2.3787 |
| Above + 16 nominal-only rows | 106 | 136.601 | 1.79127 | 0.01 (bound) | 3.4286 / 2.4270 |
| S2 + 65 hcp, literal Table 2 V0 | 90 | 125.544 | 2.06192 | 0.01 (bound) | 2.9309 / 2.2865 |
| Above + 16 nominal-only rows | 106 | 130.534 | 1.82144 | 0.01 (bound) | 2.7372 / 2.2622 |

All minimize unweighted pressure residuals with the same three free Campbell
parameters and bounds as the earlier audit. The lower q bound is explicitly
recorded; local inverse-Jacobian covariance does not justify symmetric
confidence intervals for these boundary-active solutions. At this stage, the 29-row fcc fit
was a nominal-isobar sensitivity (K0=110.432, gamma0=2.08498, q=0.01 at
its bound), not a newly recovered individual-pressure fit. Reconstruction
therefore **does not establish coefficient parity** or identify Campbell's
unpublished treatment of Seagle observations. Published EOS coefficients and
the ledger's current-study baseline fits remain unchanged.

## Pressure-scale and regression-protocol investigation

A follow-up audit separates two questions: which pressures were supplied to
Campbell's regression, and which residuals/weights that regression minimized.
Run `uv run python scripts/audit_campbell_2009_protocol.py --check` to reproduce
[the complete protocol matrix](../data/campbell-2009-protocol-audit.json).
It now contains 12 data/pressure-scale selections, each with seven fits
(**84 results**, including matched unweighted controls). The first seven
selections below preserve the original 49-result investigation; five added
selections use the recovered fcc calibration described in the follow-up.

### What the primary sources establish

Campbell (2009), Section 2.2 and Tables S2–S3, explicitly uses B1/B2 NaCl for
its current-study pressures. Section 3 adds the older Seagle observations;
Section 4.1 adopts Dewaele's hcp-Fe EOS because it describes their iron data
well. That latter statement is **not an explicit instruction to recalibrate
the Seagle inputs with Dewaele**. The article does not specify the older
rows' reduction, regression weights or residual variable.

There is also a discrepancy in later source descriptions. The
[Fischer (2011) author manuscript](https://mineralphysics.uchicago.edu/Papers/FischerEPSL2011preprint.pdf),
page 10, says the Seagle and Campbell comparison observations were recalculated
from iron volumes using Dewaele hcp Fe or Campbell fcc Fe, and describes their
earlier pressures as NaCl-based. Seagle (2008), page 657, instead describes
iron/FeO calibration routes. Fischer's statement does not recover the missing
NaCl measurements, identify exactly which prior reduction it describes, or
establish Campbell's 2009 fitting inputs. Its reported average 0.3% NaCl/Fe
agreement cannot be applied as a correction to individual Seagle rows.

The initial DOI/title, institutional-repository and Zotero searches did not
recover the complete Funamori (1996) or Boehler (1990) articles. Those access
outcomes remain in the source manifest as historical records. Both PDFs have
now been supplied by the user. An indexed Maryland Campbell author-manuscript
link returned HTTP 404. No document requests were sent to authors.

### Controlled numerical comparison

The primary sensitivity retains Seagle (2006) pressures already checked
against independent benchmarks. A separate **Dewaele hypothesis** replaces
only the 65 Seagle hcp pressures. Campbell's 25 current-study observations
retain their reported NaCl pressures; nominal-only rows remain explicit.
The independent Dewaele implementation is cross-checked against the native
`iron_dewaele_2006_vinet_thermal` record. This replacement lowers those 65
pressures by **2.2176 GPa on average**, with differences from **−3.3807 to
−1.0764 GPa**. The row-wise differences are saved in the report.

All protocols retain the published fixed coefficients and fit K0, gamma0 and
q with the same bounds, from the published starting values. The volume
objective uses numerical EOS inversion at each observed P,T, rather than
rescaling pressure residuals and calling them volume residuals. Both pressure
and volume RMS are saved, so changes of residual units do not suggest an
artificial improvement.

| FeO pressure input | Rows | Residual / weighting | K0, GPa | gamma0 | q |
|---|---:|---|---:|---:|---:|
| Seagle original hcp scale | 90 | pressure / equal | 130.320 | 2.09544 | 0.01 (bound) |
| Seagle original hcp scale | 90 | volume / equal | 144.617 | 1.53921 | 0.01 (bound) |
| Seagle original hcp scale | 86 | pressure / effective variance | 146.782 | 1.33463 | 0.01 (bound) |
| Seagle original hcp scale | 86 | pressure / pressure error | 155.869 | 1.30987 | 0.53323 |
| Dewaele hcp hypothesis | 90 | pressure / equal | 132.082 | 1.91098 | 0.33259 |
| Dewaele hcp hypothesis | 90 | volume / equal | 139.922 | 1.49501 | 0.01 (bound) |
| Dewaele hcp hypothesis | 86 | pressure / effective variance | 141.014 | 1.36490 | 0.01 (bound) |
| Dewaele hcp hypothesis | 86 | pressure / pressure error | 152.049 | 1.33820 | 0.81354 |
| Published Campbell coefficients | — | unpublished protocol | 146.900 | 1.42000 | 1.30000 |

The original seven selections additionally include 106-row FeO combinations retaining all 16
nominal-pressure rows, current-study-only Fe/FeO controls, the 29-row nominal
fcc-Fe combination, and volume/effective-variance fits.

For these initial weighted comparisons, Seagle rows 33–36 lack iron-volume errors and are
excluded, leaving 86 or 102 FeO rows. The fcc weighted subset excludes row 6,
leaving 28 rows. Every such subset has its own unweighted pressure and volume
controls; excluded IDs are explicit, and no missing error is replaced by zero.

Effective variance is fixed at the published model and includes the shared
measured temperature in hcp-derived pressure:
`var(rP)=sigmaP²+(dPmodel/dV*sigmaV)²+(dPmodel/dT*sigmaT)²−2*dPmodel/dT*cov(P,T)`,
where `cov(P,T)=dPcal/dT*sigmaT²` for the reconstructed hcp rows. The volume
weight is its square root divided by `abs(dPmodel/dV)`. Current-study and
nominal-row covariance is unavailable and assumed zero; other cross-variable
covariances and calibration-systematic errors remain unknown. These are
conditional weights, not recovered author weights or confidence intervals.

**None of the 84 tested results places all three coefficients within their
printed source error intervals.** This uses the printed intervals without
assigning an unreported confidence level; it does not redefine the common
ledger's combined-error parity test. Dewaele recalibration removes the q
boundary in one unweighted pressure fit, but that change is not robust to the
objective/weighting choice. The experiment demonstrates sensitivity, not a
unique explanation or recovery of the source fit. No protocol was selected
for promotion to a new EOS record, and no global optimum or exhaustive search
over all conceivable procedures is claimed.

## FCC follow-up with the recovered primary papers

The user supplied the full [Funamori et al. (1996)](https://doi.org/10.1029/96GL00943)
and [Boehler et al. (1990)](https://doi.org/10.1029/JB095iB13p21731) papers.
Funamori p956 specifies **K0=120 GPa, K0′=5 at 1400 K** and takes its reference
volume from [Basinski et al. (1955), *The lattice expansion of iron*](https://doi.org/10.1098/rspa.1955.0102).
That free primary paper was also recovered. Funamori does not print the
analytic compression formula here: **using BM3 is an explicit audit assumption**.

Basinski Table 2, p462, brackets 1400 K with lattice parameters 3.6524 at
1347 K and 3.6622 at 1457 K. These are **kX, not Å**. Linear interpolation
of the lattice parameter gives 3.657121818 kX at 1400 K. Applying the historical
1.00202 Å/kX convention documented in [NBS Monograph 25, Section 3, p3](https://digital.library.unt.edu/ark:/67531/metadc13211/m1/9/),
then cubing and dividing the four-atom cell by four, gives
**V0(1400 K)=7.408637775 cm³/mol**. The exact interpolation used by Funamori
is not reported. The nine fcc lattice rows and 18 Funamori Table 1 fcc
observations are retained in separate transcriptions:
[Basinski lattice data](../data/basinski-1955-fcc-lattice.csv) and
[Funamori Table 1](../data/funamori-1996-fcc-table1.csv).

Boehler Table 1 and Figure 6 provide the fcc 300 K extrapolated reference
volume 6.835 cm³/mol, ambient **mean** expansivity 7.70×10⁻⁵ K⁻¹, and
`d ln(alpha)/d ln(V)=6.5±0.5`. The mean coefficient describes a linear
isobar relative to 300 K, not an instantaneous coefficient to exponentiate.
This audit evaluates its volume dependence at the unknown 300 K volume:

```text
alpha_mean(v300) = 7.70e-5 * (v300 / 6.835)^6.5
V(T) = v300 * [1 + alpha_mean(v300) * (T - 300)]
V1400 = v300 * [1 + alpha_mean(v300) * 1100]
P = BM3(V1400 / V0_1400, K0=120 GPa, K0_prime=5)
```

The first volume equation is solved numerically from each observed volume
and temperature. Choosing the 300 K volume in the expansivity law is an
audit assumption; Seagle does not fully specify its implementation.
Boehler's ambient thermal extrapolation and Basinski's independently measured
1400 K volume are retained as separate inputs, without forcing equality.
No coefficient is adjusted to the Seagle checkpoints or Campbell targets.

All **eight unambiguous fcc matches** between Seagle Table 2 and the melting
supplement fall within their printed pressure errors. RMS difference is
**2.8429 GPa**, with **+2.7667 GPa mean bias**; every residual is positive,
so this is approximate agreement, not exact recovery. The two Funamori
1400 K points differ by −0.6013 and +0.3539 GPa. Funamori prints no row-wise
pressure error there; these differences are not classified by an invented
error threshold. Some Seagle melting checks exceed the pressure/temperature
coverage of Boehler's expansion measurements. They test extrapolated behavior.
Funamori also discusses nonhydrostatic stress and possible hydrogen effects.
These calibration systematics are not quantified by measurement-only errors.

Run `uv run python scripts/reconstruct_seagle_2008_fcc.py --check` to reproduce
[the benchmark report](../data/seagle-2008-fcc-reconstruction.json) and
[the completed row-level pressure audit](../data/seagle-2008-completed-pressure-audit.csv).
This retains all 81 source row IDs: **65 hcp plus 14 conditional fcc pressures,
with rows 61–62 still blank because no iron volume is reported**. The fcc
series have mean pressures **20.640 and 28.921 GPa**, versus nominal 20 and
27 GPa. Separate delta=6 and 7 columns expose thermal-exponent sensitivity;
they are not total uncertainty bounds. Missing volume errors remain missing.
The original 65-pressure file remains unchanged as a historical stage.

### Expanded combined fits

The added selections use 29 fcc-Fe rows (15 Campbell + 14 Seagle) and 104 FeO
rows (25 Campbell + 79 Seagle), with separate 106-row FeO sensitivities adding
only the two missing-iron rows at their explicit nominal pressures. For FeO,
both the Seagle hcp scale and the hypothetical Dewaele replacement are kept;
the fcc pressures are identical in those two comparisons.

| Input selection | Rows | K0, GPa | gamma0 | q | Refit pressure RMS, GPa |
|---|---:|---:|---:|---:|---:|
| fcc Fe, reconstructed fcc pressures | 29 | 123.958 | 2.02403 | 0.01 (bound) | 0.8792 |
| FeO, Seagle hcp + reconstructed fcc | 104 | 134.304 | 1.90635 | 0.01 (bound) | 2.3538 |
| Above + two nominal rows | 106 | 135.065 | 1.86468 | 0.01 (bound) | 2.3827 |
| FeO, Dewaele hcp hypothesis + reconstructed fcc | 104 | 132.980 | 1.83018 | 0.25041 | 2.1486 |
| Above + two nominal rows | 106 | 134.178 | 1.81947 | 0.37134 | 2.1404 |

The table shows unweighted pressure fits. All five selections also have
unweighted volume fits and the same matched-error/weighting controls as the
earlier audit. Weighted fcc fits retain 28 rows; weighted FeO fits retain
99 or 101 rows after excluding missing iron-volume errors in Seagle rows
6 and 33–36. These exclusions and all 84 results are explicit in the
[protocol report](../data/campbell-2009-protocol-audit.json).

For reconstructed fcc-Fe inputs, the same measured Fe volume supplies both
pressure and the fitted volume. Effective-variance fits therefore extend
the earlier formula with `−2*dPmodel/dV*cov(P,V)`, where
`cov(P,V)=dPcal/dV*sigmaV²`, alongside shared temperature covariance.
Here the derivatives use the same molar-volume units. This is a comparison
with an independently sourced pressure calibration, but not independent
pressure-and-volume observations. Measurement errors can largely cancel
when model and calibrant have similar slopes; unknown calibration systematics
then matter greatly. Weights remain fixed at the published model, and the
pressure-error-only fits deliberately ignore covariance as a separate
heuristic. None is promoted as the recovered author procedure.

**Adding the fcc pressures does not establish coefficient parity.** None of
the 84 variants puts all three fitted coefficients inside their printed source
intervals. The published coefficient values are retained as source evidence. The
Fe/FeO records are deferred from execution; the common ledger retains the
current-study fit numbers but classifies their overall reproduction as
`parity_not_achieved`.

### Concrete outstanding evidence

Access to the Funamori, Boehler and Basinski primary papers is now resolved.
Their exact historical interpolation/thermal-reduction implementation remains
unreported; the reconstruction above is explicitly conditional.
The decisive missing author material is Campbell's actual fit-input table or
workbook: row IDs, pressures and temperatures assigned to the Seagle points,
any NaCl observations used to re-reduce them, final row mask, residual variable,
weights, and any staged or constrained fitting instructions. Those specific
inputs would distinguish the documented alternatives without tuning to the
published answer.

## Remaining limits and recovery routes

The publication of record, both official supplements, Zotero attachments,
Seagle's primary article/supplement, the author publication list, Crossref
relations and a DOI/title correction search were checked. No correction was
located. There is no missing Campbell primary table. Exact source regression
weights, covariance, quenched-row mask and Campbell's treatment of the older
Seagle pressure scale remain unresolved. The original hcp scale is now
benchmark-checked, and a conditional fcc reconstruction now passes eight
independent pressure checkpoints. Exact historical fcc reduction details
remain unresolved. These limits affect coefficient-level reproduction but
permit diagnostic evaluation of the specified published curves, but the
Fe/FeO reproduction discrepancy prevents their acceptance into the executable
catalog. The Fei/Decker calibration route remains unregistered.

The deferred Campbell fcc-Fe record is not an executable calibration for
Fischer (2011). Its citation and reported parameters remain available, but
the automatic calibration link has been removed. Fischer’s existing fit and
historical diagnostics are unchanged.
