# Validation

Peritheos uses several complementary validation layers.

## Material-library validation levels

Structural `.eosmat` validation and scientific EOS validation are deliberately
separate. All 116 bundled material documents pass the format-3 validator. The
147 raw records transferred from Dioptas 0.10.0 were loaded with its material
implementation and construct through Dioptas's Peritheos-backed EOS wrapper.
Peritheos additionally supplies native, primary-sourced aragonite BM2,
B2-KCl P-V-T, RbCl-B2, and Correa and Benedict diamond Helmholtz records. Primary review consolidates one duplicate material,
removes two EOS reductions that their citations do not define, and excludes the
unreproducible Martinez global HT-BM3 reduction. Splitting the two distinct
phase-D reference volumes and adding the Benedict record produces 147 bundled
records; adding the Correa diamond branch produces 148. Two derived records
combine the Correa and Benedict thermal increments with the experimental
Dewaele 298 K Vinet isotherm, and adding the independently reproducible B4C
Berman public-data and Hemley neon refits plus the Campbell-Heinz RbCl record
produces 163 records in total. These checks
establish file and software interoperability only.

The primary-source audit covers all 163
bundled records. Every record is `primary_source_validated`; none remains
deferred or `pending_primary_source_check`. Promotion required a
direct trace of the equation, every stored parameter, units, reference state,
phase, published uncertainty convention, and represented data range to the
cited primary publication or official supplement. Parsing or reproducing
another library was never sufficient.

The package ships the complete record-by-record ledger as
`peritheos/data/primary-source-audit.json`. Each `.eosmat` record repeats its
audit date, source URL, DOI or stable report/article locator,
equation/table/page locations, outcome, and any
caveat. `scripts/apply_primary_source_audit.py` deterministically reapplies
the ledger after a mechanical Dioptas migration.

Primary-source traceability is complemented by the independent
[primary EOS refit campaign](primary-eos-refits.md). It attempts a Peritheos
fit for every record with sufficient direct observations and documents all 163
records, including selected columns, row count, published and refitted
coefficients, curve and refit RMSE, uncertainty comparison, and solver
diagnostics. The current campaign finds 82 uncertainty-parity matches and 33
additional numerically similar results. [Eight direct refits](primary-eos-refits.md#parity-not-achieved) do not
recover at least one published coefficient, while 39 records cannot be
directly refitted because row-level inputs or an executable source reduction
are unavailable. There are no unresolved extraction or solver failures. The
machine-readable results are in
[`docs/data/primary-eos-refits.json`](data/primary-eos-refits.json).

### Audit outcome

The validated set includes the eleven Sokolova material parameterizations;
Fei Au and Ne; the Benedict, Datchi, and Dewaele diamond, c-BN, metal, and alkali-halide
fits; and
primary EOS results for graphite, FeO B1/B8, Fe7C3, B4C, hBN, Ru, Rh, Ir, Cr,
Zr, Pd, MgO, bridgmanite/post-perovskite, SiC, high-pressure silica phases,
garnets, carbonates, hydrides, ices VI/VII, coesite, zircon, the institutional
fcc-Pb pressure-marker report, and additional well-documented solids. Stable
publisher or institutional report URLs are accepted for older primary sources
without DOIs. This expands executable coverage without treating a migrated
label as evidence.

Martinez et al.'s global aragonite P-V-T fit is intentionally not bundled. Its
printed coefficient table omits a required fitted reference volume, and an
independent reconstruction cannot recover its other coefficients from the 64
printed observations using the documented equations. The separately reported
staged BM2 P-V-T result is reproducible, validated, and executable.

Primary-source findings changed or qualified several migrated records:

- Campbell and Heinz publish the B2/B1 KCl volume ratio rather than an
  absolute B2 reference volume. The executable composite record multiplies
  `0.8483(57)` by Dewaele et al.'s experimental B1 `V0=62.36 angstrom^3`,
  yielding `52.899988(355452) angstrom^3`; the error covers the ratio only
  because Dewaele et al. publish no B1-`V0` uncertainty.
- Muñoz and Kunc fit theoretical wurtzite InN energies with the Murnaghan
  equation, not BM3. Its `V0=59.92519880888224 angstrom^3` is reconstructed
  from their theoretical `a0=3.483 angstrom` and `c0=5.7039 angstrom` in
  Table 1. Their cutoff-sensitivity estimate is documented but is not treated
  as a covariance-derived coefficient error.

- Hanfland et al. report graphite `V0 = 35.12(2) angstrom^3`; the omitted
  `0.02` uncertainty is restored with an audit correction.
- Clendenen and Drickamer's CoO coefficients belong to the Murnaghan equation
  printed as Equation 4, not Birch--Murnaghan. The corrected record uses
  Table II `a0 = 4.258 angstrom`, Table VI `B0 = 190.5 GPa` and `B0' = 3.9`,
  and the 0--30.8 GPa Table III data span. The source prints no fit errors or
  covariance and warns that these whole-range empirical constants need not be
  the true one-atmosphere derivatives. The [detailed CoO reproduction](literature-reproductions.md#coo-clendenen-1966)
  rules out decimal rounding alone. It identifies a strongly anticorrelated,
  weakly identifiable coefficient pair plus unavailable source rows, weights,
  and regression details; conditional fits recover either published
  coefficient when the other is fixed.
- Gleason et al.'s goethite record is a BM3 reference isotherm plus
  Mie--Gruneisen--Debye thermal pressure, not a standalone room-temperature
  BM3. The corrected validation now uses all 65 deposited P-V-T rows, but
  `K0=183.338 GPa` and a boundary `K0'=0` still do not reproduce the printed
  `140.3 GPa` and `4.6`. Depository row 32 is itself an isolated 100 degC,
  6.66 GPa state at `V=122.80 angstrom^3`; the published EOS expects
  `133.105 angstrom^3`. Removing it improves the curve residual but still gives
  `K0=177.642 GPa` and `K0'=1.5419`. The
  [detailed goethite reproduction](literature-reproductions.md#goethite-gleason-2008)
  documents the anomaly, a fresh Figure 3 digitization cross-check, the paper's
  malformed printed Equation 1, the 9 GPa non-hydrostatic discontinuity, and
  the unresolved source weighting/reduction details. The figure agrees with the
  primary table where points are shown, but omits the entire 100 degC series
  containing row 32. Fitting the 47 visible markers independently gives
  `K0=193.45 GPa` and `K0'=0.893` with the reconstructed thermal model, or
  `203.78 GPa` and `0.262` as plain BM3; fixing `K0=140.3 GPa` requires
  `K0'=6.112`, rather than `4.6`. This record is retained for provenance and
  published-curve reproduction but is not recommended for quantitative
  pressure-volume or thermoelastic work.
- Noguchi et al.'s NiO record is a shock-derived 300 K isotherm, not a static
  compression experiment. The official 1998 primary conference paper supplies
  the sample's pseudo-cubic `a0 = 4.177(1) angstrom`, the Mie--Gruneisen
  reduction, and the Murnaghan--Birch equation context; the final 1999 article
  reports `K0 = 191 GPa`, `K0' = 3.9`, and the 147.6 GPa range. The `V0` error
  is propagated from `a0`; neither article reports coefficient errors for
  `K0` or `K0'`.
- Somayazulu et al. call the B4C fit third-order Birch--Murnaghan in the
  abstract but second-order in the Figure 1 caption while reporting
  `K0' = 3.3(1)`. Peritheos retains BM3 because conventional BM2 fixes
  `K0' = 4`, and stores the conflict in `reported_inconsistencies`.
- Hazen and Finger's zircon value `K0 = 227(2) GPa` assumes `K0' = 6.5`.
  The inherited BM2 record was therefore corrected to BM3 with the derivative
  fixed, and its reported `V0 = 260.79(4) angstrom^3` was restored.
- Ross's unconstrained magnesite BM3 fit reports `V0 = 279.41(8)
  angstrom^3`; the migrated ambient measured volume and its uncertainty were
  replaced by the fitted reference volume and error.
- Haines et al. fit relative volumes for their Mo2C specimen. Its reference
  volume is now normalized from their measured ambient subcell rather than
  unrelated lattice constants, the represented range ends at 46 GPa, and the
  measured `V0` is fixed in the refit as required by the source's `V/V0`
  representation. The corrected result, `K0=325.874(9.500) GPa` and
  `K0'=4.909(651)`, is within combined two-sigma uncertainty of the published
  pair and is now classified `similar`; its remaining endpoint sensitivity is
  documented in the [dedicated reproduction](literature-reproductions.md#mo2c-haines-2001).
- Holmes et al.'s platinum pressure scale is Equation (11), the universal
  (Vinet) EOS, not BM3. Its reference volume, 0--550 GPa static-isotherm range,
  and Equation (12) constant thermal-pressure extension were restored. The
  separate 32--660 GPa interval describes the shock Hugoniot.
- Anderson et al.'s Au scale is Equation (29), not the inherited constant-
  expansivity reference-state approximation. It now combines the adopted
  300 K BM3 isotherm with the generic logarithmic-volume linear thermal
  pressure and reproduces the top rows of Table V. The published partial
  error on `(dKT/dT)V` is propagated; its additional unquantified contribution
  from `K0'` uncertainty remains an explicit caveat rather than an invented
  covariance.
- Frank et al.'s ice-VII record is the simultaneous three-parameter 300 K BM3
  fit from Equation (2). Its `12.4(1) cm^3/mol` reference volume is converted
  explicitly to the two-formula-unit cubic cell, and a Table 1 state is
  reproduced within its reported pressure uncertainty.
- Shim and Mao's cubic CaSiO3 fits and Richet et al.'s separate B1/B2 CaO
  fits are now traced to their primary tables. Reported one-standard-deviation
  errors are retained, including Shim's published `V0=45.58(5) angstrom^3`
  even though `V0` was fixed in that particular fit; parameters without a
  published error remain null. Molar CaO volumes are converted explicitly to
  the represented cells.
- Hazen and Finger's rutile GeO2/SnO2 reference volumes are recalculated from
  their Table 3 lattice constants. The paper reports errors for `K0` but no
  EOS-fit error for derived `V0`, so Peritheos does not manufacture one.
- Knorr et al.'s PbS record now uses their `a=5.9240(4) angstrom` sample,
  giving `V0=207.8955(421) angstrom^3`, and restores the published
  `K0=51.0(1.2) GPa` error. The unrelated migrated volume was removed.
- Katsura et al.'s wadsleyite record restores `V0=538.49(2) angstrom^3`,
  the source uncertainty on fixed `K0=169.2(9) GPa`, and the fitted
  `gamma0=1.64(2)` and `q=1.5(1)` errors. Fixed `theta0=814 K` remains
  errorless because the paper does not assign it an error.
- Zhao et al.'s jadeite 300 K reference isotherm now uses the unrounded
  `K0=124.5(4.0) GPa` and Table 1 `V0=403.32(8) angstrom^3`. Walker et al.'s
  B2-KCl record now keeps one internally consistent preferred Table 3 fit.
  Equation BE1 is represented as BM3 plus linear thermal pressure, rather than
  as a temperature-shifted reference state. Individual `V0`, `K0`, `K0'`, and
  `alpha0` errors remain null because the authors explicitly reject them as
  meaningful under the strong parameter covariance; the directly reported
  `alpha0*K0=0.0275(9) kbar/K` product and its uncertainty are retained and
  propagated as `0.00275(9) GPa/K`. The refit campaign now follows the preferred
  Table 3 staging: hold fictive `V0=53.53 angstrom^3`, fit `K0` and `K0'` to
  the eight 23-24 degC rows, and then fit the thermal-pressure product. It
  recovers `K0=23.775 GPa`, `K0'=4.416`, and
  `alpha_KT=0.002766 GPa/K`; the previous large discrepancy came from comparing
  an all-free simultaneous refit with the preferred staged row. The full
  diagnosis is in the [Walker KCl reproduction](literature-reproductions.md#kcl-walker-2002).
- Tateno et al.'s B2-KCl record now follows the final published article rather
  than the accepted manuscript: `gamma0=2.3(2)`, `q=0.8(2)`, and Equation 6's
  integrated-Gruneisen Debye-temperature law. The 39 observations are
  retranscribed from the official MSA Supplemental Table S1 workbook, which
  preserves the Pt-KCl row pairing that was lost across the split manuscript
  table. The corrected joint refit gives `K0=18.3446 GPa`, `K0'=5.60096`,
  `gamma0=2.29519`, and `q=0.82490`, restoring uncertainty parity. See the
  [Tateno KCl reproduction](literature-reproductions.md#kcl-tateno-2019).
- Chidester et al.'s B2-KCl record now declares both inputs to its published
  simultaneous fit: all 123 Dewaele et al. (2012) room-temperature B2 rows and
  all 155 new high-temperature rows. The corrected unweighted fit uses the
  thermodynamically integrated constant-`q` Debye-temperature relation and
  gives `V0=53.2036 angstrom^3`, `K0=23.9721 GPa`, `K0'=4.55798`,
  `gamma0=2.91714`, and `q=0.96524`. All five regain uncertainty parity, and
  the 1.582 GPa high-temperature RMSE matches the source's reported 1.6 GPa.
  The earlier `q=0` boundary came from fitting only the high-temperature rows
  with a different objective. See the
  [Chidester KCl reproduction](literature-reproductions.md#kcl-chidester-2021).
- Shen and Smith's ten Cu-anchored 300 K records reproduce the phase-specific
  Vinet fits in Equation (4) and Table II: Pt, Au, Ta, W, Mo, MgO, NaCl B1,
  NaCl B2, bcc Fe, and hcp Fe. The fixed reference volumes, fitted pressure
  intervals, and printed `K0`/`K0'` errors are retained. The article does not
  state the confidence level of those errors or publish their covariance, so
  Peritheos records neither and explicitly assumes independent parameters when
  propagating them. These room-temperature fits do not acquire a thermal model
  merely because the experiment was controlled at `298.5(5) K`.

The audit also restores model inputs that the interchange migration omitted:
`n` and `Z` for the eleven Sokolova compositions, `n = 3` for the two SiO2
Debye records, and `Tr = 300 K` for the Bezacier ice records. Each addition is
recorded under `audit_corrections` with its primary equation/table location.
All validated records are then instantiated in the test suite; a
scientifically validated label is never allowed to mask an incomplete model.

Martinez et al. (1996) provide two distinct aragonite reductions. Their staged
second-order route fits `V0,T` and `K0,T` along each isotherm with `K0'=4`, then
uses Equation (2) for linear `K0(T)` and Equation (3) for a directly linear
reference volume. Regressing the printed Table 6 values gives
`alpha_bar=6.484e-5 K^-1`, reproducing `6.5(1)e-5 K^-1`; unweighted and
error-weighted `K0,T` slopes are `-0.01969` and `-0.01702 GPa/K`, bracketing
the reported `-0.018(2) GPa/K`. Peritheos therefore provides this executable
BM2 P-V-T record (`V0=227.5(8) angstrom^3`, `K0=64.81(348) GPa`). The conflict
between its Table 6 `K0` error and the approximately 4.3 GPa summary error is
retained in `reported_inconsistencies`.

The paper's separate six-parameter global HT-BM3 reduction is excluded. Table
7 omits its fitted `V0(298 K)` and error, and its remaining coefficients do not
reproduce the 64 printed P--V--T observations under the documented equations
with ordinary pressure- or volume-residual least squares. No Peritheos refit is
substituted for that unreproducible published record.

### Full-text resolution of the former deferred set

Full published articles recovered during the 2026-09-01 audit resolve the
other former blockers directly from primary evidence:

| Material | Primary result represented | Important qualification |
|---|---|---|
| CsCl | Campbell and Heinz (1994) 300 K BM3, `K0=17.01(29) GPa`, `K0'=5.49(15)` | `V0` is fixed from the accepted `a0=4.123 angstrom`; all 13 new rows and nine corrected Yagi ratios are bundled. The 22-point refit achieves parity, with a qualification because the source does not publish its normalized-stress weights. |
| RbCl B2 | Campbell et al. (1994) 300 K BM3, `K0=17.9(10) GPa`, `K0'=5.23(29)` | The non-quenchable phase uses the paper's fixed hypothetical `rho02=3.3068(10) Mg/m^3`; all 24 Table 1 rows are bundled and achieve parity. |
| Fe3O4 | Mao et al. (1974) BM3, `K0=183(10) GPa` | `K0'=4.0(4)` is explicitly assumed; the `K0` error combines fit and pressure-scale contributions rather than defining covariance. |
| Li | Hanfland et al. (1999) Vinet, `K0=11.32(10) GPa`, `K0'=3.62(4)` | This is one empirical fit spanning bcc and fcc data, not a phase-specific bcc EOS; the stored conventional bcc cell contains two atoms. |
| majorite | Yagi et al. (1992) BM3, `V0=1513.1 angstrom^3`, `K0=161.2 GPa`, fixed `K0'=4` | The paper prints no fit error or covariance; the inherited `4 GPa` error was removed. |
| `(Mg0.4Fe0.6)O` | Richet et al. (1989) BM2, `K0=149(4) GPa` | The corrected formula is `Mg0.4Fe0.6O`; `V0` and its error are propagated from `a0=4.2805(4) angstrom`. |
| NiS | Campbell et al. (1993) BM3, `K0=156(10) GPa`, `K0'=4.4(12)` | The record is the metastable NiAs-type phase and uses the measured hexagonal reference cell. |
| phase D | Shieh et al. (2000) BM2, `K0=134(5) GPa` | AntA and AntB have distinct measured `V0=88.12(32)` and `87.191(97) angstrom^3`; two records preserve those reference-volume conventions instead of inventing a shared fit value. |
| cubic SnO2 | Ono et al. (2000) 300 K BM3, `V0=130.6(3) angstrom^3`, `K0=252(28) GPa`, `K0'=3.5(22)` | The paper's separate expansivity at 25 GPa is not a complete thermal EOS. The printed `130.6(3)` uncertainty is `0.3`, not `3.0 angstrom^3`. |
| SrO B1 | Liu & Bassett (1973) BM3, `K0=91.3(27) GPa`, `K0'=4.3(3)` | The fit includes the reported slight tetragonal distortion because the paper observes no volume discontinuity. |
| SrO B2 | Sato & Jeanloz (1981) BM2, `K0=160(19) GPa` | `V0=28.0224(9128) angstrom^3` is converted from the published extrapolated density; `K0'=4` is fixed. |

The same full-text pass confirms Scott et al.'s cementite result:
`V0=155.26(14) angstrom^3`, `K0T=175.4(35) GPa`, and `K0T'=5.1(3)`.
The reference volume is an adopted ambient measurement, not a simultaneously
fitted coefficient. None of these articles publishes a parameter covariance
matrix; Peritheos preserves each printed error and its stated convention, but
does not create covariance or silently reinterpret a fixed parameter as fitted.

Two cited migrated reductions remain removed rather than deferred. The
Fei-labeled FeO citation is a thermal-expansion chapter that does not define
the stored static BM3, and the Hixson tungsten shock paper publishes reduced
isotherm tables rather than the inherited standalone BM3 coefficients.

## Reference-state identities

Every isothermal model is tested for

\[
P(V_0)=0, \qquad K(V_0)=K_0.
\]

Models parameterized by derivatives are also checked numerically for their
supplied $K_0'$ and $K_0''$. Every thermal model is checked for zero thermal
pressure at its reference temperature.

## Derivative identities

Analytic isothermal bulk moduli are compared with numerical
`-V dP/dV`. Thermodynamically integrated Mie-Gruneisen characteristic
temperatures are checked against

\[
\gamma=-\frac{\partial\ln\Theta}{\partial\ln V}.
\]

The `variable_exponent` Debye-temperature law is instead checked directly
against Fei's stated $\Theta_D=\Theta_0(V/V_0)^{-\gamma(V)}$ convention. These
forms are intentionally not equated when $q\ne0$.

Caloric models are checked for

\[
C_V=\left(\frac{\partial E}{\partial T}\right)_V,
\qquad F=E-TS,
\qquad \frac{K_S}{K_T}=\frac{C_P}{C_V},
\]

and for the public-unit form

\[
C_P-C_V=10^4\alpha^2K_TVT.
\]

## Round trips and array behavior

All model families are tested for P-to-V and P,V-to-T round trips, scalar
behavior, NumPy arrays, and P-V-T broadcasting. Forward DAC-confined volume and
the two-volume DAC inversion are checked independently against their pressure
closure for every thermal model, including both reduced and complete Sokolova
configurations. Their `f_dac` domain, non-heated volume pairs, non-invertible
thermal models, and fraction sensitivity are also tested. Invalid volumes,
temperatures, non-finite parameters, singular parameter sets, and states
outside analytic domains are explicitly rejected.

## Fitting recovery

Synthetic P-V and P-V-T grids with known generating parameters are refitted.
The suite verifies parameter recovery, fixed parameters, bounds, absolute
uncertainty scaling, covariance dimensions, and input errors.

## Literature regressions

The MgO B1 Tange et al. (2009) pressure standard is regressed against values
printed in Table 5 at seven combinations of compression and temperature. These
cases exercise the Vinet reference isotherm, the publication-specific
Gruneisen law, Debye thermal energy, reference-temperature subtraction, and
the conventional-cell-to-molar-volume conversion together. The comparison
tolerance reflects the table's pressure rounding to 0.01 GPa. All catalog
standards also have reference-state, inversion round-trip, array, validity,
and uncertainty-propagation tests.

The expanded catalog is also checked against measured rows printed in Dewaele
(2019) LiF Table 4, Dewaele et al. (2008) diamond Table I, and Dewaele et al.
(2008) Ag/Ni Table II. Those are fitted experimental observations, so their
tolerances reflect the plotted/tabulated residuals rather than table rounding.
The B2 KCl/KBr regressions independently reproduce the exact linear thermal
increment in Dewaele et al. (2012), equation 2 and Table V.

The Bezacier ice VI/VII tests reproduce the temperature-dependent reference
state from equations (1)--(3), perform volume and temperature round trips, and
compare against representative measured rows in Table I within the published
fit residuals. Uncertainty tests combine the published fit errors with the
reported volume and temperature measurement errors.

The Sokolova diamond model includes numerical regression cases derived from the
accompanying Excel calculation. These intentionally validate spreadsheet
compatibility rather than a separate literal transcription of the journal
article's printed equations; the distinction is documented under
[Paper versus spreadsheet](equation-reference.md#paper-versus-spreadsheet).
The catalog additionally checks all eleven Sokolova records against the 2013
Tables 1 and 4 and regresses MgO at eight compressions against the 2016
pressure-calculation values printed in Figure 2. Fei catalog tests verify all
four Table 1 parameter sets,
the publication-specific Debye-temperature law, reference states, and thermal
round trips.
Additional equation-level cases and their source DOIs are stored in
`crates/peritheos/tests/data/literature_reference_cases.json`; keeping the
values in the packageable core crate makes changes to scientific baselines
visible in review and lets the published crate run its own tests.
The cases cover Birch-Murnaghan orders two through four, Vinet, Holzapfel,
natural strain, modified Tait, Murnaghan, Debye, Einstein, and Holland-Powell
families. Limiting forms and independent derivative identities reduce the risk
of copying the same algebraic mistake into expected values.

The fitting suite also compares a fixed-reference second-order Birch-Murnaghan
fit with its closed-form weighted least-squares solution, including the
absolute covariance. This checks the optimizer and covariance calculation
against an independently solvable statistical result.

Before publishing fitted parameters, users should still compare against the
original fitting program or publication used by their community and report the
model order, reference state, units, weighting, and covariance convention.

## External black-box comparisons

External libraries are used only through their documented public APIs. Their
source code, bundled parameters, and implementation structure are not inputs to
Peritheos. The following values were obtained with BurnMan 2.1.0 and Pytheos
0.0.2 for public classes explicitly named `Fei2007vinet`; pressure is in GPa.
The volume ratios use each Peritheos record's conventional-cell reference
volume.

| material | V/V0 | T (K) | Peritheos | Pytheos | BurnMan |
|---|---:|---:|---:|---:|---:|
| Au | 0.90 | 300 | 24.027229522 | 24.015900762 | 24.027229522 |
| Au | 0.80 | 1000 | 76.807819064 | 76.781194805 | 76.800857166 |
| Pt | 0.90 | 300 | 38.000024146 | 37.994819022 | 38.000024146 |
| Pt | 0.80 | 1000 | 112.879825670 | 112.862087358 | 112.870486265 |
| NaCl B2 | 0.75 | 300 | 15.988842533 | 15.988842533 | — |
| NaCl B2 | 0.65 | 1000 | 36.527254597 | 36.509419876 | — |
| Ne | 0.65 | 300 | 2.646283932 | 2.646283932 | — |
| Ne | 0.50 | 1000 | 14.180923257 | 14.168848666 | — |

The exact 300 K BurnMan agreement for Au/Pt and Pytheos agreement for NaCl/Ne
are useful convention checks, but are not independent validation. Pytheos's
small 300 K Au/Pt offsets and both libraries' small high-temperature offsets
are retained as observed differences: Peritheos continues to follow Fei et al.
equation 3 and Table 1 directly rather than changing equations or parameters to
match another implementation. The comparison can be rerun with
`scripts/compare_external_black_boxes.py`; neither external package is a
Peritheos dependency or test oracle.
