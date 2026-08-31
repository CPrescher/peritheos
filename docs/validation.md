# Validation

Peritheos uses several complementary validation layers.

## Material-library validation levels

Structural `.eosmat` validation and scientific EOS validation are deliberately
separate. All 116 bundled material documents pass the format-3 validator and
were loaded with the Dioptas 0.10.0 material implementation. Their 147 records
also construct through Dioptas's Peritheos-backed EOS wrapper. These checks
establish file and software interoperability only.

The 2026-08-31 primary-source audit completed the classification of all 147
migrated records. One hundred sixteen are `primary_source_validated`; 31 are
`deferred`; none remains `pending_primary_source_check`. Promotion required a
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

### Audit outcome

The validated set includes the eleven Sokolova material parameterizations;
Fei Au and Ne; the Datchi and Dewaele diamond, c-BN, metal, and alkali-halide
fits; and
primary EOS results for graphite, FeO B1/B8, Fe7C3, B4C, hBN, Ru, Rh, Ir, Cr,
Zr, Pd, MgO, bridgmanite/post-perovskite, SiC, high-pressure silica phases,
garnets, carbonates, hydrides, ices VI/VII, coesite, zircon, the institutional
fcc-Pb pressure-marker report, and additional well-documented solids. Stable
publisher or institutional report URLs are accepted for older primary sources
without DOIs. This expands executable coverage without treating a migrated
label as evidence.

The deferred set is also intentional data. It contains 20 distinct unresolved
reference groups: older papers for which the accessible primary source did not
establish every stored value or convention; and ten Shen--Smith (2026) records whose APS accepted
manuscript remains under CHORUS embargo until 24 April 2027. These records stay
available for interchange but `Material.from_eosmat()` refuses them by default.

Primary-source findings changed or qualified several migrated records:

- Hanfland et al. report graphite `V0 = 35.12(2) angstrom^3`; the omitted
  `0.02` uncertainty is restored with an audit correction.
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
  unrelated lattice constants, and the represented range ends at 46 GPa.
- Holmes et al.'s platinum pressure scale is Equation (11), the universal
  (Vinet) EOS, not BM3. Its reference volume, 0--550 GPa static-isotherm range,
  and Equation (12) constant thermal-pressure extension were restored. The
  separate 32--660 GPa interval describes the shock Hugoniot.
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
  propagated as `0.00275(9) GPa/K`.

The audit also restores model inputs that the interchange migration omitted:
`n` and `Z` for the eleven Sokolova compositions, `n = 3` for the two SiO2
Debye records, and `Tr = 300 K` for the Bezacier ice records. Each addition is
recorded under `audit_corrections` with its primary equation/table location.
All 116 validated migrated records are then instantiated in the test suite; a
scientifically validated label is never allowed to mask an incomplete model.

Some complete papers still lead to intentional deferral. For example,
Martinez et al. (1996) express aragonite expansion as
`alpha(T) = alpha0 + alpha1*T`; the migrated record only carries a constant
mean expansivity. Peritheos therefore does not substitute that approximation
for the paper's exact P--V--T model. Other remaining groups expose only an
abstract or partial preview, or do not unambiguously establish the migrated
reference volume, equation order, fit constraints, and validity range.
The audit ledger now records those blockers DOI by DOI. In particular, the
Hanfland lithium record is not promoted because the primary work identifies a
Vinet fit while the migrated record is BM3; the FeO record cites a handbook
thermal-expansion chapter rather than an unambiguous original static-EOS
source; and the Hixson tungsten data are shock-derived while the migrated
standalone 300 K BM3 reduction could not be justified from accessible primary
tables. The remaining abstract-only groups are likewise not completed using
values from secondary catalogs.

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
behavior, NumPy arrays, and P-V-T broadcasting. The two-volume DAC inversion is
checked independently against its reduced pressure equation for every thermal
model, including both reduced and complete Sokolova configurations. Its
`f_dac` domain, non-heated volume pairs, non-invertible thermal models, and
fraction sensitivity are also tested. Invalid volumes, temperatures, non-finite
parameters, singular parameter sets, and states outside analytic domains are
explicitly rejected.

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
The catalog additionally checks all eleven Sokolova Table 1 parameter sets and
regresses MgO at eight compressions against the pressure-calculation values
printed in Figure 2. Fei catalog tests verify all four Table 1 parameter sets,
the publication-specific Debye-temperature law, reference states, and thermal
round trips.
Additional equation-level cases and their source DOIs are stored in
`tests/data/literature_reference_cases.json`; keeping the values in a data file
makes changes to scientific baselines visible in review.
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
