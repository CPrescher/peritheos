# Validation

Peritheos keeps several kinds of validation separate. Passing one does not
imply another, and none of them is experimental validation of a prediction.

| Layer | What it establishes | Where to look |
|---|---|---|
| Structural | The `.eosmat` file matches its version's schema and declarations. File and software interoperability only. | `python -m peritheos validate` |
| Source check | Equation, coefficients, units, normalization, reference state, phase, uncertainty convention and represented range were traced to the cited primary publication or official supplement. Parsing or reproducing another library is never sufficient. | Each record's `check.source`; [source-audit reports](source-audits.md) |
| Reproduction | Whether a refit of the original data reproduced the published coefficients: the `check.reproduction` headline `reproduced`, `close`, `disagrees` or `untested`, plus any usage `warning`. | [Reconstruction reports](reconstruction-reports/index.md) (counts by headline and one section per record); [refit ledger](primary-eos-refits.md) |
| Numerical | The implemented equations, inversions, fits and uncertainty propagation behave as specified. | The sections below |

The words, criteria and caveats are defined once in [EOS labels](eos-labels.md).
A record that `disagrees`, is `untested` or carries a warning keeps its
published coefficients; only its use is qualified. Record-specific findings,
corrections and restrictions are stated in that record (see the
[material summaries](material-summaries/index.md)), its reconstruction report
and the linked source audits; the
[recommendation dispositions](research/eos-recommendation-dispositions-2026-10-17.md)
and [paper investigation ledger](paper-investigation-ledger.md) give the
paper-level decisions. Records removed during the primary-source review are
listed in the [material library README](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/materials/README.md).

The package ships the record-by-record source-audit ledger as
`peritheos/data/primary-source-audit.json`; each bundled record's audit sidecar
holds its audit date, source URL, DOI or stable locator, equation/table/page
locations, outcome and caveats. `scripts/apply_primary_source_audit.py`
(repository only) deterministically reapplies that ledger after a mechanical
Dioptas migration. Model inputs restored by the audit, such as Sokolova `n` and
`Z`, are recorded under `audit_corrections` with their primary location. All
validated records are instantiated in the test suite; a scientifically
validated label is never allowed to mask an incomplete model.

## Reference-state identities

Every isothermal model is tested for

\[
P(V_0)=0, \qquad K(V_0)=K_0.
\]

Models parameterized by derivatives are also checked numerically for their
supplied $K_0'$ and $K_0''$. Reference-isotherm thermal models are checked for
zero `thermal_pressure_increment` at their reference temperature. Absolute
thermal-pressure and Helmholtz models can retain nonzero zero-point or other
non-cold pressure there; their reference convention and corresponding
thermodynamic identities are tested instead.

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
[Paper versus spreadsheet](equations/thermal/multi-oscillator.md#paper-versus-spreadsheet).
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
| NaCl B2 | 0.65 | 1000 | 36.581254597 | 36.509419876 | — |
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
