# Published KCl alternatives: Walker B1, Tateno and Chidester

Five explicitly nondefault records extend the catalog without changing the
preferred B1 Dewaele, B2 Dewaele, Walker B2, Tateno MGD or Chidester BM3
coefficients. Equations, tables and footnotes were checked against final articles,
including the locally held final Walker and Tateno PDFs. The source manifest
records checksums; article PDFs are not redistributed.

Run `PYTHONPATH=. python scripts/reproduce_kcl_variants.py --check`.
The [numerical report](../data/kcl-variants-reproduction.json) contains the
independent fits, native-evaluator comparisons and derived Holmes coordinates.
The script evaluates BM3/Vinet directly and integrates Debye energy with
Gauss-Legendre quadrature; 64/96-point convergence is checked separately.
This verifies equation implementation independently of the package evaluator.
All source errors retain their printed meaning; missing covariance and error
confidence remain unavailable. No diagnostic refit replaces published values.

## Walker (2002): B1 thermal EOS

[Walker et al., American Mineralogist 87, 805-812](https://doi.org/10.2138/am-2002-0701),
Tables 1 and 3, specifies the conventional B1 cell, Z=4, the conflicting tabulated V0=249.53 A3,
K0=177 kbar=17.7 GPa and fixed K0'=5. The reference is 23 Celsius=296.15 K.
The abstract reports alpha0*K0=0.0195(5) kbar/K, represented as
0.00195 +/- 0.00005 GPa/K; using rounded Table 3 alpha0=0.00011 instead
would silently change that directly identifiable product. No individual
elastic-parameter errors are published because of correlation.

All 30 Table 1 entries are retained: 23 measured sample states, six NaCl-only
spot checks and one derived ambient-volume anchor. The latter seven entries
are explicitly excluded from the fit. The table's misleading B2 volume-column
heading is corrected only in metadata: its volumes and Z identify B1.
Temperatures and pressures remain raw Celsius and kbar in the data; the
reproduction converts them to Kelvin and GPa. The source calibrant is the
NaCl BE2 thermal equation of Birch (1986), explicitly identified on page 806.
Its exact equation is not bundled here; paired NaCl lattice values are retained.
The [Walker follow-up](walker-2002-kcl.md) verifies this ancestry for both
B1 and B2 directly, recovers the original pressure-residual objective and
Excel Solver attribution, and recovers the original Birch reference coefficients and conditional run-normalized
pressure replay, with exact author-reproduction gaps documented.
The ESDs are spectrum-fitting errors, with no additional NaCl-EOS or temperature
error. The thermally derived V0 has no fabricated uncertainty.

The selected reference now follows Figure 1: 37.50 cm3/mol, or
249.080860076 A3 per conventional Z=4 cell. The raw tabulated anchor remains
249.53 A3. At the Figure reference, the complete unweighted joint fit holds
V0 and K0' fixed and yields K0=17.68342 GPa and alpha_KT=0.001927146 GPa/K,
RMSE=0.05382736 GPa. The published curve gives RMSE=0.05443363 GPa.
K0 rounds to the published value; beta differs by 1.17%, within its printed
error width. Exact solver and uncertainty parity remain unverified.
Page 808 specifies the joint
pressure-residual objective; the B2 Table 3 staging footnote does not apply
to B1. The [follow-up](walker-2002-kcl.md) documents the figure/table reference-
volume discrepancy and remaining original Solver inputs. Table 1 also lacks enough room-temperature compression points
to determine independent elastic errors.

The printed BE1 expression has minus signs inconsistent with its printed
positive-compression strain definition. We explicitly retain the standard
compression-positive BM3 convention already tested against Walker's B2 data;
the literal inconsistent expression is not silently treated as a new model.
The published B1 curve reproduces measured off-reference states within their
coordinate-error scale. Marginal 0-1.733 GPa and 303.15-873.15 K observations
are not a rectangular phase-stability guarantee; the reference anchor is at
296.15 K.

## Tateno (2019): distinct Pt scales and thermal models

[Tateno et al., American Mineralogist 104, 718-723](https://doi.org/10.2138/am-2019-6779),
final Table 1, contains these separate rows:

| Pressure coordinate | Thermal model | K0 (GPa) | K0' | Thermal coefficients |
|---|---|---|---|---|
| Holmes Pt | MGD | 17.4(2) | 5.77(4) | gamma0=1.8(2), q=0.7(3) |
| Holmes Pt | Linear thermal pressure | 17.7(3) | 5.73(4) | alpha_KT=0.0033(2) GPa/K |
| Sokolova Pt | Linear thermal pressure | 18.3(2) | 5.60(3) | alpha_KT=0.0037(1) GPa/K |

V0=54.5 A3 is fixed in every row. The MGD variant fixes theta0=235 K,
Tr=300 K and n=2 and uses the final Equation 6 integrated-Gruneisen law.
The existing Sokolova MGD row has K0 uncertainty 0.3 GPa; the linear row's
0.2 GPa is preserved separately. Superseded accepted-manuscript thermal
coefficients are never used.

The [official MSA AM-19-56779 deposit](http://www.minsocam.org/MSA/AmMin/TOC/2019/May2019_data/AM-19-56779.zip)
contains `6779TableS1 revised.xlsx`. The workbook is bundled verbatim and a
39-row CSV preserves its full numerical precision and Excel row numbers.
Each source row's pressure, T, Pt volume and KCl volume is kept together.

**Additional transcription correction:** the older bundled CSV had correctly
paired P/T/KCl values but still permuted the Pt lattice/volume columns in
runs 3/4. This audit corrects those Pt fields by matching the official Pt volume
to the corresponding existing printed lattice/volume group. All 110 changed
fields are retained in the source manifest with before/after values. Source
P/T/KCl/stress columns are unchanged. For example, run 3 at 6.0 GPa has
Pt V=59.16 A3, rather than 55.29 A3. The raw workbook contains cell volumes,
not lattice constants; the older rounded lattice constants remain article
transcriptions. They are not described as new workbook measurements.

Holmes diagnostic pressures are calculated separately from these paired Pt
volumes using the audited [Holmes (1989) Equations 11-12](holmes-1989-platinum.md):
V0=60.4000884 A3, P_T=798.31 GPa, eta=7.2119 and
alpha*B_T=0.0069426 GPa/K. The Sokolova pressure column is never used as a
Holmes pressure target. The numerical report retains every derived coordinate,
its source Excel row and a flag for T>=2000 K. No Sokolova pressure error is
reassigned as a Holmes error. Calibration covariance and systematic error are
not available, so no complete derived-pressure uncertainty is claimed.

All Tateno fits use 39 rows with equal pressure weights. Each alternative's
refitted coefficients agree within the source error widths. However, Holmes
pressure-coordinate reproduction is **conditional**: Tateno did not deposit
its Holmes pressures or the detailed rounding/thermal-reduction recipe;
Holmes states its approximate thermal term is adequate below 2000 K, and
Tateno reaches 2560 K. An executable calibration link enables a diagnostic
recalculation; it does not establish the author's exact reduction or extend
Holmes's thermodynamic accuracy. Source weights/covariance are unspecified.
The pressure extrema of Holmes records are labeled derived-coordinate bounds;
the 233.49 GPa maximum is a room-temperature state, not heated coverage.

## Chidester (2021): Vinet + MGD alternative

[Chidester et al., Physical Review B 104, 094107](https://doi.org/10.1103/PhysRevB.104.094107),
final Table I R-V row gives V0=34.3(5) cm3/mol, K0=13(1) GPa,
K0'=6.2(1), gamma0=3.4(4), q=1.0(1). All five are fitted, including q.
The B2 Z=1 cell conversion uses the exact Avogadro constant:
V0=56.95649 A3, error=0.830270 A3. Tr=300 K, theta0=235 K and n=2 are fixed.

The complete joint fit uses all 123 Dewaele room-temperature points and all
155 [author-deposited KCl effective-temperature rows](https://knowledge.uchicago.edu/records/6t4wc-w2146/files/SuppTable_KCl.csv?download=1).
The existing raw resources and their uncertainties are retained. The script
uses ordinal author CSV headers directly, and no raw Pt volume is invented.
The Dewaele 298 K observations enter the source regression as its 300 K
reference isotherm; this fit convention does not rewrite the raw temperature.
Cold ruby and high-temperature Dorogokupets-Oganov Pt ancestry remain distinct.
The laser-heating gradient model supplies an effective KCl temperature;
this is not a directly measured homogeneous-temperature dataset.

Independent unweighted fit in molar units gives V0=34.23179 cm3/mol,
K0=13.11896 GPa, K0'=6.209818, gamma0=3.387296 and q=0.985081.
All coefficients are within their reported error widths. Joint RMSE is
1.25160 GPa; published rounded coefficients give 1.35628 GPa. The integrated
Debye-temperature relation reproduces this complete-data fit, but the paper's
Equations 3-4 leave theta(V) implicit: its identification is a reproduction
inference, not a verbatim explicit source formula. No full covariance or
confidence interpretation of the printed errors is claimed.

## Follow-up source audit

The [Tateno-Holmes/Campbell-Heinz source-gap audit](kcl-tateno-campbell-source-gaps.md)
checks the complete MSA deposit, quantifies Pt-volume precision sensitivity, and
identifies Campbell's Mao (1978) ruby scale from its recovered final methods.
Exact Tateno Holmes reduction and source regression weights remain unavailable.

## Rights and validation scope

The source manifest separates original numerical data from contributor
transcription/correction. Source data terms remain unspecified where not
provided; any CC0 dedication covers only rights held by Peritheos contributors
in transcription, normalization, metadata and arrangement. Publisher article
text and figures are not covered.

Primary-source validation means the stored published equations, parameters,
units and limitations were checked. The refit ledger separately classifies all
five additions as `similar`; it does not claim statistical parity. In particular,
Walker coefficient drift and the conditional Holmes thermal coordinate remain
visible. Existing defaults and coefficients are unchanged.

## Integration check limitation

The catalog-wide refit generator at base `8ee9153` has pre-existing stale
Wang (1996) / Shim (2002) CaSiO3 ledger entries. In particular, the recovered
Shim table has both printed and digitized pressure columns, which the generic
selector rejects as ambiguous. This change appends only the five KCl outcomes
and retains existing unrelated ledger entries. The dedicated KCl reproduction
check passes; the global `validate_primary_eos_refits.py --check` remains stale
until those separate CaSiO3 changes and their generator are reconciled.

The full Python run finishes with 2,530 passes, 90.94% coverage and one
pre-existing failure: `test_primary_table_resources_are_complete_and_unchanged`
in `test_oganov_wang_chizmeshya_eos.py` still expects the earlier Wang (1996)
CSV hash. Both that test and the current Wang CSV are byte-identical to the
base checkout. No Wang/Shim source files or tests are modified here.
