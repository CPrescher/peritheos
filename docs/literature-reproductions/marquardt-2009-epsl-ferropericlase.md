# Marquardt (2009b): EPSL ferropericlase HS and LS equations of state

Audited 2026-09-08 against Marquardt, Speziale, Reichmann, Frost and Schilling,
*Single-crystal elasticity of (Mg0.9Fe0.1)O to 81 GPa*, EPSL **287**, 345–352,
[DOI 10.1016/j.epsl.2009.08.017](https://doi.org/10.1016/j.epsl.2009.08.017).
The supplied Zotero primary PDF was inspected, including a rendered page 347.
Its filename and SHA-256 are preserved in the record's primary-source audit.

## Catalog audit and disposition

Before this addition, the entire material catalog contained no record with this
DOI. `mg090fe010o.eosmat` contained two Solomatova (2016) Table 7 coupled-fit
reference branches. Those records and their coefficients remain unchanged.
The additions are the HS record `mg090fe010o_marquardt_2009b_hs_bm3` and
12 LS records `mg090fe010o_marquardt_2009b_ls_bm3_s1_01` through `_12`,
all in that existing material. The suffix 2009b distinguishes the EPSL source from the
separate Science paper; DOI, title, and author list establish ownership.

| Source result | Disposition |
|---|---|
| Section 3.1.1 high-spin isothermal BM3 | Accepted as a published, source-owned EOS |
| Section 3.1.2 low-spin BM3 families | Accepted: 12 alternative constrained fits from the subsequently supplied EPSL Table S1 |
| Section 3.1.3 empirical spin-population mixture | Not added: numerical A, B, C, D are not supplied in the inspected main paper or supplement |
| Section 3.2 Brillouin-only adiabatic moduli | Documented as elastic constraints, not an additional isothermal P–V EOS |
| Science DOI 10.1126/science.1169365 and `marquardt.som.pdf` | Separate source; no rows or coefficients imported into this EPSL fit |
| Solomatova DOI 10.2138/am-2016-5510 | Later coupled fit; no coefficient relabeling |

The main article describes LS fits with fixed V0(LS)/V0(HS) ratios spanning
0.96–0.99 and K0' spanning 3.5–4.5, fitting K0 to observations above 63 GPa.
The initially missing EPSL Table S1 has now been supplied and verified. It
contains all 12 fitted parameter triplets plus column-wise min, max and avg
summaries. Those three summaries are preserved as source data, not promoted
to independent EOS records or substituted for a uniquely selected LS model.

## Published equation and coefficients

Sections 2.2 and 3.1.1 explicitly specify third-order Birch–Murnaghan:

P = 3 K0 f (1 + 2f)^(5/2) [1 + (3/2)(K0' − 4)f],
where f = [(V0/V)^(2/3) − 1]/2.

| Parameter | Published value | Meaning |
|---|---|---|
| V0 | 75.62 ± 0.06 Å³ | Ambient single-crystal XRD conventional B1 cell, four formula units |
| KT0 | 158.2 ± 2 GPa | Isothermal modulus with an ambient Brillouin constraint |
| KT0' | 3.98 ± 0.14 | Isothermal pressure derivative |

The catalog treats the separately measured V0 as an anchored/fixed reference
volume, retaining its measurement uncertainty; the paper does not describe
an unconstrained three-parameter fit or publish a covariance matrix. Table 2
states one-sigma volume errors; the confidence of the fitted modulus errors is
not explicitly stated. Room temperature is represented as 300 K, not as a
measured temperature precision. The volume basis is neither per atom nor per
formula unit: V0/4 = 18.905 Å³ per formula unit.

The HS domain is P < 45 GPa. Observed HS pressures span 0.0001–43.8 GPa;
45–63 GPa is mixed spin and P > 63 GPa is LS. The title's 81.2 GPa maximum
belongs to Brillouin measurements; Table 2 XRD ends at 77.4 GPa. Neither is
the validity bound of the HS EOS.

Section 3.2 separately gives Brillouin-only **adiabatic** KS0 = 161.4 ± 1.1 GPa
and KS0' = 3.98 ± 0.14, using iterative density determination up to 32 GPa.
These numbers are not substituted for KT0. Above 32 GPa, longitudinal Brillouin
signals were masked; c11 and c12 use XRD-derived bulk moduli combined with
Brillouin shear constants, so they cannot serve as independent P–V observations.
Section 2.2 describes isothermal/adiabatic conversion using MgO thermal
expansion and Grüneisen assumptions. No thermal EOS is inferred from them.

## Observations, calibration, and reproducibility

The supplied complete [Table 2 CSV](../../peritheos/data/datasets/mg090fe010o-marquardt-2009-table2-pv.csv)
is copied byte-for-byte and linked by dataset identifier and SHA-256 in the
material. All 29 pressure, volume, and volume-error triples match the primary
table: 14 HS, 9 mixed-spin, and 6 LS. Ambient pressure is 10^-4 GPa.
`spin_region` and `used_in_solomatova_2016_fit` are curation annotations, not
printed EPSL columns; the latter is preserved from the supplied transcription
and is not used to select this fit. The record uses only source orders 1–14.
No pressure errors, ruby wavelengths, density, or shear observations are invented.
Section 2.2 reports neon medium and ruby fluorescence before/after XRD, but does
not identify an exact ruby scale; calibration remains partially resolved.

Run `python scripts/reproduce_marquardt_2009_epsl.py` for an independent BM3
check and a clearly labeled **validation-only** P–V fit. Published coefficients
give pressure RMSE 0.715844 GPa, volume RMSE 0.228649 Å³, and maximum absolute
volume residual 2.67946 sigma on the 14 HS rows. A volume-error-weighted fit
with V0 anchored gives K0 ≈ 138.243 GPa and K0' ≈ 5.65047. That diagnostic
does not recover the published coefficients; it omits the source's Brillouin
constraint, whose objective weight is unavailable. This discrepancy is retained,
not hidden by replacing published values or introducing an invented prior.
The global refit ledger marks the complete constrained objective unavailable;
the separate script documents the limited P–V validation that is possible.

Tests verify unique source ownership, preserved Solomatova branches, published
coefficients, cell basis, complete dataset checksum and spin selections,
independent pressure evaluation, validation diagnostics, and aggregate audit
coverage. The pre-change Python suite passed all 1860 tests.

Verification before supplement ingestion: 139 focused Python tests and all 1863 Python tests passed;
`cargo test -p peritheos --offline` passed 93 Rust tests including doctests.
Ruff, whitespace checks, and both generated documentation-index checks passed.
No unrelated pre-existing test failures were observed.

## Follow-up: origin of 158.2 GPa and fixed-derivative hypothesis

The paper reports KT0 = 158.2 +/- 2 GPa as a **best-fit result**, not as a
fixed input. Section 3.1.1 says room-pressure Brillouin measurements tighten
the KT0 constraint, but does not specify its numerical target, objective
weight, or whether KT0' was held at the Brillouin-derived value. We therefore
cannot reconstruct exactly how the optimization produced 158.2 GPa.

Two diagnostic calculations help assess possible interpretations:

- Applying KT = KS / (1 + alpha gamma T) to the separately reported
  Brillouin-only KS0 = 161.4 GPa, with the paper's alpha = 31.2e-6 K^-1,
  gamma = 1.524 and T = 300 K, gives KT0 approximately 159.13 GPa.
  Thus 158.2 GPa is not the direct conversion of that particular KS0 value.
  The Brillouin-only result comes from the separate iterative consistency
  analysis; it is not established as the exact ambient constraint target.
- A diagnostic BM3 fit to the 14 HS Table 2 rows, fixing V0 = 75.62 A^3
  and KT0' = 3.98, and minimizing squared volume residuals divided by the
  tabulated volume variances, gives KT0 = 156.22244 GPa and volume RMSE
  0.219130 A^3. This is within the published +/- 2 GPa interval. The trial
  fixes the derivative; it does not establish that the authors did so.

The identical reported isothermal and Brillouin-only derivatives,
3.98 +/- 0.14, motivate testing a fixed/shared derivative. However, neither
that equality nor the fact that 158.2 lies between 156.22 and 159.13 proves
the original fitting protocol. A weighted combination is a plausible
explanation only. No inferred constraint weight or fixed derivative is
promoted into the published record, and coefficient parity is not claimed.

### Which pressure range applies?

| Quantity or analysis | Pressure range |
|---|---|
| Published HS EOS interpretation | P < 45 GPa |
| Actual Table 2 rows supporting the HS EOS | 0.0001-43.8 GPa, 14 rows including the ambient anchor |
| Catalog HS evaluation interval | 0-43.8 GPa |
| Mixed-spin interval identified by the authors | 45-63 GPa; Table 2 rows at 45.7-62.7 GPa |
| LS XRD observations | 63.9-77.4 GPa, 6 rows; authors assign P > 63 GPa to LS |
| Complete Table 2 XRD dataset | 0.0001-77.4 GPa, 29 rows |
| Independent Brillouin bulk-modulus/density consistency analysis | Up to approximately 32 GPa |
| Full Brillouin experiment | Up to 81.2 GPa; higher-pressure c11 and c12 depend on XRD bulk moduli |

For use of the added HS EOS, the relevant measured upper limit is **43.8 GPa**,
not the 77.4 GPa XRD maximum or the 81.2 GPa Brillouin maximum.


## Recovered EPSL supplement and low-spin validation

The user supplied `1-s2.0-S0012821X09004889-mmc1 (1).doc` and
`1-s2.0-S0012821X09004889-mmc1 (2).doc`. They are byte-identical, SHA-256
`78c3a49a41795450b2ae39e97882268fab7da08cc89630ae0eaaded026160f04`.
The article PII and Tables S1/S2 match the EPSL source; this is distinct from
`marquardt.som.pdf` for Science. Table S1's layout was checked against a native
Word thumbnail and both tables were extracted from the Word file. Source
filenames, checksums, references and table roles are attached to the datasets.

| S1 row | V0 LS / V0 HS | Printed V0 LS (Å³/cell) | Fixed KT0' | Fitted KT0 (GPa) |
|---|---|---|---|---|
| 1 | 0.99 | 74.85 | 4.5 | 144 |
| 2 | 0.98 | 74.10 | 4.5 | 153 |
| 3 | 0.97 | 73.34 | 4.5 | 163 |
| 4 | 0.96 | 72.59 | 4.5 | 175 |
| 5 | 0.99 | 74.85 | 4.0 | 156 |
| 6 | 0.98 | 74.10 | 4.0 | 164 |
| 7 | 0.97 | 73.34 | 4.0 | 175 |
| 8 | 0.96 | 72.59 | 4.0 | 186 |
| 9 | 0.99 | 74.85 | 3.5 | 167 |
| 10 | 0.98 | 74.10 | 3.5 | 177 |
| 11 | 0.97 | 73.34 | 3.5 | 187 |
| 12 | 0.96 | 72.59 | 3.5 | 199 |

Each record preserves the printed V0 rather than recalculating it from the
rounded LS/HS ratio. These are extrapolated ambient LS volumes, not measured
ambient LS states. Each fit uses the same six Table 2 observations above
63 GPa (63.9-77.4 GPa), with V0 and KT0' fixed and KT0 fitted. No member is
marked preferred. Parameter uncertainties and covariance are unreported;
the spread among assumed models is not a statistical error bar.

The complete [Table S1 transcription](../../peritheos/data/datasets/mg090fe010o-marquardt-2009-table-s1-ls-models.csv)
contains 15 rows: 12 models plus the three summary rows. Its KT columns are
derived predictions, not extra compression observations. Only the columns at
69, 70.1, 75.8 and 81.2 GPa are pure-LS BM3 checks. The 45.9, 52 and 59 GPa
columns contain the mixed-spin response and cannot be reproduced by evaluating
a pure-LS derivative alone. The 81.2 GPa checkpoint extrapolates beyond the
77.4 GPa XRD maximum. Averaging coefficient columns does not generally yield
the average of these nonlinear curves.

`reproduce_supplement()` in the reproduction script independently computes
KT = -V dP/dV and checks all 48 LS model checkpoints against their printed
integers. The maximum difference is 0.47739 GPa, within rounding. Independent
volume-error-weighted fits to the six LS rows reproduce KT0 to within
1.484 GPa (0.951%) across the 12 models. This is **numerical similarity**;
no coefficient-uncertainty parity is claimed because the source supplies no
coefficient standard errors or exact objective weights. Published values
remain unchanged. The global refit ledger likewise finds 12 similar LS fits.

The complete [Table S2 transcription](../../peritheos/data/datasets/mg090fe010o-marquardt-2009-table-s2-elasticity.csv)
preserves 26 rows, all printed values and their parenthesized uncertainties,
with parsed numerical uncertainties in separate columns. It is supporting
elasticity evidence, never a P-V regression dataset. Its ambient KS value is
**165 +/- 3 GPa**, distinct from the main text's Brillouin-only fitted intercept
**161.4 +/- 1.1 GPa**. Applying the stated conversion at 300 K gives
**162.679 +/- 2.958 GPa** from Table S2, treating alpha and gamma as fixed for
this diagnostic. This does not equal the HS published KT0 = 158.2 GPa.
The supplement does not specify the exact ambient constraint implementation;
our earlier 159.13 GPa conversion and fixed-KT0' trial remain diagnostics,
not reconstructions of that implementation. The HS full-fit parity remains
unestablished even though the LS coefficient gap is now resolved.

Post-supplement verification: all 1866 Python tests and 93 Rust tests passed.
Ruff, whitespace checks and both generated documentation-index checks passed.
The two pre-existing Solomatova records compare exactly equal to their original
JSON objects. No unrelated pre-existing failures were found.
