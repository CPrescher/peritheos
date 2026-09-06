# Sakai et al. (2025): nine mutually consistent Rydberg-Stacey scales

## Outcome

Nine production records are accepted from this paper: Cu, Re, Pt, W, Au, Mo,
MgO, B2-NaCl, and hcp-Fe.  They are source-reported room-temperature
generalized Rydberg-Stacey parameterizations, not nine refits of the
supplementary calculated grid.  The paper is also the first catalog use of the
four-parameter `RydbergStacey` model, which is implemented in both the Python
and Rust runtimes.

LitCurate exposed only the Cu row.  Direct inspection of the primary article
and official supplementary workbooks found the other eight scientifically
equivalent Table 2 rows.  All nine are retained so that the pressure-scale
network is represented without material-selection bias.

## Primary sources

Takeshi Sakai, Hirokazu Kadobayashi, Yuki Nakamoto, Haruhiko Dekura, Naoki
Ishimatsu, Saori Kawaguchi-Imada, Yusuke Seto, Oki Sekizawa, Kiyofumi Nitta,
and Katsuya Shimizu, “The equations of state of nine materials up to 0.43 TPa
for extreme pressure science,” *Communications Materials* **6**, 68 (2025),
[doi:10.1038/s43246-025-00792-5](https://doi.org/10.1038/s43246-025-00792-5).

The audit checked Equation 1 and Tables 1-2 in the article, Supplementary
Tables S1-S14, and the publisher's R-S formula workbook.  Retrieved-file
checksums are:

| Artifact | SHA-256 |
|---|---|
| Article PDF | `49ac913ab33167c50ab0da2dd654f31db1a3dd6b299ad07ad5d63b2df4ad23e3` |
| Supplementary information PDF | `63c62c71de8a332dda23b0dbe4da3bd745766fa2a6fbe2c0dfd443b64b0321cb` |
| Official R-S formula workbook | `b3d94aad7430d160d3cc2c27870d393bdda16f2c6750e3670666e6f205c4f49b` |
| Supplementary Data 1 workbook | `c559b29fea9c793831e4bc6deef4d62712b5e78821fc420969f709cbd1060ff9` |

## Equation and accepted coefficients

With `x = (V/V0)^(1/3)`, Equation 1 is

`P = 3 K0 x^(-3 Kinf') (1-x) exp[(3 K0'/2 - 3 Kinf' + 1/2)(1-x)]`.

It reduces to the conventional Vinet equation when `Kinf'=2/3`.  Sakai et al.
instead fix `Kinf'=5/3` for eight materials and fit it for Cu.  Uncertainties
below are copied as printed; the article does not state their confidence
convention or publish parameter covariance.

| Material | V0 (A3/cell) | K0 (GPa) | K0' | Kinf' | Disposition |
|---|---:|---:|---:|---:|---|
| Cu | 47.218 | 133.6 | 5.2035(6) | 2.000(2) | ACCEPT |
| Re | 29.468 | 352.6 | 4.411(8) | 5/3 fixed | ACCEPT |
| Pt | 60.409 | 274.1 | 5.128(9) | 5/3 fixed | ACCEPT |
| W | 31.724 | 296.0 | 4.415(16) | 5/3 fixed | ACCEPT |
| Au | 67.716 | 167.0 | 5.780(4) | 5/3 fixed | ACCEPT |
| Mo | 31.12 | 261.0 | 4.141(9) | 5/3 fixed | ACCEPT |
| MgO | 74.698 | 160.6 | 4.227(9) | 5/3 fixed | ACCEPT |
| B2-NaCl | 39.57(60) | 37.12(3.90) | 4.637(130) | 5/3 fixed | ACCEPT |
| hcp-Fe | 22.42(6) | 161.9(4.5) | 5.51(8) | 5/3 fixed | ACCEPT |

Most V0 and K0 values are fixed to cited ambient constraints; NaCl and Fe
also fit them.  These dependencies are represented in each record's pressure
calibration and source-lineage fields rather than presenting all coefficients
as new independent measurements.

## Data and numerical reproduction

Supplementary Tables S1-S8 contain the inter-material volume ratios used to
build the mutually consistent pressure network.  The final regression also
uses prior studies and sequential constraints, so the publication does not
provide one flat independent P-V table or the complete least-squares
covariance needed for a faithful coefficient refit.

Supplementary Table S9 is an official 26-point calculated grid for each final
EOS.  Peritheos bundles those nine grids as `derived_output_only`: they are
curve checkpoints, never mislabeled as experimental observations or used to
claim an independent refit.  Evaluating the printed Table 2 coefficients
reproduces the workbook to 0.105 GPa or better.  Cu agrees within 0.004 GPa;
the approximately 0.05-0.10 GPa residuals for the other materials are the
expected effect of the workbook retaining more internal coefficient precision
than Table 2 prints.

The deterministic check is
`scripts/reproduce_sakai_2025_nine_materials.py`.  Dataset hashes and native
model execution are covered by `tests/test_sakai_2025_nine_materials.py` and
`tests/test_rydberg_stacey.py`.

## Scope and cautions

These records are reference parameterizations at 300 K.  They should not be
interpreted as independent new measurements over every point in the nominal
range, and they do not include a thermal EOS.  Cu is the network anchor and is
reported to 1.2 TPa through the source's external constraint; the direct
simultaneous-volume experiments for the remaining materials reach roughly
0.31-0.43 TPa.  The records preserve those distinctions in their range and
calibration metadata.

No source row was rejected.  No extra curves were created from citations,
plots, rounded table values, or alternate fits not explicitly reported by the
authors.
