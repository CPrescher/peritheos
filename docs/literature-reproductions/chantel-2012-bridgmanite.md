# Chantel et al. (2012): Mg and Fe-bearing bridgmanite

Primary DOI: <https://doi.org/10.1029/2012GL053075>
Open final article: <https://zenodo.org/records/3653697>

## Outcome

One of seven source-reported LitCurate candidates was accepted: the complete pure-MgSiO3 Table 3 BM3–Mie-Grüneisen-Debye thermoelastic parameterization. The eleven pure-MgSiO3 Table 1 rows are bundled verbatim. The 300 K pressure curve follows the nine room-temperature density states with RMSE 0.628 GPa and maximum absolute residual 1.166 GPa.

The other six source rows were withheld. Two printed BM2 density fits cannot be reproduced from Table 1 within their reported coefficient errors. Three Table 2 entries are third-order finite-strain fits to *adiabatic acoustic moduli*, not volumetric BM3 pressure equations; mapping them to Peritheos BM3 would change their meaning. The Fe-bearing Table 3 model prints `V0=25.50 cm3/mol`, which implies about 3.999 g/cm3 for the stated composition, irreconcilable with the independently printed 4.161(1) g/cm3 reference density. No correction is inferred.

## Seven source-row dispositions

| LitCurate row | Source result | Disposition |
|---|---|---|
| `litcurate_31bc34e8292357ee` | Mg-Pv density BM2, KT=257(2) GPa | **Withheld:** unweighted Table 1 refit gives 243.1 GPa; published curve exceeds the high-P densities by several printed density errors. |
| `litcurate_079361f39d81f154` | Fe-Pv density BM2, KT=246(2) GPa | **Withheld:** Table 1 refit gives 229.3 GPa; coefficient parity is not achieved. |
| `litcurate_542dc7a1f73d9d73` | Mg-Pv acoustic KS=247(4), K'=4.5(2) | **Withheld from EOS catalog:** correctly reconstructed as an acoustic finite-strain modulus fit, not a volumetric pressure EOS. |
| `litcurate_259fceb6dd4123fc` | Mg-Pv combined acoustic KS=252(1), K'=4.1(1) | **Withheld:** same semantic mismatch, and the combined Li–Zhang observations are external to this paper. |
| `litcurate_59accf80318ce824` | Fe-Pv acoustic KS=236(2), K'=4.7(1) | **Withheld from EOS catalog:** acoustic finite-strain modulus fit. |
| `litcurate_4c8341508b0d1a6a` | Mg-Pv Table 3 thermoelastic model | **Accepted:** `bridgmanite_chantel_2012_bm3_mgd`. |
| `litcurate_4de85b1c154b864a` | Fe-Pv Table 3 thermoelastic model | **Withheld:** internal `V0`/density contradiction changes the pressure reference state. |

The independently reconstructed current-study Mg acoustic fit gives KS=246.07 GPa and K'=4.571, within the Table 2 uncertainties. This verifies that the Table 2 numbers are sound while also demonstrating why they must not be mislabeled as isothermal P–V BM3 records.

## Fourteen comparison rows

These are all citation-reported literature comparisons and are rejected from this paper's production count: `litcurate_646cbed7a410ad5d`, `litcurate_67cb8f3a4139db77`, `litcurate_c6e72c7743e81ed6`, `litcurate_13f91e0f265496af`, `litcurate_82fa467c469361c3`, `litcurate_32212ca6614abcc8`, `litcurate_ed458727faf26675`, `litcurate_25bcc1466760438f`, `litcurate_94fc588faf4261ea`, `litcurate_e28b7d2629a0bd84`, `litcurate_b0e9cd7bc1f10123`, `litcurate_781d39e2b88b3db0`, `litcurate_adaf5cd6f8d92989`, and `litcurate_b2e503f20b0da628`.

## Zotero-ready citation

Chantel, J., D. J. Frost, C. A. McCammon, Z. Jing, and Y. Wang (2012). “Acoustic velocities of pure and iron-bearing magnesium silicate perovskite measured to 25 GPa and 1200 K.” *Geophysical Research Letters* 39, L19307. DOI: `10.1029/2012GL053075`.
