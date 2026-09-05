# Speziale et al. (2007) ferropericlase audit

## Primary source

S. Speziale, V. E. Lee, S. M. Clark, J.-F. Lin, M. P. Pasternak, and
R. Jeanloz, “Effects of Fe spin transition on the elasticity of (Mg, Fe)O
magnesiowüstites and implications for the seismological properties of the
Earth's lower mantle,” *Journal of Geophysical Research: Solid Earth* 112,
B10212 (2007),
[doi:10.1029/2006JB004730](https://doi.org/10.1029/2006JB004730).

The audit used the complete 12-page primary PDF. Sections 2 and 3.2 establish
the polycrystalline `(Mg0.80Fe0.20)O` sample, room-temperature synchrotron
powder XRD, cubic B1 structure, and a spin-transition onset at 40 GPa.
Paragraph 9 reports a Birch-Murnaghan fit of only the source's data below that
boundary: `V0=76.03(9) A^3`, `KT0=158(3) GPa`, and `KT0'=4.4(2)`. Three free
coefficients identify the executable BM3 form. No parameter was fixed.

The source does not identify the pressure calibration used for its Fe20 rows,
instead referring to Speziale et al. (2005) for further experimental details.
The calibrated pressures are preserved, but cross-scale recalculation remains
unresolved because no row-wise calibrant observable is given.

## Primary data and reproduction

The bundled CSV transcribes all 30 `(Mg0.80Fe0.20)O` rows in Table 1, including
every pressure, volume, printed uncertainty, and compression/decompression
flag. The 18 observations from `10^-4` to 39.5 GPa are selected for the
high-spin fit. Twelve observations from 42.1 to 62 GPa are retained as explicit
spin-transition exclusions.

Run:

```bash
.venv/bin/python scripts/reproduce_speziale_2007_ferropericlase.py
```

The rounded published BM3 curve has a 1.2609 GPa pressure RMSE against the 18
selected rows. Propagating each printed volume uncertainty locally into
pressure and combining it with the printed pressure uncertainty gives a
normalized-residual RMS of 0.7034; every point lies within 1.65 combined
standard deviations of the published curve.

A naive unweighted pressure-space refit gives `V0=76.0794 A^3`, `K0=133.05
GPa`, and `K0'=6.7496`, demonstrating the strong coefficient trade-off in
these relatively uncertain data rather than reconstructing the source fit.
The paper gives no objective, weights, covariance, or fitting software. This
record therefore claims numerical **curve parity within reported observational
uncertainties**, not recovery of undocumented regression settings.

## LitCurate dispositions

All 23 same-DOI candidates were checked:

| Source row | LitCurate identifier | Disposition | Reason |
| ---: | --- | --- | --- |
| 591 | `litcurate_028cabd3df34638e` | **Accepted** | Exact Fe20 high-spin BM3 mapped to `mg080fe020o_speziale_2007_high_spin_bm3_1`; complete own Table 1 series bundled. |
| 592 | `litcurate_70bd94d16a70cf9b` | Citation trace | Fe36 coefficients are from van Westrenen et al. (2005), not this source's fit. |
| 593 | `litcurate_c7cc098836ac6280` | Citation trace | Fe27 coefficients are from Jacobsen et al. (2005), not this source's fit. |
| 594 | `litcurate_3dc0bb739beff131` | Withheld: mixed identity | The `x=0.17-0.20` HS fit combines this paper's Fe20 rows with Lin et al. (2005) Fe17 rows; it is not an exact-composition material EOS. |
| 595 | `litcurate_09abac86cb23e2b5` | Withheld: external data | The weakly constrained LS fit uses only Lin et al. (2005) data above 80 GPa and represents `x=0.17-0.20`, not this source's Fe20 measurements. |
| 596 | `litcurate_949625379761379e` | Withheld: sensitivity row | Table 2 assumes `V0LS/V0HS=0.88`; it is not an independent observation or selected EOS. |
| 597 | `litcurate_39b0f7da1d067333` | Withheld: sensitivity row | Table 2 assumes ratio 0.89. |
| 598 | `litcurate_c39f1147bb6c0f4b` | Withheld: sensitivity row | Table 2 assumes ratio 0.90. |
| 599 | `litcurate_7ef51af2b6dc3226` | Withheld: sensitivity row | Table 2 assumes ratio 0.91. |
| 600 | `litcurate_7bc7576d82c86423` | Withheld: sensitivity row | Table 2 assumes ratio 0.92. |
| 601 | `litcurate_e068444002857401` | Withheld: sensitivity row | Table 2 assumes ratio 0.93. |
| 602 | `litcurate_3070039db3de1333` | Withheld: sensitivity row | Table 2 assumes the preferred ratio 0.94, but this remains one point in a 13-model trade-off grid fitted to Lin et al. data. |
| 603 | `litcurate_fe8f3e2b3d5fdece` | Withheld: sensitivity row | Table 2 assumes ratio 0.95. |
| 604 | `litcurate_a5fcd2910d8c5a9b` | Withheld: sensitivity row | Table 2 assumes ratio 0.96. |
| 605 | `litcurate_1125192c0985d811` | Withheld: sensitivity row | Table 2 assumes ratio 0.97. |
| 606 | `litcurate_423e5c7048e61711` | Withheld: sensitivity row | Table 2 assumes ratio 0.98. |
| 607 | `litcurate_402fee0220d8dff5` | Withheld: sensitivity row | Table 2 assumes ratio 0.99. |
| 608 | `litcurate_92e42092dfcdd063` | Withheld: sensitivity row | Table 2 assumes ratio 1.00. |
| 609 | `litcurate_9c7057702f19e404` | Citation/model input | Table 3's aluminous bridgmanite row averages values from other publications for lower-mantle modeling. |
| 610 | `litcurate_542492ccabd16463` | Withheld: modeling duplicate | Table 3 reuses the mixed-composition HS fit from row 594 as an assemblage input. |
| 611 | `litcurate_918f55eef55542ec` | Withheld: modeling input | Table 3's LS row is a selected assemblage parameterization ultimately based on Lin et al. data, not an independent Fe17 source fit. |
| 612 | `litcurate_83ee405d7c5c3cf4` | Citation/model input | Table 3's CaSiO3 coefficients average Shim et al. (2000) and Lee et al. (2004). |
| 613 | `litcurate_d0481c0901afa176` | Citation trace | The MgO modulus is a literature comparison, not a new source measurement or complete EOS. |

Production result: **one** independent EOS record. The Table 2 sensitivity
grid is documented but deliberately not inflated into 13 production records.
