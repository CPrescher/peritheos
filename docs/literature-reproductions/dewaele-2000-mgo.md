# Dewaele et al. (2000): MgO model-family audit

Primary source: A. Dewaele, G. Fiquet, D. Andrault, and D. Hausermann,
“P-V-T equation of state of periclase from synchrotron radiation
measurements,” *Journal of Geophysical Research: Solid Earth* **105**,
2869–2877 (2000), <https://doi.org/10.1029/1999JB900364>. The article is
available from the [HAL author archive](https://insu.hal.science/insu-03596948/document).

## Outcome

The paper reports five source-owned EOS parameterizations. The BM3-MGD record
was already present and remains the only preferred full P-V-T model. Four
explicit Table 3 room-temperature fits are now retained as nonpreferred model
sensitivities. They use the same combined observations and must not be counted
as independent experiments.

| Source row | Peritheos record | Disposition |
|---|---|---|
| BM, K0' fixed at 4 | `mgo_dewaele_2000_bm2_sensitivity_1` | ACCEPT, nonpreferred BM2 sensitivity |
| BM, K0 fixed at 161 GPa | `mgo_dewaele_2000_bm3_mgd_5` | ACCEPT, preferred BM3-MGD P-V-T EOS |
| Murnaghan | `mgo_dewaele_2000_murnaghan_sensitivity_2` | ACCEPT, nonpreferred functional-form sensitivity |
| Vinet | `mgo_dewaele_2000_vinet_sensitivity_3` | ACCEPT, nonpreferred functional-form sensitivity |
| Logarithmic | `mgo_dewaele_2000_natural_strain3_sensitivity_4` | ACCEPT, nonpreferred natural-strain sensitivity |

The source calls the last row the logarithmic EOS and cites Poirier and
Tarantola. Peritheos therefore represents it with `NaturalStrain3`, the
third-order logarithmic/natural-strain equation, rather than inventing a new
equation family.

## Primary data and fitting scope

Table 2 contains 61 new observations: 41 heated points and 20 at 300 K, spanning
0–53 GPa and 300–2474 K. The complete table is already bundled as
`mgo_dewaele_2000_table2_pvt`. The Table 3 room-temperature regression also
uses measurements from Fiquet et al. (1996), Fei (1999), and Utsumi et al.
(1998), reaching 65 GPa. Those external tables and their heterogeneous pressure
standards are not silently merged into the Dewaele dataset, so an exact
coefficient refit from the bundled rows alone is not claimed.

The new experiment used synthetic B1 MgO mixed with Pt in an argon pressure
medium. Pressures were obtained from the Jamieson, Fritz, and Manghnani (1982)
Pt shock-Hugoniot/Debye calibration. The source estimates pressure uncertainty
as `0.03 P + 0.0004 (T - 300)` GPa.

## Numerical checks

At `V/V0 = 0.667`, Section 3 reports approximately 145 GPa for BM3, 142 GPa
for Vinet, and 138 GPa for the logarithmic EOS. Direct evaluation gives
144.947, 141.291, and 138.919 GPa, respectively. The constrained BM2 gives
145.171 GPa. The Murnaghan row gives 153.097 GPa; the source does not print an
independent high-compression checkpoint for that curve.

Run:

```bash
uv run python scripts/reproduce_dewaele_2000_mgo_table3.py
```

No digitization was needed. Parameter covariance and the confidence convention
for Table 3 parenthetical errors are not stated. The existing preferred record
documents the paper's internal K0' error discrepancy (`0.05`, `0.2`, and `0.3`
in different passages) and retains the explicitly propagated `0.2` value.
