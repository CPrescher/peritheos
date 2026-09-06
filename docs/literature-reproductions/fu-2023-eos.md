# Fu et al. (2023): bridgmanite compression and CaSiO3 thermal refit

## Primary-source audit

- Citation: S. Fu, Y. Zhang, T. Okuchi, and J.-F. Lin, “Single-crystal elasticity of (Al,Fe)-bearing bridgmanite up to 82 GPa,” *American Mineralogist* **108**, 719–730 (2023), DOI [`10.2138/am-2022-8435`](https://doi.org/10.2138/am-2022-8435).
- Author-hosted final article: <https://www.jsg.utexas.edu/lin/files/FuBgmElasticityAM2023.pdf>
- Official MSA deposit: <http://www.minsocam.org/MSA/AmMin/TOC/2023/Apr2023_data/AM-23-48435.zip>
- Article SHA-256: `15806ef52e1b89a63b8c81e133055022b877793b2adc21050c73f37ae802106a`
- Official deposit SHA-256: `12121bb769862e781e8c232c85f2f3caae83d1a6d37f6e18d76ba6fde61b9486`
- Supplement SHA-256: `f262b3c16454fbeb26a64ef45b969d955199879f1ebb5cde1ba3d43f8d6e9b05`
- Audited 2026-09-05.

Online Materials Text S1 and Table S2 resolve an ambiguity in the main article: the same 43-point Fe10-Al14-Bgm compression series has two explicitly reported sensitivity fits with fixed `V0=163.75 Å³`: BM2 `K0=256(2) GPa`, and BM3 `K0=259(4) GPa`, `K0'=3.8(2)`. These are genuine equation-order/configuration alternatives, not run splits. All 43 non-ambient P–V rows are directly transcribed.

Text S3, main-text Equations 23–28, and Table S3 additionally publish a complete source-owned CaSiO3 composite refit: BM3 `V0=45.4 Å³`, `K0=248 GPa`, fixed `K0'=4`; MGD `theta0=1000 K`, `gamma0=1.42`, and `q0=2.65`. The authors refit Gréaux et al. (2019) and Sun et al. (2016) data rather than merely copying one cited row. The pressure-volume thermal model is executable; separate shear parameters are outside Peritheos EOS scope.

## Exhaustive LitCurate disposition

| Candidate | Origin | Decision | Reason |
|---|---|---|---|
| `litcurate_5661c5f835e32b74` | source | **rejected: not a compression EOS** | `KS=326 GPa` and `K'S=3.32` are finite-strain elastic-modulus coefficients referenced to 25 GPa; no V0 exists, and mapping them to a room-pressure P–V EOS would be a category error. |
| `litcurate_a3d1728dde300a01` | source | **accepted and primary-expanded to two records** | The ledger retained only the main-text `K0=256`, fixed-4 row. Official Text S1 supplies V0 and the co-reported BM3 sensitivity fit, yielding the two independent executable source configurations. |
| `litcurate_efb817fde2eec8ae`, `litcurate_1347395c4b56d29a`, `litcurate_a70856250f14a2ee`, `litcurate_7a71575938f3f9fb`, `litcurate_a2591f9f56f627a0`, `litcurate_438c7d5970e1cd0f`, `litcurate_a6ab9ce1494d1f92`, `litcurate_b601d5b21bd6f007` | citation | **rejected here** | Table 3 comparisons to earlier elasticity publications; they are not fits produced by this DOI. |
| `litcurate_6f0c0a415873bd49`, `litcurate_dc5880a87d2a1522`, `litcurate_f150ba8e2eb3c5db`, `litcurate_7e1c43165d898ecb`, `litcurate_57faf342bab1540f`, `litcurate_c570565a82d28af4`, `litcurate_62c01575e6a31000`, `litcurate_8f536afdfadf777c`, `litcurate_b46a77f6e6173d23` | citation | **rejected here** | Main Table 3 literature-comparison compression rows; each belongs to its cited primary DOI, not Fu et al. |

The 19 ledger candidates therefore yield two production records. Primary-source expansion adds one more production record—the Table S3 CaSiO3 BM3-MGD refit—which LitCurate did not enumerate. Net production yield is three.

## Composition discrepancy

The article consistently identifies sample 5K2667 as `Mg0.88Fe0.10Al0.14Si0.90O3`, but the “this study” rows of main Table 3 print `Mg0.93Fe3+0.048Fe2+0.032Al0.10Si0.90O3`. Because Text S1/Table S2 call the compressed specimen Fe10-Al14-Bgm and match its ambient lattice constants, the production records retain the independently characterized sample formula and document the Table 3 label as an unresolved internal inconsistency.

## Reproduction

Run `uv run python scripts/reproduce_fu_2023_eos.py`. It checks both bridgmanite curves against all 43 primary Table S2 observations (pressure RMSE 0.664 and 0.669 GPa for the rounded BM2 and BM3 coefficients), verifies inverse P–V round trips, and executes the CaSiO3 thermal branch at 300, 1200, and 2200 K. Exact refitting parity is not claimed because the source does not disclose coefficient covariance or weights.
