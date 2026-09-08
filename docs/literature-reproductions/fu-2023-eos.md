# Fu et al. (2023): bridgmanite compression and CaSiO3 thermal refit

## Primary-source audit

- Citation: S. Fu, Y. Zhang, T. Okuchi, and J.-F. Lin, “Single-crystal elasticity of (Al,Fe)-bearing bridgmanite up to 82 GPa,” *American Mineralogist* **108**, 719–730 (2023), DOI [`10.2138/am-2022-8435`](https://doi.org/10.2138/am-2022-8435).
- Author-hosted final article: <https://www.jsg.utexas.edu/lin/files/FuBgmElasticityAM2023.pdf>
- Official MSA deposit: <http://www.minsocam.org/MSA/AmMin/TOC/2023/Apr2023_data/AM-23-48435.zip>
- Article SHA-256: `15806ef52e1b89a63b8c81e133055022b877793b2adc21050c73f37ae802106a`
- Official deposit SHA-256: `12121bb769862e781e8c232c85f2f3caae83d1a6d37f6e18d76ba6fde61b9486`
- Supplement SHA-256: `f262b3c16454fbeb26a64ef45b969d955199879f1ebb5cde1ba3d43f8d6e9b05`
- Audited 2026-09-08.

Online Materials Text S1 and Table S2 resolve an ambiguity in the main article: the same 43-point Fe10-Al14-Bgm compression series has two explicitly reported sensitivity fits with fixed `V0=163.75 Å³`: BM2 `K0=256(2) GPa`, and BM3 `K0=259(4) GPa`, `K0'=3.8(2)`. These are genuine equation-order/configuration alternatives, not run splits. All 43 non-ambient P–V rows are directly transcribed.

Text S3, main-text Equations 17–18 and 23–28, and Table S3 additionally publish a source-owned CaSiO3 composite coefficient set: BM3 `V0=45.4 Å³`, `K0=248 GPa`, fixed `K0'=4`; MGD `theta0=1000 K`, `gamma0=1.42`, and `q0=2.65`; and shear terms `mu0=126 GPa`, `mu0'=1.6`, and `etaS0=1.54`. The pressure-volume thermal branch is executable in Peritheos; the shear terms are retained in the audit but remain outside the EOS record's public model surface.

The official ZIP contains only `8435_supp1.pdf` and `8435fu.cif`: it has no row table, workbook, covariance, weights, or fitting code. Recovering the upstream studies changes the correct status from “no direct refit possible” to **completed refit attempt, coefficient parity not achieved**. It does not justify pretending that Fu's undisclosed regression protocol is known.

## Input-study reconstruction

| Study | Recovered numerical content | Role in Fu Figure S3 | Audit disposition |
|---|---|---|---|
| [Gréaux et al. (2019)](https://doi.org/10.1038/s41586-018-0816-5) | Official source-data XLSX, SHA-256 `cde0964709863c3b487a054597fd3dd549738aaa9dc3f0c637a996932944479a`: 34 cubic rows at 700–1700 K in `Figure 3b`, plus 13 tetragonal 300 K rows in `Figure 3a` | Solid observation symbols | The 34 cubic rows are candidate fit observations. `PNaCl` is measured; `PFS` is an independently calculated finite-strain pressure and is retained only as a diagnostic. The 13 tetragonal rows are excluded. |
| [Sun et al. (2016)](https://doi.org/10.1002/2016JB013062) | Author PDF, SHA-256 `10d8f6389b59fd7dcb94e349cbac0390e04a2ec55a815b8d238fd97797fccb2b`: all 144 Table 1 P-V-T rows with printed P and V uncertainties | Colored observation symbols/isotherms at 1200–2200 K | The 140 rows at 1200–2200 K match Fu's legend and range. Four rows at 2400–2600 K are excluded from the primary diagnostic. Temperature uncertainty is stated only as 50–100 K, not row by row. |
| [Kawai and Tsuchiya (2015)](https://doi.org/10.1002/2015GL063446) | FPMD finite-strain elastic parameterization using the authors' 2014 P-V-T EOS | Dashed calculation curves | Analytical comparison, not counted as observations. |
| [Li et al. (2006)](https://doi.org/10.1016/j.pepi.2005.12.006) | AIMD P-V-T and elastic parameterization; inspected author PDF SHA-256 `712812bdf18a83a454413e9b686f68460fd6e9536b91b166ae08f667e3f2fc73` | Dotted 2000 K calculation curves | Analytical comparison, not counted as observations. |
| [Thomson et al. (2019)](https://doi.org/10.1038/s41586-019-1483-x) | Independent experimental composite EOS; inspected accepted manuscript SHA-256 `b233f2d5da32c01e97e810a2c5a5d71667fed16d899253243188188a8af5c295` | Open circles | Explicitly excluded by Fu because of inconsistencies; never a fit observation in this audit. |

This yields 174 likely observation rows: 140 Sun P-V-T states and 34 Gréaux cubic states. The Gréaux states contribute pressure, adiabatic bulk modulus, and shear modulus, so the reconstructed joint objective has 242 scalar outputs. The 174-row selection is an inference from Figure S3, not a source-published row manifest.

## Composite model and weighting audit

The audit implements the BM3 pressure and finite-strain bulk/shear expressions in Fu Equations 17–18, followed by the volume-dependent Mie-Grüneisen-Debye corrections in Equations 23–28. Gréaux density and velocities are converted with `KS=rho*(Vp^2-4*Vs^2/3)` and `mu=rho*Vs^2`; their printed standard deviations are propagated into both moduli. Equation 27 is treated as the dimensionally standard Debye heat capacity before evaluating `Delta(CV*T)` in Equation 24.

Fu et al. do not disclose the residual definition, row or output weights, pressure-column choice, staged-versus-joint strategy, coefficient covariance, or fitting software. Four explicit joint-fit diagnostics therefore test sensitivity rather than assert a unique reproduction:

| Parameter | Fu Table S3 | Unweighted P, KS, mu | Printed sigmas | P-V propagated sigmas | Equal observable groups |
|---|---:|---:|---:|---:|---:|
| `V0` (Å³) | 45.4 | 45.1346 | 45.5064 | 45.5264 | 45.1265 |
| `K0` (GPa) | 248 | 258.271 | 237.701 | 236.596 | 257.777 |
| `mu0` (GPa) | 126 | 129.662 | 126.645 | 126.481 | 129.731 |
| `mu0'` | 1.6 | 1.48479 | 1.57804 | 1.58375 | 1.48825 |
| `gamma0` | 1.42 | 2.44196 | 1.57408 | 1.55879 | 2.39676 |
| `q` | 2.65 | 1.83946 | -1.03500 | -1.02679 | 1.80306 |
| `etaS0` | 1.54 | 1.37236 | 1.37696 | 1.37693 | 1.37206 |

The unweighted joint solution is also registered as the opt-in executable EOS record `ca_perovskite_fu_2023_candidate_data_unweighted_bm3_mgd_refit`. Its pressure-volume record stores the fitted `V0`, `K0`, `gamma0`, and `q`; the jointly fitted shear coefficients remain in the audit artifact because the EOS API does not encode a shear model. This record is a reproducible Peritheos result, not Fu et al.'s undisclosed regression.

All optimizations converge, but none reproduces Table S3; `gamma0`, `q`, and `etaS0` are especially protocol-sensitive. The propagated P-V diagnostic includes printed pressure and volume errors but cannot include temperature uncertainty because Sun reports only a 50–100 K range, not row-level values. Choosing the derived Gréaux `PFS` column instead of measured `PNaCl`, or fitting only pressure, also fails to recover the published thermal pair. The bounded conclusion is therefore non-parity caused by an underdetermined source protocol, not evidence that the published curve is numerically unusable.

## Redistribution boundary

No new Gréaux, Sun, Kawai–Tsuchiya, Li, or Thomson rows are bundled. The Gréaux article is Springer Nature copyright and neither it nor the Sun author PDF provides a reusable table-data license. The repository contains only contributor-authored code, source checksums/counts, aggregate fit results, and the already-authorized Fu Table S2 transcription. Model-generated pressure/volume round trips are explicitly labeled **analytical checkpoints** and are never registered as `fit_datasets` or called observations.

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

Run `uv run python scripts/reproduce_fu_2023_eos.py`. It checks both bridgmanite curves against all 43 primary Table S2 observations (pressure RMSE 0.664 and 0.669 GPa for the rounded BM2 and BM3 coefficients), verifies inverse P–V round trips, and reports the CaSiO3 analytical checkpoints separately from the 174-observation coefficient audit.

To recompute the CaSiO3 sensitivity fits from lawfully obtained source files:

```bash
pdftotext -layout SunLowerMantleEoSJGR2016.pdf sun-2016.txt
uv run python scripts/audit_fu_2023_casio3_refit.py \
  --sun-text sun-2016.txt \
  --sun-pdf SunLowerMantleEoSJGR2016.pdf \
  --greaux-xlsx 41586_2018_816_MOESM1_ESM.xlsx
```

The checked-in aggregate result is [`docs/data/fu-2023-casio3-refit-audit.json`](../data/fu-2023-casio3-refit-audit.json). It contains no third-party observation rows.
