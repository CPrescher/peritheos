# Miozzi (2020): comparison after the author-described three-stage replay

All runs use the official EosFit7c 7.60 console: 36 RT rows → 131 rows with the cold coefficients fixed → 131 rows with V₀, K₀, K′, γ₀ and q free. θ₀=420 K, Tr=300 K and n=1 remain fixed. All stages use equal pressure-residual weights.

The protocol is the user’s report of an in-person conversation with Miozzi on 2026-10-03, including the subsequent clarification that V₀ is free and θ₀ fixed. Original weights and exact pressure reduction remain unconfirmed.

## Final parameters against the preferred thermal EOS (Section 3.3)

| Parameter | Paper | Speziale reconstruction | Original printed P | Tange recalculation |
|---|---:|---:|---:|---:|
| V₀ (Å³/cell) | 22.81† | 22.7031 ± 0.0984 | 22.6026 ± 0.1233 | 22.5783 ± 0.0921 |
| K₀ (GPa) | 129(1) | 139.42 ± 8.37 | 146.81 ± 11.35 | 150.49 ± 8.63 |
| K′ | 6.24(4) | 5.8700 ± 0.2738 | 6.0534 ± 0.3670 | 5.8015 ± 0.2664 |
| γ₀ | 1.11(1) | 1.9307 ± 0.1143 | 2.2357 ± 0.1908 | 1.9894 ± 0.1233 |
| q | 0.3(3) | 0.4296 ± 0.2220 | 1.5142 ± 0.3284 | 0.6431 ± 0.2380 |
| θ₀ (K) | 420, fixed | 420, fixed | 420, fixed | 420, fixed |

Fit errors are conditional EosFit standard errors; the paper’s parentheses preserve its reported error widths without assigning a confidence level. †The paper prints V₀=6.87(2) cm³/mol and 22.81 Å³/cell; the molar error converts to 0.06642 Å³/cell and its central value to 22.81581 Å³/cell. The differences below use the printed cell value 22.81.

## Signed differences: replay minus paper

| Parameter | Speziale reconstruction | Original printed P | Tange recalculation |
|---|---:|---:|---:|
| V₀ (Å³/cell) | -0.1069 (-0.47%) | -0.2074 (-0.91%) | -0.2317 (-1.02%) |
| K₀ (GPa) | +10.42 (+8.08%) | +17.81 (+13.81%) | +21.49 (+16.66%) |
| K′ | -0.3700 (-5.93%) | -0.1866 (-2.99%) | -0.4385 (-7.03%) |
| γ₀ | +0.8207 (+73.94%) | +1.1257 (+101.42%) | +0.8794 (+79.23%) |
| q | +0.1296 (+43.20%) | +1.2142 (+404.74%) | +0.3431 (+114.36%) |
| θ₀ (K) | +0 (+0.00%) | +0 (+0.00%) | +0 (+0.00%) |

## Stage 1: RT-only BM3 against Table 1

| Parameter | Paper RT BM3 | Speziale reconstruction | Original printed P | Tange recalculation |
|---|---:|---:|---:|---:|
| V₀ (Å³/cell) | 22.80(2) | 22.8063 ± 0.0976 | 22.6906 ± 0.1048 | 22.6981 ± 0.1002 |
| K₀ (GPa) | 129(6) | 128.51 ± 8.12 | 136.16 ± 9.59 | 136.57 ± 9.09 |
| K′ | 6.2(2) | 6.2869 ± 0.3064 | 6.4649 ± 0.3581 | 6.3080 ± 0.3313 |

## Stage 2: cold curve fixed

| Pressure input | γ₀ | q | RMS (GPa) |
|---|---:|---:|---:|
| Speziale reconstruction | 2.0494 ± 0.0790 | 0.5708 ± 0.1451 | 0.9419 |
| Original printed P | 2.3862 ± 0.1328 | 1.6548 ± 0.2169 | 1.2887 |
| Tange recalculation | 2.1451 ± 0.0870 | 0.8236 ± 0.1560 | 1.0047 |

The paper does not give the intermediate stage-two coefficients. These errors condition on each fitted cold curve and do not propagate its uncertainty.

## Residual comparison on each pressure input

| Pressure input | Published-coefficient RMS | Final-refit RMS | Final residual range |
|---|---:|---:|---:|
| Speziale reconstruction | 5.4522 GPa | 0.9281 GPa | -2.5768 to +3.1736 GPa |
| Original printed P | 6.9826 GPa | 1.2780 GPa | -2.9054 to +9.0802 GPa |
| Tange recalculation | 7.3105 GPa | 0.9845 GPa | -2.7400 to +3.2765 GPa |

The paper reports χ²=1.79 and pressure residuals between −3 and +3 GPa; it does not supply a directly comparable aggregate pressure RMS. Our unweighted RMS cannot be compared numerically with its χ². Each final fit has 131 rows, five free parameters and 126 residual degrees of freedom.

The Speziale reconstruction is a specified variable-q Debye scenario, not confirmation of the exact author calibration. Tange is a separate calibration choice. The 15 He pressures remain unchanged in all cases. The registered Tange refit and all source coefficients remain unchanged; these replays are diagnostics.
