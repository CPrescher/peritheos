# Miozzi full-MGD replay with confirmed author n=2 — 2026-10-03

## New author information

The user relayed an in-person confirmation from Miozzi that the authors used full MGD, not q-compromise, with n=2 and V0 approximately 6.87 cm3/mol in EosFit. The previous three-stage protocol remains: cold RT fit; all data with cold coefficients fixed; all five parameters free from that solution. V0 is free in the last stage and theta0 stays fixed at 420 K. This is user-relayed author information, not a recovered original EosFit input/session file. Exact row selection, weights and pressure reduction remain unresolved.

Our previous n=1 runs used the same approximately 6.87 cm3/mol volume convention. They were internally normalized per mole of Fe, but did not reconstruct the newly confirmed author setting. The new run uses n=2 without changing Fe volumes or temperatures. It therefore doubles the vibrational energy amplitude at fixed gamma and theta. With full MGD the fitted gamma does not transform by an exact factor of two, because gamma also sets volume-dependent theta. No q-compromise interpretation is used for this replay.

## Final three-stage results

All runs use the official EosFit7c 7.60 console, unit pressure-residual weights, 131 observations, theta0=420 K, Tr=300 K and n=2. There are 36 RT observations in stage one. The final stage has five free parameters and 126 residual degrees of freedom. V0 is reported below as the original two-atom hcp cell volume; console inputs use molar volume of Fe. Errors are conditional local standard errors, excluding pressure-calibration and fixed-parameter uncertainty.

| Pressure input | V0 (Å³/cell) | K0 (GPa) | Kprime | gamma0 | q | RMS (GPa) |
|---|---:|---:|---:|---:|---:|---:|
| Published preferred coefficients | 22.81 | 129 ± 1 | 6.24 ± 0.04 | 1.11 ± 0.01 | 0.3 ± 0.3 | Not reported |
| Supplementary printed pressures | 22.60890 ± 0.12411 | 146.082 ± 11.387 | 6.07796 ± 0.37080 | 1.12654 ± 0.09609 | 1.58039 ± 0.32995 | 1.27445 |
| Tange Fit3-Vinet MgO pressures | 22.58426 ± 0.09239 | 149.777 ± 8.626 | 5.82394 ± 0.26808 | 1.00371 ± 0.06204 | 0.71483 ± 0.23857 | 0.97879 |
| Speziale variable-q Debye MgO reconstruction | 22.70980 ± 0.09880 | 138.716 ± 8.367 | 5.89337 ± 0.27576 | 0.97443 ± 0.05756 | 0.50272 ± 0.22274 | 0.92278 |

## Effect of changing n only

| Pressure input | Previous n=1 gamma0 | Author-reported n=2 gamma0 | Previous n=1 q | n=2 q |
|---|---:|---:|---:|---:|
| Supplementary printed pressures | 2.23572 | 1.12654 | 1.51423 | 1.58039 |
| Tange Fit3-Vinet MgO pressures | 1.98945 | 1.00371 | 0.64308 | 0.71483 |
| Speziale variable-q Debye MgO reconstruction | 1.93069 | 0.97443 | 0.42959 | 0.50272 |

## Published coefficients evaluated with n=1 versus n=2

These are residuals of the unchanged published coefficients, not refitted residuals.

| Pressure input | RMS with n=1 (GPa) | RMS with confirmed n=2 (GPa) |
|---|---:|---:|
| Supplementary printed pressures | 6.98261 | 2.49886 |
| Tange Fit3-Vinet MgO pressures | 7.31048 | 1.80528 |
| Speziale variable-q Debye MgO reconstruction | 5.45224 | 2.44770 |

## Normalization control and interpretation

A separate official-console replay sets n=2 and doubles every molar volume, volume error and initial V0 to the convention of moles of two-atom hcp cells (initial V0 approximately 13.74 cm3/mol). This reproduces the previous n=1 optimum within convergence and printing precision. Algebraically it leaves both gamma/V times Debye energy and V/V0 unchanged, so the two consistent normalizations describe the same full-MGD pressures. The user explicitly confirmed that the author used 6.87, not this doubled-volume control.

The reported n setting explains the dominant near-factor-of-two difference in fitted gamma. The complete published fit is still not reproduced: pressure calibration, weighting and remaining K0/q differences persist. For example, the supplied pressure column gives gamma0 close to 1.11 but q=1.58039 rather than 0.3; the Speziale reconstruction gives gamma0=0.97443 and q=0.50272. The n=2 author setup with volume per mole of Fe differs from the standard per-Fe normalization. Recovering the authors calculation and establishing physical normalization are distinct questions. No published record, default or selectable independent Tange refit was silently changed.

Evidence: [author-setting report](report.json), [independent pressure check](pressure-check.json), [doubled-volume control](../miozzi-2020-eosfit-n2-double-volume-control/report.json). All 18 fits converged. Independent pressures match the console within 0.0025 GPa; independent five-parameter full-MGD fits verify the n=2 optimum and conditional covariance.

## GUI input and replay

For the author-reported setup, use the previously exported text files unchanged and set full MGD, n=2, theta0=420 K and reference T=300 K. Initial V0 is approximately 6.87 cm3/mol. Do not double the exported volumes for this author-setting replay. The earlier n=1 file comments/README document the previous per-Fe convention.

```sh
python -m scripts.replay_miozzi_2020_staged --executable /path/to/eosfit7c --output /new/output/author-n2 --atoms 2
python -m scripts.replay_miozzi_2020_staged --executable /path/to/eosfit7c --output /new/output/n2-volume-control --atoms 2 --molar-volume-multiplier 2
```
