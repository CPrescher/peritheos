# Miozzi q-compromise MGD diagnostic — 2026-10-03

The EosFit 7.6 bundled manual (`EoSFitanddisplay.html`) says the GUI estimates MGD with q-compromise for stability and permits switching to the general model after refinement. This describes the GUI estimation workflow, not the authors' final selected model. The current documentation and source identify q-compromise as a v7.6 addition; this is not evidence that the January 2020 paper used it. The source paper reports a fitted q=0.3(3).

Q-compromise uses theta(V)=theta0 and gamma(V)/V=gamma0/V0. Therefore Pthermal depends only on T: gamma0/V0_molar times the Debye-energy difference from 300 K. It has no adjustable q. This differs from full MGD with q fixed to either zero or one.

All runs use the unmodified official EosFit7c 7.60 executable. Stage one fits V0,K0,Kprime to 36 RT rows; stage two fixes the cold curve and fits gamma0 only; stage three releases V0,K0,Kprime,gamma0. Theta0=420 K, n=1 and Tr=300 K remain fixed. All 131 rows receive equal pressure-residual weights. The alternative switches to full MGD before stage three and releases q, initially set to 0.3, retaining the preceding cold curve and gamma0 in memory. Original author weights/pressure column remain unresolved. No library coefficients or defaults were changed.

## Final q-compromise coefficients

Errors below are conditional local EosFit standard errors. They exclude pressure-scale uncertainty and fixed-parameter uncertainty. V0 refers to the two-atom hcp Fe cell. Full covariance is retained in each report.

| Pressure input | V0 (Å³/cell) | K0 (GPa) | Kprime | gamma0 | q | RMS pressure residual (GPa) |
|---|---:|---:|---:|---:|---|---:|
| Published preferred full-MGD model | 22.81 | 129 ± 1 | 6.24 ± 0.04 | 1.11 ± 0.01 | 0.3 ± 0.3 | Not reported |
| Supplementary printed pressures | 22.58353 ± 0.11879 | 151.288 ± 10.847 | 5.82850 ± 0.31650 | 1.93913 ± 0.04477 | Absent | 1.28915 |
| Tange Fit3-Vinet MgO pressures | 22.60233 ± 0.09375 | 146.697 ± 8.346 | 5.95921 ± 0.25453 | 2.14302 ± 0.03465 | Absent | 0.97958 |
| Speziale variable-q Debye MgO reconstruction | 22.74796 ± 0.10415 | 133.059 ± 8.280 | 6.16139 ± 0.27767 | 2.20012 ± 0.03422 | Absent | 0.93549 |

## Full MGD compared with q-compromise initialization

| Pressure input | Starting procedure | Final K0 (GPa) | Final gamma0 | Final q | RMS (GPa) |
|---|---|---:|---:|---:|---:|
| Supplementary printed pressures | Full MGD throughout | 146.81015 | 2.23572 | 1.51423 | 1.277982 |
| Supplementary printed pressures | Q-compromise stage 2, then full MGD | 146.80460 | 2.23582 | 1.51438 | 1.277982 |
| Tange Fit3-Vinet MgO pressures | Full MGD throughout | 150.48611 | 1.98945 | 0.64308 | 0.984468 |
| Tange Fit3-Vinet MgO pressures | Q-compromise stage 2, then full MGD | 150.48717 | 1.98942 | 0.64303 | 0.984468 |
| Speziale variable-q Debye MgO reconstruction | Full MGD throughout | 139.41940 | 1.93069 | 0.42959 | 0.928127 |
| Speziale variable-q Debye MgO reconstruction | Q-compromise stage 2, then full MGD | 139.41826 | 1.93066 | 0.42951 | 0.928128 |

## Stage-two thermal coefficients with cold curve fixed

| Pressure input | Fixed V0 (Å³/cell) | Fixed K0 (GPa) | Fixed Kprime | Q-compromise gamma0 | RMS (GPa) |
|---|---:|---:|---:|---:|---:|
| Supplementary printed pressures | 22.69057 | 136.16335 | 6.46487 | 1.98449 ± 0.02018 | 1.33861 |
| Tange Fit3-Vinet MgO pressures | 22.69807 | 136.57272 | 6.30804 | 2.19819 ± 0.01496 | 0.99182 |
| Speziale variable-q Debye MgO reconstruction | 22.80628 | 128.50888 | 6.28689 | 2.24266 ± 0.01434 | 0.94667 |

## Interpretation and verification

All 18 fits in the two new three-case replays converged. The q-compromise approximation does not recover gamma0=1.11. Switching from q-compromise to full MGD returns the previous full-MGD optimum within console convergence/printing precision. Thus this GUI starting workflow does not resolve the discrepancy for the tested data and equal weights. This does not prove the authors used these inputs, weights, model, or software version.

Independent Python pressure calculations agree with all console-calculated pressures within 0.0025 GPa, including the different thermal laws. The finite precision of printed coefficients and the gas constants (EosFit R=8.314 versus Python R=8.314462618) account for that tolerance. Independent Python four-parameter least-squares fits reproduce q-compromise RMS within 1e-6 GPa and coefficients within 1e-4 relative tolerance. Source hashes, saved model flags, constraints and conditional covariance are checked in `tests/test_miozzi_2020_staged.py`.

Evidence: [q-compromise report](report.json), [independent four-parameter fit check](independent-check.json), [q-compromise-to-full report](../miozzi-2020-eosfit-q-compromise-start/report.json), [previous full-MGD report](../miozzi-2020-eosfit-staged/report.json). Each case contains console inputs/macros/logs and all three saved EOS files.

Replay commands (require a working display/runtime for the external executable):

```sh
python -m scripts.replay_miozzi_2020_staged --executable /path/to/eosfit7c --output /new/output/q-compromise --q-compromise
python -m scripts.replay_miozzi_2020_staged --executable /path/to/eosfit7c --output /new/output/q-compromise-start --q-compromise-start
```
