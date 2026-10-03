# Whole-hcp-cell MGD fit with Peritheos

36 RT rows (298/300 K); all 131 rows with RT coefficients fixed; all 131 rows with V0,K0,K0_prime,gamma0,q free. theta0=420 K and Tr=300 K fixed throughout thermal stages.

Equal pressure-residual weights; central coordinates held fixed; linear least squares. Bounds are identical after conversion to physical cell volumes.

The n=1 fit uses volume per mole of Fe; the n=2 fit uses volume per mole of whole two-atom cells. Both start from 22.81 Å³/cell, respectively 6.86825 and 13.73650 cm³/mol (0.686825 and 1.373650 J/bar/mol). The fitted volume below is converted back to the same physical cell basis.

| Pressure input | n | V0 (Å³/cell) | K0 (GPa) | K0′ | gamma0 | q | RMS (GPa) |
|---|---:|---:|---:|---:|---:|---:|---:|
| printed | 1 | 22.602632 ± 0.123290 | 146.806468 ± 11.350304 | 6.053557 ± 0.366965 | 2.235644 ± 0.190770 | 1.514293 ± 0.328456 | 1.27798169 |
| printed | 2 | 22.602632 ± 0.123290 | 146.806469 ± 11.350304 | 6.053557 ± 0.366965 | 2.235644 ± 0.190770 | 1.514293 ± 0.328456 | 1.27798169 |
| tange_vinet | 1 | 22.578319 ± 0.092116 | 150.484856 ± 8.628765 | 5.801501 ± 0.266405 | 1.989345 ± 0.123305 | 0.643085 ± 0.237983 | 0.98446756 |
| tange_vinet | 2 | 22.578319 ± 0.092116 | 150.484853 ± 8.628764 | 5.801501 ± 0.266405 | 1.989345 ± 0.123305 | 0.643085 ± 0.237983 | 0.98446756 |
| speziale_variable_q_debye | 1 | 22.703126 ± 0.098420 | 139.419575 ± 8.365237 | 5.869960 ± 0.273780 | 1.930567 ± 0.114262 | 0.429554 ± 0.221943 | 0.92812717 |
| speziale_variable_q_debye | 2 | 22.703126 ± 0.098420 | 139.419572 ± 8.365237 | 5.869960 ± 0.273780 | 1.930567 ± 0.114262 | 0.429554 ± 0.221943 | 0.92812717 |

## Verification

| Pressure input | Independently refitted n=1 versus n=2: maximum ΔP (GPa) | n=2 versus direct cell/kB calculation: maximum ΔP (GPa) |
|---|---:|---:|
| printed | 6.35e-07 | 4.21e-10 |
| tange_vinet | 2.23e-07 | 4.37e-10 |
| speziale_variable_q_debye | 1.94e-07 | 4.38e-10 |

All 18 stages converged. The same coefficients under matching atom-count/volume scaling also give identical pressures within numerical precision, separately from optimizer stopping differences.

Conditional local standard errors and covariance scaled by RSS/degrees of freedom. Calibration, fixed-theta, predictor and inter-row covariance uncertainties excluded. Stage-two errors condition on the cold coefficients.

Independent Peritheos normalization experiment. Whole-cell n=2 also doubles the molar volumes; this differs from the user-reported author setting n=2 with V0 about 6.87 cm3/mol. It does not resolve original author inputs or physically validate a pressure scale. Library records and defaults are unchanged.

Complete stage results, residuals, covariance, normalization and input fingerprints are retained in [report.json](report.json).

Reproduce from the repository root:

```sh
.venv/bin/python -m scripts.refit_miozzi_2020_cell_basis
```
