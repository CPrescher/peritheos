# Reconstruction uncertainties and conditional numerical parity

This bounded follow-up was completed on 25 October 2026. It reconstructs
missing uncertainty metadata and reviews existing numerical qualifications.
It does not create replacement EOS coefficients.

The screen identified 51 records with missing reconstruction errors in the
common ledger, including Finger's specialized argon assessment, and 12 further
records whose existing errors or numerical flags required a scope review.
All **63 records** have been reviewed and achieve **conditional numerical
parameter parity** within their explicitly documented comparison scope. Conditional standard
errors are now available for all 51 missing tables. Ma (2024) KCl already had
diagnostic covariance; its errors were recovered and converted, not refitted.

Published coefficients, printed source errors and confidence conventions,
fixed/adopted values, observations, pressure reductions, executable equations,
default selections and recommendation holds are unchanged. Their complete
record fields and dataset declarations are protected by semantic checksums.
Earlier audit wording is retained explicitly as historical evidence; its
scientific qualifications continue to apply.

The subsequent [remaining-similar scope review](similar-parameter-review.md)
checks the other 96 similar records. Its findings do not alter this bounded
63-record reconstruction or its uncertainty values.

## Numerical criterion and statistical meaning

For 61 promoted records, every compared free central coefficient lies inside its own
positive printed source parameter error width, with finite conditional
reconstruction standard errors. The explicit comparison basis is
`within_reported_parameter_errors`. It does not interpret a published 95%
interval as one sigma, infer missing source confidence, or establish a new
formal combined-two-sigma test. Original author weighting, confidence,
calibration covariance and full workflow remain distinct questions.

The remaining two promotions, Chen (2018) tetragonal Ca-perovskite and Litasov
(2007) HTBM/T03 superhydrous phase B, meet the ledger's combined-error comparison
and numerical similarity screen. The initial follow-up unnecessarily required
them also to lie within each printed source width. No source-specific numerical
reason was identified for imposing that additional cutoff on parameter agreement.
Chen's source parameter widths are explicitly two-sigma: dividing them by two
before combining with the reconstruction standard errors still passes for both
coefficients. For Litasov, source confidence is unspecified; the retained
combined-two-width comparison is a numerical diagnostic, not a formal two-sigma
test. Neither result recovers original weighting or author workflow. Both fits
reuse source observations, and their cross-fit covariance is unknown; the
quadrature rule must not be presented as a formal test of independent estimates.

Farla (2025) RbCl meets the printed-source-width parameter criterion for all four
varied coefficients: the largest difference is 0.40 of the source error width.
Its parameter agreement is conditional on fixed V0/theta0, the documented
analyst selection, source-reported CsCl pressures and frozen diagonal effective
P/V/T weights. The unresolved CsCl calibration replay, RbCl benchmark entry and
original author-workflow findings remain separate qualifications. They do not
cap this assessment's explicitly limited `parameter_agreement` outcome.

This numerical parameter assessment does not establish experimental accuracy,
phase stability, recommended use or a complete uncertainty budget for pressure
calibration. Fixed parameters remain constraints, with their errors excluded
from these conditional calculations unless explicitly propagated below.

## Covariance reconstruction

Ordinary least-squares reconstructions use

```text
s² = sum(residual²)/(n-k)
C = s² (JᵀJ)^-1
```

The original recorded residual objective, selected observations, free mask and
central estimates are preserved. The Jacobian uses five-point central
differences in scaled parameter coordinates and is checked at half the step.
Rank-deficient parameter combinations are rejected rather than assigned
pseudoinverse errors. The row count and descriptive pressure RMS must replay
before covariance is calculated. Marginal standard errors are sqrt(diag(C));
the full ordered covariance and correlation matrices are retained.

Unweighted pressure fits condition on observed volumes. For the weighted
reconstructions, the original weights are held fixed:

- Ono (2020) argon uses the source audit's effective P/V widths evaluated at
  published coefficients.
- Chidester (2018) ThO2 uses the earlier frozen effective P/V/T widths,
  including its explicit Ar 3% cold-pressure sensitivity assumption.
- Farla (2025) RbCl uses the final frozen effective P/V/T widths at the
  existing converged reconstruction, with no pressure floor. Its unresolved
  calibration and benchmark findings remain in force.
- Diamond retains the earlier equal-volume/run-2 diagnostic on digitized
  observations, with residuals in atomic A³. Its conditional covariance does
  not estimate graphical or pressure-scale systematics.

These are diagonal reconstruction error models. Residual scaling estimates
scatter under those models; it does not recover missing calibration covariance
or turn graphical widths into measured Gaussian errors.

### Komabayashi: shared observations in two fitting stages

Post-perovskite first fits K0 to three 300 K rows, then fits alpha1 and dK/dT
to all 21 rows with that K0 held at its reconstructed value. The three cold
rows also occur in the second stage. Treating the stages as independent would
lose their covariance.

The linear influence of the cold-stage estimator is propagated into the
thermal estimator on the same pressure observations. If A maps observation
perturbations to all three fitted coefficients and J is their prediction
Jacobian, the reconstruction uses

```text
C = s² A Aᵀ
s² = sum(pressure residual²) / trace[(I-JA)(I-JA)ᵀ]
```

The effective residual degrees of freedom are approximately 18.7705. This
assumes a common independent pressure-residual variance across the original
rows and uses a local linearization of the staged procedure. V0, K0-prime
and alpha(300) remain fixed. The alpha0/alpha1 tie is maintained when taking
derivatives. The result is not the covariance of a simultaneous three-parameter
fit or a claimed reconstruction of the authors' uncertainty calculation.

### Recovered and previously available covariance

Ma's existing joint acoustic/Walker covariance is converted with
`C_cell = D C_molar Dᵀ`, where D scales V0 from cm³/mol to conventional-cell
A³ and leaves K0 and K0-prime unchanged. Both covariance axes and the standard
error receive the corresponding conversion. The same 11 acoustic and eight
Walker observations and relative-error objective are retained. This assesses
the three BM3 coefficients, not the complete thermal/calibration chain.

The 12 other reviewed records already had reconstructed marginal errors.
Those errors are retained exactly; absent off-diagonal covariance is not
invented. Reconstructing missing off-diagonal elements was not used as a
condition for accepting or rejecting their existing coefficient comparisons.

Where provided, 95% intervals use the Student-t multiplier with the recorded
residual degrees of freedom. They are conditional, marginal, local
approximations, not exact nonlinear or joint confidence regions. Source error
widths retain their original convention, including Dewaele's 95% intervals.

## Parameter agreement and separate reproduction qualifications

The final three records also meet the printed-source-width criterion for their
two varied coefficients. Their earlier similarity caps mixed numerical parameter
agreement with different reproduction questions. No fitting inputs, constraints,
central estimates or uncertainty values were changed to obtain these outcomes.

| Record | Conditional parameter result | Retained limitation |
| --- | --- | --- |
| Datchi (2007), diamond | K0 and K0-prime agree within 0.45 source error widths. | Conditional on the source-unconfirmed run-2/equal-volume hypothesis, 14/24 digitized markers and fixed source V0. Agreement helped select the hypothesis, so it does not independently validate the selection. Graphical/calibration uncertainty and author workflow remain unrecovered. |
| Fei (2004), preliminary 300 K Pt | K0 and K0-prime agree within 0.67 source error widths. | Conditional on diagnostic V0=60.38 A³ fixed and equal pressure-residual weights for seven cold rows. The source fixed/free mask is unrecovered. V0 is not independently recovered; the derivative is weakly constrained and this remains the preliminary, nonpreferred Pt scale. |
| Litasov (2013), siderite | K0 and K0-prime agree within 0.23 source error widths. | All 27 cold rows, source-fixed V0=293.4 A³, printed Au pressures and equal pressure-residual weights. Published-curve RMS is still 0.22958 GPa versus reported 0.11–0.17 GPa. Pressure reduction, source-fit statistic and thermal coefficient conflicts remain separately unresolved. |

The follow-up does not enlarge or reinterpret source errors to promote these
records. This is parameter agreement under named conditions, not a claim that
every EOS is validated, recommended, or independently reproduced under the
authors' original protocol. A finite error estimate does not remove another
reproduction problem; all source findings and use restrictions remain.

## Full results

<!-- uncertainty-followup-table -->

| Record | Outcome | Reconstructed free coefficients ± conditional standard error |
| --- | --- | --- |
| [akimotoite_wang_2004_t0133_bm3](../reconstruction-reports/akimotoite.md#record-akimotoite_wang_2004_t0133_bm3) | parity | V0: 264.32 ± 0.14; K0_prime: 4.5 ± 0.4 |
| [akimotoite_zhou_2014_300k_bm3](../reconstruction-reports/akimotoite.md#record-akimotoite_zhou_2014_300k_bm3) | parity | V0: 262.4 ± 0.4; K0: 207 ± 5 |
| [aluminum_dewaele_2019_mao_vinet](../reconstruction-reports/aluminum.md#record-aluminum_dewaele_2019_mao_vinet) | parity | V0: 66.29 ± 0.15; K0: 76.3 ± 1.4; K0_prime: 4.16 ± 0.04 |
| [aluminum_dewaele_2019_dor_vinet](../reconstruction-reports/aluminum.md#record-aluminum_dewaele_2019_dor_vinet) | parity | V0: 66.34 ± 0.16; K0: 74.2 ± 1.4; K0_prime: 4.52 ± 0.05 |
| [argon_fcc_ono_2020_vinet](../reconstruction-reports/argon_fcc.md#record-argon_fcc_ono_2020_vinet) | parity | V0: 210 ± 40; K0: 0.5 ± 0.5; K0_prime: 8.7 ± 0.9 |
| [argon_fcc_finger_1981_murnaghan2_debye](../reconstruction-reports/argon_fcc.md#record-argon_fcc_finger_1981_murnaghan2_debye) | parity | K0_prime: 6.90 ± 0.10; K0_double_prime: -0.34 ± 0.09 |
| [ca_perovskite_tetragonal_chen_2018_vinet](../reconstruction-reports/ca_perovskite_tetragonal.md#record-ca_perovskite_tetragonal_chen_2018_vinet) | parity | V0: 186.0 ± 1.1; K0: 215 ± 11 |
| [cobalt_hcp_dewaele_2019_mao_vinet](../reconstruction-reports/cobalt_hcp.md#record-cobalt_hcp_dewaele_2019_mao_vinet) | parity | V0: 22.154 ± 0.014; K0: 197 ± 3; K0_prime: 3.85 ± 0.13 |
| [cobalt_hcp_dewaele_2019_dor_vinet](../reconstruction-reports/cobalt_hcp.md#record-cobalt_hcp_dewaele_2019_dor_vinet) | parity | V0: 22.154 ± 0.015; K0: 195 ± 3; K0_prime: 4.36 ± 0.14 |
| [copper_dewaele_2019_mao_vinet](../reconstruction-reports/copper.md#record-copper_dewaele_2019_mao_vinet) | parity | V0: 47.24 ± 0.05; K0: 135.3 ± 1.8; K0_prime: 4.91 ± 0.05 |
| [copper_dewaele_2019_dor_vinet](../reconstruction-reports/copper.md#record-copper_dewaele_2019_dor_vinet) | parity | V0: 47.24 ± 0.05; K0: 133.1 ± 1.8; K0_prime: 5.38 ± 0.05 |
| [diamond_datchi_2007_vinet_1](../reconstruction-reports/diamond.md#record-diamond_datchi_2007_vinet_1) | parity | K0: 444.3 ± 1.7; K0_prime: 3.96 ± 0.04 |
| [fe09ni01_hcp_sakai_2014_p5_vinet](../reconstruction-reports/fe09ni01_hcp.md#record-fe09ni01_hcp_sakai_2014_p5_vinet) | parity | K0: 206.3 ± 2.2; K0_prime: 4.82 ± 0.06 |
| [gold_dewaele_2019_mao_vinet](../reconstruction-reports/gold.md#record-gold_dewaele_2019_mao_vinet) | parity | V0: 67.93 ± 0.09; K0: 166 ± 3; K0_prime: 5.47 ± 0.10 |
| [gold_dewaele_2019_dor_vinet](../reconstruction-reports/gold.md#record-gold_dewaele_2019_dor_vinet) | parity | V0: 67.94 ± 0.10; K0: 163 ± 4; K0_prime: 6.03 ± 0.11 |
| [gold_hirose_2008_bm3_fit2](../reconstruction-reports/gold.md#record-gold_hirose_2008_bm3_fit2) | parity | beta1: 10e-7 ± 4e-7; beta2: 4e-10 ± 4e-10; beta3: 15e-14 ± 10e-14; kprime_log_coefficient: 35.9e-5 ± 2.3e-5 |
| [gold_hirose_2008_bm3_fit1](../reconstruction-reports/gold.md#record-gold_hirose_2008_bm3_fit1) | parity | dK_dT: -0.0283 ± 0.0027; alpha0: 31.8e-6 ± 1.4e-6; alpha1: 15e-9 ± 3e-9 |
| [gold_hirose_2008_bm3_300k](../reconstruction-reports/gold.md#record-gold_hirose_2008_bm3_300k) | parity | K0_prime: 5.583 ± 0.023 |
| [iron_dewaele_2019_mao_vinet](../reconstruction-reports/iron.md#record-iron_dewaele_2019_mao_vinet) | parity | V0: 22.42 ± 0.07; K0: 165 ± 5; K0_prime: 4.96 ± 0.10 |
| [iron_dewaele_2019_dor_vinet](../reconstruction-reports/iron.md#record-iron_dewaele_2019_dor_vinet) | parity | V0: 22.35 ± 0.07; K0: 169 ± 6; K0_prime: 5.33 ± 0.10 |
| [iron_sakai_2014_p1_vinet](../reconstruction-reports/iron.md#record-iron_sakai_2014_p1_vinet) | parity | K0: 176.3 ± 1.9; K0_prime: 5.37 ± 0.05 |
| [iron_sakai_2014_p4_vinet](../reconstruction-reports/iron.md#record-iron_sakai_2014_p4_vinet) | parity | K0: 190.7 ± 2.2; K0_prime: 5.43 ± 0.06 |
| [iron_sakai_2014_p5_bm3](../reconstruction-reports/iron.md#record-iron_sakai_2014_p5_bm3) | parity | K0: 198.4 ± 2.2; K0_prime: 4.90 ± 0.06 |
| [iron_sakai_2014_p5_vinet](../reconstruction-reports/iron.md#record-iron_sakai_2014_p5_vinet) | parity | K0: 190.5 ± 1.8; K0_prime: 5.45 ± 0.05 |
| [iron_sakai_2014_p6_vinet](../reconstruction-reports/iron.md#record-iron_sakai_2014_p6_vinet) | parity | K0: 168.7 ± 1.9; K0_prime: 5.93 ± 0.06 |
| [iron_dewaele_2006_bm3](../reconstruction-reports/iron.md#record-iron_dewaele_2006_bm3) | parity | V0: 22.461 ± 0.013; K0_prime: 4.982 ± 0.021 |
| [kcl_b2_tateno_2019_holmes_vinet_mgd](../reconstruction-reports/kcl.md#record-kcl_b2_tateno_2019_holmes_vinet_mgd) | parity | K0: 17.49 ± 0.25; K0_prime: 5.76 ± 0.04; gamma0: 1.93 ± 0.21; q: 0.73 ± 0.25 |
| [kcl_b2_tateno_2019_holmes_vinet_linear_thermal](../reconstruction-reports/kcl.md#record-kcl_b2_tateno_2019_holmes_vinet_linear_thermal) | parity | K0: 17.55 ± 0.24; K0_prime: 5.75 ± 0.03; alpha_KT: 0.00321 ± 0.00014 |
| [kcl_b2_tateno_2019_sokolova_vinet_linear_thermal](../reconstruction-reports/kcl.md#record-kcl_b2_tateno_2019_sokolova_vinet_linear_thermal) | parity | K0: 18.34 ± 0.25; K0_prime: 5.60 ± 0.03; alpha_KT: 0.00368 ± 0.00014 |
| [kcl_b2_chidester_2021_vinet_mgd](../reconstruction-reports/kcl.md#record-kcl_b2_chidester_2021_vinet_mgd) | parity | V0: 56.8 ± 0.8; K0: 13.1 ± 1.2; K0_prime: 6.21 ± 0.10; gamma0: 3.4 ± 0.4; q: 0.99 ± 0.13 |
| [kcl_b2_ma_2024_bm3_mgd](../reconstruction-reports/kcl.md#record-kcl_b2_ma_2024_bm3_mgd) | parity | V0: 53.93 ± 0.09; K0: 21.3 ± 0.4; K0_prime: 4.85 ± 0.05 |
| [mgsio3_post_perovskite_komabayashi_2008_thermal_bm3](../reconstruction-reports/mgsio3_post_perovskite.md#record-mgsio3_post_perovskite_komabayashi_2008_thermal_bm3) | parity | K0: 223.2 ± 0.6; alpha1: 1.1e-9 ± 2.0e-9; dK_dT: -0.0084 ± 0.0022 |
| [molybdenum_dewaele_2019_mao_vinet](../reconstruction-reports/molybdenum.md#record-molybdenum_dewaele_2019_mao_vinet) | parity | V0: 31.138 ± 0.018; K0: 270 ± 3; K0_prime: 3.33 ± 0.07 |
| [molybdenum_dewaele_2019_dor_vinet](../reconstruction-reports/molybdenum.md#record-molybdenum_dewaele_2019_dor_vinet) | parity | V0: 31.133 ± 0.019; K0: 269 ± 3; K0_prime: 3.87 ± 0.07 |
| [nickel_dewaele_2008_vinet_1](../reconstruction-reports/nickel.md#record-nickel_dewaele_2008_vinet_1) | parity | V0: 43.83 ± 0.06; K0: 176 ± 3; K0_prime: 4.87 ± 0.08 |
| [nickel_dewaele_2019_mao_vinet](../reconstruction-reports/nickel.md#record-nickel_dewaele_2019_mao_vinet) | parity | V0: 43.83 ± 0.06; K0: 176 ± 3; K0_prime: 4.87 ± 0.08 |
| [nickel_dewaele_2019_dor_vinet](../reconstruction-reports/nickel.md#record-nickel_dewaele_2019_dor_vinet) | parity | V0: 43.82 ± 0.06; K0: 175 ± 3; K0_prime: 5.37 ± 0.08 |
| [platinum_dewaele_2019_mao_vinet](../reconstruction-reports/platinum.md#record-platinum_dewaele_2019_mao_vinet) | parity | V0: 60.396 ± 0.020; K0: 273.4 ± 2.2; K0_prime: 4.83 ± 0.07 |
| [platinum_dewaele_2019_dor_vinet](../reconstruction-reports/platinum.md#record-platinum_dewaele_2019_dor_vinet) | parity | V0: 60.394 ± 0.022; K0: 270.8 ± 2.4; K0_prime: 5.50 ± 0.07 |
| [platinum_fei_2004_bm3_300k](../reconstruction-reports/platinum.md#record-platinum_fei_2004_bm3_300k) | parity | K0: 294 ± 6; K0_prime: 2.1 ± 0.7 |
| [platinum_fei_2007_vinet_300k](../reconstruction-reports/platinum.md#record-platinum_fei_2007_vinet_300k) | parity | K0_prime: 5.08 ± 0.05 |
| [rbcl_b2_farla_2025_bm3_mgd_cscl](../reconstruction-reports/rbcl.md#record-rbcl_b2_farla_2025_bm3_mgd_cscl) | parity | K0: 19.92 ± 0.09; K0_prime: 4.995 ± 0.024; gamma0: 1.96 ± 0.04; q: 1.06 ± 0.09 |
| [rhenium_dewaele_2019_mao_vinet](../reconstruction-reports/rhenium.md#record-rhenium_dewaele_2019_mao_vinet) | parity | V0: 29.477 ± 0.014; K0: 349 ± 4; K0_prime: 4.01 ± 0.08 |
| [rhenium_dewaele_2019_dor_vinet](../reconstruction-reports/rhenium.md#record-rhenium_dewaele_2019_dor_vinet) | parity | V0: 29.471 ± 0.015; K0: 349 ± 4; K0_prime: 4.65 ± 0.08 |
| [siderite_fe095mn005_litasov_2013_bm3](../reconstruction-reports/siderite_fe095mn005.md#record-siderite_fe095mn005_litasov_2013_bm3) | parity | K0: 120.2 ± 1.3; K0_prime: 3.57 ± 0.13 |
| [silver_dewaele_2008_vinet_1](../reconstruction-reports/silver.md#record-silver_dewaele_2008_vinet_1) | parity | V0: 68.28 ± 0.17; K0: 100.4 ± 3.0; K0_prime: 5.69 ± 0.10 |
| [silver_dewaele_2019_mao_vinet](../reconstruction-reports/silver.md#record-silver_dewaele_2019_mao_vinet) | parity | V0: 68.28 ± 0.17; K0: 100.4 ± 3.0; K0_prime: 5.69 ± 0.10 |
| [silver_dewaele_2019_dor_vinet](../reconstruction-reports/silver.md#record-silver_dewaele_2019_dor_vinet) | parity | V0: 68.35 ± 0.19; K0: 97 ± 3; K0_prime: 6.21 ± 0.11 |
| [superhydrous_phase_b_lt_litasov_2007_bm3_a89](../reconstruction-reports/superhydrous_phase_b_lt.md#record-superhydrous_phase_b_lt_litasov_2007_bm3_a89) | parity | V0: 623.31 ± 0.19; K0: 138.4 ± 1.5; K0_prime: 4.95 ± 0.16 |
| [superhydrous_phase_b_lt_litasov_2007_bm3_t03](../reconstruction-reports/superhydrous_phase_b_lt.md#record-superhydrous_phase_b_lt_litasov_2007_bm3_t03) | parity | V0: 623.30 ± 0.17; K0: 134.2 ± 1.3; K0_prime: 6.04 ± 0.15 |
| [superhydrous_phase_b_lt_litasov_2007_htbm_a89](../reconstruction-reports/superhydrous_phase_b_lt.md#record-superhydrous_phase_b_lt_litasov_2007_htbm_a89) | parity | V0: 623.46 ± 0.26; K0: 135.8 ± 1.9; K0_prime: 5.25 ± 0.19; alpha0: 32.0e-6 ± 1.1e-6; alpha1: 11.5e-9 ± 2.6e-9; dK_dT: -0.0256 ± 0.0026 |
| [superhydrous_phase_b_lt_litasov_2007_htbm_t03](../reconstruction-reports/superhydrous_phase_b_lt.md#record-superhydrous_phase_b_lt_litasov_2007_htbm_t03) | parity | V0: 623.4 ± 0.3; K0: 131.9 ± 2.3; K0_prime: 6.32 ± 0.26; alpha0: 34.3e-6 ± 1.2e-6; alpha1: 14e-9 ± 3e-9; dK_dT: -0.028 ± 0.003 |
| [superhydrous_phase_b_lt_litasov_2007_thermal_pressure_a89](../reconstruction-reports/superhydrous_phase_b_lt.md#record-superhydrous_phase_b_lt_litasov_2007_thermal_pressure_a89) | parity | V0: 623.42 ± 0.26; K0: 135.7 ± 1.8; K0_prime: 5.24 ± 0.19; alpha0: 37.7e-6 ± 1.8e-6; dK_dT_V: -0.0023 ± 0.0021 |
| [superhydrous_phase_b_lt_litasov_2007_mgd_free_debye_a89](../reconstruction-reports/superhydrous_phase_b_lt.md#record-superhydrous_phase_b_lt_litasov_2007_mgd_free_debye_a89) | parity | V0: 623.44 ± 0.25; K0: 136.0 ± 1.8; K0_prime: 5.24 ± 0.18; gamma0: 1.19 ± 0.06; q: 1.8 ± 0.4; theta0: 550 ± 130 |
| [superhydrous_phase_b_lt_litasov_2007_mgd_elastic_debye_a89](../reconstruction-reports/superhydrous_phase_b_lt.md#record-superhydrous_phase_b_lt_litasov_2007_mgd_elastic_debye_a89) | parity | V0: 623.46 ± 0.26; K0: 136.4 ± 1.9; K0_prime: 5.23 ± 0.19; gamma0: 1.30 ± 0.06; q: 2.1 ± 0.5 |
| [tantalum_dewaele_2019_mao_vinet](../reconstruction-reports/tantalum.md#record-tantalum_dewaele_2019_mao_vinet) | parity | V0: 36.040 ± 0.021; K0: 197.8 ± 2.0; K0_prime: 3.17 ± 0.07 |
| [tantalum_dewaele_2019_dor_vinet](../reconstruction-reports/tantalum.md#record-tantalum_dewaele_2019_dor_vinet) | parity | V0: 36.038 ± 0.023; K0: 196.1 ± 2.2; K0_prime: 3.64 ± 0.07 |
| [tho2_cotunnite_chidester_2018_bm3_linear_thermal](../reconstruction-reports/tho2_cotunnite.md#record-tho2_cotunnite_chidester_2018_bm3_linear_thermal) | parity | V0: 164.5 ± 0.4; K0: 191.0 ± 2.9; alpha_KT: 0.0035 ± 0.0003 |
| [thorianite_chidester_2018_bm3_linear_thermal](../reconstruction-reports/thorianite.md#record-thorianite_chidester_2018_bm3_linear_thermal) | parity | K0: 203.1 ± 2.9; alpha_KT: 0.00366 ± 0.00028 |
| [tungsten_dewaele_2019_mao_vinet](../reconstruction-reports/tungsten.md#record-tungsten_dewaele_2019_mao_vinet) | parity | V0: 31.724 ± 0.024; K0: 298 ± 4; K0_prime: 3.82 ± 0.08 |
| [tungsten_dewaele_2019_dor_vinet](../reconstruction-reports/tungsten.md#record-tungsten_dewaele_2019_dor_vinet) | parity | V0: 31.715 ± 0.024; K0: 299 ± 4; K0_prime: 4.37 ± 0.08 |
| [zinc_hcp_dewaele_2019_mao_vinet](../reconstruction-reports/zinc_hcp.md#record-zinc_hcp_dewaele_2019_mao_vinet) | parity | V0: 30.29 ± 0.10; K0: 64.3 ± 1.9; K0_prime: 5.30 ± 0.07 |
| [zinc_hcp_dewaele_2019_dor_vinet](../reconstruction-reports/zinc_hcp.md#record-zinc_hcp_dewaele_2019_dor_vinet) | parity | V0: 30.31 ± 0.11; K0: 62.3 ± 2.1; K0_prime: 5.70 ± 0.08 |

<!-- /uncertainty-followup-table -->

## Evidence and reproducibility

The [historical comparison inputs](../data/uncertainty-reconstruction-inputs.json)
preserve the selected original ledger findings and their complete qualifications.
The [numerical follow-up](../data/uncertainty-reconstruction.json) records each
decision, conditional errors, ordered covariance, weighting, degrees of freedom,
derivative checks and prior assessment. Its SHA-256 manifest binds the input
evidence, linked data resources, pressure-calibration registry and the existing
Ma covariance audit. Primary references remain attached to each original EOS
and source dataset.

These are repository-only maintenance commands, not installed-package entry
points:

```sh
uv run --frozen python -m scripts.reconstruct_similar_eos_uncertainties
uv run --frozen python -m scripts.reconstruct_similar_eos_uncertainties --check
uv run --frozen python -m scripts.integrate_uncertainty_followup
```

Integration updates only the selected reconstruction assessments, their
ledger rows and aggregate counts. Material summaries, reconstruction reports,
figures and the paper ledger then render that evidence. No source coefficient,
observation, default or recommendation is replaced, and no release or
deployment is performed by these commands.

## Verification

The focused Python suite passed 108 tests, including all 63 numerical replays,
source/declaration checksums, native/legacy compatibility round trips, model
row/RMS replay and generated report/figure consistency. Independent tests cover
closed-form weighted linear covariance with mixed parameter units, conversion
of both Ma covariance axes, and 600 independently implemented simulated staged
post-perovskite refits. Singular designs and inconsistent reconstructed
objectives are rejected. The dedicated uncertainty `--check` passed.

The strict documentation build completed, including execution of the example
notebooks. Representative updated coefficient plots were visually inspected.
The Rust bundled-catalog load/round-trip test passed. The wider Rust EOSMAT
suite had 40 passes and one unrelated failure in the expected error-message
text for a synthetic unsupported double-Debye/Murnaghan combination; this
follow-up changes no Rust model or error text.

The full ledger freshness check encountered an existing fixed-parameter order
mismatch in `fe09ni01_hcp_sakai_2014_type1_thermal`. That entire row is unchanged
from the saved pre-follow-up ledger and is outside the 63 selected records.
The broader record was not regenerated or reclassified to hide the mismatch.
All unselected ledger rows and all unselected records in affected material
files were independently checked against the pre-follow-up snapshot and remain
unchanged.
