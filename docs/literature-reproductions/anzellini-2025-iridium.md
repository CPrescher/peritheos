# Anzellini et al. (2025): iridium thermal EOS

## Outcome

The catalog record is now the complete published BM3 plus Holland--Powell
thermal-pressure parameterization, rather than a 300 K BM3 curve attached to
high-temperature observations. The exact published pressure surface is
executable in Python and Rust.

The available numerical evidence supports a conditional reproduction, not the
paper's complete four-parameter regression. Holding the published 300 K BM3
coefficients fixed and fitting the thermal-expansion coefficient to all 122
rows in Supplementary Tables 1--3 gives
`alpha0 = 1.77403(33)e-5 K^-1`, within combined two-sigma uncertainty of the
published `1.87(8)e-5 K^-1`. The resulting pressure RMSE is 3.033 GPa, compared
with 3.138 GPa for the rounded published coefficients.

## Primary sources and exact model

The primary article is Anzellini et al., *Phase stability of iridium at extreme
conditions*, Communications Materials **6**, 221 (2025),
[doi:10.1038/s43246-025-00963-4](https://doi.org/10.1038/s43246-025-00963-4).
The audit used the article's Equations 1--3, Table 1, Methods, and Supplementary
Figure 5, plus Supplementary Tables 1--3. The source article and supplement
PDFs inspected for this audit had SHA-256 hashes
`3c531753df70d77a32433c6663dc54b6796650d9412a2edbd981b130ea91aa01`
and `f661efccbb391e92cb6adf3070c4a21a01f5fba0a2306075175948f344980f73`;
the PDFs are not redistributed in this repository.

The source writes

\[
P(V,T)=P(V,300\,\mathrm{K})+P_{\mathrm{Th}}(T),
\]

where the first term is third-order Birch--Murnaghan and

\[
P_{\mathrm{Th}}(T)=\frac{\alpha_0K_0\Theta_E}{\zeta_0}
\left[
\frac{1}{\exp(\Theta_E/T)-1}-
\frac{1}{\exp(\Theta_E/300)-1}
\right],
\]

\[
\zeta_0=\left(\frac{\Theta_E}{300}\right)^2
\frac{\exp(\Theta_E/300)}{[\exp(\Theta_E/300)-1]^2}.
\]

For iridium, `theta = 298 K` is fixed. Table 1 reports
`V0 = 56.62(5) A^3`, `K0 = 327(2) GPa`, `K0_prime = 5.46(3)`, and
`alpha0 = 1.87(8)e-5 K^-1`. Peritheos evaluates the algebraically equivalent
Einstein-energy form
`alpha0*K0/Cv(300) * [E(T)-E(300)]`; the atom count cancels from pressure.

The new generic `HollandPowellThermalPressure` model accepts the BM3 reference
used here. `ThermalModifiedTait` remains the backward-compatible restricted
form for the older modified-Tait endmember records.

## Temperature reduction and dataset boundary

The Methods section says that the temperature assigned to an observation is
the arithmetic mean of the simultaneous upstream and downstream spectral-
radiometry temperatures. The dataset therefore adds `mean_temperature_k`
without modifying either measured temperature or the source-printed
temperature-difference column. No synthetic uncertainty is attached to the
mean because the paper does not specify how the two side-specific uncertainty
estimates were propagated into the EOS regression.

The bundled CSV contains 122 observations from the laser-heating runs, including
two 300 K baseline rows, and has SHA-256
`2133efc7b402be8d804d098ef8badfe0e46d8a7ffc91351ca71bffaf6fc85250`.
It retains the applicable MgO or KCl calibrant, pressure, Ir lattice parameter,
both pyrometry measurements and uncertainties, their derived mean, and the
source-printed temperature difference.

The article states that its reported coefficients were fitted with EosFit7 to
a broader combined dataset:

- the 300 K helium-medium observations of Monteseguro et al. (2019),
  [doi:10.1038/s41598-019-45401-x](https://doi.org/10.1038/s41598-019-45401-x);
- the present resistive-heating data at 833 K; and
- the present laser-heating data.

The first two row sets are visible in plots but are not published numerically.
The sources also do not state the exact row selection behind the plotted
plus-or-minus-100 K temperature averages, the residual variable, regression
weights, or parameter covariance. Plot digitization would add pseudo-precision
and is not substituted for those missing inputs.

## Reproduction and identifiability diagnostic

The executable audit is `scripts/validate_primary_eos_refits.py`. Its
source-conditioned check uses all 122 published rows, their arithmetic-mean
temperatures, unweighted pressure residuals, and the published fixed values
`V0`, `K0`, `K0_prime`, `Tr`, `theta`, and `n`; only `alpha0` is refined. This
isolates the thermal term that the available hot data can independently test.

As a diagnostic, allowing `V0`, `K0`, `K0_prime`, and `alpha0` all to vary on
the 122-row subset produces `(55.8741 A^3, 475.405 GPa, 2.07378,
1.18260e-5 K^-1)` with a 2.825 GPa RMSE. Those substantially displaced
reference-isotherm coefficients demonstrate that this mostly hot subset does
not identify the published 300 K reference curve. This diagnostic is not a
replacement EOS and is not stored as a catalog record.

Accordingly, the ledger classifies the independently constrained `alpha0`
check as parity while separately recording
`complete_combined_fit_status = not_refittable_from_published_rows`. No heated
observation is collapsed onto or directly fitted as if it lay on the 300 K
curve.
