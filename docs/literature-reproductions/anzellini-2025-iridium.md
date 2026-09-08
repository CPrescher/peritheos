# Anzellini et al. (2025): iridium thermal EOS

## Outcome

The catalog record is now the complete published BM3 plus Holland--Powell
thermal-pressure parameterization, rather than a 300 K BM3 curve attached to
high-temperature observations. The exact published pressure surface is
executable in Python and Rust.

The exact tabulated evidence supports a conditional reproduction. Holding the
published 300 K BM3 coefficients fixed and fitting the thermal-expansion
coefficient to all 122 rows in Supplementary Tables 1--3 gives
`alpha0 = 1.77403(33)e-5 K^-1`, within combined two-sigma uncertainty of the
published `1.87(8)e-5 K^-1`. The resulting pressure RMSE is 3.033 GPa, compared
with 3.138 GPa for the rounded published coefficients.

Figure 4 also contains vector marker paths for the otherwise untabulated 300 K
and 830 K series. A clearly labeled plot-derived reconstruction combining those
107 markers with 65 exact laser-heating rows in the paper's stated
plus-or-minus-100 K isotherm bands gives `V0 = 56.62669 A^3`,
`K0 = 326.15546 GPa`, `K0_prime = 5.58105`, and
`alpha0 = 1.82917e-5 K^-1`, with a 1.986 GPa pressure RMSE. `V0`, `K0`, and
`alpha0` fall inside the paper's reported one-sigma intervals; `K0_prime` is
0.121 above the published value, about four times its printed one-sigma error.
This materially reconstructs the combined fit but does not establish exact
source parity.

## Primary sources and exact model

The primary article is Anzellini et al., *Phase stability of iridium at extreme
conditions*, Communications Materials **6**, 221 (2025),
[doi:10.1038/s43246-025-00963-4](https://doi.org/10.1038/s43246-025-00963-4).
The audit used the article's Equations 1--3, Table 1, Methods, and Supplementary
Figures 4 and Supplementary Figure 5, plus Supplementary Tables 1--3. The source article and supplement
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
- the present resistive-heating data represented by the 830 K Figure 4
  isotherm (the underlying run is described elsewhere as approximately 833 K);
- the present laser-heating data.

The first two row sets are visible in plots but are not published numerically.
The sources also do not state the exact row selection behind the plotted
plus-or-minus-100 K temperature averages, the residual variable, regression
weights, or parameter covariance. The plot-derived coordinates below are
therefore stored separately from exact source rows and are not substituted for
the unavailable author data.

## Figure 4 vector digitization

The official article PDF stores the Figure 4 filled circles as vector paths.
`scripts/digitize_anzellini_2025_figure4.py` selects marker circles by geometry
and series color, reads their centers directly in PDF points, and applies a
linear calibration from the vector frame and ticks:

- pressure: `(346.745 pt, 0 GPa)` to `(552.845 pt, 150 GPa)`;
- volume: `(248.247 pt, 44 A^3)` to `(50.643 pt, 59 A^3)`.

The resulting CSV contains 91 black 300 K markers attributed by the legend to
Monteseguro et al. (2019) and 16 dark-red markers on the source-labeled 830 K
isotherm. It has SHA-256
`097f35f293d7f6d93d286e86ed72166e7660817132e481b43e3301cac2d7bd8d`.
Fitted curves, the other temperature series, ambient-pressure literature
triangles, and legend symbols are excluded. Every row preserves its PDF
coordinates and carries conservative graphical bounds of 0.75 GPa and
0.08 A^3, approximately a marker radius plus calibration allowance. These are
digitization bounds, not experimental errors or regression weights.

As a cross-check, the original Monteseguro et al. Figure 2 contains the same 91
black-circle paths. After independent axis scaling, its point order correlates
with the Anzellini replot at 0.999947 in pressure and 0.999931 in volume. A
volume-residual BM3 fit to the higher-resolution Figure 4 replot gives
`(V0, K0, K0_prime) = (56.60041 A^3, 333.42839 GPa, 5.34720)` with a
0.0991 A^3 RMSE; the rounded Monteseguro coefficients give 0.1137 A^3 on those
same presentation coordinates.

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
check as parity while separately recording the reconstructed four-parameter
result and `complete_combined_fit_status = plot_derived_reconstruction_only`.
The reconstruction uses 91 Figure 4 markers at 300 K, 16 at 830 K, and the 65
exact supplementary rows within 100 K of 1700, 2000, 2200, 2600, 3300, 3800,
or 4200 K. It uses unweighted pressure residuals because the source weights are
unpublished; the graphical uncertainties remain metadata rather than invented
author weights. No heated observation is collapsed onto or directly fitted as
if it lay on the 300 K curve.
