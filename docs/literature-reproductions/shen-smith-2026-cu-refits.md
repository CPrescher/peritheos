# Shen and Smith (2026) Cu-referenced Vinet audit

## Scope and primary evidence

This audit covers the ten Vinet records reported by Shen and Smith,
*Simultaneous x-ray diffraction measurements of nine pressure calibrants to
140 GPa*, Physical Review B **113**, 144113 (2026),
[doi:10.1103/fxgq-96sg](https://doi.org/10.1103/fxgq-96sg). The inspected
primary evidence is the 17-page final article and the official Supplemental
Table S1 Excel workbook. The checked-in long-form transcription contains all
3,511 nonblank volume/error pairs from its DAC-1 and DAC-2 worksheets and
preserves experiment, run, phase, and Pt first/last labels.

The official workbook retrieved from the APS supplement page has SHA-256
`339eb0d160316ef3e32abe67c5d5e917dfc6be4b4fbf726287680cca61596565`.
A fresh two-sheet import (`DAC-1!A1:S232` and `DAC-2!A1:Y199`) was normalized
with the checked-in transcription rules: all 3,511 rows match the bundled CSV
exactly in experiment, run, phase, sequence label, volume, and standard error.

The previously unresolved dependency is Ref. 2: Fratanduono et al., *Probing
the Solid Phase of Noble Metal Copper at Terapascal Conditions*, Physical
Review Letters **124**, 015701 (2020),
[doi:10.1103/PhysRevLett.124.015701](https://doi.org/10.1103/PhysRevLett.124.015701).
Its main article Table I and official Supplemental Material Section S4,
Equation (2), completely specify the adopted 298 K Cu pressure scale.

## Exact pressure reconstruction

Fratanduono's third-order Vinet form is

\[
P(X)=3K_0\frac{1-X^{1/3}}{X^{2/3}}
\exp\left[\eta(1-X^{1/3})+\beta(1-X^{1/3})^2+
\psi(1-X^{1/3})^3\right],
\]

where \(X=\rho_0/\rho=V/V_0\). The central 298 K coefficients are
\(\rho_0=8.939\ \mathrm{g\,cm^{-3}}\), \(K_0=133.6\ \mathrm{GPa}\),
\(\eta=6.29\), \(\beta=2.06\), and \(\psi=1.65\). Using the conventional
Cu molar mass 63.546 g/mol and four atoms per fcc cell gives
\(V_0=47.21808505\ \mathrm{\mathring{A}^3}\). The supplement's Table S1 is
the reduced *isentrope*, not a 298 K lookup table, and is therefore not used as
pseudo-isotherm data.

For each selected experiment/run, the implementation evaluates this equation
at the one Cu unit-cell volume and assigns that pressure to every simultaneous
target-phase volume. Pt's explicitly labeled first and last measurements are
both retained. This is the paper's Section III.E construction, not an
interpolation or a fit to model-generated pressure points.

## Source row masks and objective

The masks follow Table II and the article's data-quality discussion:

- Pt uses DAC-2, Au uses DAC-1, and Ta uses DAC-1; the rejected complementary
  datasets remain in the bundled workbook transcription.
- MgO and both NaCl phases use DAC-2. NaCl-B1 uses reconstructed pressures
  below 30 GPa and NaCl-B2 uses pressures above 35 GPa, excluding the prose's
  approximately 30--35 GPa coexistence interval.
- Mo, W, bcc Fe, and hcp Fe combine DAC-1 and DAC-2.

Shen and Smith state that \(K_0\) and \(K'_0\) were obtained by least squares
with \(V_0\) fixed. They do not publish pressure uncertainties, regression
weights, a residual covariance, or a parameter covariance. The reproduction
therefore minimizes unweighted pressure residuals. Table S1's crystallographic
volume standard errors remain available as source observations but are not
silently converted into an errors-in-variables objective.

## Fixed-volume refit results

| Phase | n | Pressure range (GPa) | Published \(K_0,K'_0\) | Refit \(K_0,K'_0\) | RMSE (GPa) | Outcome |
|---|---:|---:|---:|---:|---:|---|
| bcc Fe | 42 | 2.786--15.735 | 162.1, 5.40 | 162.351, 5.345 | 0.215 | parity |
| Au | 228 | 3.909--137.377 | 167.5, 5.85 | 167.642, 5.844 | 0.335 | parity |
| hcp Fe | 380 | 15.361--136.622 | 168.55, 5.53 | 169.011, 5.498 | 0.727 | parity |
| MgO | 194 | 2.786--110.250 | 161.9, 4.08 | 161.827, 3.954 | 0.476 | similar |
| Mo | 283 | 2.786--110.250 | 260.6, 4.06 | 260.687, 4.054 | 0.418 | parity |
| NaCl-B1 | 55 | 2.786--29.977 | 23.5, 5.23 | 23.474, 5.231 | 0.212 | parity |
| NaCl-B2 | 125 | 35.714--109.513 | 24.92, 5.63 | 25.447, 5.583 | 0.293 | similar |
| Pt | 387 | 2.786--110.250 | 261.2, 5.75 | 261.262, 5.722 | 0.377 | parity |
| Ta | 198 | 3.909--115.733 | 195.2, 3.62 | 195.309, 3.614 | 0.366 | parity |
| W | 335 | 20.358--136.622 | 307.0, 4.02 | 304.921, 4.114 | 1.156 | parity |

Eight records meet Peritheos's combined two-standard-error parity criterion.
MgO and NaCl-B2 are defensible numerical reproductions but remain `similar`:
their point estimates meet the numerical tolerance while one coefficient lies
outside combined two-standard-error agreement.

The remaining limits are source-side and explicit. For MgO, the paper does not
publish a weighting protocol that could explain the small derivative shift.
For NaCl-B2, Table II's short qualifier says "above 33" GPa while the main text
says the approximately 30--35 GPa coexistence region was excluded; the latter
is the only unambiguous exclusion rule and is used here. Trying alternate
weights, inventing covariance, or fitting the article's smoothed V--V curves
would not be a source-faithful resolution.

The executable reconstruction and all record-level diagnostics are generated
by `scripts/validate_primary_eos_refits.py`; the exact outputs are stored in
`docs/data/primary-eos-refits.json` and summarized in
`docs/primary-eos-refits.md`.

The reconstruction delegates pressure and reference volume to
`copper_fratanduono_2020_vinet3_298k` (`Vinet3`), so the
catalog and audit share one executable equation and coefficient set. All ten
Shen–Smith records link this Cu record in their pressure-calibration metadata.
See the [Cu audit](fratanduono-2020-cu.md) for uncertainty and fit-data limits.
