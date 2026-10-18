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

**2026-10-17 disposition:** MgO and NaCl-B2 are `not_reproduced` and
`not_recommended_for_quantitative_use`, retained only for explicit diagnostic
selection. The eight other phase fits retain conditional coefficient parity.
The published coefficients, reported errors, fixed volumes, primary observations
and existing material defaults are preserved. Diagnostic refits are not
validated replacement pressure standards.

The official workbook retrieved from the APS supplement page has SHA-256
`339eb0d160316ef3e32abe67c5d5e917dfc6be4b4fbf726287680cca61596565`.
The earlier two-sheet import (`DAC-1!A1:S232` and `DAC-2!A1:Y199`) was normalized
with the checked-in transcription rules: all 3,511 rows match the bundled CSV
exactly in experiment, run, phase, sequence label, volume, and standard error.
This is the earlier workbook verification, not a fresh download on October 17;
the current APS supplement page and Table-S1.xlsx requests return HTTP 401.
The October 17 check uses the bundled transcription and locally available final
article and seven-page figures/stress supplement. The final article SHA-256 is
`fb101fe78237c91077586d8999103a10c7ee75813a13ada441ee81a880c83f72`.

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
| MgO | 194 | 2.786--110.250 | 161.9, 4.08 | 161.827, 3.954 | 0.476 | not reproduced; not recommended |
| Mo | 283 | 2.786--110.250 | 260.6, 4.06 | 260.687, 4.054 | 0.418 | parity |
| NaCl-B1 | 55 | 2.786--29.977 | 23.5, 5.23 | 23.474, 5.231 | 0.212 | parity |
| NaCl-B2 | 125 | 35.714--109.513 | 24.92, 5.63 | 25.447, 5.583 | 0.293 | not reproduced; not recommended |
| Pt | 387 | 2.786--110.250 | 261.2, 5.75 | 261.262, 5.722 | 0.377 | parity |
| Ta | 198 | 3.909--115.733 | 195.2, 3.62 | 195.309, 3.614 | 0.366 | parity |
| W | 335 | 20.358--136.622 | 307.0, 4.02 | 304.921, 4.114 | 1.156 | parity |

Eight records meet the legacy combined-two-error-width coefficient criterion.
This is conditional numerical agreement, not recovery of the author workflow;
the source does not identify the confidence convention of its printed errors.
The broad numerical similarity screen also accepts MgO and NaCl-B2, but that
screen does not resolve their published-curve discrepancies. Their ledger
outcome is now `parity_not_achieved`, with `coefficient_comparison_status:
similar` retained separately.

## Published-curve discrepancy

Residual means P(Table II EOS at measured sample volume) minus P(paired Cu).

| Phase | Published RMS (GPa) | Refit RMS (GPa) | Published mean (GPa) | Published residual range (GPa) |
|---|---:|---:|---:|---:|
| MgO | 1.161241 | 0.475879 | +0.820916 | -0.4407 to +2.8032 |
| NaCl-B2 | 0.622381 | 0.293229 | -0.547610 | -1.337712 to -0.009561 |

All 125 selected NaCl-B2 residuals are negative. At fixed K0-prime,
Vinet pressure is proportional to K0; increasing K0 therefore decreases the
least-squares objective for any fixed positive row weights on this dataset.
The printed NaCl-B2 pair cannot be a stationary solution of that objective,
even for another nonempty subset of these observations. This argument does
not exclude an unpublished processed input table, extra observations or a
different parameter-dependent statistical objective.

For MgO, the large pressure trend and K0-prime shift from 4.08 to 3.953583
remain unresolved under the direct unweighted replay. Unlike NaCl-B2, the
MgO residuals have both signs; an unspecified weighting or processing protocol
has not been mathematically excluded. The finding is failure to recover the
published fit from available inputs, not proof that the authors' unrecovered
workflow was wrong. Both refits reduce RMS by more than a factor of two, but
retain structured residuals; neither is promoted as a replacement standard.

The independent October 17 calculation transcribed both equations separately,
agreed with the executable engine within 2.85e-13 GPa, and confirmed each
two-parameter optimum with a one-dimensional profile eliminating K0 analytically.
The maintained ledger now records both published and refitted residuals.

## Could additional experimental fit inputs have been missed?

### Fitting software and original input files

The final article and inspected seven-page figures/stress supplement do not
identify the software used for the Table II EOS regressions. Section II.C
(p. 144113-2) names DIOPTAS for diffraction-image integration and background
subtraction and describes Moffat peak fitting; this does not identify the later
EOS fitting implementation. Section III.E (p. 144113-8) specifies least-squares
Vinet fitting with V0 fixed and K0/K0-prime varied, while p. 144113-12 identifies
direct P-V fits as the Table II primary results and K(V) fitting as a separate
consistency check. The EOS regression program, optimizer, residual coordinate
and weighting procedure are unspecified. Neither EosFit nor custom code can
be inferred from that description.

Targeted title/DOI searches on October 17, including EosFit, code/software,
GitHub and Zenodo queries, found no identifiable public author fitting code.
The current publisher supplement listing remains inaccessible without
authentication, so this is a bounded search finding, not proof that no fitting
files exist. The useful author request is the software name/version and the
original selected inputs with uncertainties, plus saved fitting/project files
and commands or the original script, including the Cu pressure conversion.

**Actual author fit-input completeness is unconfirmed.** The final article's
data-availability statement (p. 144113-15) distinguishes data in Table S1 from
other data available on request. That statement does not identify those other
data or establish that they entered the EOS regressions. The figures/stress
supplement contains peak-width, stress and inter-run comparisons (Figs. S1-S5),
not an additional numerical MgO or NaCl fit table.

The evidence for the input set described by the paper is:

- Section II.A, p. 144113-2: DAC-1 contains only the seven metals; MgO and
  NaCl occur only in DAC-2. The absent DAC-1 MgO/NaCl rows are therefore expected.
- Sections II.B-C, pp. 144113-2/3: 194 DAC-2 steps are reported, and the
  measured volumes/errors are assigned to Table S1. The bundled table has
  **194 MgO readings, one for every DAC-2 run 1-194**, and 131 NaCl-B2 readings.
- Section III.E and Table II caption, pp. 144113-8/9: pressures are assigned
  to the simultaneously measured volumes from the adopted Cu scale and fitted
  directly. Earlier literature parameters are presented for comparison; no
  combined fit with those earlier experiments is specified. Adopted NaCl-B2
  V0 from Fei is a fixed parameter, not evidence that Fei's observations entered.
- The paragraph following Table II, p. 144113-12, identifies direct P-V fits
  as the Table II primary results; differentiated polynomial K(V) fits are a
  separate consistency check, not additional experimental observations.

All 194 available MgO points are used. For NaCl-B2, six of the 131 rows are
excluded by the conservative P(Cu)>35 GPa mask. Table I footnote f excludes
the approximately 30-35 GPa coexistence region, whereas Table II says above
33 GPa; the exact author mask is unavailable. Including all 130 points above
33 GPa gives K0=25.631126 GPa, K0-prime=5.558143 and RMS 0.348236 GPa, moving
farther from the published K0. All 131 available B2 observations also have
negative residuals at the printed coefficients. This ambiguity alone does
not restore the printed fit.

There is a separate count inconsistency: Section II.B reports 229 DAC-1
steps, but the earlier workbook audit and bundled transcription contain 228
numbered DAC-1 runs. Its cause is unresolved. It cannot supply an identified
missing MgO/NaCl run, since those materials were absent from DAC-1.

Thus no explicitly described MgO/NaCl experimental input set has been found
missing, but unpublished inputs, a different reduction version, or undocumented
processing cannot be ruled out. Resolving this requires the exact author-used
MgO and NaCl-B2 P-V tables, row/run selections (including any added historical
measurements), Cu reduction code/constants, residual coordinate/weights and
confirmation of the Table II coefficient version. See the
[reviewed data request](../dataset-requests.md#shen-and-smith-2026-mgo-and-nacl-b2).
No author has been contacted.

Figure 7 (p. 144113-8) also reports Au P-V coefficients 168.76 GPa and 5.82,
whereas Table II gives 167.5 GPa and 5.85 at the same fixed cell volume
67.792 A^3. The bundled observations closely recover the latter. This is
evidence of inconsistent reporting, but does not establish the cause of the
MgO/NaCl-B2 discrepancies. W retains large, structured, experiment-dependent
residuals (published RMS 1.180682 GPa; refit 1.155888 GPa); coefficient parity
does not establish that these observations or any calibration are exact.

The executable reconstruction and all record-level diagnostics are generated
by `scripts/validate_primary_eos_refits.py`; the exact outputs are stored in
`docs/data/primary-eos-refits.json` and summarized in
`docs/primary-eos-refits.md`.

The reconstruction delegates pressure and reference volume to
`copper_fratanduono_2020_vinet3_298k` (`Vinet3`), so the
catalog and audit share one executable equation and coefficient set. All ten
Shen–Smith records link this Cu record in their pressure-calibration metadata.
See the [Cu audit](fratanduono-2020-cu.md) for uncertainty and fit-data limits.
