# Zha, Bassett, and Shim (2004): rhenium pressure scale

Audit date: 2026-09-15. The final publisher PDF and publisher snapshot were
read from the user's Zotero library. Equations and parameter/observation
tables were checked on rendered PDF pages. No supplementary or correction
links were present in the saved publisher page.

[Review of Scientific Instruments 75, 2409-2418](https://doi.org/10.1063/1.1765752).
Zotero item `9BAJ2FCF`, PDF attachment `JABFZJUC`.

The executable record `rhenium_zha_2004_bm3_log_thermal` represents the
**continuous Equation (6), Table III parameterization**:

\[
P(V,T)=P_{\rm BM3}(V;29.4087,360,4.5)
 +[0.00776-0.00815\ln(29.4087/V)](T-300).
\]

Pressure is in GPa, temperature in K, and volume in angstrom cubed per hcp
conventional cell (two Re atoms). The existing hcp structure remains applicable;
its reference volume differs from the JCPDS 05-0702 value fixed in this paper.
The static coefficients and 300 K reference are fixed. Neither uncertainties
nor covariance are reported for the two fitted thermal coefficients. The prose
on page 2414 reads `4.5` followed by superscript reference 3, not `4.53`.

This fits the existing `BM3` + `LogVolumeThermalPressure` implementation exactly.
No new model family or thermal interpolation is introduced.

## Observations and independent refit

All eight Table IIb Re measurements are bundled with their matching Table IIa
Au marker observations: four d spacings, four individual lattice estimates,
mean lattice parameter, its printed error and relative error, and cell volume.
Printed lattice parameters and volumes are retained separately, even where
their rounding/consistency differs. No row-level Re volume, temperature, or
pressure uncertainties are supplied. The Appendix's approximately 0.28 GPa
aggregate pressure uncertainty is not assigned as a measured error to every row.
Table VII's six diffraction spacings are also retained as strain diagnostics;
only the 100/101 pair contributes to the Table IIb fit volumes.

The pressure calibration is Anderson, Isaak, and Yamamoto (1989),
DOI `10.1063/1.342969`, with the paper's explicitly substituted Au reference
volume **67.847 angstrom cubed**. Exact recalculation must apply this reference
volume instead of using the related catalog record's default volume unchanged.
Even with that substitution, the catalog Anderson coefficients predict pressures
up to **0.318 GPa higher** than Table IIa, increasingly at high temperature.
The exact source realization of that calibration is therefore only partially
resolved. The related record is not registered as an exact executable calibration
link. The source refit uses the authors' printed gold pressures, and the
reproduction report retains all eight calibration residuals.

### Follow-up: why the Au calibration differs

The original [Anderson et al. (1989) paper](https://doi.org/10.1063/1.342969)
was checked directly in Zotero (item `6A54H6QK`, PDF `IJ9TMXKL`). Equation
(29) and the adjoining text on page 1541 confirm the catalog's final choices:
K0=166.65 GPa, K0'=5.4823, alpha KT=0.00714 GPa/K, and
(dKT/dT)V=-0.0115 GPa/K. The independent Anderson reproduction recovers all
126 entries in the upper rows of its Table V within 0.005 GPa. Thus this
disagreement does not justify changing the catalog's Anderson coefficients.

Using those coefficients with **Zha's** reference volume gives:

| Temperature (K) | Zha Table IIa (GPa) | Anderson Equation (29) (GPa) | Calculated minus printed (GPa) |
|---:|---:|---:|---:|
| 1380.3 | 7.29 | 7.2963 | +0.0063 |
| 1480.2 | 8.47 | 8.4820 | +0.0120 |
| 1506.5 | 8.00 | 8.0361 | +0.0361 |
| 1548.6 | 7.58 | 7.6281 | +0.0481 |
| 1628.5 | 7.46 | 7.5325 | +0.0725 |
| 1716.8 | 6.41 | 6.5882 | +0.1782 |
| 1801.0 | 6.48 | 6.7188 | +0.2388 |
| 1914.5 | 6.66 | 6.9777 | +0.3177 |

Three checks narrow the explanation:

- **Reference volume:** substituting 67.847 for the catalog's 67.79 angstrom
  cubed raises pressure by about 0.100-0.129 GPa across these rows. That
  substitution is required by Zha, but does not remove the discrepancy.
- **Lattice conversion and rounding:** every printed Au volume agrees with
  the cube of its printed mean lattice parameter within 0.00005 angstrom
  cubed. Replacing that mean with the arithmetic mean of the four printed
  reflection estimates changes pressure by at most 0.0193 GPa. Simultaneously
  varying the inputs and coefficients by half their last printed digit gives
  a final-row residual of +0.296 to +0.339 GPa at the interval corners,
  including rounding of the printed pressure. Rounding cannot account for
  +0.318 GPa. These are numerical precision checks, not measurement-error
  estimates or confidence intervals.
- **Thermal parameter choice:** Anderson's Section III (page 1537) explores
  four paired K0' and (dKT/dT)V choices. Substituting each pair in Equation
  (29), with the other coefficients and Zha's Va held fixed, gives:

| K0' | (dKT/dT)V (GPa/K) | Pressure RMS residual (GPa) | Largest absolute residual (GPa) |
|---:|---:|---:|---:|
| 6.39 | -0.0052 | 0.02379 | 0.05178 |
| 6.12 | -0.0071 | 0.05970 | 0.13233 |
| 5.50 | -0.0115 | 0.15840 | 0.31933 |
| 5.21 | -0.0135 | 0.20362 | 0.40360 |

These are sensitivity comparisons using published exploratory pairs, not four
equally endorsed final Anderson scales. The 6.39/-0.0052 pair is much closer
to Zha's values, but still does not reproduce every printed pressure. Zha
does not state that this pair was used, so its use cannot be inferred as fact.

The direction of the effect is understandable from the thermal term
`b * ln(Va/V) * (T-300)`: seven of the eight Au volumes exceed Va, so the
logarithm is negative. A more negative `b` produces more positive thermal
pressure, with increasing sensitivity as the gold expands. This explains why
the thermal parameter choice can produce the observed trend. It does not
identify the authors' actual implementation. These expanded states also lie
outside the nonnegative compression grid tabulated by Anderson, although the
printed formula can be evaluated there. Temperature and volume covary in
Zha's eight observations, which limits causal attribution from the residuals.

**Conclusion:** the catalog reproduces Anderson's final published equation;
Zha's tabulated realization remains unresolved. Different thermal coefficients
are a plausible explanation, but neither an exact alternative nor an author
correction is established. Printed Zha pressures and published Re coefficients
are retained. The approximately 0.28 GPa aggregate uncertainty discussed by
Zha is a separate experimental consideration and cannot specify which
calibration formula generated the table.

## Reproducing the Re thermal coefficients

Section III first constructs an isochore slope for each observation,
`s_i = [P_i - P_BM3(V_i)]/(T_i-300)`, then regresses those slopes against
`ln(Va/V_i)`. An unweighted least-squares reconstruction gives:

| Coefficient (GPa/K) | Published | Reconstructed | Relative difference |
|---|---:|---:|---:|
| alpha KT | 0.00776 | 0.0077561565 | 0.0495% |
| (dKT/dT)V | -0.00815 | -0.0081487890 | 0.0149% |

Both agree within half a unit of the last printed decimal. The pressure RMS
residual is 0.40849 GPa for the published coefficients and 0.40838 GPa for
the reconstruction; the largest published residual is 0.77234 GPa.
These residuals are disclosed rather than claiming that all observations lie
within the Appendix's aggregate uncertainty. The ledger classification is
`similar`: numerical coefficient agreement is established, while formal
uncertainty parity is unavailable.

The source does not explicitly enumerate regression weights. The reconstruction
uses equal weights on the isochore slopes, which recovers Table III. Changing
the objective to unweighted pressure residuals gives 0.007779528 and
-0.005396049 GPa/K; the second coefficient is sensitive to the objective over
this narrow volume range. Neither diagnostic fit replaces the published record.

## Distinguish the two published thermal representations

The authors subsequently refit calculated isotherms to BM3 in Table IV. They
prefer the branch with zero-pressure volume fixed from independent thermal
expansion data. Both branches (12 parameter sets total) are bundled as derived
data, together with all 147 Table V pressures.

Equation (6) reproduces the six **unconstrained** Table IV isotherms within
0.066 GPa at `Va/V = 1, 1.1, 1.2`, covering up to 3000 K. This is a comparison
to independently printed source output, not a refit to synthetic observations;
a 0.1 GPa tolerance allows the rounded coefficients and approximate isotherm
representation. Separately, the preferred **fixed-volume** Table IV branch
reproduces all Table V entries within 0.014 GPa (0.02 GPa tolerance for rounded
parameters and pressures).

The continuous Equation (6) surface differs from Table V by as much as
0.772 GPa. It is labeled as the Table III branch and must not be presented as
the authors' preferred Table V scale. The Table IV rows are not promoted as
independent static EOS records, and no unstated temperature interpolation is
assumed.

Actual observations span only **6.41-8.47 GPa and 1380.3-1914.5 K**. The
published grids extend to 20% compression and 3000 K, far beyond those
observations. The stored experimental coverage is not a rectangular stability
claim; evaluation beyond it is source extrapolation.

## Reproduction and verification

Run `python scripts/reproduce_zha_2004_xian_2022_rhenium.py` to regenerate
the [numerical reproduction report](../data/zha-2004-xian-2022-rhenium-reproduction.json).
Its `zha_*` sections contain the independent formulas, staged refit, Au
calibration diagnostics, and derived-table comparisons. These are covered by
`tests/test_zha_xian_rhenium.py`; the normal catalog tests additionally
validate resource checksums and Python/Rust loading.
