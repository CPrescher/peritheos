# Literature reproductions

This page records numerical comparisons between published equations of state
and their Peritheos implementations. It complements the material files: an
EOS record preserves the recommended literature parameterization, while this
page explains whether the primary observations are available, exactly how a
refit was performed, and whether the published parameters can be recovered.

For every future material, record the following whenever the source permits:

1. the primary publication and the exact table, supplement, or repository;
2. a local machine-readable transcription of the observations;
3. the pressure, volume, temperature, and uncertainty columns actually used;
4. the EOS form, units, reference state, fixed parameters, bounds, residual
   definition, and covariance scaling;
5. published parameters beside Peritheos refits, with goodness-of-fit and
   parameter correlations;
6. unresolved source contradictions and any reason the fit is not exactly
   reproducible; and
7. separate crystallographic provenance sufficient to calculate diffraction
   reflections: unit cell, space group, formula units per cell, and atom sites
   or a documented fallback peak list.

The source-reported parameter set remains the library record unless a correction
is explicitly justified. A refit is supporting evidence, not permission to
silently replace the publication.

## C01: boron carbide, Somayazulu et al. (2023)

### Sources and model

- Primary article: Somayazulu et al., *P-V-T equation of state of boron
  carbide*, Philosophical Transactions of the Royal Society A **381**, 20220331
  (2023), [doi:10.1098/rsta.2022.0331](https://doi.org/10.1098/rsta.2022.0331).
- Primary data: official supplementary Table 4,
  [doi:10.6084/m9.figshare.c.6751752](https://doi.org/10.6084/m9.figshare.c.6751752).
  It contains 51 rows with `P`, `V`, `T`, and a reported uncertainty for
  each coordinate. The exact transcription is embedded as dataset
  `b4c_somayazulu_2023_table4_pvt` in the bundled `b4c.eosmat` material and is
  linked to both fitted EOS records.
- EosFit method: Angel, Gonzalez-Platas, and Alvaro (2014),
  [doi:10.1515/zkri-2013-1711](https://doi.org/10.1515/zkri-2013-1711), and the
  [official open-source CrysFML `CFML_EoS` implementation](https://code.ill.fr/scientific-software/CrysFML2008/-/tree/35549dff98213088b580618a294be48475535045/Src/CFML_EoS)
  at commit `35549dff98213088b580618a294be48475535045`.
- Reference isotherm: third-order Birch--Murnaghan with conventional-cell
  `V0 = 328.4 Å³`, `K0 = 221 GPa`, and `K0′ = 3.3`, all fixed during the
  thermal fit.
- Thermal model: Mie--Gruneisen--Debye with `Tr = 298 K`, `theta0 = 1425 K`,
  `n = 5`, and the integrated constant-`q`
  relations

  \[
  \gamma(V)=\gamma_0(V/V_0)^q,
  \qquad
  \theta(V)=\theta_0\exp[(\gamma_0-\gamma(V))/q].
  \]

The paper's generalized Tange form has `a = 1` and is therefore exactly the
`integrated_gruneisen` model implemented by `MieGruneisenDebye`.

### Numerical comparison

All fits use the conventional-cell-to-molar-volume conversion shown in the
equation below, because the diffraction cell contains nine B4C
formula units. Reported sigmas are treated as independent one-standard-deviation
errors, and `absolute_sigma=True` prevents rescaling the covariance by the fit
residuals.

\[
V_m=V_{cell}N_A10^{-25}/9
\]

| Calculation | Data and objective | Result | Reduced chi-square |
|---|---|---:|---:|
| Published curve, no refit | 41 heated rows; `P` residual at printed `V,T` | `gamma0 = 0.8, q = 2.1`; RMSE 1.387 GPa; chi-square/N = 0.961 | not applicable |
| Peritheos refit, `gamma0` fixed | 41 heated rows; pressure errors only | `q = 1.024 +/- 0.478` | 0.868 |
| Peritheos errors-in-variables refit, `gamma0` fixed | 41 heated rows; errors in `P,V,T` | `q = 1.050 +/- 0.519` | 0.767 |
| Peritheos errors-in-variables refit, `gamma0` fixed | all 51 rows; errors in `P,V,T` | `q = 1.089 +/- 0.521` | 17.483 |
| Peritheos errors-in-variables refit | 41 heated rows; errors in `P,V,T` | `gamma0 = 0.952 +/- 0.106, q = 2.772 +/- 1.266` | 0.724 |

The published curve is a good direct description of the heated observations:
for example, the 289.1 Å³, 2023 K row evaluates to 40.45 GPa, compared with
40.4(2.3) GPa observed. A refit does **not**, however, independently recover the
article's `q = 2.1` when `gamma0 = 0.8` is fixed.

The poor all-row errors-in-variables statistic is caused by the ten 300 K rows,
not by the heated measurements. Their stated pressure errors are all 0.1 GPa,
while the fixed quench-data reference EOS differs systematically by 0.49--2.49
GPa. These printed sigmas therefore cannot serve as a complete independent
error model for combining those rows with the separately fitted reference
isotherm. On the heated rows, freeing both `gamma0` and `q` gives a good fit,
but the parameters are strongly correlated (`r = 0.909`) and their
uncertainties are large. That result should be treated as a diagnostic, not a
replacement catalog parameterization.

### Berman thermal-reference-state fit

The article also reports an EosFit7c “Berman” fit. This required a source-code
check because Equation 3.8 in the article uses the exponential integral of
`alpha0 + alpha1*T`, while the open-source CrysFML engine used by EosFit7
implements Berman as the truncated polynomial

\[
V_0(T)=V_{00}\left[1+\alpha_0(T-T_r)
+\frac{1}{2}\alpha_1(T-T_r)^2\right].
\]

EosFit calls the exponential-integral form “Fei,” not “Berman.” Its default
reference temperature for the Berman model is 298 K. Peritheos now represents
the actual EosFit form explicitly as `BM3 + ThermalReferenceStateEOS` with
`reference_volume_law="berman"`.

The comparison was subsequently repeated with the authors' distributed
[EosFit7c console program](https://www.rossangel.com/text_eosfit.htm), version
7.60 dated 13 May 2021. The program read the embedded Table 4 transcription in
its native `FORMAT 1 T sigT P sigP V sigV` layout. It was configured as BM3 +
Berman 1988 + linear `dK/dT`, with `Tr = 298 K`, `V0 = 328.4 A^3`,
`K0 = 221 GPa`, `K0' = 3.3`, and `alpha1 = 0.0573e-8 K^-2` fixed. Only
`alpha0` and `dK/dT` were refined, and all three reported sigmas were enabled.

The fit objectives also differ. EosFit minimizes pressure residuals using the
iteratively updated effective variance

\[
\sigma_{P,\mathrm{eff}}^2=\sigma_P^2
+\sigma_V^2(K/V)^2+\sigma_T^2(\alpha K)^2,
\]

whereas Peritheos's standard fit treats observed volumes and temperatures as
latent coordinates and minimizes their corrections together with the pressure
residual. The reproduction script evaluates both methods.

| Parameter | Published | EosFit7c 7.60, 41 heated rows | Peritheos errors-in-variables |
|---|---:|---:|---:|
| `alpha0` (K^-1) | `1.94(16)e-5` | `(1.81120 +/- 0.21554)e-5` | `1.781(211)e-5` |
| `alpha1` (K^-2) | `0.0573e-8`, fixed | fixed | fixed |
| `dK0/dT` (GPa/K) | `-0.008(3)` | `-0.01311 +/- 0.00419` | `-0.01232(411)` |

The direct software run converged after three cycles with weighted chi-square
`1.0412` and an `alpha0`--`dK/dT` correlation of `-95.8%`. Peritheos's
effective-variance emulation gives `alpha0 = 1.81113e-5 K^-1`,
`dK/dT = -0.0131165 GPa/K`, and reduced chi-square `1.0397`; the tiny
differences are consistent with EosFit's printed precision and convergence
tolerance. This is an actual software match, not merely a match to a formula
transcribed from the documentation.

The article says that the “whole dataset” was fitted. Repeating the identical
EosFit7c refinement with all 51 rows gives `alpha0 = 1.79681e-5 K^-1`,
`dK/dT = -0.01305 GPa/K`, and weighted chi-square `17.9673`. The ten 300 K
rows have 0.1 GPa pressure uncertainties but disagree with the separately
fixed reference isotherm by as much as 2.49 GPa; they therefore dominate this
statistic while barely constraining the two thermal parameters.

Thus the original Peritheos comparison was not the same calculation, but the
matched Peritheos calculation now reproduces EosFit itself. Neither the
41-row nor the literal 51-row EosFit result reproduces the two coefficients in
published Table 2. The published curve gives reduced effective chi-square
2.008 on the same 41 heated rows. The remaining publication discrepancy cannot
be assigned from the public record. Possible causes are a different internal
data revision, higher-precision values used before the supplement was rounded,
or an undocumented row/weight selection. The published parameters remain
unchanged in `b4c_somayazulu_2023_berman_2`. The directly reproduced 41-row
result is also available as the explicitly opt-in
`b4c_somayazulu_2023_berman_refit`; its metadata records EosFit version, row
selection, weighting, fit statistics, parameter correlation, and the matched
Peritheos result. It is a derived public-data record, not a claim that these
coefficients appeared in the paper.

```python
from peritheos import get_material

b4c = get_material("b4c")
published = b4c.get_eos_record("b4c_somayazulu_2023_berman_2")
refit = b4c.get_eos_record("b4c_somayazulu_2023_berman_refit")
```

Run the complete comparison from the repository root with:

```bash
uv run python scripts/reproduce_somayazulu_2023_b4c.py
```

### Source reconciliation and structure

The library keeps the final article's `gamma0 = 0.8, q = 2.1`, and separately
retains the published Berman parameterization. The matched Berman refit is a
third, explicitly derived record selected only by its `_refit` identifier. The
material record documents five source
conflicts: BM2 versus BM3 wording, 1425 versus 1450 K, different
`gamma0,q` pairs in the article and supplement, nominal 300 K versus the
explicit 298 K calculation reference, and the article's displayed exponential
thermal-volume equation versus EosFit's Berman implementation.

Diffraction provenance is tracked separately. The material uses the average
disordered rhombohedral B4C model in hexagonal `R-3m` (#166), with five
asymmetric occupied sites. Their multiplicities and split occupancies expand to
B36C9, or nine B4C formula units. The supplement's ambient `a = 5.601 Å` and
`c = 12.087 Å` give 328.383 Å³, matching the EOS reference cell volume.
The average atomic model follows Clark and Hoard (1943),
[doi:10.1021/ja01251a026](https://doi.org/10.1021/ja01251a026).

<a id="kcl-walker-2002"></a>

## KCl: Walker et al. (2002)

### Why the first validation failed

[Walker et al. (2002)](https://doi.org/10.2138/am-2002-0701) report two B2-KCl
solutions in Table 3. The preferred bold row and the italic simultaneous-fit
row are different regressions, not alternative roundings of one result. The
Table 3 footnote says that the preferred result fits `K0` and `K0'` to the
room-temperature data before fitting the thermal coefficient. The italic row
fits all parameters simultaneously.

The first Peritheos campaign released `V0`, `K0`, `K0'`, and `alpha_KT`
together against all 39 P-V-T rows, used the printed pressure and volume ESDs,
and then compared that result with the preferred staged row. It therefore
changed both the staging and the objective. The resulting
`V0 = 56.811 A^3`, `K0 = 9.482 GPa`, and `K0' = 9.629` looked alarming, but
they were not evidence that the stored preferred EOS was incorrectly evaluated.

This distinction matters especially here. B2 KCl is observed only from 3.18 to
8.14 GPa, whereas `V0` is a fictive zero-pressure reference for a high-pressure
phase. Walker et al. explicitly warn that the parameters are strongly correlated
and decline to assign meaningful individual errors. Their own simultaneous
solution (`V0 = 55.25 A^3`, `K0 = 14.8 GPa`, `K0' = 6.9`, and
`alpha0 = 0.00018 K^-1`) demonstrates the size of that covariance.

### Source-protocol reproduction

The corrected validation minimizes the unweighted squared pressure residuals
specified by the article. It uses two stages:

1. hold the preferred fictive `V0 = 53.53 A^3` fixed and fit `K0` and `K0'`
   to the eight Table 2 observations at 23-24 degC; and
2. hold that fitted reference isotherm fixed and fit `alpha_KT` to all 39
   observations from 23 to 600 degC using Equation BE1.

| Parameter | Preferred Table 3 | Peritheos staged refit | Difference |
|---|---:|---:|---:|
| `V0` (A^3) | 53.53, held | 53.53, held | 0% |
| `K0` (GPa) | 23.7 | 23.7754 | 0.32% |
| `K0'` | 4.4 | 4.41587 | 0.36% |
| `alpha_KT` (GPa/K) | 0.00275 | 0.00276625 | 0.59% |

The room-temperature stage has a pressure RMSE of 0.0261 GPa. Across all 39
rows, the printed preferred coefficients have an RMSE of 0.1046 GPa and the
staged refit has an RMSE of 0.1035 GPa. This is curve-level and coefficient-level
agreement to the precision of the printed table. The campaign labels the result
`similar`, rather than strict uncertainty `parity`, only because the source
withholds individual errors for `K0` and `K0'`.

For diagnosis, the ledger also retains the correctly unweighted all-free
four-parameter calculation. It independently reproduces the paper's italic
simultaneous solution:

| Parameter | Italic Table 3 | Peritheos simultaneous refit | Difference |
|---|---:|---:|---:|
| `V0` (A^3) | 55.25 | 55.3923 | 0.26% |
| `K0` (GPa) | 14.8 | 14.1887 | 4.13% |
| `K0'` | 6.9 | 7.15656 | 3.72% |
| `alpha_KT` (GPa/K) | 0.002664, derived from `alpha0*K0` | 0.00269692 | 1.24% |

That simultaneous fit has a 0.0816 GPa pressure RMSE. It is still not the
protocol used for the preferred row, but its agreement with the separate italic
row confirms the diagnosis from the opposite direction. The material record
remains the preferred bold Table 3 EOS, now with `V0` marked as fixed and with
both fit protocols described explicitly.

<a id="kcl-tateno-2019"></a>

## KCl: Tateno et al. (2019)

### Two source-control errors

The first validation of
[Tateno et al. (2019)](https://doi.org/10.2138/am-2019-6779) failed for two
independent reasons.

First, the record had been audited against the accepted manuscript. That
version reports `gamma0 = 0.58(5)` and `q = 0.9(2)`, but the final published
article reports `gamma0 = 2.3(2)` and `q = 0.8(2)` for the preferred Sokolova-Pt
fit. Final Equation 6 defines

\[
\theta(V)=\theta_0\exp\{[\gamma_0-\gamma(V)]/q\},
\]

which is Peritheos's `integrated_gruneisen` law, not `variable_exponent`.
The final article's approximately 10 GPa thermal pressure at 3000 K is also
consistent with `gamma0 = 2.3`, not 0.58.

Second, the accepted-manuscript table placed the Pt pressure/temperature half
and the KCl-volume half on separate pages in different row orders. Joining them
by printed position corrupted runs 3 and 4. For example, the 6.0 GPa, 300 K
observation was incorrectly paired with `V = 33.8821 A^3`; its actual KCl
volume is `45.0463 A^3`.

The authoritative observations are the 39 rows in the
[official MSA deposit AM-19-56779](http://www.minsocam.org/MSA/AmMin/TOC/2019/May2019_data/AM-19-56779.zip),
workbook `6779TableS1 revised.xlsx`. The checked-in CSV now follows that
workbook's row alignment and retains the reported Pt pressures, temperatures,
volumes, and uncertainties.

### Corrected reproduction

The corrected errors-in-variables fit uses all 39 Supplemental Table S1 rows.
It holds the source-fixed `V0 = 54.5 A^3`, `theta0 = 235 K`, `Tr = 300 K`, and
`n = 2`, and fits the remaining Vinet and MGD coefficients simultaneously.

| Parameter | Final publication | Peritheos refit | Difference |
|---|---:|---:|---:|
| `K0` (GPa) | 18.3(3) | 18.3446(362) | 0.24% |
| `K0'` | 5.60(3) | 5.60096(492) | 0.02% |
| `gamma0` | 2.3(2) | 2.29519(278) | 0.21% |
| `q` | 0.8(2) | 0.82490(283) | 3.11% |

Every fitted coefficient agrees within the combined two-standard-deviation
uncertainty and the numerical similarity threshold. The final published curve
has a 0.726 GPa pressure RMSE on the rounded deposited rows; the refit gives
0.569 GPa with reduced chi-square 0.431. The record therefore moves from
`parity_not_achieved` to `parity`. The large earlier residuals were entirely
explained by the superseded thermal parameters, wrong Debye law, and mispaired
table rows; no remaining EOS implementation discrepancy is evident.

<a id="kcl-chidester-2021"></a>

## KCl: Chidester et al. (2021)

### Why the first refit drove `q` to zero

The first validation did not reproduce the fit described by
[Chidester et al. (2021)](https://doi.org/10.1103/PhysRevB.104.094107). It used
only the 155 rows in the author-deposited high-temperature table. The paper,
however, states that those new observations were fitted simultaneously with
the room-temperature B2-KCl data of
[Dewaele et al. (2012)](https://doi.org/10.1103/PhysRevB.85.214105). Figure 3
also plots both sets and describes its residuals as belonging to the fit of all
the data.

That omission removed all observations below 26.2 GPa while leaving the
fictive zero-pressure `V0` and the other four coefficients free. Applying the
new-table temperature, pressure, and volume uncertainties as an
errors-in-variables objective further changed the regression from the
source-compatible unweighted pressure fit. In that under-anchored
calculation, `V0` fell from `53.1373` to `51.8022 A^3` and `q` reached its
lower bound. The boundary was therefore a fit-scope and weighting artifact,
not evidence that the published `q = 1.0(1)` is wrong.

The article does not state a detailed weighting formula. The Dewaele source
table has no row-wise uncertainties, and the unweighted calculation recovers
both the complete coefficient set and the reported RMSE. The validation
therefore records the unweighted objective as a reproduction-supported
inference rather than an explicitly printed method.

The corrected material record names both source inputs explicitly:

- all 123 rows of `kcl_dewaele_2012_table1_compression`, assigned to the
  300 K reference isotherm; and
- all 155 rows of `kcl_chidester_2021_supplemental_pvt`, which retain the
  authors' effective KCl temperatures and Pt-derived pressures.

The high-temperature source file is the author-deposited
[`SuppTable_KCl.csv`](https://knowledge.uchicago.edu/records/6t4wc-w2146/files/SuppTable_KCl.csv?download=1).
It contains the reduced KCl P-V-T states, but not the raw Pt cell volumes;
therefore the published pressure basis can be reproduced but cannot be
independently re-reduced observation by observation.

### Debye-temperature convention

Equations 3-4 print

\[
P_{\mathrm{thermal}}=\frac{\gamma(V)}{V}
 [E_D(V,T)-E_D(V,300\,\mathrm{K})],
\qquad
\gamma(V)=\gamma_0(V/V_0)^q,
\]

but do not separately print the volume dependence of the Debye temperature.
The corrected record uses the thermodynamically consistent integral of the
stated constant-`q` Gruneisen law,

\[
\theta(V)=\theta_0\exp\{[\gamma_0-\gamma(V)]/q\}.
\]

This convention is not chosen only by preference: with the complete data scope
it recovers the printed coefficients, including `gamma0` and `q`, whereas the
previous direct variable-exponent convention gives `gamma0 = 3.104` and
`q = 1.130` on the same rows. The reproduction therefore supplies the missing
operational detail while recording that the article itself leaves
`theta(V)` implicit.

### Corrected reproduction

The corrected calculation minimizes unweighted pressure residuals over all
278 observations. It fixes only the source-defined `Tr = 300 K`,
`theta0 = 235 K`, and `n = 2`, and fits the BM3 reference and MGD thermal
coefficients simultaneously.

| Parameter | Publication | Peritheos refit | Difference |
|---|---:|---:|---:|
| `V0` (A^3) | 53.1373 ± 0.4982, converted from 32.0(3) cm^3/mol | 53.2036 ± 0.4752 | 0.12% |
| `K0` (GPa) | 24 ± 1 | 23.9721 ± 1.1961 | 0.12% |
| `K0'` | 4.56 ± 0.05 | 4.55798 ± 0.04694 | 0.04% |
| `gamma0` | 2.9 ± 0.4 | 2.91714 ± 0.36014 | 0.59% |
| `q` | 1.0 ± 0.1 | 0.96524 ± 0.16751 | 3.48% |

Every coefficient meets both the numerical criterion and the combined
two-standard-deviation uncertainty test. Across all 278 rounded observations,
the published and refitted curves have pressure RMSE values of 1.303 and
1.263 GPa. On the 155 high-temperature rows alone, the refit gives
1.582 GPa, reproducing the article's reported 1.6 GPa. Chidester therefore
moves from `similar` with a boundary solution to full `parity`; there is no
remaining coefficient-level discrepancy after restoring the source's actual
data scope and model convention.
