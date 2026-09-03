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
