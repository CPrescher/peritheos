# Criniti et al. (2023): Fe-bearing Al-phase D

## Scope and source authority

This record represents the room-temperature, pre-crossover compression of the
single-crystal Fe-bearing Al-phase D sample reported by Criniti et al. (2023),
*American Mineralogist* **108**, 1764–1772,
[doi:10.2138/am-2022-8559](https://doi.org/10.2138/am-2022-8559).
The audit used the final article, American Mineralogist Online Materials deposit
AM-23-98559, and the authors' [deposited CIF](https://figshare.com/s/57ea586081085cf9b9ac).
LitCurate was used only to discover the candidate.

The analyzed sample composition is
Al~1.53(2)~Fe~0.22(1)~Si~0.86(1)~O~6~H~3.33(9)~. Mössbauer spectroscopy assigns
the iron to Fe^3+^. This composition is not merged with the existing Mg-rich
`phase_d` material.

## Structure and volume basis

The author-provided CIF block `data_AlPhaseD_P6322` reports:

- space group P6~3~22 (number 182);
- *a* = 4.74653(13) Å, *c* = 4.29002(19) Å;
- *V* = 83.703(6) Å³ at 293(2) K; and
- *Z* = 1.

Consequently, every EOS volume is the conventional hexagonal-cell volume in
Å³, numerically also the volume per analyzed O~6~ formula unit. The deposited
X-ray scattering model sums to Al~1.506~Fe~0.217~Si~0.847~O~6~. Hydrogen was
measured analytically but is absent from the X-ray refinement, so Peritheos does
not invent an H coordinate.

The CIF also contains a P-31m refinement. The article says Hamilton tests on
merged and unmerged data prefer different models, but interprets the diffuse
scattering and presents P6~3~22 in the abstract and conclusions. The material
therefore uses the P6~3~22 model while recording the unresolved symmetry
comparison.

## Published equations and canonical choice

The paper fits the 0–38 GPa data with EosFit7c/EosFit7-GUI using the standard
third-order Birch–Murnaghan equation

\[
P(V)=\frac{3K_0}{2}\left[\eta^7-\eta^5\right]
\left\{1+\frac{3}{4}(K'_0-4)[\eta^2-1]\right\},
\qquad \eta=(V_0/V)^{1/3},
\]

and the standard Vinet equation

\[
P(V)=3K_0\frac{1-x}{x^2}
\exp\left[\frac{3}{2}(K'_0-1)(1-x)\right],
\qquad x=(V/V_0)^{1/3}.
\]

The freely fitted Table 2 rows are preserved as separate executable records:

| Record | V~0~ (Å³) | K~0~ (GPa) | K′~0~ | Role |
|---|---:|---:|---:|---|
| `fe_bearing_al_phase_d_criniti_2023_bm3_1` | 83.68(2) | 166.3(15) | 4.46(12) | canonical |
| `fe_bearing_al_phase_d_criniti_2023_vinet_2` | 83.68(2) | 165.5(15) | 4.62(12) | published alternative |

BM3 is canonical because the abstract reports its coefficients, Figure 3 shows
only BM3 curves, and Table 3 compares the result with earlier phase-D studies
using BM equations. Vinet remains executable because the authors explicitly fit
it and state that its curve overlaps BM3 between 0.0001 and 38 GPa. The two
ambient-volume-fixed variants of each family are documented Table 2 sensitivity
fits, not additional canonical records.

Parenthetical coefficient errors are retained as published uncertainties. The
paper does not state their confidence convention, provide parameter covariance,
or print a fit statistic; none is inferred.

## Primary observations and refit

All 26 Table 1 rows are transcribed in
`fe-bearing-al-phase-d-criniti-2023-table1-compression.csv`. The table includes
run number, pressure and its uncertainty, *a*, *c*, *V*, and their printed
uncertainties. The 21 rows at or below 38 GPa are the ordinary BM3/Vinet fit
selection. The five rows from 40.75 to 52.41 GPa are retained but excluded from
these records because the article observes changed compression behavior and
models them only with a separate spin-crossover formalism.

An independent iterative effective-variance fit used

\[
\sigma_{P,\mathrm{eff}}^2=\sigma_P^2+
\left(\frac{\partial P}{\partial V}\sigma_V\right)^2
\]

and varied V~0~, K~0~, and K′~0~ simultaneously. Results are:

| Family | Refit V~0~ (Å³) | Refit K~0~ (GPa) | Refit K′~0~ | Published-curve RMSE (GPa) | Refit RMSE (GPa) |
|---|---:|---:|---:|---:|---:|
| BM3 | 83.68041718 | 166.25494346 | 4.45597358 | 0.12637766 | 0.12603947 |
| Vinet | 83.68347936 | 165.47809134 | 4.62083697 | 0.12375790 | 0.12335568 |

The BM3 refit reproduces every coefficient within its printed rounding interval.
The Vinet refit reproduces K~0~ and K′~0~ within their rounding intervals and V~0~
well inside its 0.02 Å³ reported uncertainty. The refits are validation
diagnostics, not replacement EOS records.

## Pressure calibration and validity

Both DAC runs used helium as pressure medium and the Shen et al. (2020)
IPPS-Ruby2020 fluorescence scale. Each Table 1 pressure is the mean of ruby
measurements before and after the diffraction scan; the tabulated σ~P~ is their
semi-difference. The two row-wise R1 wavelengths are not published, so exact
pressure reduction cannot be repeated even though the calibration reference is
unambiguous.

The executable records are restricted to 0–38 GPa and 293(2) K. Above 38 GPa
the data depart from the ordinary compression curve. The authors model the
broad Fe^3+^ high-to-low spin crossover over approximately 30–65 GPa, so this
range is not a rectangular high-spin phase-stability guarantee.

## Source inconsistencies retained

- Table 2 labels the 0–38 GPa BM3/Vinet group “Low-spin state EOS,” while the
  prose treats it as the behavior before the high-to-low spin crossover and
  places the distinct spin-crossover fits below it.
- The official AM-23-98559 supplemental PDF has the correct deposit number and
  phase-D Figures S1–S2, but its title line incorrectly names the adjacent
  hydrous Al-bearing silica article.

Neither inconsistency changes the reproducible room-temperature coefficients,
data selection, or volume basis.
