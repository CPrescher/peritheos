# Xu et al. (2020): Al-bearing phase D

## Scope and source authority

This record represents the 300 K compression of the polycrystalline Al-bearing
phase-D sample reported by Xu et al. (2020), *Geophysical Research Letters*
**47**, e2020GL088877,
[doi:10.1029/2020GL088877](https://doi.org/10.1029/2020GL088877). The audit used
the final free-access publisher article and the corresponding author's complete
numerical P-V-T workbook in EarthChem,
[doi:10.26022/IEDA/111629](https://doi.org/10.26022/IEDA/111629). The workbook's
published SHA-1 is `cddfb0e859c36847ec75e776495ad8ef7f95b456`; the retrieved
file has SHA-256
`1349bc1260337c1e8a3de6830e089b4ad5c32c18b3aea69f0ff6986a9b9c4d3e`.

The paper reports 20.7 wt% MgO, 44.5 wt% SiO~2~, and 18.8 wt% Al~2~O~3~ by
EDS, with about 16.0 wt% H~2~O estimated from the analytical total deficit.
Converting those values to an O~6~ basis gives the approximate formula
Mg~0.90~Al~0.64~Si~1.29~H~3.10~O~6~. It is an audit label, not a direct
site-occupancy or hydrogen analysis. This composition is kept separate from
the Mg-rich `phase_d` and Fe-bearing Al-phase-D materials.

## Structure and volume basis

The article and EarthChem Table 2 report the ambient lattice as
*a* = 4.8077(4) Å, *c* = 4.3218(6) Å, and *V* = 86.51(1) Å³. The paper assigns
phase D to P-31m. The conventional trigonal/hexagonal cell has Z = 1, so its
cell volume is also the volume of one approximate O~6~ formula unit.

There is a source inconsistency at the reference state. EarthChem Table 3 begins
with a distinct zero-pressure, 300 K measurement *a* = 4.8078(6) Å,
*c* = 4.3316(2) Å, and *V* = 86.71(2) Å³. Section 3.1 explicitly fixes this
86.71 Å³ value as the EOS V~0~. Peritheos therefore retains the 86.51 Å³
ambient lattice as the material characterization and the 86.71 Å³ in situ value
as the executable EOS reference volume instead of silently merging them.

## Published equation and fit selection

Section 3.1 says that all room-temperature observations through 20.5 GPa were
fitted by least squares with EosFit to a standard third-order Birch-Murnaghan
equation,

\[
P(V)=\frac{3K_0}{2}\left[\eta^7-\eta^5\right]
\left\{1+\frac{3}{4}(K'_0-4)[\eta^2-1]\right\},
\qquad \eta=(V_0/V)^{1/3}.
\]

The author-deposited workbook exposes exactly eight 300 K states from 0 to
20.53 GPa. The represented coefficients are:

| Parameter | Value | Status |
|---|---:|---|
| V~0~ | 86.71(2) Å³ | fixed to the in situ zero-pressure row |
| K~T0~ | 143(5) GPa | fitted |
| K′~T0~ | 5.8(7) | fitted |

The parenthetical errors have no stated confidence convention, and no covariance
matrix or scalar fit statistic is published. Table 1 prints 143(4) GPa for the
same 300 K K~T0~, whereas the Section 3.1 sentence that directly describes the
fit prints 143(5) GPa. The record stores the narrative uncertainty and preserves
both values as a reported inconsistency.

Table 1 labels the represented relation “Gold P-scale, Tsuchiya, 2003,” referring
to the Au P-V-T EOS in
[doi:10.1029/2003JB002446](https://doi.org/10.1029/2003JB002446). The methods
also say that Au and NaCl standards were used across the broader experiment,
but the Table 1 heading unambiguously selects the Tsuchiya gold scale for this
BM3. Row-wise Au volumes are not published, so the source pressures can be
refitted but cannot be independently reduced again from raw calibrant
observations.

## Complete P-V-T observations and refit

All 28 numerical states in EarthChem Table 3 are transcribed in source order in
`al-bearing-phase-d-xu-2020-earthchem-pvt.csv`. The CSV preserves pressure,
temperature, *a*, *c*, *V*, and every printed row-wise uncertainty. Eight 300 K
rows are flagged for the represented BM3; the other 20 rows remain available as
high-temperature provenance without implying an unverified thermal equation.

The article states that 33 diffraction patterns were collected. It explicitly
excludes the late-cycle measurements after partial dehydration near 20 GPa and
above 1400 K. Those five patterns are not published numerically. The 28-row
author deposit is therefore the complete numerical set retained for the source's
subsequent models, not an attempt to reconstruct the rejected data.

An independent iterative effective-variance fit combined pressure and volume
errors as

\[
\sigma_{P,\mathrm{eff}}^2=\sigma_P^2+
\left(\frac{\partial P}{\partial V}\sigma_V\right)^2,
\]

fixed V~0~ = 86.71 Å³, and varied K~0~ and K′~0~. It gives
K~0~ = 143.07336748 GPa and K′~0~ = 5.86138397. The published and refitted
curves have pressure RMSE values of 0.12191813 and 0.11079143 GPa,
respectively. The refit differs from the published values by only 0.08 GPa and
0.07, far inside the printed uncertainties.

## Why the broader models are not encoded

The 300 K BM3 is exactly representable and is the only executable record added.
The paper's other parameterizations fail the same exactness test:

- The prose says a high-temperature fit through 1300 K gives K~0~ = 143(4) GPa,
  K′~0~ = 5.8(6), dK~0~/dT = -0.027(12) GPa/K, and
  α~0~ = 3.8(4) × 10^-5^ K^-1^. Table 1 instead presents only 300–900 K
  variants with different coefficients, fixed/free choices, and a nonzero
  quadratic expansivity coefficient. The source does not choose one exact
  thermal form.
- The pressure-scale-free acoustic finite-strain section prints several
  alternative fixed/free coefficient sets and delegates its equations to Text
  S2. The EarthChem archive contains P-V-T lattice observations, not the full
  acoustic Table S1, and the article does not designate a unique acoustic set as
  an absolute P-V-T relation.

Encoding either broader model would require choosing among conflicting source
variants or inventing a mapping. They are documented but intentionally left
unrepresented.
