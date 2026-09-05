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

<a id="casio3-caracas-2005"></a>

## CaSiO3 perovskite: Caracas et al. (2005)

### Source identity and equation convention

The primary source is Caracas, Wentzcovitch, Price, and Brodholt,
*CaSiO3 perovskite at lower mantle pressures*, Geophysical Research Letters
**32**, L06306 (2005),
[doi:10.1029/2004GL022144](https://doi.org/10.1029/2004GL022144). The audit used
the publisher HTML and the UCL author-deposited five-page article. Neither page
lists supporting information or a data attachment.

This is a static-lattice, zero-Kelvin first-principles study, not an experimental
pressure calibration. Section 4 defines the Eulerian strain

\[
f=\frac{1}{2}\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]
\]

and identifies the Table 2 fits as third- and fourth-order Birch--Murnaghan.
In Peritheos notation their pressure forms are

\[
P_{\mathrm{BM3}}=3K_0f(1+2f)^{5/2}
\left[1+\frac{3}{2}(K'_0-4)f\right]
\]

and

\[
P_{\mathrm{BM4}}=3K_0f(1+2f)^{5/2}
\left[1+\frac{3}{2}(K'_0-4)f
+\frac{3}{2}\left(K_0K''_0+(K'_0-4)(K'_0-3)+\frac{35}{9}\right)f^2\right].
\]

Thus a BM3 row does not have a missing free `K0_double_prime`: truncation fixes
its implied value to

\[
K''_0=-\frac{(K'_0-4)(K'_0-3)+35/9}{K_0}.
\]

That derived value must not be confused with the independently fitted and much
less negative `K0_double_prime` printed in the BM4 block.

### Resolution of all Table 2 rows

LitCurate's repeated rows correspond to 18 distinct source parameterizations:
nine crystallographic distortions, each fitted once as BM3 and once as BM4.
Volumes are the paper's `V0/Z`, in A^3 per CaSiO3 formula unit. The final BM3
column below is derived from the truncation identity above; every BM4 `K0''`
is directly source-reported in GPa^-1.

| Structure | Glazer tilt | BM3 `V0`, `K0`, `K0'` | BM3 implied `K0''` | BM4 `V0`, `K0`, `K0'`, `K0''` |
|---|---|---|---:|---|
| Pm-3m | a0a0a0 | 44.579, 250, 4.098 | -0.015985972 | 44.588, 248, 4.206, -0.002 |
| I4/mcm | a0a0c- | 44.537, 249, 4.090 | -0.016012004 | 44.547, 247, 4.213, -0.002 |
| Imma | a0b-b- | 44.567, 249, 4.094 | -0.016031024 | 44.576, 247, 4.218, -0.002 |
| R-3c | a-a-a- | 44.821, 247, 4.100 | -0.016189834 | 44.832, 244, 4.236, -0.002 |
| P4/mbm | a0a0c+ | 44.629, 247, 4.124 | -0.016308765 | 44.641, 244, 4.261, -0.002 |
| I4/mmm | a0b+b+ | 44.599, 250, 4.103 | -0.016009992 | 44.609, 248, 4.219, -0.002 |
| Im-3 | a+a+a+ | 44.600, 250, 4.104 | -0.016014820 | 44.610, 247, 4.229, -0.002 |
| P42/nmc | a+a+c- | 44.576, 248, 4.092 | -0.016086100 | 44.566, 251, 3.977, -0.001 |
| Pnma | a-b+a- | 44.576, 249, 4.104 | -0.016079136 | 44.588, 246, 4.248, -0.002 |

These are not interchangeable composition variants. Pm-3m is the ideal cubic
parent; the other eight rows represent distinct tilted structures. The paper
finds Pm-3m dynamically unstable throughout the investigated static pressure
range and identifies I4/mcm as the lowest-energy static configuration. It does
not publish the optimized atomic coordinates required for new diffraction-ready
low-symmetry material records.

### Selected record and numerical reproduction

The one executable addition is `ca_perovskite_caracas_2005_bm3_3`, the Pm-3m
BM3 fit (`V0 = 44.579 A^3`, `K0 = 250 GPa`, `K0' = 4.098`). It belongs in the
existing cubic `ca_perovskite` material because Pm-3m has one formula unit in
its primitive conventional cell, so the source's `V0/Z` is directly the public
cell volume. The low-symmetry rows do not belong in that material. The cubic
BM4 row is retained as explicit alternative-fit metadata rather than being
misrepresented as a duplicate or added as a second record in this focused
change.

As an independent check, inversion of the published BM3 equation at 130 GPa
gives `V = 33.436317 A^3`. Using the CaSiO3 formula mass and one formula unit per
cell gives `rho = 5.76885 g/cm^3`, reproducing the paper's separately stated
`5.77 g/cm^3` cubic density. The zero-pressure conversion gives
`4.32690 g/cm^3`, consistent with the paper's truncated `4.32 g/cm^3` statement.

The underlying first-principles energy-volume points, coefficient uncertainties,
fit weights, residuals, and covariance are not published. Table 2's stated
calculation uncertainty of about 20 meV per molecule describes the energy scale,
not an uncertainty on `V0`, `K0`, `K0'`, or `K0''`; the coefficient errors remain
explicitly unavailable. Consequently no source-faithful coefficient refit is
possible; the refit ledger classifies this record as
`theoretical_parameterization_only` and `not_refittable` rather than digitizing
or synthesizing observations.
## Stishovite, Wang et al. (2012)

### Source, observations, and calibration

- Primary article: F. Wang, Y. Tange, T. Irifune, and K. Funakoshi,
  *P-V-T equation of state of stishovite up to mid-lower mantle conditions*,
  *Journal of Geophysical Research: Solid Earth* **117**, B06209 (2012),
  [doi:10.1029/2011JB009100](https://doi.org/10.1029/2011JB009100).
- Primary data: all 56 rows of Table 1, spanning 16.85--54.5 GPa and
  300--1700 K. The bundled CSV preserves the stishovite lattice parameters and
  conventional-cell volumes, simultaneous Au volumes, and all printed
  parenthetical estimated standard deviations. The publisher also lists
  tab-delimited supplements for Tables 2--4; these contain fitted parameters
  and calculated P-V-T/thermoelastic grids rather than additional observations.
- Pressure calibration: each pressure was calculated from simultaneous Au
  volume and temperature with the Tsuchiya (2003) Au thermal EOS,
  [doi:10.1029/2003JB002446](https://doi.org/10.1029/2003JB002446). The paired
  calibrant observations are preserved, but that Au EOS is not yet executable
  in Peritheos, so exact pressure re-reduction is recorded as pending rather
  than approximated with a different gold scale.

The sample is pure, Al-free SiO2 stishovite. The diffraction cell is the
existing rutile-type tetragonal `P42/mnm`, `Z = 2` structure in
`sio2_stv_andr.eosmat`; therefore the Table 1 conventional-cell volume is the
model volume without a crystallographic conversion.

### Model choice and numerical reproduction

The paper fits one Mie--Gruneisen--Debye thermal formulation with two alternative
300 K reference curves. The implemented record uses the Vinet curve because
Figure 2 plots it and Tables 3--4 list it first. The authors do not claim that
it is statistically preferred: they say both Vinet and BM3 reproduce the data
and agree within analytical uncertainty. To keep this audit to one executable
literature record, the complete BM3 coefficients are retained in the Vinet
record's notes.

For `x = V/V0`, the stored thermal law is exactly equations (4)--(8), with

\[
\gamma=\gamma_0\{1+a(x^b-1)\},\qquad
\Theta=\Theta_0\exp\left[-\gamma_0\left((1-a)\ln x+
\frac{a}{b}(x^b-1)\right)\right],
\]

and a Debye thermal-pressure increment referenced to 300 K. The atom count is
`n = 3` per SiO2 formula unit. Table 2 fixes `V0 = 46.55 Å³` and `a = 1`, and
reports the Vinet solution `K0 = 292(2) GPa`, `K0' = 5.01(12)`,
`theta0 = 1130(110) K`, `gamma0 = 1.67(7)`, and `b = 3.0(4)`. The alternative
BM3 solution is `K0 = 294(2) GPa`, `K0' = 4.85(12)`,
`theta0 = 1130(100) K`, `gamma0 = 1.66(7)`, and `b = 2.9(4)`, with the same
fixed `V0` and `a`.

The fixed `V0` is the rounded form of the separately reported recovered-sample
measurement `46.553(19) Å³`; it has no fit uncertainty because it was held
fixed, although the underlying ambient measurement does have an uncertainty.

The published Vinet coefficients reproduce the Table 3 values: at
`V/V0 = 0.98` and 1000 K Peritheos gives 11.103 GPa versus 11.11 GPa, and at
`V/V0 = 1` and 3000 K it gives 22.619 GPa versus 22.62 GPa. Across all 56
observations the printed curve has a pressure RMSE of 0.213 GPa. An independent
errors-in-variables refit using the printed pressure and volume uncertainties,
with `V0`, `Tr`, `a`, and `n` fixed, gives `K0 = 294.68 GPa`, `K0' = 4.879`,
`theta0 = 1134.7 K`, `gamma0 = 1.6555`, and `b = 3.0469`; every coefficient is
within the combined two-standard-deviation interval of the published value.

Several measurements lie beyond the cited stishovite--CaCl2-type boundary, but
no transition was observed on the experimental timescale. The authors report
that excluding those points leaves the fit unchanged. Consequently the stored
P-T bounds describe the observations, not a rectangular phase-stability field.

## Corundum and Rh2O3(II)-type Al2O3: Shi et al. (2022)

### Source classification and equation

Shi et al., *Thermal Equations of State of Corundum and Rh2O3 (II)-Type
Al2O3 up to 153 GPa and 3400 K*, Journal of Geophysical Research: Solid Earth
**127**, e2021JB023805 (2022),
[doi:10.1029/2021JB023805](https://doi.org/10.1029/2021JB023805), is the
primary EOS source. Its [official Zenodo deposit](https://doi.org/10.5281/zenodo.5771198)
contains the Supporting Information under CC BY 4.0.

The LitCurate discovery result is not a scientifically correct record model.
The paper does not report six interchangeable, isothermal BM3 equations. It
uses the following full Mie--Gruneisen--Debye equation, with a 300 K BM3 cold
term:

\[
P(V,T)=P_{300}(V)+\frac{\gamma(V)}{V}
\left[E_D(V,T)-E_D(V,300\ \mathrm{K})\right],
\]

\[
\gamma(V)=\gamma_0(V/V_0)^q,\qquad
\theta(V)=\theta_0\exp[-(\gamma(V)-\gamma_0)/q].
\]

Equation 2 is algebraically third-order Birch--Murnaghan even though the prose
calls it Murnaghan. Every fit fixes `K0' = 4`, so the cold curve numerically
reduces to BM2. Equation 4 uses `n = 5` atoms per Al2O3 formula unit. This maps
exactly to `BM3 + MieGruneisenDebye` with
`debye_temperature_law: integrated_gruneisen`.

The apparent multiple records are sensitivity trials. Main Tables 2 and 3
contain three corundum trial-Debye-temperature fits and five Rh2O3(II)
trial-`q` fits:

| Phase/model | `K0` (GPa) | `K0'` | `V0` (A^3/cell) | `theta0` (K) | `gamma0` | `q` |
|---|---:|---:|---:|---:|---:|---:|
| corundum M1 | 245(1) | 4 fixed | 255.1 fixed | 500 fixed | 1.21(6) | 0.7(3) |
| corundum M2 | 245(1) | 4 fixed | 255.1 fixed | 800 fixed | 1.26(6) | 0.7(3) |
| corundum M3, selected | 246(1) | 4 fixed | 255.1 fixed | 1100 fixed | 1.32(7) | 0.8(4) |
| Rh2O3(II) M1 | 253(6) | 4 fixed | 165.5(7) | 600(200) | 1.33(5) | 0.6 fixed |
| Rh2O3(II) M2, selected | 256(6) | 4 fixed | 165.2(7) | 600(200) | 1.47(5) | 1.0 fixed |
| Rh2O3(II) M3 | 259(6) | 4 fixed | 164.8(7) | 600(200) | 1.62(6) | 1.4 fixed |
| Rh2O3(II) M4 | 262(6) | 4 fixed | 164.5(6) | 600(200) | 1.79(7) | 1.8 fixed |
| Rh2O3(II) M5 | 262(6) | 4 fixed | 164.2(6) | 500(200) | 1.97(7) | 2.2 fixed |

Main Table 1 and Supporting Information Table S3 select corundum M3 and
Rh2O3(II) M2 for downstream calculations. Table 1 prints a corundum `K0`
error of 2 GPa while Table 2 prints 1 GPa for the same central value; this
internal discrepancy is one reason not to materialize the sensitivity table
as a set of equivalent EOS records. Peritheos adds exactly two records from
this article: selected corundum M3 in `alumina.eosmat` and selected
Rh2O3(II) M2 in its phase-specific material. The other six rows remain
documented sensitivity trials rather than executable alternatives.

### Identity, data, and calibration

Rh2O3(II)-type alumina is a distinct orthorhombic `Pbcn` (#60) polymorph, so it
does not belong in the existing corundum `alumina.eosmat` material. The
separate `alumina_rh2o3_ii.eosmat` structure follows Lin et al. (2004) at
113 GPa and 300 K and contains four Al2O3 formula units. Consequently the EOS
`V0 = 165.2(7) A^3` is a conventional-cell value, equivalent to
41.30 A^3/formula unit and 24.871 cm^3/mol. It is a high-pressure fit
extrapolated to zero pressure, not an observed ambient cell.

All 75 corundum P--T--V rows from Supporting Information Table S1 are bundled
in `alumina-shi-2022-table-s1-pvt.csv`, including every printed pressure,
temperature, and volume uncertainty. The table groups 59 rows under Pt
pressures, including the fixed ambient reference row, and 16 under NaCl
pressures. The dataset keeps this calibrant identity. The paper reports that
NaCl and Pt pressures agree within 1 GPa from 35 to 67 GPa up to 3000 K, but
neither row-wise calibrant volumes nor a separate numerical NaCl
parameterization are published. Calibrant-level recalculation is therefore
unavailable even though the reduced sample pressures are complete.

All 75 Pt-calibrated Rh2O3(II) P--T--V rows from Supporting Information Table
S2 are bundled in `alumina-rh2o3-ii-shi-2022-table-s2-pvt.csv`, including every
printed pressure, temperature, and volume uncertainty. The source does not
state an uncertainty confidence convention or publish a coefficient
covariance matrix; errors from both tables are therefore retained as generic
uncertainties. No row-wise Pt lattice parameters are published. Pressures use
the self-consistent Fei et al. (2007) Pt thermal scale, so observation-level
recalculation is not possible from the supplement alone.

The article describes the EOS as extending to 153 GPa and 3400 K. Table S2
also prints two identical 156.7(4.0) GPa, 3670(370) K Rh2O3(II) volumes, while
the main text identifies that state as the onset of CaIrO3-type Al2O3. The
rows are retained and flagged, but the executable record conservatively keeps
the authors' 153 GPa and 3400 K scope.

### Numerical reproduction and refit

For selected corundum M3, Peritheos fixes the conventional six-formula-unit
cell `V0 = 255.1 A^3`, `K0' = 4`, `theta0 = 1100 K`, `Tr = 300 K`, and `n = 5`
exactly as reported, and varies `K0`, `gamma0`, and `q`. The source does not
state its weighting objective. An ordinary nonlinear least-squares fit of
pressure residuals at all 75 printed volume-temperature states independently
recovers the selected coefficients:

| Parameter | Published | Peritheos unweighted refit |
|---|---:|---:|
| `K0` (GPa) | 246 +/- 2 | 246.308 +/- 1.333 |
| `gamma0` | 1.32 +/- 0.07 | 1.35931 +/- 0.06725 |
| `q` | 0.8 +/- 0.4 | 0.800868 +/- 0.35270 |

The published corundum parameterization has a 1.23305 GPa pressure RMSE on
Table S1; the refit RMSE is 1.19139 GPa. All three free coefficients agree
within the published uncertainties. An errors-in-variables fit using the
printed pressure and volume uncertainties drives `q` close to zero, so the
successful comparison is recorded specifically as unweighted
pressure-residual parity and is not presented as knowledge of an unpublished
weighting scheme. The executable record retains the published coefficients.

At 300 K the thermal increment vanishes. The selected cold curve gives
`P(129.8 A^3, 300 K) = 100.1123 GPa`, reproducing the first Table S2 value
`100.0(1.0) GPa`. At the independent heated state `V = 125.5 A^3` and
`T = 1560 K`, it gives 130.7277 GPa versus `129.6(2.0) GPa`. Across all 75
printed rows, the published parameterization has a pressure RMSE of
1.1434 GPa.

The generic Peritheos errors-in-variables refit uses all 75 rows, pressure and
volume uncertainties, and the published fixed values `K0' = 4`, `q = 1`,
`Tr = 300 K`, and `n = 5`. Temperature uncertainties are not used because the
source omits them for every 300 K row and for the 710 K row.

| Parameter | Published | Peritheos refit |
|---|---:|---:|
| `V0` (A^3/cell) | 165.2 +/- 0.7 | 167.194 +/- 1.813 |
| `K0` (GPa) | 256 +/- 6 | 239.415 +/- 14.703 |
| `theta0` (K) | 600 +/- 200 | 766.258 +/- 572.847 |
| `gamma0` | 1.47 +/- 0.05 | 1.5502 +/- 0.1622 |

The refit pressure RMSE is 0.8660 GPa and reduced chi-square is 0.152. Every
coefficient agrees within the combined two-standard-deviation uncertainty,
but `theta0` differs by 27.7%, so the campaign classifies the result as
`similar` rather than strict parity. Unpublished weighting, covariance,
temperature-error handling, or additional point selection can explain the
remaining coefficient tradeoff; the published parameterization remains the
executable record.

## C02: cubic boron nitride, Datchi et al. (2007)

### Primary evidence and equation convention

The source is Datchi, Dewaele, Le Godec, and Loubeyre, *Equation of state of
cubic boron nitride at high pressures and temperatures*, Physical Review B
**75**, 214104 (2007),
[doi:10.1103/PhysRevB.75.214104](https://doi.org/10.1103/PhysRevB.75.214104).
The equation and data were audited in the author-posted
[arXiv manuscript](https://arxiv.org/abs/cond-mat/0702656): Section II,
Equations (1)--(4), Tables IV--V, and Figure 4.

The source explicitly defines

\[
P(V,T)=P_0(V)+P_{th}(V,T), \qquad P_{th}(V,0)=0.
\]

Table V's Vinet parameters are consequently a 0 K cold curve, not the Table I
295 K isotherm. On Peritheos's conventional-cell public volume basis they are
`V0 = 5.9026 * 8 = 47.2208 A^3`, `K0 = 397 GPa`, and `K0_prime = 3.62`.
The thermal term uses `theta0 = 1700 K`, `gamma0 = 1.04`, `q = 4`, and `n = 2`
atoms per BN formula unit. The source definitions
`gamma=-d ln(theta)/d ln(V)` and `gamma=gamma0(V/V0)^q` require the integrated
characteristic-temperature law

\[
\theta(V)=\theta_0\exp[(\gamma_0-\gamma(V))/q].
\]

The production record therefore uses `thermal_pressure_reference =
"absolute_zero"`. Its positive `Tr = 295 K` is only the record's default and
the baseline for explicitly requested thermal-pressure increments; it is not
subtracted from total pressure. Composing the old 295 K Vinet record with a
reference-subtracted MGD wrapper would be a different parameterization.

### Data, refit, and numerical checks

All 66 Table IV rows are bundled as
`cubic_boron_nitride_datchi_2007_table4_pvt`: 38 rows at 295 K, 21
simultaneous high-pressure/high-temperature rows, and seven ambient-pressure
thermal-expansion rows. The reproduction uses `V=a^3/8`, rather than the
table's more coarsely rounded printed atomic-volume column. The paper reports
an average lattice-parameter uncertainty of `5e-4 A` and temperature
uncertainty of `+/-5 K`, but no row-wise errors or regression weights.

| Calculation | Result |
|---|---:|
| Published Table V curve, all 66 rows | pressure RMSE `0.5838 GPa` |
| Paper's reported Figure 4/Section V statistic | pressure rms `0.6 GPa` |
| Published curve, `P=0`, 300 K | `5.90550 A^3/atom` |
| Independent Table VI value, `P=0`, 300 K | `5.9055 A^3/atom` |
| Peritheos unweighted pressure-residual refit, q only | `q=7.35 +/- 1.89`; RMSE `0.5654 GPa` |

The refit holds every quantity that the paper says was fixed: the cold Vinet
coefficients, `theta0`, `gamma0`, and `n`. Its q discrepancy is less than two
combined standard deviations from the published `4 +/- 1.5`. Because the
primary table is rounded and the objective, weights, covariance, and unrounded
inputs are not published, the diagnostic refit is not promoted as an EOS
record and does not replace Table V. Run it with:

```bash
python scripts/reproduce_datchi_2007_cbn.py
```

### Scope and pressure calibration

The phase is zinc-blende c-BN, `F-43m` (#216), with four BN formula units per
conventional cubic cell; the article defines atomic volume as `a^3/8` and
reports no transition during 295 K compression. The experiment reaches 162.5
GPa at 295 K, 84.2 GPa at 600 K, 54.2 GPa at 900 K, and 948 K only at ambient
pressure, so the marginal bounds are not a rectangular validated domain.

Reported pressures are on Holzapfel's 2005 ruby scale. Heated points use a
SrB4O7:Sm2+ calibration adjusted to that scale; 295 K points above 100 GPa use
the Dewaele et al. 4He EOS with an H2005 correction. The printed source omits
the row-wise marker shifts and helium volumes, preventing exact recalculation.
Parenthetical coefficient errors are fit standard deviations, not absolute
uncertainties; the authors identify pressure calibration as dominant and
estimate a 3% spread among plausible ruby scales at 160 GPa. No coefficient
covariance or explicit data-reuse license is stated. Section IV prints
`gamma_th0=1.04(1)`, while the final parameter set in Table V prints
`1.04(2)`; the production record conservatively follows final Table V.

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

<a id="campbell-heinz-1994-cscl-and-rbcl"></a>

## CsCl and RbCl: Campbell and Heinz (1994)

### Source-table correction

[Campbell and Heinz (1994)](https://doi.org/10.1029/94JB00127) report separate
room-temperature BM3 equations of state for CsCl and the high-pressure B2 phase
of RbCl. Table 1 on page 11767 is divided into two consecutive blocks: 24 RbCl
observations from 1.11 to 32.3 GPa, followed by 13 CsCl observations from 6.97
to 28.7 GPa.

The original Peritheos resource labeled as CsCl instead contained the first 21
rows of the RbCl block. It also omitted the final RbCl observations at 28.0,
30.1, and 32.3 GPa. Evaluating those RbCl lattice parameters with CsCl's much
larger reference volume produced the spurious validation result
`K0 = 2.4298 GPa`, `K0' = 17.1164`, and a 12.260 GPa published-curve RMSE.
This was a data-to-material assignment error, not an EOS implementation or
optimizer failure.

The corrected library now contains:

- `cscl_campbell_1994_table1_compression`, the 13-row CsCl block; and
- `rbcl_campbell_1994_table1_compression`, the complete 24-row RbCl block.

The paper used gold as its internal pressure standard, following the Heinz and
Jeanloz (1984) gold EOS, and corrected gold strain anisotropy with the procedure
of Meng et al. (1993). The published table already contains the reduced
pressures. Raw gold diffraction values and row-wise correction terms are not
provided, so the pressure reduction cannot be independently repeated.

### RbCl complete reproduction

The B2 phase of RbCl is not quenchable to ambient pressure. Campbell and Heinz
therefore fixed a hypothetical zero-pressure density
`rho02 = 3.3068(10) Mg/m^3`, chosen to satisfy the measured B1-B2 volume change
and the other experimental constraints. Using `M(RbCl) = 120.9208 g/mol`,
Peritheos converts this to the one-formula-unit reference
`V0 = 60.72146(184) A^3`; `V0` is fixed while `K0` and `K0'` are fitted.

The validation treats the parenthesized pressure and lattice-parameter errors
as one-standard-deviation uncertainties and propagates
`sigma_V = 3 a^2 sigma_a`. It then performs an errors-in-variables BM3 fit to
all 24 rows.

| Parameter | Publication | Peritheos refit | Difference |
|---|---:|---:|---:|
| `K0` (GPa) | 17.9 +/- 1.0 | 17.8808 +/- 1.0582 | 0.11% |
| `K0'` | 5.23 +/- 0.29 | 5.23815 +/- 0.35259 | 0.16% |

The published and refitted curves have pressure RMSE values of 0.290 and
0.279 GPa, respectively. Both coefficients satisfy the numerical and combined
two-standard-deviation criteria, so the new
`rbcl_b2_campbell_1994_bm3_1` record achieves full parity.

### CsCl complete reproduction with fit-protocol qualification

For CsCl, the paper fixes `a0 = 4.123 A`, hence
`V0 = 70.087408867 A^3`. The published `K0 = 17.01(29) GPa` and
`K0' = 5.49(15)` come from a joint least-squares fit to the 13 new Table 1
observations and earlier Yagi (1978) data. Yagi's Table 1 supplies nine 25 degC
`V/V0` values at 10--90 kbar. Campbell and Heinz state that those values used
the older `a0 = 4.118 A` reference and corrected them to 4.123 A. Peritheos
therefore applies

`(V/V0)_corrected = (V/V0)_Yagi * (4.118/4.123)^3`,

where the correction factor is `0.9963662826302614`. Yagi states that the
room-temperature volume uncertainty is below 0.15%; this is retained as an
accuracy bound rather than silently interpreted as a one-standard-deviation
weight.

Because Campbell and Heinz do not publish the numerical weights used for their
straight-line normalized-stress regression, the reproducible ledger uses an
unweighted pressure-residual BM3 fit to all 22 observations:

| Parameter | Publication, Campbell + Yagi | Peritheos, all 22 rows | Difference |
|---|---:|---:|---:|
| `K0` (GPa) | 17.01 +/- 0.29 | 17.6967 +/- 0.5706 | 4.04% |
| `K0'` | 5.49 +/- 0.15 | 5.20260 +/- 0.18002 | 5.23% |

The published and refitted curves have pressure RMSE values of 0.390 and
0.344 GPa, respectively, on the complete combined dataset. Both coefficients
meet the numerical limits and overlap within combined two-standard-deviation
uncertainty, so `cscl_campbell_1994_bm3_1` remains classified as `parity`.
The residual central-value difference is documented as a fit-protocol
qualification, not missing-data parity.

<a id="phase-egg-mookherjee-2019"></a>

## Phase Egg: Mookherjee et al. (2019)

### Publication identity and scientific scope

The publication of record is Mookherjee, Panero, Wunder, and Jahn,
*Anomalous elastic behavior of phase egg, AlSiO3(OH), at high pressures*,
*American Mineralogist* **104**, 130--139 (2019),
[doi:10.2138/am-2019-6694](https://doi.org/10.2138/am-2019-6694). The separate
[10.2138/am-2018-6694](https://doi.org/10.2138/am-2018-6694) registration
resolves to the accepted manuscript. That manuscript labels itself a preprint
and instructs readers to cite the final 2019 DOI. Peritheos therefore creates
one literature EOS record under the 2019 DOI and retains the 2018 identifier
only in `source_lineage` as a discovery alias.

The implemented record is the ordered low-pressure proton configuration,
phase Egg (LP), with composition `AlSiO3(OH)` (`AlSiO4H`). Its diffraction
structure is independently sourced to Schmidt et al. (1998),
[doi:10.2138/am-1998-7-820](https://doi.org/10.2138/am-1998-7-820), Table 4:
monoclinic `P21/n` (the alternate setting of `P21/c`, number 14),
`a=7.14409 A`, `b=4.33462 A`, `c=6.95253 A`, `beta=98.396 deg`,
`V=212.99 A^3`, and `Z=4`. Seven fully occupied general `4e` sites reproduce
four `AlSiO4H` formula units. This measured structural cell is retained
separately from the static-DFT EOS reference cell.

Mookherjee et al. calculated the primitive cell with GGA-PAW/VASP, an 800 eV
cutoff, a 6 by 9 by 6 k-point mesh, and the stated van der Waals correction.
For primitive monoclinic `P21/c`, that calculation cell is also the
conventional `Z=4` cell used by Peritheos. The source fits total energy to

\[
E=E_0+\frac{9}{2}K_0V_0\left[f_V^2+(K'_0-4)f_V^3\right],\qquad
f_V=\frac{1}{2}\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right],
\]

whose negative volume derivative is the standard third-order
Birch--Murnaghan pressure equation. The LP fit reports
`V0=210.21+/-0.14 A^3`, `K0=164.4+/-1.8 GPa`, and
`K0'=7.14+/-0.24`; all three were fitted. The source does not define the
confidence convention or publish covariance. It describes these calculations
as static conditions (0 K), so the record uses a 0 K reference isotherm rather
than relabeling it as a room-temperature EOS.

### Supplement and numerical reproduction

The [official MSA data deposit](https://www.minsocam.org/MSA/AmMin/TOC/2019/Jan2019_data/AM-19-16694.zip)
contains `6694_supp.xlsx`. Supplementary Table 1 reports 11 LP volumes from
220 to 180 `A^3` and corresponding static pressures from -6.4 to 43.7 GPa;
no row uncertainties are given. The separate five-row HP block and elastic
tensor are outside this one-record EOS dataset.

At the paper's proton-transfer volume `V=196 A^3`, the published LP
coefficients calculate `P=14.725955 GPa`, reproducing the workbook's
`14.7 GPa` to its 0.1 GPa printed precision. Across all 11 LP rows, the
published curve has a pressure RMSE of `0.03011 GPa` and a maximum absolute
residual of `0.04836 GPa`.

An independent unweighted BM3 pressure-residual fit gives:

| Parameter | Publication | Diagnostic P-V refit | Difference |
|---|---:|---:|---:|
| `V0` (A^3) | 210.21 +/- 0.14 | 210.21760 | +0.00760 |
| `K0` (GPa) | 164.4 +/- 1.8 | 164.83396 | +0.43396 |
| `K0'` | 7.14 +/- 0.24 | 7.08137 | -0.05863 |

The refit RMSE is `0.02339 GPa`; every coefficient agrees within its published
error. This is a strong curve-level check, but not an exact reconstruction of
the authors' regression: the source fitted total energies, while the deposited
workbook omits those energies, row uncertainties, objective, weights,
covariance, fitting-software details, and fit statistic. The refit ledger
therefore records that limitation explicitly.

The article defines LP below the approximately 15 GPa proton transfer. Its
supplement continues the LP starting configuration metastably to 43.7 GPa.
Peritheos preserves every row but limits the preferred LP validity branch to
-6.4--15 GPa at 0 K. The separately published HP parameterization is neither
merged with LP nor added as a second record in this publication-scoped change.

<a id="coo-clendenen-1966"></a>

## CoO: Clendenen and Drickamer (1966)

### What the publication actually provides

[Clendenen and Drickamer (1966)](https://doi.org/10.1063/1.1726610)
measured the compression of nine oxides and sulfides. For rocksalt CoO,
Table II gives `a0 = 4.258 A`, Equation (4) is the Murnaghan relation, and
Table VI gives the whole-range empirical coefficients `B0 = 1905 kbar`
(`190.5 GPa`) and `B0' = 3.9`. The authors caution that coefficients obtained
over the complete pressure interval need not equal the true one-atmosphere
derivatives.

Table III is not a row-level diffraction dataset. It is explicitly described
as smoothed data and reports only nine pairs:
`a/a0 = 1.000, 0.995, ..., 0.960` with corresponding pressures from 0 to
308 kbar. The normalized lattice values have three decimal places and the
pressures are rounded to whole kilobars. No row uncertainties, parameter
errors, covariance, or regression weights are printed. Figure 2 contains the
individual observations, but only graphically; the paper does not provide the
numerical rows or say exactly how they were weighted in the Table VI fit.

### Objective-sensitivity audit

Peritheos fixes `V0 = 4.258^3 = 77.199941512 A^3` and evaluates the same
Murnaghan equation under three ordinary least-squares conventions. This tests
whether the failed derivative comparison is merely a residual-variable choice.

| Fit to the nine printed Table III points | `K0` (GPa) | `K0'` |
|---|---:|---:|
| Publication (Table VI; exact rows/weights unreported) | 190.5 | 3.9 |
| Pressure residuals, `P(a)` | 176.797 | 5.14810 |
| Lattice-ratio residuals, `a(P)/a0` | 180.765 | 4.74012 |
| Volume-ratio residuals, `V(P)/V0` | 181.290 | 4.68300 |

Changing the residual variable narrows the disagreement but does not recover
the Table VI derivative when both coefficients are free. The published curve
nevertheless follows the printed table closely: its pressure RMSE is
`0.519 GPa`, compared with `0.423 GPa` for the unconstrained pressure-residual
minimum.

### Why rounding alone is not an explanation

The displayed precision permits at most `+/-0.0005` in `a/a0` and
`+/-0.5 kbar` in pressure. Even allowing both coordinates to move anywhere
inside those rounding intervals, four rows (`a/a0 = 0.975`, `0.970`, `0.965`,
and `0.960`) cannot lie on the published Murnaghan curve. Their minimum
remaining pressure gaps are `1.98`, `4.44`, `0.63`, and `3.64 kbar`,
respectively. The discrepancy therefore cannot be described as ordinary
decimal rounding of points generated by the published curve.

The stronger diagnostic is the coefficient covariance. For the nine printed
points, the two-parameter pressure-residual refit has a `K0`--`K0'`
correlation of `-0.971`: a higher modulus is compensated by a lower pressure
derivative along a shallow objective valley. Two conditional fits expose this
directly:

| Constraint | Best remaining coefficient | Pressure RMSE (GPa) |
|---|---:|---:|
| `K0 = 190.5 GPa` fixed | `K0' = 3.8391` | 0.516 |
| `K0' = 3.9` fixed | `K0 = 189.090 GPa` | 0.504 |

Thus each published coefficient is recovered when its strongly correlated
partner is held fixed. The unconstrained optimizer moves both together to
`176.797 GPa` and `5.1481` for only a `0.0965 GPa` reduction in pressure RMSE.
Using the residual variance estimated by that fit, the published parameter
pair also remains inside the approximate joint 95% least-squares confidence
region of the printed table.

This behavior agrees with the paper's own methodological warning immediately
after Table VI: the calculated pressures are relatively insensitive to the
individual values of `B0` and `B0'` when the pair is balanced. Consequently,
the large relative change in `K0'` is visually alarming but is not evidence of
a comparably large difference between the two pressure-volume curves.

### Parity conclusion

The stored `K0` and `K0'` remain the publication's authoritative Table VI
coefficients. The result remains `parity_not_achieved` under the strict
central-value criterion, but the earlier shorthand explanation of a rounding
error is rejected. The defensible cause is a strongly non-identifiable
two-parameter fit combined with unavailable numerical observations, weights,
and source regression details; Table III is only a smoothed surrogate. Exact
coefficient parity would require the numerical Figure 2 observations and the
authors' weighting/objective. Changing the EOS family or optimizer is not
supported by the evidence.

<a id="goethite-gleason-2008"></a>

## Goethite: Gleason et al. (2008)

### Correction of the represented equation and fit scope

[Gleason et al. (2008)](https://doi.org/10.2138/am.2008.2942) do not report
a standalone room-temperature BM3 fit. Their Equation 1 combines a BM3
reference curve with Mie-Gruneisen thermal pressure,

\[
P(V,T)=P_{\mathrm{BM3}}(V)+\frac{\gamma(V)}{V}
\left[E_D(V,T)-E_D(V,300\,\mathrm{K})\right].
\]

Table 1 gives `V0 = 138.75(2) A^3`, `K0 = 140.3(37) GPa`,
`gamma0 = 0.91(7)`, and `theta0 = 740(5) K`; the text gives
`K0' = 4.6(4)` and assumes `d ln(gamma)/d ln(V) = q = 1`. Peritheos uses
`n = 4` atoms per FeOOH formula unit and converts the four-formula-unit
conventional cell volume to molar volume for the Debye energy calculation.
The measured `V0` is fixed because the paper says specifically that all P-V-T
data determine `K0` and `K0'`. The thermally derived `gamma0` and `theta0`, the
assumed `q`, the reference temperature, and the stoichiometric atom count are
also held fixed in the reproduction.

The first validation had represented this record as plain BM3 and therefore
discarded all but the 27 room-temperature observations. The material record
now executes the combined thermal equation and the refit uses all 65 rows in
MSA depository Table 1, matching the source's stated fit scope. The deposited
temperatures span 23--310 degC; this is retained even though the abstract
summarizes the experimental range as 23--250 degC.

### An anomalous row in the official deposit

Depository row 32 reports `P = 6.66 GPa`, `T = 100 degC`,
`a = 4.159(2) A`, and `V = 122.80(8) A^3`. The two crystallographic values
are internally consistent with the reported `b` and `c`, but the state is
physically isolated: the neighboring 5.55 and 5.65 GPa rows have volumes
133.84 and 134.0 A^3. The published EOS predicts 133.105 A^3 at the row-32
pressure and temperature, or 22.953 GPa at its reported volume. Its pressure
residual is therefore 16.293 GPa.

A transposition from `a = 4.159` to `4.519 A` would imply approximately
133.44 A^3 and would agree with the neighboring rows, but there is no primary
source authority for making that correction. Peritheos consequently retains
the official values verbatim, flags the row, and reports sensitivity both with
and without it.

### Primary tables compared with Figure 3 digitization

The deposited numerical tables were checked against a fresh raster-axis
digitization of Figures 3a and 3b. This figure check is diagnostic only; the
primary tables, rather than digitized markers, remain the fit inputs.

Figure 3a does not show all 65 alpha-FeOOH rows. Its legend and markers contain
the 27 observations at 23 degC, 10 at 200 degC, and 10 at 250 degC: 47 table
positions in total. The 11 observations at 100 degC, six at 150 degC, and the
single 310 degC observation are omitted. In particular, anomalous depository
row 32 is a 100 degC point and therefore has no marker in Figure 3a.

After calibrating the pressure and volume axes from their labeled ticks, all
47 plotted alpha-FeOOH table coordinates fall on the corresponding marker or
an overlapping marker cluster. A local marker-center extraction differs from
the tabulated values by no more than approximately 0.07 GPa and 0.09 A^3,
both well inside the plotted symbol width. The figure therefore supports the
table transcription but supplies no alternate value for row 32. It also
visibly confirms that several high-pressure room-temperature points lie above
the printed EOS curve.

Figure 3b similarly agrees with all 49 epsilon-FeOOH Table 2 coordinates.
The exact 200, 250, 300, 350, and 400 degC rows use the five legend symbols;
the isolated 310, 324, 385, 390, and 395 degC measurements are grouped into
the neighboring nominal-temperature symbol sets. Because several markers
overlap, not every row is individually distinguishable, but every table
coordinate lands on a visible marker or marker cluster. The local calibrated
differences are at most approximately 0.05 GPa and 0.03 A^3.

This comparison also confirms the phase separation: Figure 3a and depository
Table 1 contain alpha-FeOOH with volumes near 119--139 A^3, whereas Figure 3b
and Table 2 contain epsilon-FeOOH with volumes near 60--65 A^3. No
epsilon-FeOOH coordinate entered the alpha-FeOOH refit.

### Figure-marker-only refit and equation audit

Fitting the Figure 3a marker centers themselves does not recover the published
coefficients. The pressure and volume axes were calibrated from the labeled
ticks, marker centers were extracted locally, and the marker shape was used to
assign the plotted temperature. Because overlapping symbols cannot always be
separated blindly, the primary table positions were used only to identify the
corresponding plotted marker cluster. The fits below use unweighted pressure
residuals; `V0 = 138.75 A^3` is fixed except in the explicitly free row.

| Figure 3a fit | Points | `V0` (A^3) | `K0` (GPa) | `K0'` | Pressure RMSE (GPa) |
|---|---:|---:|---:|---:|---:|
| Publication | -- | 138.75 | 140.3 | 4.6 | -- |
| BM3 + Mie-Gruneisen-Debye, all visible markers | 47 | 138.75 fixed | 193.45 | 0.893 | 0.702 |
| BM3 + Mie-Gruneisen-Debye, all visible markers, `V0` free | 47 | 138.668 | 196.96 | 0.728 | 0.700 |
| Plain BM3, all visible markers | 47 | 138.75 fixed | 203.78 | 0.262 | 0.840 |
| Plain BM3, 23 degC markers only | 27 | 138.75 fixed | 189.35 | 1.158 | 0.863 |

The same fits to the primary table values for the visible rows give,
respectively, `(K0, K0') = (193.77 GPa, 0.849)`, `(197.62 GPa, 0.671)` with
`V0 = 138.660 A^3`, `(204.03 GPa, 0.229)`, and `(189.77 GPa, 1.103)`.
The small digitization-to-table shifts are negligible compared with the
approximately 50--64 GPa discrepancy from the publication. Letting `V0` vary
also does not resolve it.

Conditional fits rule out a simple coefficient-correlation explanation. With
`K0 = 140.3 GPa` fixed, the digitized 47-marker thermal fit requires
`K0' = 6.112`; plain BM3 requires `6.209`, and the 23 degC markers alone
require `5.862`. The corresponding primary-table results are `6.057`, `6.153`,
and `5.806`. Conversely, fixing `K0' = 4.6` in the visible primary rows gives
`K0 = 157.02 GPa` with the thermal model (`158.93 GPa` as plain BM3). Thus
neither reported coefficient conditionally recovers the other.

The low derivative is already present without any high-temperature reduction.
An unweighted BM3 fit to the 27 primary 23 degC rows gives
`K0 = 189.775 GPa` and `K0' = 1.103`; the corresponding digitized markers give
`189.349 GPa` and `1.158`. Restricting the primary data to the 23 room-
temperature compression points, rather than mixing compression and
decompression, gives `194.207 GPa` and `1.041`. For the all-room-temperature
fit, the calculated static bulk modulus rises only from `189.8 GPa` at `V0` to
`199.5 GPa` at `V = 130 A^3`, then decreases to `196.6 GPa` at
`V = 120 A^3`. The fitted `K0'` must therefore not be interpreted as credible
evidence that goethite scarcely stiffens: it is the compensation chosen when a
single smooth BM3 curve is forced through a discontinuous, non-hydrostatic
compression path.

The source text and figure captions do not support an external compression-
dataset fit. Nagai et al. (2003) are shown as a comparison in Figure 2 and
their separately fitted `K0 = 111(2) GPa`, with `K0' = 4` fixed, is contrasted
with the new result in the discussion. Figure 3 contains only Gleason et al.'s
compression and decompression symbols, and its 47 plotted positions map to the
new depository table. The only explicitly imported numerical thermodynamic
input is the heat capacity from Majzlan et al. (2003a), used in deriving
`gamma0` and `theta0`; it is not another P-V dataset.

As a direct sensitivity test, the 12 pressure-volume rows in
[Nagai et al. (2003)](https://doi.org/10.2138/am-2003-1005) were normalized
from their `V0 = 140.45 A^3` to Gleason et al.'s `138.75 A^3` and appended to
the Gleason rows. A joint fit to the 27 room-temperature Gleason rows plus
Nagai gives `(K0, K0') = (180.99 GPa, 0.060)`; appending Nagai to all 47
visible Gleason rows gives `(188.98 GPa, -0.332)` with the thermal correction.
Thus even an undeclared Gleason--Nagai combination drives the derivative lower
and does not explain `(140.3 GPa, 4.6)`.

The paper's intended static relation is the standard third-order
Birch-Murnaghan equation used by Peritheos: the curve drawn in Figure 3a is
consistent with the standard sign and would leave the plotted volume range at
high pressure if the opposite sign were used. Equation 1 as printed is not a
usable literal specification, however. With the paper's positive-compression
strain definition, its static bracket contains a minus sign where standard BM3
requires a plus sign, and the closing delimiter/addition around the thermal
term is missing. The printed thermal term is consequently dimensionally
ambiguous. The sign is best treated as a typographical error because the
Figure 3 caption explicitly calls the curve third-order Birch-Murnaghan.

Peritheos and the paper otherwise share the same intended physical ingredients:
an additive Mie-Gruneisen thermal pressure referenced to 300 K, Debye internal
energy, `gamma = -d ln(theta)/d ln(V)`, and the stated
`d ln(gamma)/d ln(V) = 1`. Peritheos makes the implied constant-`q` volume laws
for `gamma(V)` and `theta(V)` explicit and thermodynamically consistent. The
paper does not disclose the complete regression implementation, weighting, or
whether the derived thermal quantities were actually applied while optimizing
`K0` and `K0'`; Figure 3 shows only one plain-BM3 curve and says thermal
expansivity at pressure could not be resolved. Thus the intended formalism is
close, but an exact implementation of the source fit is not recoverable from
the publication.

**Use recommendation:** retain `goethite_gleason_2008_bm3_1` for provenance
and for reproducing the published curve only. Do not use it as a recommended
quantitative goethite pressure-volume or thermoelastic EOS. The primary table,
the independent figure-marker fit, and every tested reasonable interpretation
of Equation 1 fail to reproduce the reported elastic coefficients.

### Corrected refit and sensitivity

The primary ledger uses all 65 observations and the normal
errors-in-variables policy. Parenthesized cell and temperature errors are the
deposit's 95% Rietveld confidence limits; the paper does not disclose whether
or how they entered its EOS regression.

| Fit | `V0` (A^3, fixed) | `K0` (GPa) | `K0'` |
|---|---:|---:|---:|
| Publication | 138.75 | 140.3 +/- 3.7 | 4.6 +/- 0.4 |
| Peritheos errors-in-variables, all 65 rows | 138.75 | 183.338 +/- 4.974 | 0.000 +/- 0.413 |
| Peritheos errors-in-variables, excluding row 32 | 138.75 | 177.642 | 1.54188 |
| Unweighted pressure residuals, all 65 rows | 138.75 | 185.505 | 0.000 |
| Unweighted pressure residuals, excluding row 32 | 138.75 | 177.800 | 1.59886 |

The published curve's pressure RMSE is 2.540 GPa on all rows and 1.551 GPa
after removing row 32. Thus the anomaly accounts for much of the curve-level
error, but its removal does not recover the published coefficients under
either tested objective. The failed parity is not an optimizer artifact or
solely the consequence of the initial room-temperature-only validation.

The pressure reduction also cannot be audited row by row. The experimental
section identifies ruby (Mao et al. 1978), NaCl (Birch 1986), and gold
(Shim et al. 2002) as room-temperature standards, with gold alone used at high
temperature. The deposit contains only reduced sample pressures: it supplies
neither raw calibrant observations nor the calibrant assignment for each
room-temperature row. A pressure-scale or mixed-calibrant contribution to the
remaining mismatch therefore cannot be tested from the published table.

The residual structure provides a second clue. The paper itself reports a
room-temperature discontinuity near 9 GPa and attributes it to solidification
of the 4:1 methanol-ethanol medium and resulting non-hydrostatic stress. If the
reproduction is restricted diagnostically to the first 11 room-temperature
rows, ending below that discontinuity, it gives `K0 = 144.775 GPa` and
`K0' = 6.4865` with a 0.212 GPa RMSE. That subset is statistically compatible
with the printed pair but is not adopted because it contradicts the explicit
statement that all P-V-T data were used.

### Parity conclusion

The record remains `parity_not_achieved`. Two concrete problems are now
separated:

- the official deposited data contain one extreme but internally propagated
  lattice/volume anomaly; and
- much of the remaining disagreement begins where the authors themselves say
  the experiment becomes non-hydrostatic.

Neither deleting row 32 nor changing from errors-in-variables to unweighted
pressure residuals produces the publication's coefficients. Exact parity
would require the authors' row weights, regression residual definition,
treatment of non-hydrostatic observations, and any preprocessing not described
in the article. The library retains the published parameterization as an
explicitly warned archival record while exposing the complete source-faithful
refit and all sensitivity results; it is not a recommended quantitative EOS.

<a id="mo2c-haines-2001"></a>

## Mo2C: Haines et al. (2001)

### Source, data, and fit scope

- Primary article: Haines et al., *Experimental and theoretical investigation
  of Mo2C at high pressure*, Journal of Physics: Condensed Matter **13**,
  2447--2454 (2001),
  [doi:10.1088/0953-8984/13/11/303](https://doi.org/10.1088/0953-8984/13/11/303).
- Primary observations: Figure 2 contains 16 plotted `P`--`V/V0` markers. The
  first nine, through 14.87 GPa, are measurements in 21:4:1
  methanol--ethanol--water; the remaining seven, from 19.71 to 45.50 GPa, are
  measurements with NaCl after laser heating. The article does not print the
  numerical pressure--volume array or the regression weights.
- Digitization: vector marker centers were extracted directly from the
  publisher PDF and affinely calibrated from the printed axes. The marker
  shapes and therefore the two experimental regimes are retained in dataset
  `molybdenum_carbide_haines_2001_figure2_digitized`.
- Model: standard third-order Birch--Murnaghan. The article plots relative
  volume and publishes `K0 = 307(5) GPa` and `K0' = 6.2(3)`, but does not report
  `V0` as a fitted coefficient. The reproduction therefore fixes `V0` at
  `148.9071 A^3`, calculated from the same specimen's measured ambient
  hexagonal subcell and expressed as the equivalent four-formula-unit
  orthorhombic conventional cell.

The first validation incorrectly varied `V0`. That added a third correlated
degree of freedom that is absent from the source's relative-volume fit and
produced `V0 = 148.713 A^3`, `K0 = 343.734 GPa`, and `K0' = 4.147`.

### Corrected refit

| Calculation | Points | `V0` (A^3) | `K0` (GPa) | `K0'` | Pressure RMSE (GPa) |
|---|---:|---:|---:|---:|---:|
| Publication | 16 plotted | measured reference | 307 +/- 5 | 6.2 +/- 0.3 | 0.704 on digitized markers |
| Peritheos errors in variables | 16 | 148.9071 fixed | 325.874 +/- 9.500 | 4.909 +/- 0.651 | 0.375 |
| Unweighted pressure residuals | 16 | 148.9071 fixed | 326.197 | 4.859 | 0.618 |
| Unweighted normalized-stress regression | 15 nonzero-pressure points | 148.9071 fixed | 320.955 | 5.158 | -- |

The corrected errors-in-variables result is statistically compatible with
both published coefficients at combined two-standard-deviation uncertainty.
It is classified `similar`, rather than strict `parity`, because the `K0'`
point estimate differs by 1.291, or 20.82%, narrowly outside the library's
numerical limit of 1.0 absolute or 20% relative. The curve-level discrepancy is
small: the refit reduces the pressure RMSE on the marker centers from 0.704 to
0.375 GPa.

The all-marker result is also exposed as the opt-in EOS record
`molybenum_carbide_mo2c_haines_2001_bm3_refit`, with explicit lineage to the
published record `molybenum_carbide_mo2c_haines_2001_bm3_1`. This makes the
reproducible Peritheos reduction usable without presenting it as a coefficient
set reported by Haines et al. The published record remains the authoritative
literature parameterization.

The normalized-stress calculation tests whether the difference is merely due
to a conventional finite-strain linearization. It shifts the coefficients only
to `320.955 GPa` and `5.158`, so that convention alone does not recover the
published pair.

### Pressure-regime and endpoint sensitivity

The coefficient tradeoff is driven by curvature across the two experimental
regimes rather than by a uniform offset:

| Diagnostic errors-in-variables fit | Points | `K0` (GPa) | `K0'` |
|---|---:|---:|---:|
| Methanol--ethanol--water, at or below 15 GPa | 9 | about 284.2 | about 10.8 |
| NaCl and laser-heated regime, above 15 GPa | 7 | about 352.6 | about 3.31 |
| All markers at or below 40 GPa | 14 | 309.277 | 6.621 |
| All Figure 2 markers | 16 | 325.874 | 4.909 |

The under-40-GPa diagnostic reproduces the published pair closely, while the
two final markers at 43.10 and 45.50 GPa pull the fitted curvature down. This
does not justify deleting them: they are explicitly plotted as observations,
and the source does not identify a row exclusion or alternative weighting for
them. Their influence is reported to localize the discrepancy, not to create
an opportunistic parity result. Consequently, the under-40-GPa calculation is
documented only as a sensitivity test and is not stored as another EOS record.

### Parity conclusion

The original large failure was primarily a validation-scope error and is now
resolved to `similar`: `V0` is fixed as required by the source representation,
and both remaining coefficients overlap the published two-sigma intervals.
The small residual point-estimate mismatch cannot be assigned to digitization
alone, because it persists across pressure-residual and normalized-stress
fits. The most bounded explanation is sensitivity to the two measurement
regimes, the last two high-pressure observations, and unpublished row weights
or preprocessing. Exact coefficient parity requires the authors' numerical
P--V array, row mask, and fitting weights; the published parameterization
remains the default library EOS, while the separately labeled refit record
provides the direct all-marker Peritheos result.

## Luo et al. (2023) MgO

### Source, identity, and equation

The primary source is Luo et al., *Equation of state of MgO up to 345 GPa and
8500 K*, *Physical Review B* **107**, 134116 (2023),
[doi:10.1103/PhysRevB.107.134116](https://doi.org/10.1103/PhysRevB.107.134116).
The sample is stoichiometric B1 (NaCl-structure) MgO/periclase, represented by
the conventional cubic cell with four MgO formula units. The paper gives an
initial density of `3.590(4) g/cm3` and globally optimizes a zero-pressure,
zero-temperature specific volume `V0K = 0.2767(1) cm3/g`, bulk modulus
`B0 = 169.8(5) GPa`, and pressure derivative `B' = 4.501(7)`. Using the
paper's MgO molar mass, the cold specific volume is
`74.0741025123 A3/conventional cell`; the ambient initial density corresponds
to `74.5697677586 A3/conventional cell`.

The catalog record is not a static Vinet fit. It implements Appendix B's
executable pressure scale,

\[
P(V,T)=P_{\mathrm{Vinet},0K}(V)+c_0+c_1\delta\eta+c_2\delta T
+\tfrac12c_3\delta\eta^2+\tfrac12c_4\delta T^2
+\tfrac12c_5\delta\eta\delta T,
\]

where `eta = 1 - V/V0K`, `delta eta = eta - 0.02`, and
`delta T = T - 300 K`. The printed coefficients are `c0=0.5096 GPa`,
`c1=-13.4246 GPa`, `c2=6.3295e-3 GPa/K`, `c3=36.2194 GPa`,
`c4=5.4705e-8 GPa/K2`, and `c5=3.2238e-3 GPa/K`. No uncertainties are
reported for these six coefficients. The three cold-curve uncertainties have
an unstated confidence convention and no published covariance.

The paper's Equations (1)--(6) give the underlying quasi-Debye Helmholtz
construction, but not the numerical longitudinal- and shear-wave velocity
fits needed to reconstruct its effective sound velocity and Debye-temperature
law. The Appendix-B polynomial is therefore the complete independently
executable form supported by the primary article. The APS article record has
no linked official supplement or machine-readable deposit.

### Data reproduction and refit boundary

Peritheos bundles all five new shock states from Table I and all 576 derived
P-V-T values in Tables II--III. Re-evaluating the printed Appendix-B equation
at the table coordinates gives `0.484096 GPa` RMS and `1.435303 GPa` maximum
absolute pressure residual. For example, the state at ambient compression
`0.42` and `8500 K` is reproduced within that documented sub-GPa consistency
at the table's printed `342.01 GPa`. The table was generated from the fuller
quasi-Debye calculation, so exact equality after decimal rounding is neither
claimed nor expected.

A diagnostic unweighted least-squares fit of only `c0`--`c5` to that derived
grid gives `(1.423975, -17.564509, 0.006092604, 45.696538,
8.377914e-8, 0.004084916)` and lowers the RMS only to `0.445789 GPa`. This is a
fit to model output, not to primary observations, and is not stored as another
EOS record. An observation-level refit is not possible from the publication:
the global optimization selected data from multiple earlier studies, while
the complete row set, unpublished sound-velocity fits, objective weights, and
covariance are unavailable. The refit ledger records that precise limitation
rather than circularly fitting the authors' derived table.

Table I pressure is independently checkable from the Rankine--Hugoniot
relation `P=rho0*D*up`; all five printed pressures are recovered within
`0.6 GPa`. Shots 2, 3, and 5 use impedance matching to the stated Ta or Pt
flyer Hugoniots, while shots 1 and 4 use measured particle velocity. The
published `0--345 GPa` and `300--8500 K` limits are retained as the stated
pressure-scale domain, not as a rectangular observation envelope or a B1
phase-stability claim.
<a id="mgo-dewaele-2000"></a>

## MgO: Dewaele et al. (2000)

### Authority, identity, and represented equation

The discovery lead attributed DOI
[10.1029/1999JB900364](https://doi.org/10.1029/1999JB900364) to Fei et al., but
the version of record is Dewaele, Fiquet, Andrault, and Hausermann (2000),
*P-V-T equation of state of periclase from synchrotron radiation
measurements*. Fei (1999) is an input to the paper's combined room-temperature
comparison, not the author of this DOI.

The single implemented record is the paper's preferred complete P-V-T model:

\[
P(V,T)=P_{BM3}(V,300\,\mathrm{K})+
\frac{10^{-4}\gamma(V)}{V}\left[E_D(V,T)-E_D(V,300\,\mathrm{K})\right],
\]

with `V` in cubic angstroms per mole-equivalent formula-unit basis inside the
engine, `gamma(V) = gamma0 (V/V0)^q`, and
`theta(V) = theta0 exp[(gamma0-gamma(V))/q]`. This is Peritheos's explicit
`integrated_gruneisen` convention. The conventional crystallographic MgO
volumes in the material document contain four formula units and are converted
by the material wrapper.

The stored source parameters are `V0 = 74.71(3) A^3`, `K0 = 161(1) GPa`,
`K0' = 3.94 +/- 0.2`, `theta0 = 800(50) K`, `gamma0 = 1.45(10)`, and
`q = 0.8(5)`, with `n = 2`. The provenance is staged: `V0` and `K0` are
adopted ambient/acoustic constraints; `K0'` is obtained from the combined 300 K
data and then fixed in the thermal analysis; `theta0` is obtained from ambient
heat-capacity data; `gamma0` is constrained by the zero-pressure thermal-pressure
intercept and then fixed; and `q` is constrained by the pressure dependence.
Table 3's Murnaghan, Vinet, logarithmic, and alternate constrained BM3 rows are
extrapolation-sensitivity comparisons and are deliberately not represented as
additional full P-V-T records.

The source is internally inconsistent about the uncertainty on `K0'`: the
fit discussion gives a formal `+/-0.06` (rounded to `(5)` in Table 3) and an
expanded `+/-0.2` after propagating pressure errors; the abstract repeats
`+/-0.2`, whereas the conclusion prints `+/-0.3` without explanation. The
record stores `+/-0.2` because it is tied explicitly to the paper's pressure
uncertainty analysis and is repeated in the abstract; all three printed values
remain documented here and in the record provenance.

### Primary observations and pressure basis

Dataset `mgo_dewaele_2000_table2_pvt` contains every Table 2 observation: 41
laser-heated rows and 20 rows at 300 K, spanning 0--53 GPa and 300--2474 K.
Each row retains MgO conventional-cell volume and the simultaneously measured
Pt lattice parameter. Parenthesized pressure and volume errors are transcribed
as printed; the paper assigns +/-200 K to heated temperatures. No official
supplement or machine-readable deposit accompanies the article.

Pressure was reduced from Pt using the Jamieson, Fritz, and Manghnani (1982)
shock-Hugoniot/Debye P-V-T scale. The paper gives the row-wise Pt lattice
parameters and the pressure-error law

\[
\Delta P = 0.03P + 0.0004(T-300),
\]

for pressure in GPa and temperature in kelvin. The exact Jamieson Pt
implementation is not presently bundled as an executable reference EOS, so
the source-scale pressure values are preserved and the calibration is marked
`reference_eos_not_bundled`; no substitute scale is inferred.

### Independent value and refit

The room-temperature discussion reports that the BM3 curve reaches 145 GPa at
`V/V0 = 0.667`. Direct evaluation of the stored record gives 144.947 GPa, a
0.053 GPa difference attributable to the source's integer rounding.

For an observation-level check, the automated campaign fits `q` to the 41
heated Table 2 rows while holding `V0`, `K0`, `K0'`, `theta0`, `gamma0`, `Tr`,
and `n` at their source-staged values. It uses the complete positive row-wise
P, V, and T uncertainties in an errors-in-variables objective.

| Quantity | Published | Current-study refit |
|---|---:|---:|
| `gamma0` | 1.45 +/- 0.10 | 1.45 fixed, following the source's staged procedure |
| `q` | 0.80 +/- 0.50 | 0.847 +/- 0.122 |
| Pressure RMSE | 0.952 GPa for the published curve | 0.950 GPa |
| Reduced chi-square | -- | 0.925 |

The refitted `q` agrees within combined two-standard-deviation uncertainty and
the numerical similarity limit, so the ledger classifies the result as
`parity`. It is a
conditional current-study reproduction: the source's thermal analysis also
uses Fei (1999) observations that are not reprinted in Dewaele et al. Exact
coefficient identity from the 41 new rows alone is therefore neither expected
nor claimed. The published parameters, not the refit, remain the executable
literature record.
## CaSiO3 perovskite: Kawai and Tsuchiya (2014)

The primary source is Kawai and Tsuchiya, *P-V-T equation of state of cubic
CaSiO3 perovskite from first-principles computation*, JGR Solid Earth **119**,
2801--2809, [doi:10.1002/2013JB010905](https://doi.org/10.1002/2013JB010905).
It reports LDA first-principles molecular-dynamics calculations for ideal
cubic `Pm-3m` CaSiO3 and fits both Vinet and third-order Birch--Murnaghan
reference isotherms. Peritheos stores only the preferred Vinet branch, which
is the curve plotted in Figure 2 and the parameterization summarized in the
abstract.

The complete stored thermal EOS is a Vinet reference isotherm at 1000 K with
`V0 = 46.17 A^3/formula unit`, `K0 = 203.95 GPa`, and `K0' = 4.76`, plus a
constant-`q` Mie--Gruneisen--Debye increment with fixed `theta0 = 1100 K`,
`gamma0 = 1.576`, `q = 0.96`, and `n = 5`. The article says it follows the
Tange et al. (2009) procedure; their equations 4--10 define
`gamma(V) = gamma0 (V/V0)^q`, the thermodynamically integrated Debye
temperature, and thermal pressure relative to the reference temperature.
Because cubic `Pm-3m` CaSiO3 has one formula unit per conventional cell, the
reported formula-unit volume is also the public conventional-cell volume.

The underlying calculations used LDA in PWSCF, Vanderbilt pseudopotentials
for Ca and O and a Troullier--Martins pseudopotential for Si, a 50 Ry cutoff,
and canonical-ensemble FPMD with 1 fs steps. The 80-atom `2x2x1`
tetragonal-conventional-cell supercell was constrained to cubic metrics and
sampled on a `2x2x2` Monkhorst--Pack grid. Each state was equilibrated for
1 ps and normally averaged for the next 5 ps. The methods sentence prints a
300--1500 K simulation range, but Figure 2 explicitly shows calculated points
at 1000, 2000, 3000, and 4000 K. With no raw state table available, this
internal temperature-range inconsistency cannot be resolved; the executable
range follows the unambiguous plotted and tabulated 1000--4000 K EOS output.
Figure 1 supplies one convergence benchmark near 15 GPa and 2000 K: the
`2x2x2` sampling gives diagonal stresses of 15.3, 15.5, and 15.2 GPa, whereas
Gamma-only sampling gives 16.0, 16.1, and 14.6 GPa. Because the corresponding
volume is not printed, these values test k-point convergence but cannot serve
as a row in an EOS refit.

There is a source discrepancy in `K0`. The abstract prints 203.95 GPa, while
the PDF text layer and later comparison tables render Table 2 as 203.5 GPa.
The article's 60 one-decimal Table 1 isochors resolve the executable value:
203.95 GPa gives an RMSE of 0.03734 GPa and maximum absolute difference of
0.09071 GPa. At `V/V0 = 0.70` and 4000 K, it gives 184.4177 GPa, which rounds
to the published 184.5 GPa. Using 203.5 GPa increases the grid RMSE to
0.20047 GPa and gives 184.0565 GPa at that endpoint. The benchmark-consistent
abstract value is therefore retained, with the conflicting Table 2 extraction
recorded in provenance.

The publisher exposes no supporting data file, and the numerical FPMD stress
averages appear only as Figure 2 markers. The bundled dataset
`ca_perovskite_kawai_2014_table1_isochors` contains every printed Table 1
benchmark but is explicitly classified as a fitted-model grid, not as primary
observations. A direct independent refit is consequently not possible and no
synthetic refit record is created. The record's 0--150 GPa, 1000--4000 K
numerical envelope follows the stated pressure limit and published grid; the
paper's warning that the low-pressure/high-temperature corner may be outside
the cubic stability field or near melting remains attached to the validity
metadata.

## C04: CaSiO3 perovskite, Noguchi et al. (2013)

The primary source is Noguchi, Komabayashi, Hirose, and Ohishi,
*High-temperature compression experiments of CaSiO3 perovskite to lowermost
mantle conditions and its thermal equation of state*, *Physics and Chemistry
of Minerals* **40**, 81--91,
[doi:10.1007/s00269-012-0549-1](https://doi.org/10.1007/s00269-012-0549-1).
The complete article was checked through institutional Springer access.

### Identity, reference state, and equation

The starting material was pure CaSiO3 glass, mixed with Pt at 9:1 by weight;
Pt was the pressure standard and laser absorber. The source assigns every
accepted observation to cubic Ca-perovskite. It reports no peak splitting and
finds essentially identical 200-peak shapes at 700 and 1600 K. Three 700 K
observations with anomalously broad peaks are explicitly rejected. This is not
a room-temperature EOS: the authors deliberately set `Tr = 700 K`, above the
cubic transition, to avoid the tetragonal/cubic mixing in earlier work.

Equations (1)--(5) are represented directly as a second-order
Birch--Murnaghan reference isotherm plus referenced Mie--Gruneisen--Debye
thermal pressure:

`P(V,T) = P_BM2(V; V700, K700) + gamma(V)/V * [E_D(V,T)-E_D(V,700 K)]`,

with `gamma(V) = gamma700*(V/V700)^q`. The paper writes the equivalent
integrated relation
`theta = theta700*exp[(gamma700-gamma(V))/q]`. Equation (4) defines `n` as the
number of atoms per formula unit, hence `n = 5`. Model 1 in Table 2 gives
`V700 = 46.5(1) A^3`, `K700 = 207(4) GPa`, implicit `K0' = 4`,
`theta700 = 1300(500) K`, `gamma700 = 2.7(3)`, and `q = 1.2(8)`.

The paper calls `V` the cubic unit-cell volume. The established `Pm-3m`
CaSiO3 cell contains one formula unit, so the source value is simultaneously
`46.5 A^3/cell` and `46.5 A^3/formula unit`. Peritheos converts that public
cell volume to molar volume only inside the Debye-energy calculation. The
article does not refine atomic positions; the material's existing ideal
`Pm-3m` diffraction model is reused without attributing it to this experiment.

### Pressure scale, scope, and alternatives

Preferred model 1 uses Pt 111 and 200 spacings with Fei et al. (2004),
[doi:10.1016/j.pepi.2003.09.018](https://doi.org/10.1016/j.pepi.2003.09.018).
Table 1 also prints the same states reduced with Holmes et al. (1989), but that
scale produces the distinct model 3 (`V700 = 45.8 A^3`, `K700 = 238 GPa`,
`gamma700 = 2.8`, `q = 2.1`) and is not mixed into the production record.
Model 2 fixes `theta700 = 1000 K`; model 4 uses a different thermodynamic
thermal-pressure equation. Neither is promoted as a duplicate C04 record.

The measured table spans 51.3--127.2 GPa. The source describes the nominal
thermal scope as 700--2300 K; printed temperatures span 697--2278 K because
the 697--704 K external-heating series is treated as a single 700 K isotherm.
Executable validity follows the reported 700--2300 K cubic scope and the
observed pressure extrema. It is a marginal experimental envelope, not a
rectangular phase-stability claim. The paper does not define the confidence
level of its printed parameter errors and publishes no covariance matrix.

### Independent refit and redistribution limit

A complete local audit transcription of Table 1 contained 54 rows. The staged
source selection leaves nine externally heated rows for the 700 K BM2 fit and
42 laser-heated rows for the thermal fit; the three broad-peak rows are not
used. With no source-stated weighting rule, unweighted pressure-residual fits
through Peritheos provide the least-assumptive comparison:

| Parameter | Published model 1 | Peritheos staged refit |
|---|---:|---:|
| `V700` (A^3/cell) | 46.5 +/- 0.1 | 46.50049 +/- 0.13399 |
| `K700` (GPa) | 207 +/- 4 | 207.36877 +/- 3.53803 |
| `theta700` (K) | 1300 +/- 500 | 1292.42481 +/- 528.87816 |
| `gamma700` | 2.7 +/- 0.3 | 2.72149 +/- 0.29823 |
| `q` | 1.2 +/- 0.8 | 1.26461 +/- 0.77855 |

The BM2 pressure RMSE is 0.47281 GPa. The thermal-stage RMSE is 1.22099 GPa;
the rounded published parameters give 1.22277 GPa over the same 42 rows. This
is parameter and error-magnitude parity. The published coefficients remain
the executable record; no refit record is created.

The subscription article gives no reusable data license. Following the
repository's redistribution policy, the complete numeric table and local
audit transcription are not committed. The reproduction script therefore
checks the shipped parameterization itself, including the exactly zero
thermal increment at 700 K and the article's independently reported 300 K
extrapolation (`V300 = 45.8 A^3`, `K300 = 225 GPa`). Rounded model-1
coefficients give `45.7896 A^3` and `226.47 GPa`; the small bulk-modulus
difference is consistent with propagating coefficients printed to only two or
three significant digits.

## Phase Egg: Schulze et al. (2018)

### Source identity, structure, and the two printed BM3 sets

The record follows [Schulze et al. (2018)](https://doi.org/10.2138/am-2018-6562)
and the official MSA deposit `AM-18-126562`. The analyzed S5050 crystal is
printed as `Al0.98(1)Si0.92(1)O3OH1.39(5)`. The deposit writes the same
composition as `Al0.98(1) Si0.92(1) H1.39(5) O4`, so the canonical record
formula is `Al0.98Si0.92H1.39O4`; the hydroxyl notation does not mean O4.39.
The ideal phase name remains `AlSiO3OH`.

The deposited ambient structure is the conventional monoclinic `P21/n` cell,
space-group number 14, with `Z = 4`, `a = 7.1835(2) A`, `b = 4.3287(2) A`,
`c = 6.9672(2) A`, `beta = 98.202(2) deg`, and `V = 214.431(13) A^3`.
All Al, Si, and O coordinates come from that single-crystal refinement. The
authors report that fixed and refined Al/Si occupancies were indistinguishable
and adopt the fully occupied model. Its site multiplicities therefore expand
to ideal Al4Si4O16 per cell even though the analyzed central composition is
Al3.92Si3.68H5.56O16. The 2018 CIF contains no hydrogen atom; the H coordinate
used in the article's structural drawing and in Peritheos is explicitly traced
to Schmidt et al. (1998), Table 4. It is not represented as a 2018 refinement.

Table 2 contains two BM3 parameter sets because it compares different data,
not because the 2018 single-crystal observations have two alternative fits:

| Table 2 column | Data represented | `V0` (A^3/cell) | `K0` (GPa) | `K0'` |
|---|---|---:|---:|---:|
| This study | 15 synchrotron single-crystal rows | 214.08(17) | 153(8) | 8.6(1.2) |
| Vanpeteghem refitted | Vanpeteghem et al. (2003) powder rows | 211.41(11) | 155(5) | 6.7(5) |

Only the first is stored as `phase_egg_schulze_2018_bm3_1`. The comparison
reduction is source-lineage evidence, not a second Phase Egg EOS contributed by
this audit. The source's abstract and EOS prose give the first set's `K0'`
error as 1.3 while Table 2 gives 1.2; Peritheos retains the tabulated 1.2 and
records the contradiction. All three coefficients were fitted in EosFit7c;
none was fixed or adopted. The paper does not state a confidence level or
publish a covariance matrix.

### Numerical reproduction and calibration scope

The complete 16-row Table 1 transcription is embedded. The authors exclude
the in-house ambient point to avoid an inter-technique bias, leaving the 15
synchrotron observations from 1.09 to 23.33 GPa for the EOS fit. With the
standard Eulerian third-order Birch--Murnaghan equation, the rounded published
coefficients predict 23.0656 GPa at the final `V = 193.75 A^3` state, compared
with 23.33 GPa observed. Their pressure RMSE over the 15 fitted rows is 0.2270
GPa and the maximum absolute residual is 0.4535 GPa.

An independent errors-in-variables refit uses the printed pressure and volume
uncertainties and `absolute_sigma=True`:

| Parameter | Published | Peritheos refit |
|---|---:|---:|
| `V0` (A^3/cell) | 214.08 +/- 0.17 | 214.0811 +/- 0.0747 |
| `K0` (GPa) | 153 +/- 8 | 152.8459 +/- 3.4520 |
| `K0'` | 8.6 +/- 1.2 | 8.5975 +/- 0.5104 |

The refit has a 0.2117 GPa pressure RMSE and reduced chi-square 5.486; every
coefficient agrees within combined two-standard-deviation uncertainty, so the
campaign classifies it as `parity`. The fitted `V0` is a zero-pressure
extrapolation and need not equal the separately measured ambient structure
volume.

Pressures were assigned from ruby R1 fluorescence using Dewaele et al. (2008),
with neon as the pressure medium. The corresponding `A = 1920 GPa`, `B = 9.61`
power-law calibration is bundled as `ruby_dewaele_2008`. Observation-level
pressure re-reduction is not possible because neither Table 1 nor the deposit
provides the row-wise ruby wavelengths.
<a id="palladium-baty-2024"></a>

## Palladium: Baty et al. (2024)

### Source and exact data check

- Primary article: Baty et al., *Palladium at high pressure and high
  temperature: A combined experimental and theoretical study*, Journal of
  Applied Physics **135**, 075103 (2024),
  [doi:10.1063/5.0179469](https://doi.org/10.1063/5.0179469).
- Primary data: all 78 room-temperature `P`--`V` observations in official
  supplementary Table S1,
  [doi:10.60893/figshare.jap.c.7021680](https://doi.org/10.60893/figshare.jap.c.7021680).
  The official supplementary PDF and its deposited LaTeX source were checked
  independently. Every bundled pressure and atomic-volume value agrees with
  the source; the library's conventional-cell volume is exactly four times
  the printed fcc atomic volume.
- Model: Equation (1) is the standard third-order Birch--Murnaghan relation
  written using `eta = rho/rho0 = V0/V`. It is algebraically identical to the
  Peritheos `BM3` implementation.
- Published experimental fit: `V0 = 14.72(4) A^3/atom`, `K0 = 190(3) GPa`,
  and `K0' = 5.3(2)`, described as a fit to the authors' 0--80 GPa data. In
  conventional-cell units, `V0 = 58.88(16) A^3`.

These checks exclude transcription order, atomic-to-cell conversion, and EOS
formalism as explanations for the failed validation. The supplementary table
ends at 73.7963 GPa; “80 GPa” is the paper's rounded experimental-range
description, not a missing 80 GPa table row.

### What supplementary Table S3 actually compares

Table S3 is easy to misread because it places six EOS curves side by side at a
common set of volumes. It is **not** a table of experimental or calculated fit
observations. Each pressure column is generated by evaluating an already fitted
EOS, so fitting Table S3 would merely recover the parameters used to construct
it. The represented curves are now separate Palladium records:

| Source | Physical basis | EOS and conventional-cell parameters | Data available for an independent refit |
|---|---|---|---|
| Guigue et al. (2020) | 300 K experiment, Pd in Ne to 65 GPa | Vinet: `V0=59.00(16) A^3`, `K0=162(2) GPa`, `K0'=6.2(9)` | Plot only |
| Fedotenko et al. (2020) | 300 K experiment, Pd in Ne | BM3: `V0=58.868(2) A^3`, `K0=157(3) GPa`, `K0'=9.9(4)`; also BM2 with `K0=203(3) GPa` and `K0'=4` fixed | Complete 15-row Table 1 bundled |
| Frost et al. (2023) | 300 K experiment, three Pd-in-Ne runs to 182 GPa | Vinet: `V0=58.678(73) A^3`, `K0=189.3(3.0) GPa`, `K0'=5.473(63)`; BM3: `V0=58.588(76) A^3`, `K0=197.5(3.3) GPa`, `K0'=4.996(66)` | Complete 93-row supplementary Tables I-III bundled |
| Baty et al. (2024), this work in Table S3 | **0 K DFT cold curve** | BM3: `V0=58.4668 A^3`, `K0=195.0(2.5) GPa`, `K0'=5.1(1)` | Calculated input grid not tabulated |

Baty's separate **experimental** BM3 (`V0=58.88(16) A^3`, `K0=190(3)
GPa`, `K0'=5.3(2)`) is not one of the Table S3 columns. It is the record tested
against supplementary Table S1 below. Frost's Vinet and BM3 alternatives are
now both independent library records because Table I of the primary paper
reports both complete parameter triplets and uncertainties.

### What the Frost paper and supplement resolve

The complete Frost article corrects two earlier uncertainties. Its pure-Pd
measurements extend to **182 GPa**, not approximately 62 GPa. Figure 1 combines
three neon-loaded, hydrogen-excluded DAC experiments: Runs 1 and 2 use a
tungsten marker with the Dewaele et al. (2004) EOS, while Run 3 uses a gold
marker with the Dorfman et al. (2012) EOS. The source assigns `0.1 GPa`
pressure uncertainty to all runs and derives volume uncertainties from the Le
Bail refinements.

Both Frost parameter sets are regressions of the same experimental data. The
authors used nonlinear Marquardt--Levenberg least squares that accounts for
both pressure and volume uncertainties; the Vinet form is plotted as the black
curve in Figure 1, while Table I additionally gives the BM3 alternative. This
is experimental fitting, not fitting of calculated volumes.

The supplement supplies the full experimental input: 26 BX18 rows, 20 BX32
rows, and 47 high-pressure rows, for 93 Pd observations from 2.99 to 181.79
GPa. The first two loadings use W and the high-pressure loading uses Au. The
library bundles the printed pressure and lattice parameter for every row and
calculates the conventional-cell volume and its first-order uncertainty as
`V=a^3` and `sigma(V)=3a^2 sigma(a)`.

The independent ambient powder measurement is not silently added to the fit.
The supplement introduces `a=3.89088(2) A` explicitly "in addition to Tables
I, II and III," while those three tables contain the compression dataset. Its
equivalent `V=58.9038 A^3/cell` is larger than both fitted extrapolations,
`58.678(73) A^3` for Vinet and `58.588(76) A^3` for BM3. The unusual fitted
`V0` values therefore come from the source and are not an atomic-versus-cell
conversion error.

### Frost direct refits

Peritheos fits all 93 rows with the source-stated `sigma(P)=0.1 GPa` and the
propagated volume errors. Both models are **similar**, but not strict
uncertainty parity because `K0'` is just outside the combined two-sigma
interval:

| Model | Published | Peritheos errors-in-variables refit | Published/refit pressure RMSE |
|---|---|---|---:|
| Vinet | `V0=58.678(73) A^3`, `K0=189.3(3.0) GPa`, `K0'=5.473(63)` | `V0=58.7665(145) A^3`, `K0=184.191(613) GPa`, `K0'=5.60371(1381)` | `0.772 / 0.276 GPa` |
| BM3 | `V0=58.588(76) A^3`, `K0=197.5(3.3) GPa`, `K0'=4.996(66)` | `V0=58.6846(141) A^3`, `K0=191.477(622) GPa`, `K0'=5.14738(1393)` | `0.845 / 0.293 GPa` |

This difference is not a Peritheos optimizer or EOS-formalism issue. An
independent SciPy orthogonal-distance regression gives the same coefficients
as the library to the displayed precision for both Vinet and BM3. An
unweighted pressure-residual fit instead gives Vinet
`(58.5924, 193.437, 5.38154)` and BM3 `(58.4774, 202.948, 4.88068)`, placing
the published coefficients between two ordinary objective conventions rather
than reproducing either. Adding the separate ambient point moves the fit away
from the publication, confirming that its omission is not the cause.

The stated point uncertainties also do not describe the complete scatter. The
standard errors-in-variables fits have reduced chi-square values of 19.0
(Vinet) and 22.1 (BM3). Separate Vinet fits range from
`(V0, K0, K0')=(58.730, 190.83, 5.103)` for BX18 through
`(59.032, 169.07, 6.129)` for BX32 to `(57.860, 226.12, 4.813)` for the Au
high-pressure loading. Those individual triplets are strongly correlated
because the runs cover different pressure intervals, but their spread exposes
run-level/systematic variation that is far larger than the printed point
errors.

The remaining mismatch is therefore best attributed to an undocumented detail
of the source regression: weighting/error transformation, treatment of
run-level scatter, or use of higher-precision pre-publication values. As a
sensitivity check, scaling the propagated volume errors to roughly 40--45% of
their tabulated values moves both refits close to the published triplets, but
the source gives no justification for that scaling and Peritheos does not adopt
it. Table rounding alone is too small to explain the result. These records are
accordingly classified as `similar`, not `parity` and not a failed EOS record.

### Fedotenko direct refits

The Fedotenko paper is the one comparison study that publishes its numerical
pure-Pd input. Peritheos fits all 15 Table 1 rows. All three BM3 coefficients
are free because the paper does not report a fixed-`V0` constraint; the BM2
calculation fixes its source-defined `K0'=4` and uses the reported reference
volume:

| Model | Published | Peritheos refit | Pressure RMSE |
|---|---|---|---:|
| BM3 | `V0=58.868(2) A^3`, `K0=157(3) GPa`, `K0'=9.9(4)` | `V0=58.8537 A^3`, `K0=158.668 GPa`, `K0'=9.7835` | `0.0750 GPa` |
| BM2 (`K0'=4`) | `K0=203(3) GPa` | `K0=202.667 GPa` | `0.9335 GPa` |

Both achieve parity. The much poorer BM2 residual is expected: Table 2 presents
it as a constrained comparison, whereas Figure 1 and the main discussion use
BM3. The primary source itself contains a small internal inconsistency: Table 2
prints `V0=58.863(2) A^3`, while the ambient row of Table 1, the fit paragraph,
and the Figure 1 caption all print `58.868(2) A^3`. The repeated value is stored,
and changing it to the isolated Table 2 value does not alter the scientific
conclusion.

### Refit and objective sensitivity

The source does not state the residual variable, regression weights, whether
`V0` was fixed, or how the quoted parameter uncertainties were calculated.
The following diagnostics therefore test the ordinary plausible conventions
without treating any one of them as the undocumented source protocol:

| Calculation | Points | `V0` (A^3/cell) | `K0` (GPa) | `K0'` | Pressure RMSE (GPa) |
|---|---:|---:|---:|---:|---:|
| Publication evaluated on Table S1 | 78 | 58.88 | 190 | 5.3 | 1.273 |
| Unweighted pressure residuals | 78 | 59.960 | 152.057 | 6.099 | 0.704 |
| Unweighted volume residuals | 78 | 59.470 | 178.824 | 5.010 | 0.777 |
| Maximum-error EIV diagnostic | 78 | 59.639 | 169.063 | 5.383 | 0.734 |
| Pressure residuals, published `V0` fixed | 78 | 58.88 fixed | 207.477 | 4.296 | 0.972 |
| Published `V0` fixed, `P >= 40 GPa` | 28 | 58.88 fixed | 193.664 | 5.024 | 0.754 |

The maximum-error errors-in-variables calculation assigns the source's largest
stated errors to every row: `sigma(P) = 0.12 GPa` and
`sigma(V) = 0.01 A^3/atom`. This is deliberately more generous than the paper,
which says both are upper limits and gives a smaller pressure error at lower
pressure. Even so, its reduced chi-square is 8.11 and its coefficients do not
recover the published triplet.

Changing the residual convention matters because all three BM3 coefficients
are correlated, but it does not explain the source result. Fixing the published
ambient volume likewise exchanges `K0` against `K0'` rather than recovering
both. The high-pressure-only calculation comes close, but the paper explicitly
describes a 0--80 GPa fit and provides no row mask or justification for
discarding the lower-pressure observations. It is therefore only a localization
test, not a parity result.

### Cross-check against Frost et al. (2023)

Figure 3 does not plot the individual observations of Frost et al.; its black
line is their published room-temperature Vinet EOS with
`V0 = 58.678 A^3/cell` (`14.6695 A^3/atom`), `K0 = 189.3 GPa`, and
`K0' = 5.473`. Re-evaluating exactly that curve at every Table S1 pressure
confirms the offset visible in the figure:

| Comparison with Baty Table S1 | Result |
|---|---:|
| Observations above the Frost curve | 77 of 78 |
| Mean `V_observed - V_Frost` | `0.07687 A^3/atom` |
| RMS `V_observed - V_Frost` | `0.09355 A^3/atom` |
| Maximum `V_observed - V_Frost` | `0.20003 A^3/atom` |
| Pressure-equivalent RMSE of the Frost curve | `2.039 GPa` |

The mean volume offset is 7.7 times the `<0.01 A^3/atom` uncertainty quoted
for the Baty data, and the maximum is 20 times that bound. The corresponding
pressure residual is also far larger than the stated `0.06--0.12 GPa`
pressure uncertainty. The discrepancy is concentrated most strongly at low to
intermediate pressure: from 10 to 20 GPa, the Baty observations exceed the
Frost curve by `0.16543 A^3/atom` on average.

This cannot be assigned to a 0 K versus 300 K comparison because both the
Frost curve and the Baty observations are room-temperature results. Baty's
published BM3 curve is not numerically identical to Frost's curve: its larger
`V0` shifts it upward by about `0.04455 A^3/atom` over the measured pressures.
Nevertheless, the closeness of the two papers' reported `K0` and `K0'`, while
the new Table S1 points sit systematically above both curves, makes the
source-fit discrepancy more—not less—significant. It does not prove that the
Frost result was used as a constraint; the source gives no such fitting detail.

### Internal scatter and likely origin

After sorting Table S1 by pressure, the measured volume increases at 12 adjacent
steps even though a static compression curve must be monotonic. The largest
increase is `0.04197 A^3/atom` (`0.16787 A^3/cell`), more than four times the
stated atomic-volume upper uncertainty before accounting for the expected
decrease over the accompanying pressure increment. This indicates that the
printed pointwise uncertainty bounds do not capture the complete run-to-run or
systematic scatter relevant to a regression.

The source fit may therefore have used unreported weighting, constraints,
averaging, or row selection; a pressure-scale/reduction problem or source-side
fitting or reporting inconsistency also remains possible. Rounding cannot
account for the result because Table S1 prints five pressure decimals and eight
volume decimals, far beyond the stated measurement precision. The exact origin
cannot be selected among these options without the authors' fit script, row
weights, covariance, row mask, and underlying W-marker observations.

### Parity conclusion

Parity is not achieved. The all-row Peritheos pressure refit moves `K0` by
19.97%, from `190` to `152.057 GPa`, and fits the checked table substantially
better than the published coefficients. Unlike the earlier Mo2C problem, there
is no single source-mandated correction that resolves the discrepancy.

The library therefore retains the published parameterization as a literature
record and documents the failure rather than silently replacing it. A separate
refit record is not added at this stage because several defensible objectives
produce materially different coefficients and the paper does not identify
which objective its reported values represent.

<a id="ringwoodite-katsura-2004"></a>

## Mg2SiO4 ringwoodite: Katsura et al. (2004) blocked MGD audit

### Primary source, identity, and observations

[Katsura et al. (2004)](https://doi.org/10.1029/2004JB003094) report pure
Mg2SiO4 ringwoodite measured in a Kawai-type multianvil apparatus. The official
free-access Wiley version of record identifies the cubic conventional cell,
`a0 = 8.0663(7) A` and `V0 = 524.8(1) A^3` at 300 K and 0 GPa. Ringwoodite has
eight Mg2SiO4 formula units per conventional cell, so each formula unit has
seven atoms. The composition, phase, and cell basis match the existing
`ringwoodite` material, but the Katsura parameterization is independent of the
existing Meng et al. static record.

All 127 P-T-V/V0 observations in the article's Table 2 are preserved in
`docs/data/ringwoodite-katsura-2004-table2-pvt.csv`. They cover four runs,
-0.01 to 23.18 GPa, 300 to 2000 K, and V/V0 from 0.9106 to 1.0001. The
parenthetical pressure and normalized-volume errors are split into separate
columns exactly as printed. The CSV SHA-256 is
`aaae5542dcc46eee4fe626ca17e8b24ba1f791028f70cfb4d502036a18211649`.
The article calls the parenthetical values errors but does not state a
confidence level, so they are not relabeled as standard deviations. Wiley
marks the article free access, but the page states Copyright 2004 American
Geophysical Union and does not identify a reuse license; Peritheos therefore
makes no broader licensing claim for the transcription.

Pressures were calculated from simultaneous MgO cell volumes with the
[Matsui, Parker, and Leslie (2000)](https://doi.org/10.2138/am-2000-2-308)
high-temperature MgO EOS. The article reports typical pressure uncertainty up
to 0.04 GPa at high temperature, and the table prints row-wise errors. It does
not publish the corresponding MgO volumes, so observation-level pressure
recalculation is impossible even though the calibration reference is
unambiguous.

### Published equations and parameter status

Equation (4) is the standard 300 K third-order Birch-Murnaghan isotherm. The
article fixes `K0 = 182 GPa` to the Meng et al. (1994) result and fits the 300 K
data for `K0' = 4.6(2)`. Its MGD section defines

`P(V,T) = P(V,300) + Pth(V,T) - Pth(V,300)`,

with `Pth = gamma E_th / V`, Debye thermal energy, and
`gamma = gamma0 (V/V0)^q`. It reports fitted `theta = 846(26) K`,
`gamma0 = 1.93(3)`, and `q = 3.5(3)`. `V0` is measured, `K0` is adopted and
fixed, `K0'` is fitted separately, and the three thermal coefficients are
fitted. No covariance matrix, fit statistic, objective, or weighting rule is
given. The paper defines `n` as the number of atoms per formula unit but does
not print the numeric value used, and it does not print a volume law for the
Debye temperature.

### Numerical audit and blocker

The deterministic audit is
`scripts/reproduce_katsura_2004_ringwoodite.py`. It converts the conventional
cell to Peritheos' internal molar-energy volume with `Z = 8` and evaluates both
available Debye-temperature laws. With the chemically required `n = 7`, the
published coefficients give a 1.889936 GPa pressure RMSE over all 127 rows and
a 2.013420 GPa RMSE over the 110 rows above 304 K; the maximum residual is
3.263616 GPa. The alternate direct variable-exponent Debye-temperature law is
no better (1.904889 GPa overall). These discrepancies are far larger than the
printed row errors.

An unweighted Peritheos refit with `n = 7`, the reference BM3 fixed, and
`theta0`, `gamma0`, and `q` free gives `1035.4123 K`, `1.399632`, and
`2.826559`, with a 0.296233 GPa pressure RMSE. In particular, `gamma0` moves
far outside the paper's 0.03 reported error. The separate 300-304 K BM3 check
gives `K0' = 4.83790 +/- 0.22972`, consistent with the reported 4.6(2), so the
failure is localized to the thermal normalization rather than the cell basis
or BM3 convention.

A diagnostic substitution of `n = 5` with all published coefficients nearly
restores the curve: 0.333829 GPa overall and 0.242193 GPa for heated rows.
However, `n = 5` is neither reported nor chemically possible for Mg2SiO4. It is
evidence of an unresolved hidden normalization or source-side implementation
detail, not permission to invent a coefficient.

Candidate C05 is therefore blocked. No production Katsura EOS record and no
refit record are added. Resolution requires the authors' fitting code or an
authoritative statement of the numeric `n`, Debye-temperature volume law, and
thermal-energy/volume normalization used to obtain the published coefficients.
