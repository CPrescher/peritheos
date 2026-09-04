# Qin et al. (2023) calcium-ferrite NaAlSiO4 EOS audit

This note documents the source audit and numerical reproduction for two
composition-specific room-temperature EOS records from Qin et al. (2023),
"Crystal structure of calcium-ferrite type NaAlSiO4 up to 45 GPa,"
*American Mineralogist* 108, 2331-2337,
[doi:10.2138/am-2022-8432](https://doi.org/10.2138/am-2022-8432).

The final article was checked together with official MSA deposit
[AM-23-128432](https://msaweb.org/MSA/AmMin/TOC/2023/Dec2023_data/AM-23-128432.zip).
The downloaded archive has SHA-256
`465ff9b139fba3c60539f4bf1e0a86d877a716ab4522cefca682dd1f4fa1f07d`;
its `8432_supp.pdf` and `8432qin.cif` files have SHA-256 values
`f78e80d89fb11bb070cc3757c554e01adab3d4dcbf1ace93c2fc166aa5f3d1c7`
and `18b3a58bcefba6797dfee976025604f8165ada6678a8113e0918e37eb7fc03cc`,
respectively. LitCurate was not used as scientific authority.

## Material and volume identity

The specimens are synthesis run 5K2124, with electron-microprobe composition
Na0.93Al1.02Si1.00O4, and synthesis run 5K2681, with composition
Na0.88Al0.99Fe0.13Si0.94O4. The latter contains 90% ferric iron according to
the Mossbauer result cited from Wu et al. (2017). They are kept as separate
materials rather than folded into the existing nominal NaAlSiO4 record.

Both samples remain orthorhombic Pbnm (space-group number 62) with Z=4 over
the reported compression series. Every Table S3-S5 volume is therefore one
conventional Pbnm unit-cell volume in A^3, not a molar or formula-unit volume.

The official CIF contains idealized positional refinements for the Fe-free
crystal. Its 1.7 GPa, 293(2) K block is the Fe-free diffraction model. The
Fe-bearing article and Table S4 establish Pbnm, Z=4, and the lattice, but the
deposit contains no Fe-bearing fractional-coordinate refinement. Its material
record consequently uses the deposited isostructural Fe-free coordinates as
an explicitly qualified topology proxy. The ideal occupancies are not changed
to emulate the microprobe formula: no Fe/Al/Si site distribution is published,
and none is inferred. This structural limitation does not alter the
composition-specific EOS coefficients or measured conventional-cell volumes.

## Published parameterizations

The article states that the P-V data were fit by error-weighted least squares
in EosFit7c using a third-order Birch-Murnaghan EOS. Peritheos uses the standard
Eulerian finite-strain form

\[
P(V)=\frac{3K_0}{2}\left[\left(\frac{V_0}{V}\right)^{7/3}
-\left(\frac{V_0}{V}\right)^{5/3}\right]
\left\{1+\frac{3}{4}(K'_0-4)
\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]\right\}.
\]

Online Materials Table S5 and the article report:

| Record | Composition | V0 (A^3, Z=4 cell) | K0 (GPa) | K0 prime | Observed range |
|---|---|---:|---:|---:|---:|
| `na093al102si100o4_calcium_ferrite_qin_2023_bm3_1` | Na0.93Al1.02Si1.00O4 | 241.6(1) | 220(4) | 2.6(3) | 0.0001-41.0 GPa |
| `na088al099fe013si094o4_calcium_ferrite_qin_2023_bm3_1` | Na0.88Al0.99Fe0.13Si0.94O4 | 244.2(2) | 211(6) | 2.6(3) | 0.0001-44.0 GPa |

All three coefficients were fitted; none was fixed. The source gives neither a
confidence convention for the parenthetical coefficient errors nor a parameter
covariance matrix. Peritheos stores the errors as generic published
uncertainties, stores covariance as missing, and does not assume independence.
The compression experiments are described as room temperature. The explicit
temperature in the official CIF is 293(2) K, which is used as the record
reference; no row-wise temperatures are available.

## Primary observations and numerical reproduction

All 22 Table S3 rows and all 10 Table S4 rows are bundled as direct CSV
transcriptions. They preserve pressure, a, b, c, conventional-cell volume, and
every printed parenthetical lattice/volume uncertainty. The source publishes no
row-wise pressure uncertainty. The unexplained superscript `a` on the two
1.7 GPa rows is retained as a flag without interpretation.

Evaluating the authors' rounded coefficients at the printed volume rows gives:

| Sample | Rows | Pressure RMSE (GPa) | Maximum absolute residual (GPa) | High-pressure calculated / observed (GPa) |
|---|---:|---:|---:|---:|
| Fe-free | 22 | 0.326278 | 0.643858 | 40.609612 / 41.0 at 207.3 A^3 |
| Fe-bearing | 10 | 0.388733 | 0.790621 | 43.209379 / 44.0 at 206.6 A^3 |

The independent diagnostic refits use every printed row:

| Sample and objective | V0 (A^3) | K0 (GPa) | K0 prime | Fit statistic |
|---|---:|---:|---:|---:|
| Fe-free, unweighted pressure residual | 242.013801 | 209.585089 | 3.103588 | pressure RMSE 0.274888 GPa |
| Fe-free, volume residual / printed volume uncertainty | 242.027302 | 208.876126 | 3.136837 | reduced chi-square 1.118947 |
| Fe-bearing, unweighted pressure residual | 244.326475 | 205.775568 | 2.916648 | pressure RMSE 0.330565 GPa |
| Fe-bearing, volume residual / printed volume uncertainty | 244.057323 | 201.299582 | 3.107154 | reduced chi-square 0.536114 |

These results are not replacement EOS records. The authors specify
"error-weighted" EosFit7c regression but do not publish pressure uncertainties,
parameter covariance, or the exact weight settings. Coefficient-level parity
therefore cannot be demanded from the printed tables. The published curves do,
however, reproduce all primary observations within the stated residual bounds,
and the complete equation, coefficients, units, volume basis, and reference
state are executable. The deterministic calculation is in
`scripts/reproduce_qin_2023_naalsio4.py`.

## Pressure scale and source inconsistencies

Both diffraction samples were loaded in neon with Pt foil, and the article
identifies the internally consistent Pt EOS of Fei et al. (2007),
[doi:10.1073/pnas.0609013104](https://doi.org/10.1073/pnas.0609013104), as the
pressure calibration. Tables S3-S4 provide only reduced pressures. They do not
provide row-wise Pt volumes or lattice parameters, so observation-level
pressure recalculation is not possible. The exact Fei Pt scale is also not yet
a shared-library EOS record; the calibration is marked partially resolved.

Three transcription-relevant inconsistencies are retained rather than repaired:

- At 38.7 GPa, Fe-free Table S3 prints `a=9.689(2) A`, which is non-monotonic.
  The printed a, b, and c values multiply to 210.1499 A^3, not the printed
  `V=208.8(2) A^3`. EOS reproduction uses the reported V column.
- Fe-free Table S3 prints the ambient volume as 241.5(3) A^3, whereas the
  article prose and Table S5 print the measured value as 241.5(1) A^3. The
  observation dataset preserves the Table S3 uncertainty.
- The superscript `a` attached to both 1.7 GPa rows is not defined in the
  Online Materials.

The declared validity envelopes are the actual Table S3-S4 pressure and volume
ranges at 293 K. They are experimental coverage, not rectangular phase-stability
guarantees or permission to extrapolate beyond the data.
