# Thompson et al. (2017) hydrogen-centered epsilon-FeOOH EOS

## Scope and source

This record audits Thompson, Campbell, and Tsuchiya (2017), *Elasticity of
epsilon-FeOOH: Seismic implications for Earth's lower mantle*, Journal of
Geophysical Research: Solid Earth 122, 5038–5047,
[doi:10.1002/2017JB014168](https://doi.org/10.1002/2017JB014168). The final
publisher article and its official five-page Supporting Information S1 were
checked. The supplement PDF used for transcription has SHA-256
`55653c46d1d48762fcf7a9ec2582ae91f4319ce8827075ba0a0ca3721be080f5`.

The executable record is
`e_feooh_hc_low_spin_thompson_2017_bm3_1` in the phase-specific material
`e_feooh_hc_low_spin`. It is a static 0 K, low-spin ferromagnetic,
hydrogen-centered Pnnm branch. It is deliberately separate from the existing
hydrogen-off-center P21nm `e_feooh` material.

## Exact volumetric model

Section 3.2 prints the standard third-order Birch–Murnaghan pressure equation:

\[
P(V)=\frac{3K_0}{2}\left[\left(\frac{V_0}{V}\right)^{7/3}
-\left(\frac{V_0}{V}\right)^{5/3}\right]
\left\{1+\frac{3}{4}(K'_0-4)
\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]\right\}.
\]

Table 2 gives the preferred self-consistent-U hydrogen-centered fit over
30–140 GPa:

| Parameter | Published value | Peritheos value |
|---|---:|---:|
| Conventional-cell \(V_0\) | 58.62(5) Å³ | 58.62 Å³ |
| \(K_0\) | 223(2) GPa | 223 GPa |
| \(K'_0\) | 4.07(3) | 4.07 |

The article explicitly states that \(V_0\) was fitted freely. No coefficient
was fixed. Parentheses are defined only as uncertainty on the last digit; the
confidence convention, weights, covariance matrix, and fit statistic are not
reported. The volume is a conventional Pnnm unit cell with two FeOOH formula
units. The 0 GPa reference volume is an extrapolated static-lattice reference;
the supported calculation interval begins at 30 GPa.

There is no experimental pressure scale. Pressures are static first-principles
stresses from PBE GGA calculations with a self-consistent Hubbard U, so
pressure calibration and observation-level recalibration are not applicable.

## Primary observations and numerical reproduction

Official Supporting Information Table S1 contains 17 calculated structural
states from 0.01 to 139.99 GPa. All columns and rows are preserved in
`epsilon-feooh-thompson-2017-supplement-table-s1.csv`. The 12 rows from 30.01
to 139.99 GPa are the hydrogen-centered fit interval. Table S1 reports no
row-level uncertainties, and its values are rounded to two decimal places.

Evaluating the published HC BM3 curve at those 12 printed volumes gives a
pressure RMSE of 0.10819 GPa and a maximum absolute residual of 0.16968 GPa. A
diagnostic unweighted pressure-residual refit to the rounded rows gives
\(V_0=58.61201868\) Å³, \(K_0=223.41827131\) GPa, and
\(K'_0=4.06653530\), with a 0.09518 GPa RMSE. The source coefficients are
retained: the authors did not publish their weights or unrounded calculation
outputs, so this diagnostic does not justify a separate refit record.

Table 3 supplies independent density benchmarks. Inverting the stored BM3 and
using two FeOOH formula units per cell gives 6.05438 g cm⁻³ at 60 GPa and
6.42998 g cm⁻³ at 90 GPa, compared with the published 6.055 and 6.428 g cm⁻³.

## Excluded low-pressure parameterization

The paper also reports a separate HOC P21nm fit over 0–20 GPa:
\(V_0=57.43(7)\) Å³, \(K_0=188(4)\) GPa, and \(K'_0=5.19(12)\). It is not
bundled. With the paper's own BM3 equation, the printed \(V_0\) predicts
negative pressure at the supplement's near-zero-pressure state
(0.01 GPa, 59.88 Å³), and its pressure RMSE across the five printed
0.01–20.01 GPa rows is 7.765 GPa. That RMSE is only an inconsistency
diagnostic: Table S1 already shows a centered hydrogen bond from 10 GPa onward,
despite Table 2 labeling the HOC fit interval as 0–20 GPa, so the exact
metastable HOC row selection is not recoverable. The same V0 is printed in the
publisher article and the author's dissertation, so there is no defensible
source-based correction. Peritheos does not infer a replacement coefficient.

## Volumetric versus elastic results

Only the unit-cell BM3 fit is an EOS record. The paper's single-crystal elastic
constants are obtained from stress–strain calculations, and its aggregate
bulk and shear moduli and sound velocities are Voigt–Reuss–Hill averages.
Those elastic-density results are not substituted for the volumetric \(K_0\)
or \(K'_0\), are not fitted as pressure–volume observations, and are not
exposed as additional EOS records. The two Table 3 densities are used only as
independent numerical checks of the volumetric curve.

## Diffraction structure

Thompson et al. publish the Pnnm identity, lattice metrics, and bond geometry,
but not fractional atomic coordinates. The phase-specific diffraction model
therefore uses the complete 50 GPa static low-spin Pnnm structure in Table 1
of Insixiengmay and Stixrude (2023),
[doi:10.2138/am-2022-8839](https://doi.org/10.2138/am-2022-8839). A documented
origin shift maps the source coordinates onto standard Pnnm sites Fe 2a, O 4g,
and H 2d. Their multiplicities give Fe₂O₄H₂, or \(Z=2\). This structure is a
diffraction reference from a different calculation and is not rescaled or
treated as an observation in the Thompson EOS fit.
