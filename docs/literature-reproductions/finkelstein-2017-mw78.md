# Finkelstein et al. (2017) Mw78 compression audit

## Scope and sources

This audit covers DOI `10.2138/am-2017-5966`, “Single-crystal equations of
state of magnesiowüstite at high pressures,” American Mineralogist 102,
1709-1717. LitCurate was used only to discover the candidate. Scientific
authority comes from:

- the final typeset article hosted by the authors:
  `https://web.gps.caltech.edu/~jackson/pdf/Finkelstein2017_AmMin.pdf`
  (SHA-256 `b4caf3b4ee71869be0c7d69dd3f557745f0ec93b2a03d39f8bc10ed87a3fb6e3`);
- the official American Mineralogist deposit AM-17-85966:
  `https://msaweb.org/MSA/AmMin/TOC/2017/Aug2017_data/AM-17-85966.zip`
  (SHA-256 `80f895e66372fe27ba8104a7b65960fea67846b3aea10225c644df77ffc845d2`);
- the deposit's `HarvestTable 2.xlsx`, SHA-256
  `2b7320f0f6fce96045edf19d3cd9ea583032446173e558bb204f82654c755264`;
  and
- the deposit's supplemental PDF, SHA-256
  `f8c853ff3ef7a66788f6b225cde3fc48ab7dd4dc5654085056dd7af6043bd6f5`.

The final PDF, accepted manuscript, official supplement, typeset Table 2 PDF,
and Table 2 workbook were cross-checked. The official workbook has a missing
closing parenthesis in two `16.199(65)` cells; the typeset Table 2 PDF resolves
the intended value.

## Material identity, spin, and structure

The DAC fragments came from the same synthetic crystal initially described as
nominal `(Mg0.22Fe0.78)O`, or Mw78. Room-pressure synchrotron Mössbauer
spectroscopy separated the iron population into 71.6 mol% high-spin Fe2+ and
4.6 mol% tetrahedrally coordinated high-spin Fe3+. Charge balance gives the
defect-normalized composition
`Mg0.215Fe2+0.716Fe3+0.046□0.023O`, or
`Mg0.215Fe0.762□0.023O` when valence is not shown.

The room-pressure Bragg peaks can be indexed in the B1 cubic cell, with
`a=4.2898(4) A` and `V=78.94(2) A^3`. The conventional rocksalt cell has
`Z=4`. Weak satellites and diffuse scattering show that this average cell does
not fully describe the defect clusters. The paper does not publish an ordered
defect-cluster refinement. The bundled `Fm-3m` atom sites are therefore an
explicit average scattering proxy, not a claim that tetrahedral Fe3+ occupies
the average octahedral cation site.

The helium-loaded crystal remains metrically cubic through 55.5 GPa. The
neon-loaded crystal remains cubic-indexable through about 20 GPa, then becomes
increasingly hexagonally distorted as the ruby pressure dispersion rises. The
authors attribute the different compression behavior and distortion to
nonhydrostatic stress. Their “hexagonal” description is an indexing cell only:
no space group or atomic model is assigned.

The only spin measurement in this study is the room-pressure spectrum. No spin
state is assigned to the compressed cubic or distorted branches.

## Equation and published parameterizations

Page 1710 prints the exact third-order Birch-Murnaghan form:

\[
P(V)=\frac{3K_{0T}}{2}
\left[\left(\frac{V_0}{V}\right)^{7/3}
-\left(\frac{V_0}{V}\right)^{5/3}\right]
\left\{1+\frac{3}{4}(K'_{0T}-4)
\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]\right\}.
\]

Table 3 labels these as 300 K fits. No row-wise temperatures were measured, so
300 K is a nominal room-temperature reference.

| Branch | Fit rows | V0 (A3) | K0T (GPa) | K'0T | Treatment |
|---|---:|---:|---:|---:|---|
| Helium, cubic | 17, 1.823-55.542 GPa | 78.87(6) | 148(3) | 4.09(12) | executable default |
| Neon, cubic | 5, 1.314-19.244 GPa | 78.742(14) | 163.0(1.0) | 4.02(10) | executable non-default |
| Neon, hexagonal indexing cell | 7, 24.084-53.28 GPa | 58.7(4) | 176.8(1.1) | 4.00(10) | excluded |

All parameter errors are reported 1-sigma fit uncertainties. The authors
calculated pressure-dependent parameter covariance matrices and plotted 68%
and 95% joint-confidence ellipses, but neither the numerical matrices nor a fit
statistic are printed or deposited. Stored covariance is therefore `null`.

For both executable records, `V0` was constrained by the measured ambient
volume prior `78.94 +/- 0.1 A^3`; it was not fixed. The helium `K0` and
`K0_prime` were free. The neon cubic fit additionally used a `4.0 +/- 0.1`
prior on `K0_prime`; it too was not fixed. `fixed_parameters` is empty for both
records. The same `K0_prime` prior was used for the excluded hexagonal fit,
whose `V0` refers to a three-formula-unit hexagonal indexing cell. At the
undistorted cubic limit, `V_hex=3 V_cubic/4`.

## Primary observations and pressure calibration

The official Table 2 provides all numerical observations needed to examine the
published fits:

- one ambient cubic measurement;
- 17 helium pressure steps;
- 12 neon pressure steps;
- cubic and hexagonal lattice parameters and volumes at every pressure step;
- all printed lattice, volume, c/a, and pressure uncertainties; and
- the two ruby R1 wavelengths used at each pressure.

The CSV keeps all 30 rows and explicitly flags the ambient prior, the 17-row
helium fit, the five-row neon cubic fit, and the seven-row excluded neon
hexagonal fit. Alternate cell indexings are retained rather than cleaned away.
At 31.820 GPa in helium, ruby 1 was not measured. The authors estimate it from
ruby 2 with `P_ruby1 = 0.9973 P_ruby2 + 1.8804 GPa`; the row is flagged as an
estimate. The 19.244 GPa neon row has no printed pressure uncertainty and is
left blank rather than set to zero.

Pressures use the Dewaele, Loubeyre, and Mezouar (2004) ruby scale, DOI
`10.1103/PhysRevB.70.094112`, which is already executable in Peritheos as
`ruby_dewaele_2004`. Reported pressure is the mean from the rubies in a chamber;
the pressure uncertainty is their standard deviation. Absolute wavelengths are
printed, but the specimen-specific zero-pressure R1 reference wavelength is
not. Applying the calibration's generic `694.24 nm` reference does not exactly
recover the listed pressures (for example, `694.95 nm` gives 1.956 GPa instead
of 1.823 GPa). Exact raw-wavelength reduction is therefore not claimed.

## Numerical checks

Direct evaluation of the published curves on the rounded fit rows gives:

| Record | Pressure RMSE (GPa) | Maximum absolute residual (GPa) |
|---|---:|---:|
| Helium cubic | 0.2779224099 | 0.5288116855 |
| Neon cubic | 0.1182512894 | 0.1542757551 |

A diagnostic effective-variance fit used
`sigma_eff^2 = sigma_P^2 + (dP/dV sigma_V)^2` and the published Gaussian
priors. It gives `(V0,K0,K0')=(78.88229097,146.69244887,4.13085052)` for
helium, within the reported 1-sigma parameter errors. For neon it gives
`(78.76818584,161.64662941,3.99361735)`, within two reported 1-sigma errors.
The exact MINUTI 2.0 objective, weights, and P-V covariance treatment are not
specified, so these diagnostics are not published as refit records.

The article independently extrapolates the helium and neon cubic curves to
135.8 GPa and reports densities of 8.17 and 8.03 g/cm3. The curves reproduce
8.16674 and 8.03710 g/cm3 when the nominal `(Mg0.22Fe0.78)O`, `Z=4` mass is
used. The defect-normalized composition instead gives 8.02497 g/cm3 for the
helium curve. This source-level mass-basis inconsistency is recorded and is not
used to alter the measured material composition.

## Inclusion decision

Two published BM3 records are executable:

- `mg0215fe0762vac0023o_finkelstein_2017_bm3_helium_1`;
- `mg0215fe0762vac0023o_finkelstein_2017_bm3_neon_cubic_2`.

The helium record is the default equilibrium-oriented baseline. The neon cubic
record is retained as a non-default stress-sensitive experimental alternative.
The high-pressure hexagonal-indexing-cell parameterization is excluded because
the primary source does not establish a space group or atomic structure and
ties the distortion to nonhydrostatic stress. Its parameters and all seven
observations remain in validation metadata and the dataset.
