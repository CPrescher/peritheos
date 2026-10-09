# Yokoo et al. (2009): Pt pressure-scale branches and source audit

The supplied primary publication, [Yokoo et al., Physical Review B 80,
104114](https://doi.org/10.1103/PhysRevB.80.104114), was inspected on
2026-10-09. All nine pages were rendered; Equations (1)-(18), Tables I-V,
footnotes and Section IV were checked visually. Its checksum, complete Pt
Table IV transcription, source equations and missing dependencies are in
`peritheos/data/datasets/platinum-yokoo-2009-source.json`.

## Accepted 300 K Vinet branch

`platinum_yokoo_2009_vinet_300k` preserves the **published fitted 300 K
Vinet branch**, Table IV and footnote b: K0=276.4 GPa and K0'=5.48.
The pressure is

`P = 3 K0 (1-x) x^-2 exp[1.5 (K0'-1) (1-x)]`, `x=(V/V0)^(1/3)`.

Table IV prints ambient density 21.40 Mg/m³, not a cell volume. Converting
with four atoms per fcc cell, Pt molar mass 195.084 g/mol and exact N_A gives
V0=60.5503931552 Å³/cell, rounded to **60.55 Å³/cell** at the density's
precision. Molar volume is 9.11607 cm³/mol of atoms. The conversion constant
and rounding are explicit conventions; no uncertainty or hidden unrounded
author V0 is recovered. The commonly used 60.38 Å³ belongs to other Pt
scales and must not be substituted here. The source's 101.9 a.u.³ footnote
refers to **Holmes's** density column, not the new ambient density.

The 300 K Vinet coefficients are a fit to the full EOS's isotherm. They are
neither its 0 K BM3 coefficients nor an exact evaluation of the complete
thermal EOS. The source does not give uncertainty widths or covariance for
this Vinet triplet. The ±0.10 error on the **0 K BM3** pressure derivative
must not be attached to K0'=5.48. The record is isothermal, without a partial
phonon-only thermal extension. Compression V/V0=0.6-1 corresponds to
0-522.631 GPa in this Vinet fit; this is modeled coverage, not static
experimental coverage or a new phase-stability guarantee.

## Full EOS and branch distinctions

Equation (1) is the absolute decomposition

`P(V,T)=Pc(V,0)+Pph(V,T)+Pel(V,T)`.

The cold branch is BM3 with **Vc at 0 K**, B0=288.4 GPa and B0'=5.05±0.10.
Vc differs from ambient V0 at 300 K and is not tabulated. The separately
fitted 300 K BM3 uses B0=276.4 GPa and B0'=5.12. The Vinet fit uses 5.48.
Combining 288.4, 5.05 and ambient V0 would change the source reference state.

Equations (4)-(8), (10) and (12) give the phonon terms:

`gamma(V)=2.63 [1+0.39((V/V0)^5.2-1)]`,

`theta(V)=230 (V/V0)^[-(1-0.39)2.63] exp[-(gamma(V)-2.63)/5.2] K`,

`Eph=3 n kB T D3(theta/T)`, `Pph=gamma Eph/V`.

With molar energies, replace `n kB` by R per mole of Pt atoms and use
`Vm=Vcell N_A/4`. The source energy expression has no explicit zero-point
term. Phonon pressure is absolute; a 300 K-subtracted term added to the 0 K
cold branch would not implement Equation (1). a=0.39±0.08 and b=5.2±1.1 are
printed errors with no confidence convention or covariance given. The
explicit 2σ statement on the **Au shock-velocity regression** does not
establish 2σ for all Pt EOS parameters. gamma0 and theta0 errors are absent.

The electronic pressure and energy are taken from the first-principles
calculations of [Tsuchiya and Kawamura (2002), PRB 66,
094115](https://doi.org/10.1103/PhysRevB.66.094115). Yokoo gives no explicit
coefficient set or numerical electronic grid. That paper also has an
[erratum, PRB 67, 019902](https://doi.org/10.1103/PhysRevB.67.019902), whose
existence is confirmed by the [APS volume-67 author
index](https://journals.aps.org/prb/issues/67/24/deliverables/miscellaneous/print).
The subsequently supplied article and erratum have now been inspected.
All 51 Pt pressure nodes are bundled; the erratum changes Au only. Pt's
0 K-referenced 300 K electronic pressure is 0.04 GPa and is not subtracted
in Yokoo's absolute decomposition. Linear interpolation and volume
independence are explicit reconstruction assumptions. Numeric electronic
energy remains unavailable for original shock-temperature replay, but is
not a PVT blocker. The printed phonon coefficients plus these electronic
nodes leave a maximum 0.3806656 GPa pressure mismatch using the tabulated
cold outputs. Its cause is not established.

A separate **derived** full PVT record,
`platinum_yokoo_2009_pvt_reconstruction`, is now registered. It retains an
explicit empirical residual-pressure correction alongside the unchanged
electronic table. The [PVT library audit](yokoo-2009-pvt-library.md) documents
the pressure equation, derivation, 166-state reconstruction and 79-state
holdout. This does not promote diagnostic coefficients to published values.

## Recovered outputs and experimental constraints

All **168 Table V Pt states** (21 volume ratios, eight temperatures) are
bundled in `platinum-yokoo-2009-table5-isochores.csv`. They are derived
full-model output, spanning V/V0=0.6-1.0 and 0-3000 K, up to 550.69 GPa.
They are not measurements or independent fitting targets. Two 3000 K entries
at V/V0=1.00 and 0.98 are parenthesized. The preceding Table III identifies
parentheses as the first liquid state; Table V does not repeat that caption
convention, so this phase interpretation is recorded explicitly. The grid
does not establish a solid-state rectangular validity domain through melting.

Section IV reuses Holmes et al. (1989) Pt shock data to 660 GPa. The existing
seven new Holmes Table III experiments are retained with source density,
shock/particle velocities, impedance-match provenance and 2σ errors. They
are recovered upstream constraints, not a proven complete Yokoo fit selection.
Shock pressures are conservation-law reductions, rather than Pt EOS output;
the upstream Ta impactor/calibration assumptions remain relevant. Ambient
thermodynamic inputs are cited to Touloukian, Hultgren, Macfarlane and
Collard/McLellan (references 19, 22, 38, 39); Yokoo supplies no numerical
Pt thermodynamic input table.

The already bundled 36 Dewaele (2004) Au/Pt volume pairs underpin the
Figure 6 comparison. Au has one additional table row without a Pt partner;
the audit joins the 36 pairs by both original ruby-pressure coordinates and
retains both CSV row indices. Their ruby-pressure columns are source-calibration
metadata, not Yokoo fit targets. This compares two pressure models using
simultaneous measured volumes; it does not independently measure pressure.
The Au/Pt comparison remains separate from the new derived PVT records;
its measured ruby-pressure columns are not replaced by reconstructed values.

Section II.D describes simultaneous steepest-descent least squares with
respect to total pressures, with B0'(0 K), a and b adjustable and gamma0
obtained secondarily. It delegates procedural detail to Tange (2009), and
omits exact Pt selections, numerical ambient inputs, row weights and
covariance. Thus the published parameterization can be implemented while
its original joint fit remains **not_refittable**.

## Independent numerical audit

Run `.venv/bin/python -m scripts.reproduce_yokoo_2009_platinum --check`.
The script directly evaluates Vinet, BM3 and the Debye integral using SI
constants, independently of the library EOS evaluator.

- The canonical Vinet record agrees below 10^-9 GPa over six states spanning
  V/V0=1.0-0.6. Reconstruction from `.eosmat` and pressure-volume inversion
  are tested.
- Comparing Vinet to the full model's 21 Table V 300 K values finds a maximum
  difference of **6.009257 GPa**, at V/V0=0.6. Their distinction is material
  at high compression; this is not a failure of Vinet equation parity.
- A conditional one-parameter fit of Vc/V0 to the **21 derived 0 K grid
  values**, fixing the printed BM3 coefficients, gives **0.99382001**, with
  maximum residual **0.005269 GPa**. This is an effective cold-volume
  checkpoint diagnostic, not an independent refit or recovery of the
  unrounded author reference volume.
- The integrated theta law satisfies `d ln theta/d ln V=-gamma` numerically.
  At V/V0=0.6 and 3000 K, the table's total pressure minus its 0 K value minus
  the printed-coefficient phonon term is **1.993696 GPa**. This residual
  exposes the missing term/normalization information; it is not an extracted
  electronic datum or a validated electronic fit.
- Holmes shock momentum pressures reproduce the reported reductions within
  last-digit rounding intervals (maximum raw pressure difference 0.14396 GPa),
  separately from their experimental 2σ errors. The Table IV linear Us=3.635+1.543 up relation is compared
  to those seven shots as a constraint check, without declaring complete
  source-regression parity.

The saved JSON retains input checksums, every thermal residual, branch
comparison and upstream shock check. Original errors remain separate from
conditional diagnostics.

## Sakai (2018) Re calibration

The existing Re record explicitly names the **Yokoo Pt Vinet scale**.
That equation reference resolves to `platinum_yokoo_2009_vinet_300k`, rather
than a BM3 or incomplete thermal curve. Original paired Pt/Re volumes and
pressure uncertainties remain untabulated, and Sakai does not document an
unrounded Pt normalization separately. Resolving the published branch is
therefore distinct from reproducing the historical pressure reduction:
`recalculation.status=missing_calibrant_observations`. The recovered Re
Figure 10 fit and its qualified numerical parity are preserved unchanged.
