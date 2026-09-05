# Lee and Wan (2008): matched LDA and GGA MgO EOS

## Outcome

Both LitCurate candidates are accepted as distinct, source-computed 0 K BM3
fits. They use matched pseudopotential setups but different exchange-correlation
functionals, and the difference between the two EOS is the central object of
the paper. Neither curve is made a default pressure standard.

## Primary source and method

Shun Hang Lee and Jones Tsz Kai Wan, “Pressure correction in
density-functional calculations,” *Physical Review B* **78**(22), 224103
(2008),
doi:[10.1103/PhysRevB.78.224103](https://doi.org/10.1103/PhysRevB.78.224103).
The audit used the
[official APS full text](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.78.224103/fulltext).

Section III describes a 64-atom, 32-MgO-unit rocksalt supercell. The Mg and O
pseudopotentials for the two branches have the same cutoff radii and reference
states; only the exchange-correlation choice differs. The plane-wave and charge
density cutoffs are 30 and 240 Ry. Section IV.A states that pressures are fitted
to a third-order Birch-Murnaghan EOS and gives all three coefficients directly.

| Candidate | Branch | Source V0 (bohr3, 32 MgO) | K0 (GPa) | K0' | Disposition |
|---|---|---:|---:|---:|---|
| `litcurate_b4138675da02eb45` | LDA, 0 K | 3782.30 | 177.486 | 4.026 | ACCEPT: independent source-computed BM3. |
| `litcurate_8a7a5b23526434e1` | GGA, 0 K | 4046.43 | 149.320 | 4.080 | ACCEPT: independent source-computed BM3. |

The source also displays Karki et al. and Isaak et al. curves in Figures 1-2,
but labels them as external comparisons and draws them from their published EOS
parameters. They are not additional Lee-Wan records.

## Volume conversion and numerical reproduction

The article's V0 values are in bohr3 for the complete 32-formula-unit
simulation cell. Peritheos stores MgO volumes on the conventional rocksalt Z=4
basis. Using the exact CODATA/SI Bohr-radius conversion represented by
`1 bohr = 0.529177210903 A`, the transformation is

`V0(A3, Z=4) = V0(bohr3, 32 MgO) * 0.529177210903^3 * 4/32`.

This gives 70.059879275145 A3 for LDA and 74.952382755288 A3 for GGA. The
scaling changes only the volume basis; volume ratios and therefore BM3
pressures are invariant.

The source does not tabulate the pressure-volume simulation grid or the fit
weights. Figure 1 is explicitly drawn from the published parameters, so
digitizing that rendered curve and refitting it would add graphical rounding
without independent information. Instead, the reproduction script evaluates
both source BM3 expressions at the five volume locations visible in Figures 1,
3, and 4 and compares `P_GGA-P_LDA` with the independent Figure-4 polynomial

`delta P = 1.53e4/V + 1.06e8/V^2`,

where V is in bohr3. Across 2740-3800 bohr3, the maximum disagreement is 0.141
GPa, well below the paper's stated roughly 1 GPa temperature-to-temperature
scatter and visually consistent with Figure 4. Each stored EOS also returns
zero pressure at its converted V0.

No coefficient uncertainties, covariance, or raw calculation table are
published. These fields therefore remain null rather than being inferred from
plot resolution.
