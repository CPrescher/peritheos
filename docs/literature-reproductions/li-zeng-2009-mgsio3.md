# Li and Zeng (2009) MgSiO3 EOS audit

## Scope and primary authority

- Primary paper: Yanling Li and Zhi Zeng, “First-Principles Study of the Structural, Electronic and Optical Properties of MgSiO3 at High Pressure,” *International Journal of Modern Physics C* **20** (2009), 1093–1101, DOI [10.1142/S0129183109014242](https://doi.org/10.1142/S0129183109014242).
- Primary full-text copy inspected: [author-uploaded article](https://www.researchgate.net/publication/228511351_First-Principles_Study_of_the_Structural_Electronic_and_Optical_Properties_of_MgSiO3_at_high_pressure).
- Authority locations: computational method in Section 2, energy-volume curves in Figure 1, and the complete EOS comparison in Table 1.
- Scope: all ten LitCurate rows carrying this DOI were audited. Six are this paper's source-generated fits. Four are comparison values explicitly attributed by Table 1 to references 4 or 7 and are not records from this paper.

The authors calculated static PBE-GGA energies with the LAPW+lo method for Pnma perovskite and Cmcm post-perovskite. They fitted each phase's same energy-volume calculation independently with third-order Birch-Murnaghan, Vinet, and third-order natural-strain equations. Those are genuine source-reported model-sensitivity branches: they are not duplicated runs, but they are also not independent physical observations. The BM3 branch is the default and the Vinet/natural-strain branches are retained as non-default sensitivities.

Volumes in Table 1 are conventional four-formula-unit MgSiO3 cells. All six records therefore use `formula_units: 4` and the MgSiO3 molar mass 100.387 g mol-1. The paper publishes the coefficients and energy curves, but not the numerical energy grid, fitting weights, covariance, or coefficient uncertainties.

## Candidate-by-candidate disposition

| LitCurate row | Phase and Table 1 entry | Disposition | Reason |
|---|---|---|---|
| `litcurate_7aff77eec160b0a9` (694) | Pnma, this study, BM3: 168.12, 234.8, 4.199 | accepted | Complete author-generated BM3 fit. |
| `litcurate_bfcc5153b94a67c9` (695) | Pnma, ref. 4, BM3: 162.3, 259.5, 3.69 | rejected citation row | Table 1 explicitly attributes it to Fiquet et al.; accepting it under this DOI would corrupt provenance. |
| `litcurate_efaa92c0d651d646` (696) | Pnma, this study, Vinet: 168.08, 235.7, 4.294 | accepted sensitivity | Complete independent fit of the source calculation with another supported EOS family. |
| `litcurate_bfbff5de9c3f67b5` (697) | Pnma, ref. 7, Vinet: 167.42, 230.05, 4.142 | rejected citation row | Explicitly attributed to Oganov et al. |
| `litcurate_45b95058dca79981` (698) | Pnma, this study, natural strain 3: 168.04, 237.0, 4.384 | accepted sensitivity | Primary Table 1 resolves LitCurate's “other named EOS” to third-order natural strain. |
| `litcurate_be47411b3936fcff` (699) | Cmcm, this study, BM3: 168.29, 223.73, 4.152 | accepted | Complete author-generated BM3 fit. |
| `litcurate_429a5c1961d33dda` (700) | Cmcm, ref. 4, BM3: 167.64, 199.96, 4.541 | rejected citation row | Explicitly attributed to reference 4. |
| `litcurate_0d351dc013942d75` (701) | Cmcm, this study, Vinet: 168.05, 224.2, 4.406 | accepted sensitivity | Complete independent fit of the source calculation. |
| `litcurate_a6db893a40f4491a` (702) | Cmcm, ref. 7, Vinet: 168.26, 223.7, 4.152 | rejected citation row | Explicitly attributed to Oganov et al. |
| `litcurate_ea7e1400594f5b03` (703) | Cmcm, this study, natural strain 3: 167.75, 226.7, 4.684 | accepted sensitivity | Primary Table 1 resolves the equation identity and supplies all coefficients. |

The four comparison rows are not silently reassigned to their underlying papers in this batch. They remain citation-trace work because their own primary publications must be audited independently for equation convention, basis, and duplication against existing Peritheos records.

## Numerical reproduction

`scripts/reproduce_li_zeng_2009_mgsio3.py` loads all six serialized records through the public material API and evaluates their analytical pressure functions. It verifies zero pressure at each fitted `V0`, checks the source coefficient triplets, and reports pressure at `V/V0 = 0.8`. The three fits agree closely for Pnma (83.69–83.90 GPa), while the Cmcm sensitivity spread is 79.51–82.17 GPa. This is a useful executable check on equation identity and parameter wiring; it is not a refit because the source energy grid was not published.

## Production outcome

- Accepted: **6** net-new records.
- Rejected as comparison/citation rows: **4**.
- Primary numerical dataset bundled: none available; coefficients and plotted curves only.
- Zotero action: not performed in this tranche. Citation metadata are complete above and ready for import.
