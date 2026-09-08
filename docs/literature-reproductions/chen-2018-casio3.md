# Chen et al. (2018): tetragonal CaSiO3 perovskite

Primary source: H. Chen, S.-H. Shim, K. Leinenweber, V. Prakapenka, Y. Meng, and C. Prescher, “Crystal structure of CaSiO3 perovskite at 28–62 GPa and 300 K under quasi-hydrostatic stress conditions,” *American Mineralogist* **103**, 462–468 (2018), <https://doi.org/10.2138/am-2018-6087>.

## Fundamental correction

The previous Peritheos record called this fit BM2. That is not supported by the primary source. The EOS paragraph says the `I4/mcm` pressure-volume results were fitted to the **Vinet equation**, Figure 5 calls its solid line the Vinet fit, and Table 2 defines `V` as Vinet while separately defining `BM` as third-order Birch-Murnaghan. The corrected production record is therefore `ca_perovskite_tetragonal_chen_2018_vinet`.

The source reports `V0 = 46.3(1) A3` per formula unit, `K0 = 223(6) GPa`, and `K0' = 4.0` fixed. Table 2 defines the parenthetical parameter errors as 2-sigma. Peritheos stores the volume on the conventional `I4/mcm`, `Z=4` basis: `V0 = 185.2(4) A3`. Because the phase is unstable at one bar, `V0` is a fitted extrapolation rather than an ambient measurement.

## Authoritative table and redistribution

Table 1 was visually checked in two lossless manifestations: the [version-of-record PDF](https://pubs.geoscienceworld.org/msa/ammin/article-pdf/103/3/462/4085107/am-2018-6087.pdf) and the [author manuscript deposited by S.-H. Shim](https://www.researchgate.net/profile/Sang-Heon-Shim/publication/323517011_Crystal_structure_of_CaSiO3_perovskite_at_28-62_GPa_and_300_K_under_quasi-hydrostatic_stress_conditions/links/5e02cd4b4585159aa49853cf/Crystal-structure-of-CaSiO3-perovskite-at-28-62-GPa-and-300-K-under-quasi-hydrostatic-stress-conditions.pdf). Both independently confirm every table cell and the model attribution.

The PDFs are not bundled. The [MSA publishing policy](https://msaweb.org/publishingpolicy/) reserves copyright in the version of record, while the author-hosted manuscript is marked as author content subject to copyright. The repository instead bundles a checksum-tracked CSV containing only the table's numerical facts, row identifiers, and transparently derived values. This preserves all seven source rows and both structural refinements at printed precision without redistributing article text or layout.

## Fit reconstruction

The source protocol can be recovered as follows:

- Seven runs, in Table 1 order: `71013`, `81021`, `81030`, `71050`, `81066`, `81074`, and `71088`.
- No exclusions are reported. The EOS text explicitly selects the `I4/mcm` refinements; the parallel `P4/mmm` refinements are retained only for audit diagnostics.
- The conventional-cell fit volumes are derived without intermediate rounding as `a^2 c`; division by four gives the source's pseudo-cubic/formula-unit convention.
- Pressure is calculated from the Pt scale of Ye et al. (2017), DOI `10.1002/2016JB013811`. Pt lattice parameters or volumes are not printed, so the row pressures cannot be recalculated.
- `K0' = 4.0` is fixed; `V0` and `K0` are fitted.
- Table 1 lattice and angle parentheses are 2-sigma. Pressure uncertainties, residual direction, numerical weights, fit software, unrounded lattice values, and parameter covariance are not published.

The complete rounded table therefore supports a production Vinet parameterization and a direct source-scope fit attempt, but not exact reconstruction of the authors' optimizer or covariance.

## Numerical result

All values below use only the seven published `I4/mcm` rows and hold `K0' = 4`:

| Objective | `V0` (A3/f.u.) | `K0` (GPa) | pressure RMSE (GPa) |
|---|---:|---:|---:|
| Published curve | 46.3000 | 223.000 | 0.6885 |
| Unweighted pressure residual | 46.4945 | 215.291 | 0.6345 |
| Unweighted volume residual | 46.4685 | 216.362 | 0.6356 |
| Volume residual / propagated 1-sigma volume error | 46.4062 | 218.269 | 0.6523 |

The largest published-curve pressure residual is 1.627 GPa. The spread among reasonable objectives is larger than the tabulated lattice precision, showing that the omitted pressure errors and weighting convention matter. An unsupported diagnostic that treats both structural refinements as duplicate observations gives `V0 = 46.3160 A3/f.u.` and `K0 = 221.833 GPa`, much closer to the printed coefficients, but it contradicts the source's explicit `I4/mcm` fit statement and is not used for production.

The paper also reports a rejected sensitivity fit with `V0 = 45.58 A3/f.u.` fixed and both `K0` and `K0'` free (`290 GPa`, `2.3`). Applying unweighted pressure residuals to the rounded `I4/mcm` table gives `298.6 GPa` and `1.82`, further confirming that the exact source objective cannot be recovered from the printed rows alone.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 1008 | `litcurate_cb91837589707f8c` | ACCEPT, CORRECTED | The primary paper supplies a production Vinet record and complete seven-row table; the discovery ledger's BM2 identity was wrong. |

Result: **1 corrected production Vinet record, 1 complete seven-row dataset, and a qualified source-scope refit**.
