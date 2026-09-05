# Metsue and Tsuchiya (2012): Fe2+-bearing bridgmanite

## Outcome

Seven static BM3 records are accepted: six source-computed fits for three
high-spin and three low-spin configurations of
`(Mg0.9375Fe0.0625)SiO3`, plus the study's pure-MgSiO3 control. Four
comparison rows are citation traces and are rejected under this DOI.

The three configurations within each spin state are retained as explicitly
reported, non-identical computational sensitivity records. They are not
represented as different materials or independent experiments. The source
selects high-spin model 1 for subsequent thermodynamic calculations, so it is
the only default record. All low-spin records carry the source's instability
warning.

## Primary source and method

Arnaud Metsue and Taku Tsuchiya, “Thermodynamic properties of
(Mg,Fe2+)SiO3 perovskite at the lower-mantle pressures and temperatures: an
internally consistent LSDA+U study,” *Geophysical Journal International*
**190**, 310-322 (2012),
doi:[10.1111/j.1365-246X.2012.05511.x](https://doi.org/10.1111/j.1365-246X.2012.05511.x).
The audit used the [official full-text article](https://academic.oup.com/gji/article/190/1/310/598134).

The source uses PWSCF, internally consistent LSDA+U, a 70 Ry cutoff, and three
80-atom supercell shapes sampled with 2x2x2, 2x3x1, and 3x2x1 k-point meshes.
Cells were relaxed at 0, 30, 60, 90, 120, and 150 GPa. Section 3.2 states that
third-order Birch-Murnaghan equations were fitted and that `K0'=3.94` was
fixed for all seven Table 1 parameterizations.

## Volume conversion and checks

Table 1 reports molar volume in cm3/mol. Peritheos uses a conventional-cell
volume, so values are converted to the standard Pbnm `Z=4` basis with

`V(A3, Z=4) = V(cm3/mol) * 4 * 10^24 / N_A`,

using the exact SI Avogadro constant. No coefficient uncertainty or covariance
is published.

The source does not tabulate its six P-V points per model, so a numerical
refit is impossible without digitizing highly overlapping curves. Instead,
the deterministic reproduction verifies the exact unit conversion, zero-
pressure BM3 anchor, source statement that same-spin configuration parameters
differ by less than 0.07%, and the reported composition slopes. Maximum
configuration spreads are below 0.01% in `V0` and below 0.07% in `K0`. Model-1
bulk-modulus slopes are 0.0343 (HS) and 0.1244 (LS), consistent with the
source's rounded 0.035 and 0.123.

The paper finds low-spin Fe dynamically unstable at 0 GPa and energetically
unfavourable throughout lower-mantle pressures. These curves remain useful as
documented metastable comparisons but are never defaults.

## Candidate dispositions

| LitCurate identifier | Disposition | Reason |
|---|---|---|
| `litcurate_598c417d4fed2c69` | ACCEPT | HS model 1; source-selected representative static BM3. |
| `litcurate_0bfa6ed709638405` | ACCEPT | HS model 2 configuration sensitivity BM3. |
| `litcurate_c918eff9cd068980` | ACCEPT | HS model 3 configuration sensitivity BM3. |
| `litcurate_db9128e2f40c19d5` | ACCEPT | LS model 1 metastable comparison BM3. |
| `litcurate_32c57baac8c9a39d` | ACCEPT | LS model 2 metastable configuration BM3. |
| `litcurate_14060ec60b6b2a23` | ACCEPT | LS model 3 metastable configuration BM3. |
| `litcurate_2602c226a511a4a9` | ACCEPT | Source-generated pure-MgSiO3 static control BM3. |
| `litcurate_c0f543d032abfa5d` | REJECT / CITATION TRACE | Experimental 253 GPa comparison lacks V0 and belongs to cited experiments. |
| `litcurate_37bb4ec9bd673079` | REJECT / CITATION TRACE | Karki et al. (2001) LDA result; audit under its primary DOI. |
| `litcurate_83016da834a34cd7` | REJECT / CITATION TRACE | Hsu et al. ferric-iron comparison, not generated here. |
| `litcurate_5d9aa6b0f56950f5` | REJECT / CITATION TRACE | Li et al. Fe-Al comparison, not generated here. |

The finite-temperature QHA curves in Figure 3 are not additional Table 1 BM3
records. They require the source's phonon free-energy model and are not
silently approximated as static isotherms.
