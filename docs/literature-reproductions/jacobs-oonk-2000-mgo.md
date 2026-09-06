# Jacobs and Oonk (2000): GGK equation for MgO

Primary source: M. H. G. Jacobs and H. A. J. Oonk, “A new equation of state based on Grover, Getting and Kennedy's empirical relation between volume and bulk modulus. The high-pressure thermodynamics of MgO,” *Physical Chemistry Chemical Physics* **2**, 2641–2651 (2000), <https://doi.org/10.1039/A910247G>.

The primary article's equation development was checked directly. It assumes a linear relation between molar volume and `ln(K)` at pressure, integrates `K = -V(dP/dV)`, and evaluates the resulting power-series GGK equation for MgO over 100–3100 K and 0–225 GPa. This is a distinct analytical EOS family, not Birch–Murnaghan. The tabulated ambient MgO values (`V0 = 11.25 cm3/mol`, `K0 = 161.5 GPa`, and `Kp = 4.769`) do not by themselves make the GGK pressure law executable in Peritheos. Mapping them to BM3 would change the source model.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 263 | `litcurate_312cfb82b6e5ea3b` | HOLD | Source-owned and complete as an ambient derivative triple, but the source uses the unsupported GGK pressure law. Requires a dedicated GGK implementation and its complete coefficients. |
| 264 | `litcurate_a729ebfc380f380b` | REJECT | Table 2 literature comparison attributed to Mao and Bell; derivative only. |
| 265 | `litcurate_a0681128211bd534` | REJECT | Table 2 literature comparison attributed to Richet et al.; derivative only. |
| 266 | `litcurate_1629693d225db4dc` | REJECT | Table 2 literature comparison attributed to Inbar and Cohen; derivative only. |
| 267 | `litcurate_6cf386736b0c0c50` | REJECT | Table 2 literature comparison attributed to Carter et al.; derivative only. |
| 268 | `litcurate_2fe1f2dbd66413f2` | REJECT | Table 2 literature comparison attributed to Chang and Barsch; derivative only. |
| 269 | `litcurate_25e6a9d1b8bbf5f3` | REJECT | Table 2 literature comparison attributed to Jackson and Niesler; incomplete and not source-owned. |

Result: **0 production records, 1 held source row, and 6 rejected citation rows**. The source is also a concrete future model-family candidate for Peritheos.
