# Taniguchi et al. (1995): Ca-silicate calculation models

Primary source: T. Taniguchi, M. Okuno, and T. Matsumoto, “Calculation models of structures and physical properties of CaSiO3, CaMgSi2O6, and Ca2BaSi3O9 crystals,” *Mineralogical Journal* **17**, 290–300 (1995), <https://doi.org/10.2465/minerj.17.290>.

The free J-STAGE scan was checked directly. Table 6 reports this study's calculated cubic CaSiO3-perovskite cell volume (`47.34 Å3`), bulk modulus (`0.247 TPa`), and pressure derivative (`5.09`), alongside three cited comparisons. The methods describe energy minimization and pressure-dependent calculations, but the paper does not identify an analytical pressure–volume EOS family. Table 7 is a separate summary of bulk moduli and bond compressibilities. A BM3 mapping would therefore be an unsupported inference.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 118 | `litcurate_b971d10e61a3bf6a` | HOLD | Complete local derivative triple for this study's CaSiO3 model, but the analytical EOS identity and reference temperature are unreported. |
| 119 | `litcurate_f9e11fce6913dbcb` | REJECT | Cited Wolf and Jeanloz (1985) lattice-dynamics comparison. |
| 120 | `litcurate_66e25cf31d173526` | REJECT | Cited Hemley et al. (1987) lattice-dynamics comparison. |
| 121 | `litcurate_b206623029639c22` | REJECT | Cited Tamai and Yagi (1989) experimental comparison. |
| 122 | `litcurate_c3b4ba06d01aa598` | HOLD | Source Table 7 bulk modulus only; no `V0`, derivative, or EOS family. |
| 123 | `litcurate_1dae652537766448` | REJECT | Cited Matsui et al. (1987) MgSiO3 bulk modulus only. |

Result: **0 production records, 2 held source rows, and 4 rejected citation rows**.
