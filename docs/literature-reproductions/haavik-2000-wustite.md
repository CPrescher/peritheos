# Haavik et al. (2000): defect-clustered wüstite

Primary publication: C. Haavik, S. Stølen, M. Hanfland, and C. R. A. Catlow, “Effect of defect clustering on the high-pressure behaviour of wüstite. High-pressure X-ray diffraction and lattice energy simulations,” *Physical Chemistry Chemical Physics* **2**, 5333–5340 (2000), <https://doi.org/10.1039/B006026G>.

The RSC publication record and abstract verify the paper identity, compositions, experimental/simulation scope, and its conclusion that the bulk modulus is approximately composition-independent. The LitCurate rows preserve several fitted bulk moduli, but none contains the reference volume needed for an executable pressure–volume EOS. Two Fe0.99O rows also omit `K0'`. The publisher's full PDF could not be retrieved during this audit, so an ambient lattice parameter or adopted `V0` could not be verified from the primary text. No volume was inferred from another composition and no incomplete row was padded.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 270 | `litcurate_aefa9f295df29546` | HOLD | Source Fe0.925O BM3 (`K0=150 GPa`, `K0'=4`), but `V0` is missing. |
| 271 | `litcurate_902a71b5a2a8e2d4` | HOLD | Distinct Fe0.925O multi-reflection fit (`K0=167 GPa`, `K0'=4`), but `V0` is missing. |
| 272 | `litcurate_4a977ec84c750861` | HOLD | Source Fe0.99O reduction lacks both `V0` and `K0'`. |
| 273 | `litcurate_126ee62862574030` | HOLD | Source Fe0.99O restricted-data reduction lacks both `V0` and `K0'`. |
| 274 | `litcurate_9e6f3a7361abfa17` | REJECT | Citation-reported Zhang comparison and incomplete coefficients. |

Result: **0 production records, 4 held incomplete source rows, and 1 rejected citation row**. The held rows can be revisited if the primary article's ambient cell parameter and exact fit conventions become available.
