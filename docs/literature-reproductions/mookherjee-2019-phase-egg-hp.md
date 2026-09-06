# Mookherjee et al. (2019): Phase Egg HP completion

Primary source: M. Mookherjee, W. R. Panero, B. Wunder, and S. Jahn, “Anomalous elastic behavior of phase egg, AlSiO3(OH), at high pressures,” *American Mineralogist* **104**, 130–139 (2019), <https://doi.org/10.2138/am-2019-6694>.

The source reports separate static 0 K third-order Birch–Murnaghan fits for the low-pressure and proton-transferred high-pressure configurations. The LP record was already bundled. This audit adds only the scientifically distinct HP branch: `V0=207.74(1) Å3`, `K0=222.8(2) GPa`, and `K0'=4.44(2)`. The official MSA spreadsheet's five HP P–V rows reproduce the rounded published curve with 0.0397 GPa pressure RMSE. The publication's 2018 DOI is an accepted-manuscript discovery alias; the final 2019 DOI is canonical.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 1012 | `litcurate_dfe076f60c19b6f9` | REJECT | Duplicate accepted-manuscript identifier for the already-bundled LP record. |
| 1013 | `litcurate_2d449f80a6565739` | REJECT | Duplicate accepted-manuscript identifier for the HP record implemented under the final DOI. |
| 1014 | `litcurate_2d02ffd4982dbbbc` | REJECT | Citation-reported Vanpeteghem experimental EOS, not source-owned here. |
| 1015 | `litcurate_8e22bda29924b666` | REJECT | Canonical-DOI duplicate of the already-bundled LP record. |
| 1016 | `litcurate_2ecdcd601464227e` | ACCEPT | Distinct source-owned static HP proton-transferred BM3. |
| 1017 | `litcurate_49fc1ef85c2217b1` | REJECT | Citation-reported Vanpeteghem experimental EOS, not source-owned here. |

Result: **1 net-new production record and 5 rejected duplicate/citation rows**.
