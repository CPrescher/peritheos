# Liu et al. (2010): MgSiO3 post-perovskite

Primary source: Z.-J. Liu, X.-W. Sun, C.-R. Zhang, L.-N. Tian, and Y. Guo, “Thermodynamic properties of MgSiO3 post-perovskite,” *Modern Physics Letters B* **24**, 315–324 (2010), <https://doi.org/10.1142/S0217984910022391>.

The author-provided primary article was checked directly. Sections 2–3 and Table 1 report a plane-wave pseudopotential LDA calculation whose static `E-V` states were fitted to third-order Birch–Murnaghan: `V0 = 163.3 A3` per `Z=4` Cmcm unit cell, `K0 = 219.3 GPa`, and `K0' = 4.4`. The article's Debye thermodynamic calculations span 0–150 GPa and 0–2000 K, but those outputs are not additional isothermal BM3 fits.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 731 | `litcurate_cf1303ecf47a5001` | ACCEPT | Complete source-owned static LDA BM3 for the Cmcm post-perovskite phase. |
| 732 | `litcurate_fbd06dbe37c7c804` | REJECT | Table 1 literature comparison, Oganov and Ono LDA. |
| 733 | `litcurate_599b59bb893eb3e1` | REJECT | Table 1 literature comparison, Oganov and Ono GGA. |
| 734 | `litcurate_6898acc217bb615d` | REJECT | Table 1 literature comparison, Tsuchiya et al. |
| 735 | `litcurate_581cdfb9a8544c1e` | REJECT | Table 1 experimental comparison, Guignot et al. |
| 736 | `litcurate_0f1b27e4ae5f1bb0` | REJECT | Table 1 experimental comparison, Ono et al. |
| 737 | `litcurate_acc36ce5b9b6c41a` | REJECT | Table 1 experimental comparison, Shieh et al. |

The reproduction script evaluates the exact BM3 curve at five compression anchors; the focused test verifies the bundled coefficients, volume basis, and inverse pressure–volume round trip.

Result: **1 production record and 6 rejected citation rows**.
