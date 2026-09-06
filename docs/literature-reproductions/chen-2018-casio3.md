# Chen et al. (2018): tetragonal CaSiO3 perovskite

Primary source: H. Chen, S.-H. Shim, K. Leinenweber, V. Prakapenka, Y. Meng, and C. Prescher, “Crystal structure of CaSiO3 perovskite at 28–62 GPa and 300 K under quasi-hydrostatic stress conditions,” *American Mineralogist* **103**, 462–468 (2018), <https://doi.org/10.2138/am-2018-6087>.

The full author manuscript indexed with the primary article was checked directly. Thermally stress-annealed CaSiO3 in a neon pressure medium refines as tetragonal `I4/mcm`. The EOS discussion and Figure 5 report a second-order Birch–Murnaghan fit (`K0' = 4` fixed) with `V0 = 46.3(1) A3` on the one-formula-unit normalized basis and `K0 = 223(6) GPa`. The record converts only the volume basis to the conventional `Z=4` cell (`185.2(4) A3`). The source explicitly fits `V0` because the phase is unstable at one bar.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 1008 | `litcurate_cb91837589707f8c` | ACCEPT | Complete source-owned 300 K BM2 after resolving the ledger's omitted `V0` and model identity from the primary manuscript. |

The reproduction script evaluates five compression anchors through the exact BM2 equation; the focused test checks the volume normalization, coefficients, printed uncertainties, fixed derivative, and inverse pressure–volume round trip.

Result: **1 production record**.
