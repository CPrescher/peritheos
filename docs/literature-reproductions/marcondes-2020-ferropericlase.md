# Marcondes et al. (2020): dilute ferropericlase spin configurations

Primary source: M. L. Marcondes, F. Zheng, and R. M. Wentzcovitch, “Phonon
dispersion throughout the iron spin crossover in ferropericlase,” *Physical
Review B* **102**, 104112 (2020),
<https://doi.org/10.1103/PhysRevB.102.104112>.

The accepted-manuscript/arXiv copy (`arXiv:2003.12348`) was checked against the
computational methods and Table I. Static LSDA+Usc energy-volume curves were
fitted to a third-order finite-strain EOS. The 64-atom supercells represent
`XFe=0.03125` (one Fe) and `XFe=0.0625` (two Fe). At 6.25% Fe, 11th-neighbor and
2nd-neighbor arrangements are distinct atomistic configurations, and HS, LS,
and one-HS/one-LS mixed-spin states are distinct calculated branches. They are
therefore not arbitrary duplicate splits.

No numerical energy-volume grid was deposited. All eight complete Table I
parameterizations are retained, and an independent BM3 implementation checks
their curves at fixed compression ratios. The records are static theoretical
parameterizations, not room-temperature experimental EOSs.

## LitCurate disposition

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 1070 | `litcurate_8b61325b3373dcc9` | ACCEPT | 3Fp high-spin branch. |
| 1071 | `litcurate_b5a171ea0852a37e` | ACCEPT | 3Fp low-spin branch. |
| 1072 | `litcurate_77c20002800f4789` | ACCEPT | 6Fp 11nn high-spin branch. |
| 1073 | `litcurate_ec75255e7d6eae38` | ACCEPT | 6Fp 11nn low-spin branch. |
| 1074 | `litcurate_2116ff41406a0b93` | ACCEPT | 6Fp 11nn mixed-spin branch. |
| 1075 | `litcurate_cd07018138aa391c` | ACCEPT | 6Fp 2nn high-spin branch. |
| 1076 | `litcurate_c871ad45b4ff997f` | ACCEPT | 6Fp 2nn low-spin branch. |
| 1077 | `litcurate_442cdb585530380d` | ACCEPT | 6Fp 2nn mixed-spin branch. |

Result: **8 accepted production records**.
