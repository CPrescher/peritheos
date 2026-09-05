# Leonov et al. (2017): DFT+DMFT magnesiowustite

Primary source: I. Leonov, A. V. Ponomareva, R. Nazarov, and I. A. Abrikosov,
“Pressure-induced spin-state transition of iron in magnesiowüstite (Fe,Mg)O,”
*Physical Review B* **96**, 075136 (2017),
<https://doi.org/10.1103/PhysRevB.96.075136>.

The author-repository copy was checked against the methods, Table I, and Figures
1–4. The calculation is fully charge-self-consistent DFT+DMFT for paramagnetic
B1 `(Fe1-x,Mgx)O` at an electronic temperature of 1160 K. The eight-composition
energy-volume calculations use third-order Birch–Murnaghan fits with `K0'=4.1`.
Table I reports one ambient-pressure `V0` per composition in bohr³/formula unit,
plus separate high-spin and low-spin bulk moduli. The tabulated `V0` is the
ambient, high-spin equilibrium volume; no low-spin zero-pressure volume is
reported. The source energy-volume grid is plotted but not deposited.

The accepted records convert `V0` using the exact Bohr definition and multiply
by four for a conventional B1 cell. Their executable BM3 curves were checked
independently at fixed `V/V0` values. These checks establish equation and
coefficient transcription, not a reconstruction of the unpublished DFT+DMFT
energy regression.

## LitCurate disposition

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 964 | `litcurate_5b0bfea83d12c0fb` | ACCEPT | FeO high-spin BM3; complete Table I coefficient set after correcting the volume unit/basis. |
| 965 | `litcurate_e92cab581af0e56b` | REJECT | Low-spin `K0` is reported, but the low-spin `V0` required by BM3 is not. |
| 966 | `litcurate_47f21997eaa8b634` | ACCEPT | Mg=0.125 high-spin BM3; complete Table I coefficient set. |
| 967 | `litcurate_2f0698ae65c320b2` | REJECT | Low-spin `V0` is not reported. |
| 968 | `litcurate_5fa3c1713ba7d985` | ACCEPT | Mg=0.25 high-spin BM3; complete Table I coefficient set. |
| 969 | `litcurate_9bd8d5b53f395ddf` | REJECT | Low-spin `V0` is not reported. |
| 970 | `litcurate_6776cec20630c908` | ACCEPT | Mg=0.375 high-spin BM3; complete Table I coefficient set. |
| 971 | `litcurate_b9d8bb0bcbf290ea` | REJECT | Low-spin `V0` is not reported. |
| 972 | `litcurate_00fe1f0d9c6653b4` | ACCEPT | Mg=0.5 high-spin BM3; complete Table I coefficient set. |
| 973 | `litcurate_183858ab1b656cf9` | REJECT | Low-spin `V0` is not reported. |
| 974 | `litcurate_e2fee3875c64e323` | ACCEPT | Mg=0.625 high-spin BM3; complete Table I coefficient set. |
| 975 | `litcurate_edb317e72af3b863` | REJECT | Low-spin `V0` is not reported. |
| 976 | `litcurate_68b3a06dda40c4fe` | ACCEPT | Mg=0.75 high-spin BM3; complete Table I coefficient set. |
| 977 | `litcurate_edeba2970d95ba8d` | REJECT | Low-spin `V0` is not reported. |
| 978 | `litcurate_2c3aa7dfbb7a0b84` | ACCEPT | Mg=0.875 high-spin BM3; complete Table I coefficient set. |
| 979 | `litcurate_19730defa86587a4` | REJECT | Low-spin `V0` is not reported. |

Result: **8 accepted production records and 8 rejected incomplete rows**. The
rejections prevent the reported low-spin moduli from being turned into arbitrary
BM3 curves by silently borrowing the high-spin reference volume.
