# Vanpeteghem et al. (2002): delta-AlOOH

Primary source: C. B. Vanpeteghem, E. Ohtani, and T. Kondo, “Equation of state
of the hydrous phase δ-AlOOH at room temperature up to 22.5 GPa,”
*Geophysical Research Letters* **29**, 23-1–23-3 (2002),
<https://doi.org/10.1029/2001GL014224>.

All nine Table 1 decompression rows are transcribed. The paper identifies
orthorhombic `P21nm` δ-AlOOH and publishes two third-order Birch–Murnaghan
sensitivity fits on the same `Z=2`, `V0=56.54(9) Å3` cell: `K0=252(3) GPa`
with `K0'=4` fixed, and `K0=228(7) GPa`, `K0'=7(1)` with both elastic
coefficients fitted. Table 2 explicitly lists both “present work” rows, so they
are genuine alternate parameterizations rather than a phase split.

Direct evaluation against the rounded Table 1 volumes gives pressure RMSEs of
1.43 and 1.36 GPa. This is a transparent curve check, not a claim to recover
the authors’ unreported regression weights/covariance from rounded printed
values.

## LitCurate disposition

| Row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 382 | `litcurate_b8a1879cafff7c10` | ACCEPT | Complete fixed-derivative BM3 parameterization in the abstract and Table 2. |
| 383 | `litcurate_1fcd2fb80fadaf43` | ACCEPT | Complete unconstrained BM3 sensitivity fit in Table 2. |
| 384 | `litcurate_0e1318b28316f7df` | REJECT | Stishovite comparison is citation-reported. |
| 385 | `litcurate_a0fe3d01aac6951a` | REJECT | Phase-B comparison is citation-reported. |
| 386 | `litcurate_d106fd5343db2b68` | REJECT | Diaspore comparison is citation-reported. |

Result: **2 production records, 3 rejected citation traces**.
