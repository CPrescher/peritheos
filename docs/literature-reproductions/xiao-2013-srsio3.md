# Xiao et al. (2013): SrSiO3 audit

Primary source: W. Xiao, D. Tan, W. Zhou, J. Liu, and J. Xu, “Cubic perovskite
polymorph of strontium metasilicate at high pressures,” *American Mineralogist*
**98**, 2096–2104 (2013), <https://doi.org/10.2138/am.2013.4470>.

The paper was fully audited as the fourth replacement candidate. It reports
three defensible EOS records: experimental cubic `Pm-3m` SrSiO3
(`V0=49.18(5) Å3`, `K0=211(3) GPa`, `K0'=4` fixed), static-GGA cubic SrSiO3
(`49.97(1)`, `207.7(4)`, `4` fixed), and static-GGA 6H `P63/mmc` SrSiO3
(`311.23(10)`, `183.8(6)`, `4` fixed). The experimental Table 1 has 17 rows,
with the 6.2 GPa amorphizing point explicitly excluded from the fit.

These three valid rows are **HOLD (batch boundary)** rather than partially or
silently implemented: the preceding natural three-record Mao family exactly
closed the settled +100 target. They remain a ready complete family for a
future batch.

## LitCurate disposition

| Row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 826 | `litcurate_ea2fac3fd93bdc94` | HOLD | Valid experimental cubic fit; deferred only by the exact batch boundary. |
| 827 | `litcurate_968c1305823124d8` | HOLD | Valid static-GGA cubic fit; deferred with its complete family. |
| 828 | `litcurate_066df1f097ee9c37` | HOLD | Valid static-GGA 6H fit; deferred with its complete family. |
| 829 | `litcurate_54d86988e43e1882` | REJECT | MgSiO3 comparison belongs to its cited primary source. |
| 830 | `litcurate_49c328eec3263a4f` | REJECT | CaSiO3 comparison belongs to its cited primary source. |

Result: **0 production records, 3 held complete records, 2 rejected citation traces**.
