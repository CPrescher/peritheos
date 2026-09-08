# Ye (2017) 300 K fits and Zhu (2025) thermal P-V-T standards

## Attribution and record split

The 300 K Vinet fits are now attributed directly to Ye et al. (2017):

- `gold_ye_2017_vinet_300k`
- `platinum_ye_2017_vinet_300k`
- `mgo_ye_2017_vinet_300k`

The former `*_zhu_2025_vinet_300k` identifiers remain aliases, so existing
callers continue to resolve, but they are not canonical record names. Zhu et
al. fix Ye's 300 K coefficients and fit the broader temperature-dependent
models. Those are separate production records:

- `gold_zhu_2025_pvt`
- `platinum_zhu_2025_pvt`
- `mgo_zhu_2025_pvt`

The scientific sources are Ye et al., *JGR Solid Earth* 122 (2017),
DOI [`10.1002/2016JB013811`](https://doi.org/10.1002/2016JB013811), and Zhu et
al., ESSOAr version 1 (2025),
DOI [`10.22541/essoar.176236186.65259830/v1`](https://doi.org/10.22541/essoar.176236186.65259830/v1).

## Ye 300 K reproduction

Ye Section 3.3 fixes `V0` and `K0` to Dorogokupets and Dewaele (2007), fits
only `K0'`, and reports approximately 1 GPa random scatter. The corrected
publisher Data Sets S1-S2 are bundled as lossless numerical normalizations.

| Ye record | Selected pathway | Refit K0' | Refit error | Published | Rows | RMSE (GPa) |
|---|---|---:|---:|---:|---:|---:|
| Au | MgO-D07 -> Au | 5.896768 | 0.022221 | 5.897(22) | 28 | 0.5712 |
| MgO | Au-D07 -> MgO | 4.182461 | 0.018958 | 4.182(19) | 28 | 0.5714 |
| Pt | Au-D07 -> MgO -> Pt | 5.230473 | 0.032655 | 5.230(33) | 39 | 0.4766 |

All point estimates and uncertainties reproduce at published precision. The
publisher does not state an open license for the corrected supplements; the
adjacent dataset notice applies CC0 only to Peritheos's factual normalization
and arrangement.

## Zhu thermal model

Zhu's pressure is the Ye 300 K Vinet isotherm plus a 300 K-referenced Debye
thermal pressure and a volume-dependent excess term:

```text
gamma(V) = gamma0 {1 + a[(V/V0)^b - 1]}
F_ex(V,T) = -1/2 beta0 (V/V0)^m T^2
P_ex(V,T) = 1/2 beta0 m/V0 (V/V0)^(m-1) T^2
P(V,T) = P_Vinet(V) + [P_D(V,T)-P_D(V,300)]
                       + [P_ex(V,T)-P_ex(V,300)]
```

`beta0` is stored in J mol^-1 K^-2. The executable records use the revised
coefficients in the released Mendeley Data v3 optimizers and downstream
property scripts:

| Record | theta0 (K) | gamma0 | a | b | n | beta0 (J mol^-1 K^-2) | m |
|---|---:|---:|---:|---:|---:|---:|---:|
| Au | 180 | 2.93 | 0.75 | 2.7 | 1 | 0 | 1 |
| Pt | 240 | 2.75 | 0.39 | 5.1 | 1 | 0.002145 | 0.65 |
| MgO | 761 | 1.53 | 1.0 | 1.43 | 2 | -0.0008061 | 4.8 |

The v3 release is authoritative for the executable records, but it is not
internally identical in every file. The Au and Pt optimizer files use
`theta0=180` and `240 K`, while their pressure calculators use `170` and
`230 K`; the MgO optimizer gives `gamma0=1.531...` and its property script
rounds that to `1.53`, while the calculator hard-codes `1.52`. The preprint's
version-1 table also contains the earlier
Pt (`gamma0=2.79`, `a=0.60`, `b=3.5`, `m=0.1`) and MgO (`b=1.5`) values.
Peritheos therefore treats the optimizer/property set as the thermal-fit
target. The standalone calculators are retained only as an inconsistency
diagnostic.

## Released v3 data and thermal-refit parity

Mendeley Data version 3,
DOI [`10.17632/6kxnhc2g73.3`](https://doi.org/10.17632/6kxnhc2g73.3), is CC BY
4.0. Every active numeric row used by the released fits is bundled:

| Material | Input | Rows | Bundled file |
|---|---|---:|---|
| Au | shock states | 12 | `peritheos/data/datasets/zhu-2025-au-shock.csv` |
| Au | zero-pressure thermal expansion | 10 | `peritheos/data/datasets/zhu-2025-au-zero-pressure-thermal-expansion.csv` |
| Pt | shock states | 68 | `peritheos/data/datasets/zhu-2025-pt-shock.csv` |
| Pt | zero-pressure thermal expansion | 17 | `peritheos/data/datasets/zhu-2025-pt-zero-pressure-thermal-expansion.csv` |
| MgO | P-V-T states | 213 | `peritheos/data/datasets/zhu-2025-mgo-pvt.csv` |

The reproduction independently translates the source optimizer rather than
calling MATLAB. It follows the iterative energy-balance transformation,
recovers shock temperatures where required, applies the published stress and
excess-energy corrections, and performs leverage-adjusted Tukey-bisquare
nonlinear least squares. The optimizer fixes `a`, `theta0`, `beta0`, and `m`
and refits `gamma0` and `b`.

| Record | Rows | Refit gamma0 | Refit b | Stored target | Result |
|---|---:|---:|---:|---|---|
| Au | 22 | 2.93599996 | 2.62559912 | 2.93, 2.7 | parity at source precision; both differences are well inside the reported Au uncertainties |
| Pt | 85 | 2.75224510 | 5.10989873 | 2.75, 5.1 | parity at source precision |
| MgO | 213 | 1.53118931 | 1.43014746 | 1.53, 1.43 | parity at source precision |

All 320 released observations participate in the corresponding refits. The
standalone calculator constants produce detectable differences from the
fit-derived records (up to about 0.015 GPa for Au, 0.029 GPa for Pt, and
0.113 GPa for MgO over the diagnostic grid), confirming that calculator
agreement must not be substituted for thermal-refit parity.

Run:

```console
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run --frozen python scripts/reproduce_zhu_2025_pressure_standards.py
```

## Exhaustive LitCurate disposition

| Candidate | Decision | Reason |
|---|---|---|
| `litcurate_b1a3d7a39cb0fae8` | **accepted after primary correction** | The candidate transcribed MgO `V0=74.71` and `K0=160.3`, but column-shifted Pt's `K0'=5.230` into MgO. Ye Table 1 gives MgO `K0'=4.182(19)`. |
| `litcurate_fdde71245d253100` | **rejected here** | `K0'=4.367` is a comparison value without a complete same-source parameterization. |
