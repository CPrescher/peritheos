# Fiquet et al. (2000) MgSiO3 bridgmanite audit

## Primary source and experiment

G. Fiquet, A. Dewaele, D. Andrault, M. Kunz, and T. Le Bihan,
“Thermoelastic properties and crystal structure of MgSiO3 perovskite at lower
mantle pressure and temperature conditions,” *Geophysical Research Letters*
27, 21-24 (2000),
[doi:10.1029/1999GL008397](https://doi.org/10.1029/1999GL008397).

The sample is pure orthorhombic `Pbnm` MgSiO3 perovskite synthesized from
synthetic enstatite or glass. Angle-dispersive synchrotron XRD was collected in
a laser-heated DAC at ESRF ID30. Platinum was the internal pressure calibrant,
using the Jamieson et al. (1982) P-V-T EOS. The source reports pressures to 94
GPa and temperatures near 3000 K.

## Accepted room-temperature BM3

Section 3.3 separately fits the volume data recorded at room temperature after
laser heating with third-order Birch-Murnaghan:

- `V0=162.27(1) A^3` per conventional `Z=4` cell;
- `K0=253(9) GPa`;
- `K0'=3.9(2)`;
- all three parameters fitted.

The bundled dataset transcribes all 25 Table 1 rows at exactly 298 K from both
side-by-side printed blocks, including every pressure and volume uncertainty.
The 13 high-temperature rows are excluded because they feed the separate
global thermal inversion.

Run:

```bash
.venv/bin/python scripts/reproduce_fiquet_2000_bridgmanite.py
```

The rounded published curve has 0.7070 GPa pressure RMSE. Holding the precisely
reported `V0` only for a deterministic audit and refitting the other two
coefficients gives `K0=253.3244 GPa` and `K0'=3.91619`, both well inside the
published uncertainties. An all-free unweighted diagnostic is also emitted;
it is not treated as the authors' generalized-least-squares fit because the
full covariance and numerical input configuration are unavailable.

## LitCurate dispositions

| Row | Identifier | Disposition | Reason |
| ---: | --- | --- | --- |
| 238 | `litcurate_33cdcf0b45faeb53` | **Accepted** | Exact post-heating 298 K BM3 mapped to `bridgmanite_fiquet_2000_bm3_1`. |
| 239 | `litcurate_6f7c6fbdde4dc27b` | Withheld: incomplete thermal model | Table 3's Mie-Grüneisen alternative combines this study with Fiquet et al. (1998); LitCurate does not supply the complete executable thermal parameter set and the row is not an independent isotherm. |
| 240 | `litcurate_77110b144bbc0b4f` | Withheld: composite thermal inversion | Table 3 fit (1) combines this paper and Fiquet et al. (1998), 63 P-V-T points, with `V0=162.3 A^3` fixed; it is distinct from the accepted 25-point room-temperature BM3. |
| 241 | `litcurate_3ec948abda241bcd` | Citation/comparison composite | Table 3 fit (2) combines 213 points from Wang et al. (1994), Funamori et al. (1996), and Saxena et al. (1999), not measurements made in this paper. |
| 242 | `litcurate_5137740ead6facf4` | Citation trace | `K0=254(13) GPa`, fixed `K0'=4`, is from Ross and Hazen (1990). |
| 243 | `litcurate_3c41e6fd61cbd9e9` | Citation trace | `K0=256(7) GPa` is from Fiquet et al. (1998). |
| 244 | `litcurate_2128003b5d814531` | Citation trace | `K0=261 GPa`, fixed `K0'=4`, is from Mao et al. (1991). |

Production result: **one net-new EOS record**. The thermal composites and
literature comparisons are documented without inflating the production count.
