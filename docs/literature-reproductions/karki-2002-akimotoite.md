# Karki and Wentzcovitch (2002): MgSiO3 ilmenite

## Outcome

Four source-generated fourth-order Birch-Murnaghan (BM4) isotherms are accepted as non-default theoretical records: the static lattice and the 300, 1000, and 2000 K quasiharmonic results. One experimental comparison row extracted by LitCurate is a citation trace and is not attributed to this paper.

Primary publication:

- B. B. Karki and R. M. Wentzcovitch, “First-principles lattice dynamics and thermoelasticity of MgSiO3 ilmenite at high pressure,” *Journal of Geophysical Research: Solid Earth* **107**(B11), 2267 (2002), DOI [10.1029/2001JB000702](https://doi.org/10.1029/2001JB000702).
- Full-text archival copy used for the audit: [CiteSeerX](https://citeseerx.ist.psu.edu/document?doi=61abb95b1fabd25290e81ac8133068de9d9b1f25&repid=rep1&type=pdf).
- Publisher version: [Wiley/AGU](https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2001JB000702).

LitCurate is discovery evidence only. The equation family, complete coefficient sets, temperature labels, volume basis, and method were checked in the primary article.

## Source equation and calculation

Section 2 describes density-functional perturbation theory within the local-density approximation, with a 70 Ry plane-wave cutoff and two special Brillouin-zone k points. Structures were fully optimized at each pressure/volume before the phonon calculation. The calculation regime is stated as approximately -11 to 33 GPa.

Section 3.2 defines the quasiharmonic Helmholtz free energy as the static internal energy plus zero-point and thermal phonon contributions. It then states explicitly that a series of thermal EOS isotherms was obtained by fitting a **fourth-order finite-strain equation** to calculated free energy versus volume at each temperature. This is the Eulerian BM4 convention implemented by `birch_murnaghan_4` in Peritheos. Figure 3 plots the static, 300, 1000, and 2000 K curves.

Table 3 reports the volume **per MgSiO3 formula unit**. The existing `akimotoite` material document uses the conventional hexagonal ilmenite cell with Z=6, so each source volume is multiplied by exactly six. No pressure calibration applies to these first-principles calculations.

## Accepted parameterizations

| Stored record | T (K) | Source V0 (A3/formula) | Stored V0 (A3, Z=6) | KT0 (GPa) | K' | K'' (GPa^-1) |
|---|---:|---:|---:|---:|---:|---:|
| `akimotoite_karki_2002_static_bm4_1` | 0 | 43.61 | 261.66 | 210 | 4.57 | -0.041 |
| `akimotoite_karki_2002_300k_bm4_2` | 300 | 44.20 | 265.20 | 201 | 4.64 | -0.042 |
| `akimotoite_karki_2002_1000k_bm4_3` | 1000 | 45.03 | 270.18 | 182 | 4.86 | -0.051 |
| `akimotoite_karki_2002_2000k_bm4_4` | 2000 | 48.59 | 291.54 | 153 | 5.20 | -0.067 |

The source prints no parameter uncertainties or covariance. None is invented. These theoretical isotherms do not replace the existing experimental equilibrium default. The 2000 K record carries the paper’s explicit warning that the quasiharmonic approximation generally becomes inadequate at high temperature and zero pressure.

## Numerical verification and primary-data limitation

The deterministic reproduction script performs the exact per-formula-to-Z=6 conversion and instantiates every BM4 curve. Each curve gives P(V0)=0 and K(V0)=KT0. Compression to 0.90 V0 produces a positive pressure for every isotherm. The source’s expected temperature trends are also preserved: zero-pressure volume rises monotonically and zero-pressure bulk modulus falls monotonically from the static curve to 2000 K.

The fitted free-energy/volume grid is not tabulated or deposited. It is represented by the curves in Figure 3, while the complete fitted coefficient sets are printed in Table 3. Accordingly these records are classified as primary-source-validated parameterizations with plot-only underlying observations, not independently refitted datasets. Digitizing a smooth fitted curve would not supply independent information and would introduce false precision.

Run:

```bash
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run --frozen python scripts/reproduce_karki_2002_akimotoite.py
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run --frozen pytest -q tests/test_karki_2002_akimotoite.py
```

## LitCurate candidate disposition

All five rows associated with this publication DOI are disposed exactly once:

| Candidate | Disposition | Reason |
|---|---|---|
| `litcurate_608b6f9d8d5880b7` | ACCEPT | Source Table 3 static BM4. LitCurate omitted K'', which the primary table supplies as -0.041 GPa^-1. |
| `litcurate_c807f442c329918a` | ACCEPT | Source Table 3 300 K QHA BM4; complete K''=-0.042 GPa^-1 recovered from the primary source. |
| `litcurate_59ae4c4e75952ce6` | ACCEPT | Source Table 3 1000 K QHA BM4; complete K''=-0.051 GPa^-1 recovered. |
| `litcurate_4c132bee48eebb52` | ACCEPT | Source Table 3 2000 K QHA BM4; complete K''=-0.067 GPa^-1 recovered and QHA limitation retained. |
| `litcurate_7f29af4420bee96c` | REJECT / CITATION TRACE | Experimental 300 K comparison values are attributed by Table 3 to Horiuchi et al., Weidner and Ito, and Reynard et al.; they are not a Karki-Wentzcovitch source fit. The relevant experimental akimotoite EOS is already represented under its own primary lineage. |

## Files

- Material: `peritheos/data/materials/akimotoite.eosmat`
- Reproduction: `scripts/reproduce_karki_2002_akimotoite.py`
- Tests: `tests/test_karki_2002_akimotoite.py`
