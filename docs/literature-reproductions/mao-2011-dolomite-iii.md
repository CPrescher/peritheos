# Mao et al. (2011): Fe-dolomite and dolomite-III

## Scope and primary source

This audit covers all three LitCurate candidates attached to Mao et al.,
*Dolomite III: A new candidate lower mantle carbonate*, *Geophysical Research
Letters* **38**, L22303 (2011),
[doi:10.1029/2011GL049519](https://doi.org/10.1029/2011GL049519). The
free-access publisher article is the scientific authority. LitCurate was used
only for discovery and candidate enumeration.

The measured natural Windham, Vermont sample is
`Ca0.988Mg0.918Fe0.078Mn0.016(CO3)2`. It was mixed with 5 wt.% Pt, loaded in Ne,
and measured by synchrotron XRD at GSECARS. Pt and the pressure scale of Fei et
al. (2007) provided pressure. The phase-III observations are room-temperature
measurements after laser annealing near 1500 K, not high-temperature P-V points.

## Candidate disposition

| Candidate | Source quantity | Decision |
|---|---|---|
| `litcurate_ca7038ed5d3c92bc` | Trigonal Fe-dolomite BM3: `V0=321.77(6) A^3`, `K0=94.1(4) GPa`, fixed `K0'=4` | **Accepted** as `ca0988mg0918fe0078mn0016c2o6_dolomite_mao_2011_bm3_1`. |
| `litcurate_883e64429604ae2d` | Low-spin monoclinic dolomite-III BM3: `V0=231.8(7) A^3`, `K0=184(4) GPa`, fixed `K0'=4` | **Accepted** as `ca0988mg0918fe0078mn0016c2o6_dolomite_iii_mao_2011_bm3_low_spin_2`. |
| `litcurate_9fea9bd5e6a1abd0` | High-spin monoclinic dolomite-III BM3: `V0=239.2(15) A^3`, `K0=164(8) GPa`, fixed `K0'=4` | **Accepted** as `ca0988mg0918fe0078mn0016c2o6_dolomite_iii_mao_2011_bm3_high_spin_1`. |

The three records are not combined: the first is a different crystal phase,
and the two phase-III parameterizations bound different volume/spin branches
across the approximately 47 GPa discontinuity.

## Crystal and volume basis

The starting Fe-dolomite is trigonal R-3 with conventional-hexagonal
`a0=4.8126(8) A`, `c0=16.0413(18) A`, and `Z=3`. Those lattice constants give
`321.7583 A^3`, independently confirming that the printed `V0=321.77 A^3` is a
conventional-cell volume.

For dolomite-III, most peaks could be indexed by a monoclinic cell with
`beta=108.3 degrees`, but the paper explicitly allows alternative indexings and
says that future work is needed to solve the structure. The phase-transition
volume collapses in Figure 3 put all three cells on a directly comparable
`Z=3` basis. Peritheos therefore preserves the monoclinic indexing-cell volume
but does not invent a space group, lattice metric at zero pressure, or atomic
coordinates. Both phase-III `V0` values are extrapolated EOS parameters, not
recovered ambient structures.

## Figure 3 digitization and fit masks

Neither the article nor its auxiliary material includes a numerical P-V table.
The publisher's 370 x 512 pixel Figure 3 asset was retained by SHA-256
`cb20211a4551d2ec74a1ece62ec9c65754e63f1331cac5c57498d43bb6b388fe`.
Linear tick calibration gives:

- panel a: 3.225 pixels/GPa and 51 pixels per 40 A^3; seven separable
  open-diamond Fe-dolomite centers are bundled;
- panel b: 4.125 pixels/GPa and 43 pixels per 10 A^3; seven separable
  larger-volume and 18 separable lower-volume phase-III centers are bundled.

One-pixel center-picking bounds are stored separately and are not experimental
uncertainties. Overlapped symbols and error-bar strokes were excluded.

The visible spin mask is the larger-volume decompression branch for high spin
and the lower-volume compression/decompression branch for low spin. There is a
source-level ambiguity: the prose says both compression and decompression data
were used for the high-spin fit, whereas the plotted high-spin branch is
represented by open decompression symbols and the filled compression symbols
follow the lower-volume curve. The repository therefore records the visible
branch mask and this contradiction; it does not claim recovery of the authors'
private row selection or regression weights.

## Numerical reproduction

The deterministic script `scripts/reproduce_mao_2011_dolomite.py` evaluates
the standard BM3 pressure curve and performs unweighted pressure- and
volume-residual refits with `K0'=4` fixed.

| Record | Points | Published-curve pressure RMSE | Published-curve volume RMSE | Diagnostic unweighted pressure refit |
|---|---:|---:|---:|---|
| Fe-dolomite | 7 | 0.3992 GPa | 0.9515 A^3 | `V0=320.6463 A^3`, `K0=98.6354 GPa` |
| phase III, high spin | 7 | 2.9551 GPa | 1.9050 A^3 | `V0=239.6659 A^3`, `K0=172.3440 GPa` |
| phase III, low spin | 18 | 2.9152 GPa | 1.2268 A^3 | `V0=228.7388 A^3`, `K0=209.4044 GPa` |

The low-pressure curve follows every separable marker within 0.70 GPa. The
phase-III raster refits are not substitutes for the published parameters: they
move strongly along the `V0-K0` tradeoff because all observations are 26--84
GPa while `V0` is extrapolated to zero pressure. The authors explicitly warn
that `V0` is unknown and that changing fixed `K0'` by 0.5 changes `K0` by up to
15% and `V0` by 1.3%. The published curves still follow their respective
plotted branches within 1.95 A^3, which is sufficient to validate equation and
branch identity without pretending that the unreported weighting has been
recovered.

## Pressure-scale and validity limits

The exact Pt calibration reference is resolved, but Figure 3 does not expose
row-wise Pt volumes, so pressure recalculation is impossible. The trigonal
record is bounded at the approximately 17 GPa transition to dolomite-II. The
phase-III records are bounded to their plotted branch envelopes and 300 K;
they are not equilibrium phase-boundary models. The paper could not determine
an EOS for dolomite-II from its limited pressure range, so no phase-II record
was created.

An exact DOI search of the local Zotero library on 2026-09-05 returned no item.
No Zotero import was performed in this scoped audit.
