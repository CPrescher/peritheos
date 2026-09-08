# Dorogokupets et al. (2015) MgSiO3 reconstruction

## Outcome

The three generalized Rydberg-Stacey records remain faithful transcriptions of
Dorogokupets et al. (2015), but they are no longer described as having no usable
observations. The three supplied primary papers add the complete Wang (2004),
Komabayashi (2008), and Zhou (2014) numerical tables. This audit applies the cited
Au or MgO pressure scale where the required calibrant information survives and
performs executable room-temperature fits with `V0` and `k=5` fixed.

None of those fits reproduces the published `K0`, `K0_prime` pair. The source
coefficients are therefore retained. The result is **partial reconstruction;
parity not achieved**, not a substitute regression.

## What the source optimized

Dorogokupets et al. describe a simultaneous thermoelastic optimization using heat
capacity, thermal expansion, adiabatic bulk modulus, and P-V-T measurements. The
paper says that `V0` is usually taken from reference data, `K0` is determined from
ultrasonic measurements after conversion from `KS`, and `k` is fixed at 5. It does
not state the phase-specific residual vector, a point-weighting rule, fitted-parameter
covariance, or unrounded coefficients. Unweighted use of the published points is
therefore the primary reconstruction; an uncertainty-weighted regression is retained
only as a sensitivity check. A 298 K P-V regression can test the reference isotherm
but cannot reproduce the joint thermal/acoustic objective.

## Reconstructed observations

| Phase / branch | Recoverable rows | Pressure coordinate | Explicit limitation |
| --- | ---: | --- | --- |
| Bridgmanite | 9 Katsura (300-308 K) + 10 Tange (300 K) | Sokolova MgO | Tange retains absolute MgO volumes; Katsura prints rounded MgO `V/V0`, so the absolute-volume mapping is approximate. |
| Akimotoite, Wang | 23 Wang 298 K rows (134-row source table complete) | T0133 recalculated to Sokolova Au; T0150 original Decker NaCl | The Au-marker volume is recovered by inverting Wang's stated Anderson (1989) Au EOS before evaluating the 2015 Au scale. |
| Akimotoite, ruby | 16 Reynard rows | Published ruby pressure | Same experiment as the ice-VII branch; ruby calibration is not identified in the primary paper. |
| Akimotoite, ice VII | 14 Reynard rows | Published preferred ice-VII pressure | Bracketed pressure values are maximum gradients, not one-sigma uncertainties. |
| Post-perovskite | 24 Guignot + 3 Komabayashi 300 K rows | Sokolova MgO from measured lattice parameters | The full Komabayashi table has 22 rows: 21 usable post-perovskite volumes and one explicitly missing volume. Ono's offset Au-scale comparison branch remains excluded. |

The combined reconstruction preserves source DOI, source row, temperature,
volume and uncertainty, recalculated pressure and uncertainty proxy, pressure
coordinate, and the inclusion decision. The source transcriptions preserve all
printed table columns, parenthetical uncertainties, calibrant lattice parameters,
and explicitly missing cells. The Zhou table contributes 58 akimotoite rows, of
which 55 contain `KS(P,T)`.

## Numerical diagnostics

The primary result minimizes unweighted pressure residuals. The second result uses
the source-reported pressure error or spread proxy combined in quadrature with
volume uncertainty projected through the published curve. That second objective is
a sensitivity diagnostic only; no inspected paper says Dorogokupets used it.

| Phase / branch | Published `K0`, `K0_prime` | Unweighted fit | Diagnostic weighted fit | Published-curve RMSE | Refit RMSE |
| --- | --- | --- | --- | ---: | ---: |
| Bridgmanite | 252.0, 4.38 | 256.443, 3.850 | 257.510, 3.814 | 1.578 GPa | 0.302 GPa |
| Akimotoite, ruby | 215.3, 4.91 | 213.574, 6.545 | 212.073, 6.734 | 0.970 GPa | 0.609 GPa |
| Akimotoite, ice VII | 215.3, 4.91 | 229.174, 3.180 | 206.573, 5.856 | 0.691 GPa | 0.656 GPa |
| Akimotoite, Wang only | 215.3, 4.91 | 270.911, 1.000 (lower bound) | 270.911, 1.000 (lower bound) | 1.593 GPa | 0.736 GPa |
| Akimotoite, Wang + Reynard ruby | 215.3, 4.91 | 254.928, 2.296 | 240.924, 3.889 | 1.372 GPa | 0.844 GPa |
| Post-perovskite, Guignot + Komabayashi | 253.7, 4.03 | 246.783, 4.016 | 248.974, 3.953 | 3.844 GPa | 0.755 GPa |

For akimotoite, the 2015 caption says that Reynard pressures are retained on their
“original” scale, but Reynard publishes both ruby and preferred ice-VII assignments.
Choosing either branch materially changes the fitted elastic coefficients. The
ambiguity cannot be resolved from the cited source.

The Wang-only result is not evidence of bad data. Wang's two runs publish ambient
cell volumes from 263.22 to 264.46 A3, whereas the 2015 model fixes `V0` at
262.531 A3. Forcing those observations into a two-parameter static slice therefore
pushes `K0_prime` to the diagnostic bound. The full 2015 optimization also contains
thermal expansion and the Zhou acoustic constraints and is not equivalent to this
static projection.

## Direct weighting check from Zhou

The original-scale [Wang, Zhou and Komabayashi refits](wang-zhou-komabayashi-mgsio3.md)
now support four separate published EOS records. Their free coefficients are
recovered within printed errors, without claiming recovery of this later
Dorogokupets joint objective.

An unweighted linear regression of all 55 printed `KS(P,T)` observations gives
`K0S=219.867 GPa`, `dKS/dP=4.6070`, and `dKS/dT=-0.02322 GPa/K`. Zhou reports
`219.4 GPa`, `4.62`, and `-0.0228 GPa/K`. Weighting by the printed `KS`
uncertainties gives `219.703 GPa`, `4.6197`, and `-0.02339 GPa/K`. Both are close;
the unweighted result has the slightly smaller ordinary RMSE. This supports treating
the rows as unweighted unless a source explicitly says otherwise, and shows that
weighting cannot explain the remaining Dorogokupets coefficient mismatch.

## Missing inputs and source-scale boundary

- The 2015 paper states no residual normalization, point-weighting rule, covariance,
  or phase-specific fixed/free parameter mask.
- Wang's Table 1 prints derived Au pressures rather than measured Au volumes. Its
  pressure-scale conversion is nevertheless reproducible by inversion of the stated
  Anderson (1989) Au EOS; it is one step removed from a directly printed calibrant volume.
- Ashida et al. (1988) calorimetric and thermal-expansion rows remain outside this
  audit. The newly supplied Wang, Zhou, and Komabayashi tables are complete.
- The implemented MgO scale follows Peritheos' corrected Sokolova workbook path.
  Its distinction from the literal 2013 high-temperature expression is negligible
  at 298-308 K but would matter for a future full P-V-T reconstruction.

## Reproduction

```console
UV_CACHE_DIR=/tmp/peritheos-uv-cache uv run --frozen python scripts/reproduce_dorogokupets_2015_mgsio3.py --write-dataset --write-report
```

Artifacts:

- `peritheos/data/datasets/dorogokupets-2015-298k-literature-reconstruction.csv`
- `peritheos/data/datasets/dorogokupets-2015-298k-reconstruction-source.json`
- `peritheos/data/datasets/mgsio3-post-perovskite-guignot-2007-table1-300k-compression.csv`
- `peritheos/data/datasets/akimotoite-wang-2004-tables1-2-pvt.csv`
- `peritheos/data/datasets/akimotoite-zhou-2014-table1-elasticity.csv`
- `peritheos/data/datasets/mgsio3-post-perovskite-komabayashi-2008-table1-pvt.csv`
- `docs/data/dorogokupets-2015-298k-refit.json`

## Primary references

- Dorogokupets et al. (2015), DOI `10.1016/j.rgg.2015.01.011`.
- Katsura et al. (2009), DOI `10.1029/2009GL039318`.
- Tange et al. (2012), DOI `10.1029/2011JB008988`.
- Reynard et al. (1996), DOI `10.2138/am-1996-1-206`.
- Guignot et al. (2007), DOI `10.1016/j.epsl.2007.01.025`.
- Wang et al. (2004), DOI `10.1016/j.pepi.2003.08.007`.
- Komabayashi et al. (2008), DOI `10.1016/j.epsl.2007.10.036`.
- Zhou et al. (2014), DOI `10.1016/j.pepi.2013.06.005`.
