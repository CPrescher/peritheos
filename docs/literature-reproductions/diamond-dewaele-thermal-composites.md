# Dewaele-anchored Correa and Benedict diamond composites

## Disposition

`diamond_correa_2008_dewaele_anchored` and
`diamond_benedict_2014_dewaele_anchored` are **source-equation
reconstructions, not independent coefficient refits**. Neither Correa et al.
(2008), Benedict et al. (2014), nor Dewaele et al. (2008) published or fitted
these combined parameterizations. Peritheos deliberately constructs

\[
P_{\rm comp}(V,T)=P_{\rm Dewaele}(V,298\ {\rm K})
 +P_{\rm theory}(V,T)-P_{\rm theory}(V,298\ {\rm K}).
\]

This replaces the theoretical cold/reference pressure with the experimental
298 K Vinet isotherm while retaining the theory model's thermal pressure. The
equivalent Helmholtz construction subtracts the theory free energy at 298 K.
It preserves entropy and fixed-volume thermal internal-energy increments, but
its absolute energy zero is conventional.

## Primary sources and available data

- [Dewaele et al. (2008)](https://doi.org/10.1103/PhysRevB.77.094106),
  Table I, prints the diamond P-V-T observations; Table III gives the 298 K
  Vinet parameters. Section II explicitly places the pressures on Holzapfel's
  2005 H05 ruby scale. The previous Dorogokupets--Oganov attribution in the
  diamond cards was incorrect and is now fixed.
- [Correa et al. (2008)](https://doi.org/10.1103/PhysRevB.78.024101) gives the
  cold Vinet, logarithmic-moment double-Debye, and anharmonic equations and
  their fitted coefficients. It does not tabulate the underlying cold-energy
  grid, phonon DOS/moments, regression weights, or covariance. All 57 diamond
  DFT-MD pressure markers in vector Figure 8 are checked in as a conservative
  plot digitization.
- [Benedict et al. (2014)](https://doi.org/10.1103/PhysRevB.89.224109) gives
  the later first-moment double-Debye equations and coefficients. Its
  [author-supplied source](https://arxiv.org/abs/1311.4577) includes an exact
  96-row solid DFT-MD table of volume, temperature, pressure, and internal
  energy. Those rows validate the finished model; they are not the upstream
  motionless-ion cold-energy and phonon fitting inputs.

The machine-readable resources are
`diamond-dewaele-2008-table1-pvt.csv`,
`diamond-correa-2008-figure8-dft-md-vector-digitized.csv`, and
`diamond-benedict-2014-supplement-solid-dft-md.csv` under
`peritheos/data/datasets/`. Their source locations, extraction methods,
artifact hashes, column roles, and exact versus plot-only status are recorded
in `diamond.eosmat`.

## Numerical reconstruction

The Dewaele 298 K component is independently refitted to the 23 observations
nearest 298 K with `K0=444.5 GPa` fixed, matching the source protocol. It gives
`V0=45.35739774 A^3/conventional cell` and `K0'=4.09473672`, versus published
`45.3544` and `4.18`; both are within combined two-sigma uncertainty. Pressure
RMSE is `0.231284 GPa`.

The complete source thermal branches—not merely their Vinet anchors—are then
evaluated against the available calculations:

| Thermal branch | Source checkpoints | Pressure RMSE | Maximum absolute pressure residual | Caloric check |
|---|---:|---:|---:|---:|
| Correa 2008 | 57 vector-digitized Figure 8 states, 3 isochores | 1.79089 GPa | 5.19024 GPa | source does not tabulate energy rows |
| Benedict 2014 | 96 exact supplementary states, 14 isochores | 3.17793 GPa | 6.64492 GPa | 0.0612094 eV/atom fixed-volume energy-increment RMSE over 82 comparisons |

For every source checkpoint, each derived Peritheos model is also compared
with the explicit equation above. Maximum pressure-composition error is
`1.14e-13 GPa`, and the 298 K reference-isotherm error is exactly zero. Maximum
thermal internal-energy-increment error is `2.39e-9 J/mol` for Correa and
`1.82e-10 J/mol` for Benedict (floating-point roundoff).

These checks answer two different questions and must remain separate:

1. the Dewaele rows support an independent refit of the experimental anchor;
2. the theory rows/markers test the published complete thermal models;
3. the exact identity test reconstructs the Peritheos-derived composites.

The available observations cannot independently optimize all coefficients of
either composite. Doing so would conflate downstream validation states with
the absent electronic-structure and phonon fitting grids.

## Reproduction

Run:

```text
uv run python scripts/reproduce_diamond_thermal_composites.py
uv run python scripts/validate_primary_eos_refits.py
```

The first command regenerates
`docs/data/diamond-thermal-composite-reconstruction.json`. The record-level
ledger classifies both derived records as `reconstructed`, retains the two
standalone theory records as `not_refittable` at coefficient level, and embeds
the independent anchor refit, source-model diagnostics, and complete-model
identity checks.
