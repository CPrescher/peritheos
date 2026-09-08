# Wang and Weidner (1994): room-temperature CaSiO3 perovskite BM2

Audit date: 2026-09-08.

## Outcome

The primary publication supports one room-temperature second-order
Birch--Murnaghan fit for metastable cubic CaSiO3 perovskite:
`V0 = 45.83(7) A3`, `K0 = 280(23) GPa`, with `K0' = 4` fixed. The source says
explicitly that only four room-temperature decompression points above 2.0 GPa
were fitted. All four marker centers are visible in Figure 3 and are now bundled
as a `plot_only` dataset. An unweighted pressure-residual BM2 refit gives
`V0 = 45.83216 A3` and `K0 = 279.2707 GPa`, so both coefficients recover the
published values within their printed uncertainties.

This is complete plot-scope recovery, not recovery of an authoritative numerical
P--V table. The article contains no numerical table of the room-temperature
observations and does not disclose the regression weights, residual coordinate,
optimizer, covariance, or confidence convention.

## Primary evidence and archive search

The controlling source is Y. Wang and D. J. Weidner, “Thermoelasticity of
CaSiO3 perovskite and implications for the lower mantle,” *Geophysical Research
Letters* **21**, 895--898 (1994),
[doi:10.1029/94GL00976](https://doi.org/10.1029/94GL00976). The audit checked the
publisher record and the complete author-uploaded paper.

The numerical-data search also checked the following plausible routes:

- Table 1 of the paper is not the missing P--V table. It is a single diffraction
  refinement at 11.2 GPa and 1172 K, ending in `a = 3.4998(5) A` and
  `V = 44.73(2) A3`.
- The 11.7 GPa value in the paper belongs to one of the two thermal-expansion
  isobars. It is not the upper bound of the room-temperature decompression fit.
- The official 1994 NSLS activity report, BNL-52455,
  [doi:10.2172/93726](https://doi.org/10.2172/93726), contains a later
  “Thermal equation of state of CaSiO3 perovskite” beamline contribution. It
  describes a subsequent, broader P--V--T experiment but prints no numerical
  room-temperature table for the 1994 GRL fit.
- Wang's 1991 dissertation predates this experiment. Wang, Weidner, and Guyot's
  later JGR paper, [doi:10.1029/95JB03254](https://doi.org/10.1029/95JB03254),
  does publish a 14-row run-13 table, but it is a distinct subsequent experiment
  and cannot be substituted for the four GRL observations.
- No supplement, publisher data attachment, thesis table, or repository deposit
  exposing the four numerical rows was found. Later reviews repeat the fitted
  coefficients rather than the observations.

The exact no-table limitation is therefore narrow: the four fit inputs survive
only as open-circle markers in Figure 3. The bundled CSV does not claim otherwise.

## Experimental reconstruction

The experiments used the DIA-6-type SAM-85 large-volume apparatus with
4 x 4 mm sintered-diamond anvil truncations at NSLS beamline X17. The starting
material was synthetic high-purity wollastonite containing 1 wt% Cr2O3 and
1 wt% NiO. A 1 mm by 2 mm chamber was filled half with wollastonite and half
with an NaCl--BN mixture; both were dried at 570 K for more than 10 hours.

The authors first raised pressure at room temperature, heated at constant load
to synthesize perovskite, and then measured during decompression. Heating was
used to reduce nonhydrostatic stress. Energy-dispersive diffraction supplied as
many as nine Ca-perovskite peaks over 2.6--0.9 A, and the paper states a unit-cell
volume precision of 0.02 A3.

Pressure was calculated from the NaCl thermal EOS of Decker (1971),
[doi:10.1063/1.1660714](https://doi.org/10.1063/1.1660714), with an estimated
pressure uncertainty of about 0.1--0.2 GPa. This resolves the identity of the
calibration. It does not make row-wise pressure recalculation executable because
the four paired NaCl lattice parameters or volumes are not published.

## Fit scope, exclusions, and weighting

The room-temperature paragraph says: four observations above 2.0 GPa, BM2,
`K0' = 4` assumed. Figure 3 and its caption establish the remaining disposition:

- the four open circles above 2.0 GPa are the complete fit set;
- lower-pressure open circles are excluded because amorphization produces
  anomalous volume behavior; and
- the solid 1-bar Kanzaki et al. point is external comparison data, not a Wang
  and Weidner measurement.

The source describes the near-linear stable branch as extending between about
2 and 9 GPa. The catalog range is therefore corrected from `[2.0, 11.7]` to the
qualitative source range `[2.0, 9.0]`; the digitized marker range is
2.74--9.24 GPa.

No source statement supports a specific weighting scheme. The catalog therefore
does not promote the reported 0.1--0.2 GPa pressure uncertainty or 0.02 A3 volume
precision into row-specific regression sigmas. The reproducible audit fit uses
unweighted pressure residuals. As sensitivity checks, unweighted volume
residuals give `V0 = 45.81954 A3`, `K0 = 283.5404 GPa`, and a diagnostic
errors-in-variables objective using only the stored digitization bounds gives
`V0 = 45.82867 A3`, `K0 = 280.4416 GPa`. All three reasonable objectives remain
inside the broad published coefficient uncertainties, but none identifies the
authors' undisclosed regression protocol.

## Bundled data and numerical result

The machine-readable resource is
[`ca-perovskite-wang-weidner-1994-figure3-digitized.csv`](../../peritheos/data/datasets/ca-perovskite-wang-weidner-1994-figure3-digitized.csv).
It stores the four marker centers, rendered-plot coordinates, conservative
digitization bounds (`0.10 GPa`, `0.01 A3`), and an explicit source-fit flag.
Those bounds describe digitization only; they are not the source's statistical
weights.

| Quantity | Published | Unweighted pressure refit | Difference |
|---|---:|---:|---:|
| `V0` (A3) | 45.83 +/- 0.07 | 45.83216 +/- 0.06979 | +0.00216 |
| `K0` (GPa) | 280 +/- 23 | 279.2707 +/- 21.1323 | -0.7293 |

The published curve has a pressure RMSE of 0.24175 GPa against the four
digitized markers; the refit RMSE is 0.24168 GPa. This negligible improvement
and the recovery under alternate objectives support the printed BM2 without
overstating the precision of plot-derived observations.
