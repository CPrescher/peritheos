# Matsui et al. (2012) high-spin ferropericlase EOS

## Source and record scope

Matsui et al., *American Mineralogist* **97**, 176-183 (2012),
[doi:10.2138/am.2012.3937](https://doi.org/10.2138/am.2012.3937), report
static compression of stoichiometric `(Mg0.83,Fe0.17)O` and
`(Mg0.75,Fe0.25)O` at 300, 700, and 1100 K. The final article was audited in
full. The official January 2012 issue page has no deposit-item link for this
paper, but the primary article itself contains the complete EOS equations,
parameters, and observation tables.

Peritheos adds two composition-specific published records:

- `ferropericlase_mg83fe17_matsui_2012_bm3_mgd_1`
- `ferropericlase_mg75fe25_matsui_2012_bm3_mgd_1`

Both records represent only the high-spin B1 (NaCl-type) phase below the
authors' 47 GPa fit cutoff. They do not model the mixed-spin or low-spin
regime.

## Material and volume identity

The authors synthesized stoichiometric ferropericlase and identified both
samples as B1-NaCl. The conventional crystallographic representation is cubic
`Fm-3m` with four formula units per cell: mixed Mg/Fe occupies `4a` and O
occupies `4b`. The reported ambient lattice parameters and fitted reference
volumes are distinct source quantities:

| Material | Ambient `a` (A) | Ambient table volume (A^3) | Fitted `V0` (A^3/cell) |
|---|---:|---:|---:|
| `(Mg0.83,Fe0.17)O` | 4.2331(2) | 75.855(11) | 75.849(11) |
| `(Mg0.75,Fe0.25)O` | 4.2427(2) | 76.368(13) | 76.372(16) |

The small differences are preserved rather than forcing the structural cell
to equal the fitted zero-pressure reference volume.

## Equation mapping

Equations 1-6 define pressure relative to the 300 K isotherm:

\[
P(V,T)=P_{\mathrm{BM3}}(V,300\,\mathrm{K})+
\frac{\gamma(V)}{V}\left[E_D(V,T)-E_D(V,300\,\mathrm{K})\right],
\]

with

\[
\gamma(V)=\gamma_0(V/V_0)^q,
\qquad
\Theta(V)=\Theta_0\exp[(\gamma_0-\gamma(V))/q].
\]

The Debye energy uses `n=2` atoms per `(Mg,Fe)O` formula unit. Because the
source volumes are conventional four-formula-unit cells, Peritheos uses
`formula_units_per_cell: 4` when converting cell volume to molar volume. This
maps exactly to `BM3 + MieGruneisenDebye` with
`debye_temperature_law: integrated_gruneisen`.

| Parameter | `(Mg0.83,Fe0.17)O` | `(Mg0.75,Fe0.25)O` | Fit status |
|---|---:|---:|---|
| `V0` (A^3/cell) | 75.849(11) | 76.372(16) | fitted at 300 K |
| `K0T` (GPa) | 160 | 160 | fixed from prior MgO data |
| `K0T'` | 4.08(2) | 4.22(3) | fitted at 300 K |
| `theta0` (K) | 500 | 500 | fixed from prior ferropericlase estimates |
| `gamma0` | 1.53(4) | 1.64(4) | fitted at 700 and 1100 K |
| `q` | 0.7(2) | 0.7(2) | fitted at 700 and 1100 K |
| `Tr` (K) | 300 | 300 | equation reference |
| `n` | 2 | 2 | formula stoichiometry |

The source does not define the confidence convention for parenthetical
uncertainties and gives no coefficient covariance matrix. Missing errors are
stored as `null`, not zero.

## Primary observations and pressure scale

Primary-article Tables 1 and 2 are bundled verbatim as CSV resources. They
contain 34 and 39 rows, respectively, including sample cell volume, sample
volume uncertainty, paired Au `V/V0`, Au-ratio uncertainty, reported pressure,
pressure uncertainty, authors' calculated pressure, and printed residual.
Fit-selection flags retain the authors' rule of using observations below
47 GPa: 23 rows for the 17% Fe sample and 30 rows for the 25% Fe sample.

Pressure was determined from simultaneously measured Au using Matsui's 2010
thermal Au EOS, [doi:10.1088/1742-6596/215/1/012197](https://doi.org/10.1088/1742-6596/215/1/012197).
The row-wise Au compression ratios and temperatures are available, but that Au
EOS is not yet bundled as an executable Peritheos reference record. Exact
observation-level pressure recalculation is therefore deferred; published
source-scale pressures are preserved unchanged.

## Numerical reproduction

The executable models reproduce the authors' printed calculated pressures for
all fitted rows within 0.014 GPa for `(Mg0.83,Fe0.17)O` and 0.006 GPa for
`(Mg0.75,Fe0.25)O`. Differences reflect printed pressure precision and the gas
constant used by the implementation. Against observed pressure, the published
parameterizations give RMSE values of 0.1233 and 0.1426 GPa, respectively.

An independent staged unweighted least-squares reproduction follows the source
sequence: fit `V0` and `K0T'` to the 300 K rows with `K0T=160 GPa`, then fit
`gamma0` and `q` to the 700 and 1100 K rows using the published cold curve with
`theta0=500 K` and `n=2`.

| Material | `V0` refit | `K0T'` refit | `gamma0` refit | `q` refit |
|---|---:|---:|---:|---:|
| `(Mg0.83,Fe0.17)O` | 75.80705 | 4.13454 | 1.52863 | 0.68230 |
| `(Mg0.75,Fe0.25)O` | 76.37279 | 4.22184 | 1.64550 | 0.74501 |

All 25% Fe refit coefficients agree with the published values within the
printed parameter errors. For 17% Fe, the thermal coefficients agree within
error, but the independently refitted room-temperature pair does not. The
paper does not specify the room-temperature objective, weights, or covariance,
so this non-parity is retained as a qualification rather than used to replace
the published record.

## Validity and spin transition

The executable ranges end at the largest observation below the stated 47 GPa
fit cutoff: 46.39 GPa for 17% Fe and 46.21 GPa for 25% Fe, over 300-1100 K.
Higher-pressure rows remain in the source datasets because they document the
onset of the spin transition. The reported onset intervals are:

| Material | 300 K (GPa) | 700 K (GPa) | 1100 K (GPa) |
|---|---:|---:|---:|
| `(Mg0.83,Fe0.17)O` | 48.9-51.3 | 50.4-52.7 | 51.6-53.5 |
| `(Mg0.75,Fe0.25)O` | 48.1-52.3 | 50.3-54.3 | 47.4-51.4 |

These marginal bounds are not a rectangular phase-stability guarantee, and
the models must not be presented as mixed-spin or low-spin EOSs.
