# Tange et al. (2009) pressure-scale-free MgO audit

## Outcome

The Fit 3 Vinet--Mie--Grüneisen--Debye implementation is partially validated,
but the paper's full global least-squares fit is not reproducible from published
row-level data. The audit recovers 102 redistributable observations from four
of the eight input sources. At the published coefficients these rows reproduce
the paper's independently reported group residuals: the LASL shock-pressure
RMSE is `1.3723 GPa` (published `1.4 GPa`) and the Li et al. adiabatic-modulus
RMSE is `1.5022 GPa` (published `1.5 GPa`). The ambient constraints also give
`K_T0 = 160.6273 GPa` and `theta0 = 761.476 K`, reproducing the printed
`160.63 GPa` and `761 K`.

This evidence does **not** establish coefficient parity for the global fit.
Zha et al. publish polynomial elastic summaries and figures, but not the 27
row-level volume--`K_S` values used by Tange. Exact weights are also unavailable
for some legacy elastic data, and COD exposes 36 Fiquet MgO states where Tange
Table 1 states 37. The production record therefore remains source-authoritative.
The qualified reconstruction below is classified `similar`, while exact
coefficient parity remains unavailable.

A second, explicitly approximate route reconstructs 164 of the stated 165 rows
without redistributing the restricted legacy tables. It reaches the repository's
numerical `similar` criterion for all four free coefficients. This is strong
evidence that the implemented coupled objective is close to the source analysis,
but the curve-based Zha surrogate and reconstructed elastic weights rule out an
upgrade to strict parity.

## Primary model and fixed conventions

The scientific authority is Tange, Nishihara, and Tsuchiya, “Unified analyses
for P-V-T equation of state of MgO: A solution for pressure-scale problems in
high P-T experiments,” *Journal of Geophysical Research* **114**, B03208
(2009), [doi:10.1029/2008JB005813](https://doi.org/10.1029/2008JB005813).
The 2010 correction changes Figure 11 only.

For `x = V/V0`, Fit 3 uses the Vinet 300 K reference curve and

\[
\gamma(V)=\gamma_0\{1+a(x^b-1)\},
\]

\[
\Theta(V)=\Theta_0\exp\left[-\gamma_0\left((1-a)\ln x+
\frac{a}{b}(x^b-1)\right)\right].
\]

The Debye thermal pressure is referenced to 300 K. The independent fit variables
are only `K0_prime`, `gamma0`, `a`, and `b`. The following quantities are fixed:

- `V0 = 74.698 A^3` per conventional `Z=4` cell;
- `K_S0 = 162.83 GPa`;
- `alpha0 = 3.17e-5 K^-1`;
- `C_P0 = 37.4 J mol^-1 K^-1`;
- `T0 = 300 K`; and
- `n = 2` atoms per MgO formula unit.

At every trial `gamma0`, the implementation recomputes the two dependent
parameters rather than fitting them independently:

\[
K_{T0}=\frac{K_{S0}}{1+\alpha_0\gamma_0T_0},\qquad
C_{V0}=\frac{C_{P0}}{1+\alpha_0\gamma_0T_0}.
\]

`theta0` is then the root of the Debye heat-capacity equation
`C_V(theta0,T0)=C_V0`. The public
`Tange2009Debye.constrained_reference_parameters` method makes these source
constraints executable and testable.

## Row-level source inventory

Tange Table 1 identifies eight earlier datasets. The audit traced each one
instead of treating Tange's derived Tables 4--5 as observations.

| Source class | Source | Rows stated | Rows located | Disposition |
|---|---|---:|---:|---|
| thermal expansion | Dubrovinsky & Saxena (1997) | 25 | 25 | bundled from COD records 9006456--9006480 (CC0) |
| thermal expansion | Fiquet et al. (1999) | 37 | 36 | all COD MgO records bundled (CC0); one stated row unresolved |
| ambient-pressure `K_S(T)` | Isaak et al. (1989) | 16 | 16 | located in NIST SRD 30; not redistributed |
| ambient-pressure `K_S(T)` | Sinogeikin & Bass (2000) | 15 | 15 | located in article Table I; not redistributed |
| 300 K `K_S(V)` | Zha et al. (2000) | 27 | 0 exact; 27 local surrogates | paper gives finite-strain summaries and plots only |
| 300 K `K_S(V)` | Li et al. (2006) | 17 | 17 | official Table 1 already bundled |
| shock Hugoniot | Marsh (1980), LASL | 24 | 24 | machine-readable public-domain block bundled |
| shock Hugoniot | Duffy & Ahrens (1993) | 4 | 4 | located in article Table 1; not redistributed |

The [Crystallography Open Database](https://www.crystallography.net/cod/)
places its database content under CC0. The bundled COD rows retain their stable
record identifiers, temperature, lattice parameter, and reported volume. For
Dubrovinsky, the selection stops at 2986 K, matching Tange's 2990 K upper limit;
the five hotter COD records are not silently included. For Fiquet, all 21
IR-wire and 15 Re-wire MgO records are retained and the 36-versus-37 discrepancy
is not filled by interpolation.

The Marsh CSV contains only the source quantities consumed by the archived
LASL program: shot identifier, initial density, shock velocity `U_s`, and
particle velocity `u_p`. Pressure and compression are recalculated, not copied
from an EOS-derived table:

\[
P_H=\rho_0U_su_p,\qquad V_H/V_0=1-u_p/U_s.
\]

The scoped [redistribution notice](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/mgo-tange-2009-redistributable-inputs.LICENSE.md)
covers only rights held in the CSV transcriptions. It does not license the
cited journal articles. No Isaak, Sinogeikin, or Duffy table is copied into the
repository because the inspected access points do not grant a compatible
redistribution right.

## Coupled partial objective

The executable subset uses the same residual variables and fixed conventions
that Tange describes:

- thermal-expansion rows contribute calculated pressure divided by the stated
  `0.3 GPa` uncertainty at nominal ambient pressure;
- Li elasticity rows contribute `(K_S,calc-K_S,obs)/sigma_KS`. Tange maps both
  the modulus difference and its uncertainty to pressure with the same local
  derivative, so the derivative cancels in the standardized residual. The
  pressures listed parenthetically for the high-pressure static-elasticity
  sources are derived coordinates and are not treated as fit observations;
- LASL shock rows contribute `(P_H,calc-P_H,obs)/3 GPa`; and
- `K0` and `theta0` are recomputed from the ambient anchors on every evaluation.

For each observed shock compression, Rankine--Hugoniot energy is

\[
\Delta E_H=\tfrac12P_H(V_0-V_H).
\]

The 300 K compression-energy change is evaluated from the Vinet Helmholtz
integral and the MGD entropy change. The remaining energy produces thermal
pressure `gamma(V) Delta E_th/V`. This implements the energy reduction without
using a measured shock temperature; Tange likewise lists the Duffy temperatures
as observations not used in the optimization.

## Numerical result and interpretation

| Check at published coefficients | Rows | RMSE |
|---|---:|---:|
| ambient-pressure thermal residual | 61 | 0.2633 GPa |
| Li `K_S` residual | 17 | 1.5022 GPa |
| LASL Hugoniot-pressure residual | 24 | 1.3723 GPa |

The exact dependent-parameter reconstruction and agreement with the two
source-quoted group residuals validate the implemented thermodynamic and shock
conventions unusually strongly for a partial dataset.

As a deliberate identifiability test, the script also optimizes the four Fit 3
variables against only these 102 rows. It converges to `K0_prime=4.3312`,
`gamma0=1.5105`, `a=1.0000`, and `b=0.5033`; `a` hits its allowed upper bound,
far from the published `0.138`. This is not a replacement MgO EOS. It shows
why a subset fit must not be reported as reproduction of the eight-source
global coefficients.

The machine-readable result is
[`tange-2009-mgo-partial-validation.json`](../data/tange-2009-mgo-partial-validation.json).
It records every source disposition, objective convention, group metric,
dependent-parameter check, and the subset-fit diagnostic.

## Local approximate global reconstruction

The second route adds four local CSVs that are deliberately not placed in the
package: the 16 Isaak states, 15 Sinogeikin states, a 27-row Zha surrogate, and
four Duffy shock states. The generated report retains only row counts, input
hashes, assumptions, residual summaries, and fit results.

The Isaak `K_S` values are calculated from the NIST elastic tensor as
`(C11+2*C12)/3`; the separate NIST page labelled “Bulk Modulus” is inconsistent
with that tensor and was not substituted. The Sinogeikin rows use the article's
1% stated modulus-accuracy bound. Zha reports 28 calculated moduli including the
ambient constraint, so the local reconstruction uses 27 non-ambient compression
locations and evaluates the published third-order finite-strain relation with
`K_S0=162.5 GPa` and `K_S0'=3.99`. A 2% Zha modulus uncertainty is the central
assumption. Duffy's four `P-rho_H` states use Tange's `3 GPa` shock weighting.

| Free coefficient | Published | Approximate refit | Relative difference |
|---|---:|---:|---:|
| `K0_prime` | 4.367 | 4.28527 | 1.87% |
| `gamma0` | 1.442 | 1.43678 | 0.36% |
| `a` | 0.138 | 0.14715 | 6.63% |
| `b` | 5.4 | 5.38459 | 0.29% |

All four pass the standard parameter-specific similarity limits. Repeating the
fit with Zha uncertainties of 1%, 2%, and 3% also leaves all four coefficients
inside those limits, so the conclusion is not confined to one tuned weight.
The result is recorded in
[`tange-2009-mgo-approximate-refit.json`](../data/tange-2009-mgo-approximate-refit.json).
It therefore enters the canonical ledger as `similar`, with an explicit
qualification because the external rows are not redistributed, one Fiquet state
is missing, and neither the exact Zha observations nor every original row weight
is published. Those limitations prevent `parity`; they do not negate the
completed numerical-similarity result.

## Reproduction

```console
uv run python scripts/reproduce_tange_2009_mgo.py
uv run python scripts/reproduce_tange_2009_mgo.py --check
uv run python scripts/reproduce_tange_2009_mgo.py \
  --approximate-inputs /path/to/tange-2009-local-inputs
uv run python scripts/reproduce_tange_2009_mgo.py \
  --approximate-inputs /path/to/tange-2009-local-inputs --check
```

The remaining step for exact parity is source recovery rather than numerical
optimization: obtain the exact 27-row Zha volume--`K_S` dataset, resolve the
missing Fiquet state, and obtain the precise uncertainties or final residual
weights used for the legacy elastic series. Until then, the source's Fit 3
coefficients remain authoritative.
