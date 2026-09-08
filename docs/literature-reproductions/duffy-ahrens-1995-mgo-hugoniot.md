# Duffy and Ahrens (1995): MgO B1 principal Hugoniot

## Source and data status

Duffy and Ahrens report four equation-of-state shots for low-porosity
polycrystalline MgO in Table 3 and the phase-specific relation

`Us = 6.87(10) + 1.24(4) up`,

with velocities in km/s and parenthesized errors defined as one standard
deviation. The authoritative source is the [author-hosted article PDF](https://duffy.princeton.edu/sites/g/files/toruqf616/files/duffy_ahrens_seismo_1995.pdf),
checked at SHA-256
`daaf0a41440bb125dd0a062ca0582c358ceabf04f39c08c40cda2ea73dc67b2b`.
The numerical table is on journal page 533; Equation (7) and the plotted fit are
on pages 533-534.

The PDF states copyright 1995 by the American Geophysical Union and contains no
open data license. Peritheos therefore does not redistribute the article, page
images, captions, or table typography. The bundled CSV is an independently
created transcription of factual measurements in a new machine-readable
arrangement. Its CC0 dedication is limited to contributors' rights in that
transcription, normalization, column names, and arrangement; it does not
relicense the source article or any third-party material.

## Selection and model boundary

All four final shock states in Table 3 are included. The fit excludes the
reverse- and forward-impact sound-velocity experiments in Tables 1-2, the one
resolved elastic-precursor velocity, the Marsh single-crystal comparison data,
and the separately calculated 300 K hydrostat. The selected states are the
untransformed polycrystalline B1/periclase principal-Hugoniot branch from 14 to
133 GPa. The paper assigns no transition in this range. This record is a
one-dimensional shock path; it is not an equilibrium isotherm or a general
`P(V,T)` model.

## Numerical reconstruction

`scripts/reproduce_duffy_ahrens_1995_mgo_hugoniot.py` performs three transparent
diagnostics on the four `Us-up` pairs:

| Objective | `c0` (km/s) | `s` | Standard errors | Result |
|---|---:|---:|---|---|
| Unweighted vertical OLS | 6.873124 | 1.235534 | residual-scaled 0.02115, 0.01077 | point estimate only |
| Vertical WLS using printed `sigma(Us)` | 6.870296 | 1.238037 | absolute 0.09962, 0.04069 | rounds to source |
| Errors in variables using printed `sigma(Us)` and `sigma(up)` | 6.870440 | 1.237986 | absolute 0.10035, 0.04120 | rounds to source |

The paper says only “least squares fit” and does not publish its residual
definition, parameter covariance, or code. Both uncertainty-aware fits recover
the printed coefficients and their one-sigma errors at stated precision; the
errors-in-variables result is the primary independent reconstruction because it
uses every printed velocity uncertainty. Its covariance is diagnostic and is
not promoted to a source-reported covariance on the published record.

The script also recomputes stress as `rho0 Us up` and shocked density as
`rho0 Us/(Us-up)`. Agreement with Table 3 is limited by the paper's rounded
values, providing an independent transcription check without treating stress or
density as extra regression observations.
