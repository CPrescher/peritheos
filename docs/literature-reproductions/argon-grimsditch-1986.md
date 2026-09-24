# Grimsditch, Loubeyre and Polian (1986): acoustic constraints and a density refit

[Physical Review B 33, 7192-7200](https://doi.org/10.1103/PhysRevB.33.7192),
published 15 May 1986. The complete nine-page publisher PDF was obtained from
[APS Harvest](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.33.7192/fulltext).
Title, all three authors, journal, volume, date, pages and Tables I-II were
visually verified. Its SHA256 is
`3bd0e12705ebafd084ce72b1f2debb7ac47b5be37017f7e26be4c7bffc8c3729`.
The copyrighted PDF is outside git at
`/Users/clemens/Documents/Peritheos/literature/argon/grimsditch-1986/grimsditch-1986.pdf`.

The new record `argon_fcc_grimsditch_1986_density_polynomial` implements the
**published refit of earlier diffraction data**. It is not a new independent
Brillouin P-V determination. `scientific_validation=primary_source_validated`
means the equation and source transcription were checked;
`reproduction.fit_status=not_reproduced` explicitly records that the original
regression has not been recovered. No replacement BM/Vinet fit is introduced.

## Equation, normalization and domain

Equations (6)-(7), pp.7193-7194, give

```
P(rho) = 12.65 - 11.43*rho + 1.5*rho^2 + 0.68*rho^3
B_T(rho) = rho*(-11.43 + 3*rho + 2.04*rho^2)
```

Pressure and bulk modulus are in GPa; mass density is in g/cm³. The coefficient
of the cubic term is **0.68**, checked on the page image (text extraction can
misread it). Coefficient errors, covariance and original regression weights are
not supplied. The fitted input was X-ray density data from references 17-19:
Hazen et al., Carnegie Yearb.79,348 (1980); Zou et al., Yearb.81,392 (1982);
Xu, Mao and Bell, High Temp. High Pressures16,495 (1984). The paper states that
these data extend to 77 GPa, but does not state their precise low-P cutoff.
The paper's own Brillouin solid table spans 1.31-33.58 GPa; these ranges have
different meanings and are stored separately.

`DensityPolynomial3` evaluates the literal polynomial in Rust and Python using
`rho = rho0*V0/V`. The public volume is Å³ per conventional four-atom fcc cell.
With Ar molar mass 39.948 g/mol and Avogadro constant 6.02214076e23 mol^-1,
`rho*V = 265.3408586218433`. The chosen computational anchor is `rho0=2` and
`V0=132.67042931092166 Å³`, where **P=1.23 GPa**, not zero. These are conversion
anchors, not fitted published coefficients. The polynomial has no physical
ambient solid reference. Evaluation and inversion use its positive-compressibility
high-density branch and reject unstable states. Generic reference-state thermal
shifts are unsupported. The model is mathematical beyond the source range;
that is not evidence for physical extrapolation.

The experiment is described as room temperature, without a numerical value.
The record's 298 K follows the library's explicit `assumed_room_temperature`
convention. It is not a measured temperature and no thermal correction is made.
The solid structure is fcc (explicit in the calculation description, p.7197,
and the cubic elastic analysis). Liquid rows are retained as liquid. Wittlinger's
separate hcp record is untouched.

## Observations versus derived quantities

`argon-grimsditch-1986-table1.csv` preserves all **98** printed rows, in source
order: 49 left-column rows then 49 right-column rows. The first 23 are liquid,
the remaining 75 solid. Phase follows the printed groups, not a pressure cutoff:
the 1.36 GPa row is liquid and the subsequent 1.31 GPa row is solid. Repeated
pressures and anisotropic shifts are retained.

Only pressure and Brillouin shift are direct experimental quantities in this
comparison. The product nv is calculated from the shift (514.5 nm laser), the
liquid density is adopted from ref.13, and solid density is taken from the
prior diffraction EOS. Liquid refractive indices come from ref.16; solid indices
are calculated with Eqs.(3)-(4). The effective longitudinal modulus C=rho*v²
is derived. It is **not** a bulk modulus in the solid. The Table I caption's
cross-referenced equation numbers for refractive index are inconsistent with
the body; the body and actual equations resolve this without changing values.

`argon-fcc-grimsditch-1986-table2.csv` retains eight rows, 2-30 GPa, and all
printed errors. C11 is an **upper bound**, C12 a **lower bound**, and C44 is
derived using the adopted bulk modulus. Cstar is the acoustic envelope value;
the prose calls it correct on p.7195 and a lower bound on p.7199. Treat it as
an envelope constraint, not an unqualified independent elastic measurement.
Table III repeats rounded experimental bounds alongside HFD-C theory values;
it is not a second independent experimental dataset and is not duplicated.
No confidence level is invented for the printed errors.

The source explicitly approximates B_S=B_T except in a small region near
melting. The library computes **B_T** from the polynomial and does not silently
apply this approximation to convert solid C into B_T. The paper provides no
complete heat-capacity/expansivity correction for an independent acoustic EOS.
The source estimates B uncertainty from different fits as 15% near 2 GPa,
2% at 4 GPa, and 1% at higher P below 60 GPa; this is separate from its Table II
printed errors. Figure 9 triangles are evaluations of Eq.(6); open triangles
are extrapolations. Neither is a set of independent measured P-V points.
BFW/Ross/HFD-C/HFD-D self-consistent phonon curves are theory comparisons;
this new record does not relabel them as the experimental polynomial.

## Pressure calibration

Equation (2) explicitly uses the nonlinear Mao et al. (1978) ruby scale,
`P=380.8*((lambda/lambda0)^5-1)` GPa, for the acoustic experiments. Raw ruby
wavelengths are absent. Pressure differences across two rubies reached 1.5 GPa
around 35 GPa (p.7196). This scale identification does not resolve the separate
pressure calibration histories of the three earlier X-ray sources. The record
therefore marks calibration `partially_resolved` and recalculation unavailable.

## Reproduction diagnostics

Run `python -m scripts.reproduce_grimsditch_1986_argon` to regenerate
[the report](../data/argon-grimsditch-1986-reproduction.json). It checks the
native EOS against an independent polynomial and analytic derivative, inversion,
and printed Table II bulk moduli without fitting the adopted Table I densities.
Native maximum errors over the 75 solid rows are about 2.2e-14 GPa for pressure,
7.2e-14 GPa for bulk modulus and 6.4e-13 Å³ for inversion. All eight Table II
B values agree within their printed errors.

The polynomial versus the printed **derived** Table I densities has 0.175943 GPa
RMS discrepancy and 1.216051 GPa maximum at source row 98: the printed
33.58 GPa, rho=3.936 g/cm³ gives P=32.363949 GPa. Both are visually verified;
the discrepancy is preserved, not repaired by refitting. The arithmetic checks
of printed nv and C differ by at most 0.01287 km/s and 0.18659 GPa, respectively.
Those checks concern rounded derived columns, not measurement uncertainties.
No independent coefficient refit, measured P-V residual or covariance is claimed.

## Integration requirements

The native eosmat dispatcher supports `DensityPolynomial3` /
`density_polynomial_3`; rebuild the Studio engine to expose it. Display the
record as a published density refit, keeping the 298 K assumption and nonzero
anchor pressure visible in source details. Acoustic and elastic datasets must
use their own axes; do not feed them to the ordinary measured P-V overlay or
fit pipeline. `fit_datasets` is deliberately absent. The dataset links mean
supporting evidence, not EOS fitting input. Keep the source-regression status
`not_reproduced` visible alongside the checked equation.
