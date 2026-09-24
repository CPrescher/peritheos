# Argon: fcc additions and the Wittlinger hcp check

The `argon_fcc` material adds the thermal EOS of
[Dewaele et al. (2021)](https://doi.org/10.1038/s41598-021-93995-y)
and the experimental 300 K Vinet EOS of
[Ono (2020)](https://doi.org/10.1038/s41598-020-58252-8).
Wittlinger remains a separate historical **hcp** material.

Run `python -m scripts.reproduce_argon` to regenerate
[the numerical report](../data/argon-reproduction.json).
The independent equations use conventional-cell volumes and 96-point
Gauss-Legendre quadrature, rather than calling the native EOS to generate
expected results. Native pressures agree within 2e-11 GPa over the checked
room-temperature and cryogenic observations. Public API and native tests
also cover inversion, pressure increments, and material interchange.

## Dewaele (2021)

Table 2 and equations (2)-(5) give a **0 K** Vinet curve with atomic
V0=38.0 A^3, K0=2.65 GPa (fixed), K0'=7.423, theta0=93.3 K,
and gamma(V)=0.5+2.2(V/V0). The fcc conventional cell has four atoms:
the stored V0 is 152 A^3, and the thermal molar normalization uses one
atom per Ar formula unit. No coefficient uncertainties or covariance are printed.

The existing `Dewaele2006` thermal model represents this equation exactly
with gamma0=2.7, gamma_inf=0.5, beta=1, and zero anharmonic/electronic
coefficients. Its new `thermal_pressure_reference="absolute_zero"` option
adds the Debye pressure without subtracting a room-temperature baseline.
Zero-point pressure is omitted as in the paper. Tr=296 K specifies the
anchor for **increments**, so the increment at 296 K is zero while the
absolute thermal contribution at 296 K is positive. State temperatures must
remain positive. The existing reference-temperature mode is unchanged.

All **288 supplementary rows** are transcribed: 88/34/35/62/69 from runs
1/2/3/4/5. Original Au and ArNe2 lattice readings, missing entries and
repeated rows are retained. `ArNe12_017` has 4.8579998 printed under aAu;
it remains in that column with a note, without silently moving it to ArNe2.
It has no pressure and cannot enter a fit. No coordinate errors are supplied.

The published selection gives **95 usable rows**: runs 2-3 and pure-Ar
runs 1 and 4 at P<=5 GPa. The printed curve's RMS pressure residual is
0.341 GPa. Equal-pressure-weight refits with K0 and thermal coefficients
fixed give V0=152.384 A^3 and K0'=7.41438 at 296 K, versus the source
152 and 7.423. At 300 K they give 152.315 and 7.41826. Both temperatures
are checked because the fitting paragraph says 300 K while the Methods
specify 296 K. The source objective and weights are unspecified; the
refits are diagnostic and do not replace the published coefficients.

The 62 run-5 rows with pressure provide a separate thermal check (RMS
0.473 GPa); the supplementary temperatures include 5.5 and 300 K even
though the headline experimental range is 10-296 K. Low-temperature
coverage is below about 20 GPa, and the 114 GPa pure-Ar measurements are
comparison data rather than the quasi-hydrostatic fit. Marginal ranges
are not a rectangular phase-stability guarantee.

Benchmarks include low, intermediate and maximum fitted pressure. The
0.15 GPa floor and 2.5% pressure tolerance are empirical curve-reproduction
tolerances covering source scatter, not invented measurement uncertainties.
Au pressures cite Takemura-Dewaele (2008). Run 5 uses Dewaele (2008) ruby
with the Syassen (2008) temperature correction; raw ruby readings are absent.
The thermal model is classified hybrid because its experimental constraints
are combined with an assumed infinite-compression Gruneisen limit.

## Ono (2020)

Supplementary Table 2 and equation (4) give the conventional-cell
Vinet coefficients V0=184.5(38.5) A^3, K0=1.07(1.33) GPa and
K0'=8.02(0.95), all fitted. This extrapolated V0 is not an ambient
solid volume. All **19 Table 1 observations** and their one-sigma P,V
errors are bundled, spanning 4.44-136.7 GPa at 300 K.

The published pressure RMS is 0.847 GPa. Equal-pressure weights recover
184.330 A^3, 1.07506 GPa and 8.01417, close to the printed coefficients.
A sensitivity fit with fixed propagated P,V errors gives 211.745 A^3,
0.45844 GPa and 8.66468. This demonstrates the strong weighting and
extrapolation sensitivity; covariance and original regression settings
are unavailable. The selected source checkpoints reproduce within 2 GPa,
a tolerance based on table scatter and the large published parameter
errors, not a claim that every observation lies within its coordinate sigma.

The gold calibration is Dorogokupets-Dewaele (2007); no Au observations
are printed with Table 1. This record excludes the separate AIMD thermal
extension: its quadratic coefficient is mislabeled in Supplementary Table 2
and the AIMD state table is absent. The room-temperature fit is independent
of that ambiguity.

## Wittlinger (1997): published fit not reproduced

**NOT REPRODUCED.** The available digitization does not establish a
source-faithful reproduction of the published fit. The card carries an explicit
`reproduction_status: not_reproduced` flag. The refit ledger classifies it as
`not_refittable`: original observations, regression settings and covariance
are unavailable, so the diagnostic fits below do not constitute an independent
reproduction. It no longer counts as a reproduced paper. This is a limitation
of the available evidence, not proof that the published EOS is false.

The earlier primary-source audit transcribed BM2, V0=78(3) A^3 and
K0=6.5(1.3) GPa, with K0'=4 implicit, for hcp Ar (two atoms per cell),
over **1.2-8.5 GPa**. These coefficients remain unchanged. A fresh full-text
attempt on 2026-09-24 reached institutional-access requirements; this check
uses the previously audited coefficients and the stored Figure 3 digitization,
not a newly verified reading of the entire paper.

On the nine digitized experimental points the source curve has **0.772 GPa
RMS** and **1.308 GPa maximum** pressure residual. An equal-pressure refit
gives V0=74.612 A^3 and K0=7.09488 GPa with 0.440 GPa RMS. An
errors-in-variables fit, including plotted and digitization uncertainties,
gives about 72.782 A^3 and 8.0603 GPa with **0.451 GPa RMS at the original
observations**.

The former ledger's approximately 0.013 GPa refit RMS was evaluated at
**adjusted coordinates**. It is not comparable to the published curve's
0.772 GPa residual at observed coordinates. The dedicated ledger entry now
reports the original-coordinate metric and retains adjusted-coordinate RMS
separately. This corrects the misleading apparent 60-fold improvement.

The plotted volume errors are large (roughly 5%); a simple propagated-error
diagnostic puts the published reduced chi-square near 0.48. Thus the plot
does not establish a transcription error or statistically reject the source
fit. The relative volumes also share its fitted V0 normalization, so treating
those errors as independent cannot establish an accurate replacement EOS.
No diagnostic refit is promoted to the executable catalog.

At the same **atomic** volume where Dewaele's fcc curve gives 5 GPa,
Wittlinger's hcp curve gives about 6.50 GPa. At 50 and 100 GPa fcc states,
extrapolated Wittlinger values are about 41.5 and 72.2 GPa. The latter
comparisons are far outside Wittlinger's fitted range and involve different
phases; they show why this record is unsuitable as a general argon pressure
standard, not a direct proof that one phase's measurements are wrong.

The card now recommends historical hcp use only, with no high-pressure
extrapolation. For fcc argon, use the separate Dewaele thermal record.

## Structure and data provenance

The diffraction reference cell is a=5.0868(20) A at 1.28(6) GPa and
293 K from [Finger et al. (1981), page 892, Table I](https://doi.org/10.1063/1.92597),
which identifies Fm-3m and Z=4. The monatomic fcc basis is fully occupied
4a (0,0,0). This measured reference lattice intentionally differs from
both extrapolated EOS reference volumes. It supplies a diffraction-ready
structure without claiming an ambient solid cell.

The two official Scientific Reports supplements are CC BY 4.0. CSV
transcriptions retain source numerical precision and include source-PDF
and CSV SHA-256 hashes. Missing values and uncertainties are not fabricated.
