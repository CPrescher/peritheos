# Fischer et al. (2011): B1 and B8 FeO thermal EOS

## Source, identity, and temperature correction

The primary source is Fischer et al., *Equation of state and phase diagram of
FeO*, *Earth and Planetary Science Letters* **304**, 496–502 (2011),
[doi:10.1016/j.epsl.2011.02.025](https://doi.org/10.1016/j.epsl.2011.02.025).
The audit used the [author manuscript](https://mineralphysics.uchicago.edu/Papers/FischerEPSL2011preprint.pdf)
and [official Table S1 workbook](https://ars.els-cdn.com/content/image/1-s2.0-S0012821X11000963-mmc1.xls),
with SHA-256 hashes `c69717d039d80421dcdac45ff0207844dd734ad284039b0a573e6a30cb153d24`
and `9a8807dcdb779840b4ad25811309cc7eeb179c42ffd53361eae36b95ebbaf561`.

The Fe0.94O starting powder was mixed with Fe; the authors treat the
equilibrated high-P-T oxide as iron-saturated, approximately stoichiometric
FeO. Surface temperatures were corrected downward by 3% before publication.
Table S1 contains those corrected temperatures. The refits use them directly:
they neither treat hot states as 300 K nor apply a second correction.

## Equation and source parameters

Equation (1) is a 300 K BM3 reference isotherm plus Mie–Grüneisen–Debye
thermal pressure,

\[
P(V,T)=P_{300}(V)+\gamma(V)[E(\theta_D(V),T)-E(\theta_D(V),300)]/V,
\]

where `gamma=gamma0*(V/V0)^q` and
`thetaD=theta0*exp[(gamma0/q)*(1-(V/V0)^q)]`. This is Peritheos'
`integrated_gruneisen` law. The paper explicitly omits anharmonic and
electronic thermal-pressure terms; `n=2` is the FeO atom count.

| Phase | Fixed | Fitted |
|---|---|---|
| B1 | `V0=12.256 cm3/mol`, `theta0=417 K`, `q=0.5`, `Tr=300 K`, `n=2` | `K0=149.4+/-1.0 GPa`, `K0'=3.60+/-0.04`, `gamma0=1.41+/-0.05` |
| B8 | `K0'=4`, `theta0=417 K`, `q=1`, `Tr=300 K`, `n=2` | `V0=11.997+/-0.018 cm3/mol`, `K0=137.8+/-0.9 GPa`, `gamma0=1.73+/-0.12` |

The molar volumes are stored as conventional-cell volumes with `Z=4` for B1
and `Z=2` for B8. Table 1's B1 `K0'` error is 0.04; an abstract rendering that
drops the zero does not override the parameter table.

## Selection, exclusions, pressure basis, and weights

B1 uses the 42 Table S1 rows with numerical B1 volumes. B8 uses all 21 rows
with numerical B8 volumes, including coexistence states. One-peak B8
detections and the room-temperature rhombohedral B1 detection have no
reportable phase volume and are excluded. Three two-peak B8 volumes have blank
workbook errors; their raw errors stay blank, while a separate audit-weight
column applies Figure 3's explicit `+/-0.1 cm3/mol` fallback. Kondo et al.
(2004) and Fei and Mao (1994) are comparison points, not B8 fit inputs.

Current-study pressure comes from simultaneous hcp-Fe volume and temperature
using Dewaele et al. (2006); NaCl is secondary. The Fe observations are
bundled, and the exact Dewaele et al. (2006) thermal/electronic reference EOS is
now executable as `iron_dewaele_2006_vinet_thermal`. Current-study pressures can
therefore be recalculated directly.

The source's global B1 fit additionally uses pressure-recalculated Campbell,
Ozawa, and Seagle rows; B8 additionally uses Ozawa. The source recovery is now:

| Study | Available row-level source | Recovered scope | Pressure status |
|---|---|---:|---|
| Seagle et al. (2008) | Official supplementary Table 1 | 81 B1 volume-temperature rows: 65 hcp-Fe, 14 fcc-Fe, 2 without Fe volume | 65 ready with Dewaele; 14 require Campbell fcc-Fe; 2 retain the source's 50 +/- 4 GPa series pressure |
| Ozawa et al. (2010) | Printed article Table 2; no supplement exists | 12 B1 and 8 B8 P-V-T rows | Reported pressures and hcp-Fe volumes bundled |
| Campbell et al. (2009) | Official supplementary Table S2 | 25 B1 P-V-T rows: 10 hcp-Fe and 15 fcc-Fe | hcp rows ready; fcc-Fe calibration still to register |

Thus Ozawa's printed table is sufficient at the publication's own precision,
and neither Seagle nor Campbell is now a missing-data source. What remains
undocumented by Fischer et al. is the exact combined row selection, regression weights,
residual variable, covariance, and optimizer. The available refits therefore
remain explicitly conditional rather than claiming exact global-fit parity.

## Reproducible refits

The B1 current-study diagnostic uses the published pressure, volume, and
temperature uncertainties. The combined B8 reproduction is unweighted because
Fischer does not publish numerical regression weights and Ozawa does not print
a complete set of row-wise uncertainties. These are auditable choices, not
claims about the undocumented source objective.

| Phase | Rows | Published free coefficients | Conditional refit | Result |
|---|---:|---|---|---|
| B1 | 42 | `149.4, 3.60, 1.41` | `144.508, 3.71562, 1.63020` | `parity` |
| B8 | 29 (21 Fischer + 8 Ozawa) | `39.843, 137.8, 1.73` | `39.8745, 137.863, 1.66616` | `parity` |

Weighting sensitivity confirms the remaining B1 limitation. B1 gives
`162.018, 3.43046, 0.97854` unweighted and `147.054, 3.69873, 1.33440` with
pressure-error weights for the 42-row current-study subset. For B8, the
current-study-only unweighted result is `37.9814, 168.125, 2.06333`; adding all
eight printed Ozawa B8 rows recovers the published coefficients within their
reported errors. The combined fit is unweighted because Fischer does not
publish the numerical weights and Ozawa does not provide complete row-wise
uncertainties. This is numerical parameter parity, not exact procedural
identity.

As a broader B1 diagnostic, all 160 recoverable observations can be retained
when the two Fe-missing Seagle rows use their source-series assignment of
`50 +/- 4 GPa`. The unweighted result is `151.360, 3.57859, 1.44606`, compared
with Fischer's published `149.4 +/- 1.0, 3.60 +/- 0.04, 1.41 +/- 0.05`.
`K0_prime` and `gamma0` lie within the published one-sigma intervals, while
`K0` lies within two sigma (1.96 sigma). The full-source B1 diagnostic is
therefore classified as numerically `similar`; the two series-pressure rows
have negligible leverage on the result.

The source parameterizations remain authoritative; the audit makes both
thermal models executable without fabricating missing rows or temperatures.
