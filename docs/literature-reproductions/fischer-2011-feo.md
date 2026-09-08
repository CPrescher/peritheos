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
bundled, but that exact thermal/electronic reference EOS is not executable, so
recalculation is `reference_eos_not_bundled` rather than unresolved.

The source's global B1 fit additionally uses pressure-recalculated Campbell,
Ozawa, and Seagle rows; B8 additionally uses Ozawa. Fischer et al. do not print
those final numerical rows, exact regression weights, residual variable,
covariance, or optimizer. Digitizing plotted temperature bins would invent
missing temperatures, so the available refits are explicitly conditional.

## Conditional refits

The ledger uses published errors in variables. B1 uses pressure, volume, and
temperature uncertainties. B8 uses pressure and volume uncertainty because
its 300 K row reports zero temperature uncertainty; the actual 300 K value is
retained. This is an auditable choice, not a claim about undocumented weights.

| Phase | Rows | Published free coefficients | Conditional refit | Result |
|---|---:|---|---|---|
| B1 | 42 | `149.4, 3.60, 1.41` | `144.508, 3.71562, 1.63020` | `parity` |
| B8 | 21 | `39.843, 137.8, 1.73` | `37.2815, 183.822, 1.99615` | `parity_not_achieved` |

Weighting sensitivity confirms the limitation. B1 gives
`162.018, 3.43046, 0.97854` unweighted and `147.054, 3.69873, 1.33440` with
pressure-error weights. B8 gives `37.9814, 168.125, 2.06333` unweighted and
`39.4343, 141.131, 2.33926` with pressure-error weights. The B8 EIV covariance
is nearly singular because its ambient `V0` is extrapolated from only
131–156 GPa. Huge formal errors are not accepted as numerical parity.

The source parameterizations remain authoritative; the audit makes both
thermal models executable without fabricating missing rows or temperatures.
