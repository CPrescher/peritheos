# Anderson and Swenson (1975): five argon isotherms

M. S. Anderson and C. A. Swenson, *Experimental equations of state for the
rare gas solids*, Journal of Physics and Chemistry of Solids **36**(3),
145–162 (1975), DOI
[10.1016/0022-3697(75)90004-9](https://doi.org/10.1016/0022-3697%2875%2990004-9).

## Outcome and source identity

The user-supplied primary PDF resolves the earlier source-access blocker.
Equation (2) and all five Table 1 argon isotherms are implemented directly in
Python and Rust, with EOSMAT loading and serialization support. The status is
**equation_reproduced**. An independent fit of the original observations is
**not_reproduced**, for the normalization and temperature reasons below.
Wittlinger's separate hcp entry is unchanged.

The 18-page journal PDF was checked against its title, authors, volume and
pages. Its stable local path is
`/Users/clemens/Documents/Peritheos Sources/argon/anderson-swenson-1975/anderson-swenson-1975.pdf`;
SHA256 `904c5a9a67a29c48c5b93027b3f8a002db43829683c98cd71e61d9305007fddf`.
The user supplied it as `1-s2.0-0022369775900049-main.pdf`; its actual download
URL is unknown. No copyrighted PDF is committed.

The paper points to Anderson's 1972 thesis, *Experimental equations of state
for solid neon, argon, krypton, and xenon*, IS-T-537, for the observations.
The [OSTI thesis PDF](https://www.osti.gov/servlets/purl/4658659) is stored beside
the paper as `anderson-1972-thesis-IS-T-537.pdf`;
SHA256 `12fb53f4855419ca03b76900240a60bab4ada9523d128c876c07df260e5040f8`.

## Equation, units and reference state

The original equation is

\[
P_{\rm kbar}(V_m)=\sum_{n=3,5,7,9} A_n V_m^{-n},
\]

with molar volume in cm³/mol. Table 1's printed coefficient multipliers are
undone in `argon-anderson-swenson-1975-table1.csv`. The `OddInversePower`
model (`odd_inverse_power`) evaluates the algebraically identical form

\[
P_{\rm GPa}(V)=\sum_n C_n(V_0/V)^n,
\qquad K_T(V)=\sum_n n C_n(V_0/V)^n.
\]

Let \(v_*\) be the zero-pressure root of the rounded published polynomial.
The catalog uses \(V_0=v_*4\times10^{24}/N_A\) in Å³ per conventional
four-atom fcc cell and \(C_n=A_n/(10v_*^n)\) in GPa, with exact
\(N_A=6.02214076\times10^{23}\) mol⁻¹. This normalization preserves the
published curve: no pressure baseline is subtracted and no fit is performed.
The rounded printed zero-pressure volumes remain in the summary dataset.
`V0` is fixed because simultaneous fitting of the normalization and all
coefficients would introduce a redundant degree of freedom.

| EOS record suffix | Temperature (K) | Published volume RMS (cm³/mol) |
|---|---:|---:|
| `4p2k` | 4.2 | 0.010 |
| `20p0k` | 20 | 0.009 |
| `40p0k` | 40 | 0.011 |
| `60p0k` | 60 | 0.020 |
| `77p0k` | 77 | 0.012 |

Every identifier starts with `argon_fcc_anderson_swenson_1975_`.
Table 1 prints 4 K; the experiment used 4.2 K. These are discrete isotherms,
not a continuous thermal EOS. Do not interpolate temperature or extrapolate
above 20 kbar (2 GPa). The low-pressure bounds depend on holder and
measurement temperature. Catalog bounds use usable thesis observations near
the nominal isotherm; they are not an independently recovered final fit mask.
The zero-pressure state is an extrapolated reference. The sample normalization
at 4.2 K and 4 kbar is not the EOS reference pressure.

The fcc association reflects the low-pressure solid, not a diffraction phase
determination by this piston-displacement experiment. Absolute molar volumes
are tied to Peterson et al. (1966). The overall accuracy of approximately
0.001 of the zero-temperature, zero-pressure volume is not a per-point sigma
or coefficient uncertainty. No coefficient covariance is available.

## Observations and unresolved fit reconstruction

The journal states that its plotted markers for this work were calculated
from Table 1. They are not imported as independent measurements. The table
CSV is explicitly a `fit_parameter_summary`, not a P–V dataset.

Thesis Appendix A, printed pp. 93–97 (PDF pp. 101–105), supplies 394 argon
table cells: 379 usable observations and 15 printed zero sentinels marked
unreliable by the source. All are retained in
`argon-anderson-1972-appendix-a.csv`, with holder diameter, sample, run,
actual temperature and original reported molar volume. Unreliable entries
have `usable=0` and a blank converted cell volume. No common-temperature or
sample-length adjustment is silently applied. Appendix A p. 89 explains
that the tabulated values already incorporate apparatus corrections but use
unadjusted sample lengths.

Appendix B describes common-temperature isobaric reduction and additive
sample-length corrections before fitting. Table 11 assigns +0.253 cm³/mol to
0.250-inch sample I and +0.449 cm³/mol to sample II. Applying these literal
assignments does not reproduce the journal RMS. The diagnostic intentionally
preserves these assignments; it does not silently exchange samples or infer
new corrections.

For comparisons only, observations within 1.1 K of each nominal isotherm are
selected. The 67.9 K run is excluded, not interpolated. With literal Table 11
corrections, volume RMS values remain about 0.16–0.18 cm³/mol, considerably
larger than the published 0.009–0.020 cm³/mol. Larger-holder observations alone
give about 0.009–0.029 cm³/mol. These diagnostics expose the unresolved
small-holder normalization and common-temperature reduction; they are not a
reproduction of the source least-squares fit. Records link the raw dataset as
`comparison_datasets`, never as exact `fit_datasets`.

## Pressure calibration

The hydraulic ram force was set with a deadweight pressure balance and
converted using piston area. Increasing and decreasing force readings at
the same displacement were averaged to account for friction. Indium samples
were used to calibrate holder deformation. The method is documented, but
individual calibration readings needed for recalculation are unavailable.
There is no ruby pressure scale. Temperature was stabilized to approximately
±0.2 K with quoted accuracy ±0.5 K.

## Verification and integration

`scripts/reproduce_anderson_swenson_1975.py` produces the committed JSON
report in `docs/data`. It independently evaluates the original molar
polynomial and compares pressure, inversion and source bulk moduli with the
catalog models. Pressure differences are below 4e-14 GPa and inversion
errors below 1e-9 cm³/mol. Computed zero-pressure bulk moduli agree with
Table 1 within 0.051 kbar. The separately rounded reduced Eq. (9) differs
from Table 1 by less than 0.0006 GPa; it is diagnostic, not another EOS.
Python and Rust tests also check the analytic bulk modulus by finite
differences. Rust evaluates the family natively; Python uses an independent
NumPy implementation and generic inversion, without a new PyO3 wrapper.

Focused checks cover the six new Python tests, existing argon and EOSMAT
tests, Rust library and EOSMAT tests, and the new independent Rust source
coefficient test. See the handoff JSON for exact results.

Studio should expose the five discrete isotherms, preserve the raw-data
qualification and actual observation temperatures, filter `usable=0`, and
keep fit summaries distinct from measured points. The integration task owns
Studio and Zotero changes, including attachment of the verified primary PDF.
