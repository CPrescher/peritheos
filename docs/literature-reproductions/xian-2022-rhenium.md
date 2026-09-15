# Xian et al. (2022): rhenium pressure scale

Audit date: 2026-09-15. The final publisher PDF and publisher snapshot were
read from the user's Zotero library. Equations and parameter/observation
tables were checked on rendered PDF pages. No supplementary or correction
links were present in the saved publisher page.

An additional public correction search did not locate an erratum.

[AIP Advances 12, 055313](https://doi.org/10.1063/5.0089292).
Zotero item `ZNQIGGGF`, PDF attachment `EPQGNQ83`.

**Disposition: investigated, deferred from the executable catalog.** The cold
Vinet coefficients and literal temperature polynomials are retained in
[`xian-2022-rhenium-source-audit.json`](../data/xian-2022-rhenium-source-audit.json).

The source gives a zero-kelvin Vinet curve with specific volume
0.047893 +/- 0.00031 cm3/g, K0=319.01 +/- 0.86 GPa, and K0'=4.13 +/- 0.11.
Errors are described as Monte Carlo results, without confidence convention or
covariance. These are cold-curve coefficients, not a room-temperature EOS.

Two independent routes remain unresolved:

1. **Full quasi-Debye/electronic model:** Equations (13)-(16) require linear
   longitudinal and shear velocity fits versus density. Their four coefficients
   are not published. Regressing upstream Ikuta observations anew would be a
   separately qualified reconstruction. Additional printed notation issues
   include `ln(1-exp(chi))` in Equation (6) for positive `chi=Theta/T`, and `h`
   instead of the reduced Planck constant in the angular-wavevector convention
   of Equation (7). The electronic heat-capacity coefficient's printed units
   also disagree with `Cv=beta*T`. The existing sound-velocity model implements
   a physical interpretation of this mechanism, but does not supply the missing
   source-specific coefficients or the electronic contribution for Re.
2. **Pressure surrogate, Equations (37)-(40):** The final PDF visibly prints
   `VT = 0.04789 - 3.635e-7*T - 1.651e-11*T^2` cm3/g. Thus its zero-pressure
   volume decreases from 0.04777946 at 300 K to 0.04655774 at 3200 K, with
   negative thermal expansivity throughout. At the reported cold reference
   volume the literal surrogate predicts -0.743 GPa at 300 K and -7.526 GPa
   at 3200 K. This is inconsistent with the intended thermal-pressure scale.
   No sign correction is assumed without author evidence. The prose calls the
   fits cubic, although the printed expressions stop at quadratic order.

The title specifies 130 GPa, while the abstract/body say 140 GPa, both with
3200 K; the metal comparison extrapolates to 500 GPa. The Au comparison reports
only K0=133.96 GPa and K0'=6.05, omitting a complete cold-volume and thermal
parameter set, so it also yields no executable Au record.

In particular, Xian's Au numbers should not be substituted into Anderson's
300 K BM3 expression. Section III B obtains Au through a new reduction of
shock data using its own thermodynamic model; Section II defines its fitted
`B0` and `B'` for the cold Vinet component. The text introduces the resulting
room-temperature scale but does not separately identify a complete 300 K Au
parameter set. Its 133.96 GPa value is substantially below Anderson's 166.65
GPa, but the missing Au inputs prevent checking how the full model translates
that fit into room-temperature pressure or bulk modulus. A reference-temperature
distinction alone therefore does not establish that this numerical difference
is physically correct. The printed unit for B' in that paragraph is also
misleading: a pressure derivative of bulk modulus is dimensionless.

The paper contains no new observation table or source code and declares data
sharing inapplicable. It cites Marsh shock data and Ikuta acoustic measurements.
A source-author correction and numerical velocity fits are needed for an exact
Re thermal record. Neither sign-repaired polynomials nor guessed reference
volumes are included in the executable bundle.

## Reproduction and verification

Run `python scripts/reproduce_zha_2004_xian_2022_rhenium.py` to regenerate
the [numerical reproduction report](../data/zha-2004-xian-2022-rhenium-reproduction.json).
Its `xian_literal_equations_37_40` section evaluates the printed temperature
polynomials and records their negative expansivity and pressure predictions.
The literal coefficients and the decision to withhold an executable thermal
record are covered by `tests/test_zha_xian_rhenium.py`.
