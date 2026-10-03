# Miozzi et al. (2020): hcp iron equations of state

Miozzi, F., Matas, J., Guignot, N., Badro, J., Siebert, J., and Fiquet, G.
(2020), *A New Reference for the Thermal Equation of State of Iron*,
Minerals **10**, 100, [doi:10.3390/min10020100](https://doi.org/10.3390/min10020100).

Three published current-study parameterizations are retained in `iron.eosmat`, all with
`scientific_validation.status: not_reproduced`, `default: false`, and
`catalog_access: explicit_selection`. These preserve source coefficients for
inspection; they are not validated pressure standards or core predictions.
The preferred thermal record uses **n=2 with volume per mole Fe**, as used in
the paper according to personal communication with Miozzi et al., relayed by
the user and recorded on 2026-10-03. Its printed coefficients are preserved.
Separate, nondefault **n=1** Peritheos refits are available on the reconstructed
Speziale (2001) variable-q Debye and Tange (2009) Fit3-Vinet MgO scales.
The existing iron default is retained.

## Sources and observations

The user supplied the publication-of-record `minerals-10-00100-v2.pdf` and the
official `minerals-10-00100-s001.zip`. The ZIP contains two numerical PDFs and
a supplementary discussion/figure PDF. All four PDFs were read, and parameter,
equation, and numerical-table pages were visually checked on 2026-10-02.
The original input names, byte lengths, and SHA-256 hashes are recorded in
`peritheos/data/datasets/miozzi_2020_sources/manifest.json`. Publisher files
are checksummed rather than redistributed. The article states CC BY 4.0;
normalized observations retain Miozzi et al. attribution. The publisher landing
page returned HTTP 429 during the live check; no claim is made that every
correction channel was checked.

All **131 observations** are bundled:

- `iron-miozzi-2020-he-pvt.csv`: 15 rows at the printed 298 K, 10.7–30.8 GPa.
- `iron-miozzi-2020-mgo-pvt.csv`: 116 rows at 300–3399 K,
  41.5691–127.9736 GPa, including the paired MgO cell volumes and errors.

The original pressure, temperature, volume and coordinate-error precision is
retained. Source page and observation order are explicit. Blank separator lines
with an isolated pressure-error zero are not measurements. Blanks remain missing,
and printed zeros remain zero. MgO row 1 lacks a pressure error, row 27 lacks an
Fe volume error, and row 76 lacks a temperature error. Confidence levels and
coordinate/parameter covariance are not supplied. The cold selection contains
**36** observations at 298/300 K; the other **95** are heated measurements.
Actual marginal ranges are recorded, not a rectangular stability guarantee.
The paper's rounded ranges (140 GPa/3500 K in the abstract, 45–135 GPa and
1300–3300 K in Section 3.3) are distinct from the supplied rows.

## Published coefficients and source inventory

| Record | V0 (A³/cell) | K0 (GPa) | K0′ | Thermal parameters | Source |
|---|---:|---:|---:|---|---|
| `iron_miozzi_2020_bm3` | 22.80(2) | 129(6) | 6.2(2) | none | Table 1, current study |
| `iron_miozzi_2020_vinet` | 22.81 | 125(5) | 6.5(2) | none | Table 1, current study |
| `iron_miozzi_2020_bm3_mgd` | 22.81 | 129(1) | 6.24(4) | theta0=420 K fixed; gamma0=1.11(1); q=0.3(3); n=2 | Section 3.3 coefficients; n from personal communication |

Parentheses retain the authors' error notation without assigning a confidence
level. The Vinet reference-volume error is missing, not zero or implicitly fixed.
The thermal reference volume retains the independent reported molar-volume
error: 0.02 cm³/mol converts to 0.06642156 A³ for a two-atom cell. Section 3.3
prints 6.87(2) cm³/mol and 22.81 A³; converting the central molar value gives
22.81580678 A³. Section 3.3 says V0 is fixed to the room-temperature result,
where Table 1 and the Figure 4 inset print 22.80(2) A³. The thermal record follows the explicitly
printed thermal cell value, with the discrepancy documented; it does not borrow
the tighter Table 1 error. Personal communication clarifies that V0 was free
in the final fitting stage, despite the article's fixed-V0 wording; the
record retains that discrepancy in its parameter provenance.

The source also prints two combined current-study + Dewaele (2006) cold fits:
(22.75(1), 134(5), 6.1(2)) with pressure errors and
(22.57(1), 153(3), 5.3(1)) without them. A combined thermal alternative uses
V0=22.57 A³ (6.79(1) cm³/mol), K0=153(1) GPa, K0′=5.39(4),
theta0=420 K, gamma0=1.09(1), q=0.26(4).
These nonpreferred combined-study alternatives are documented here rather than
represented as fits to the two current-study tables alone; their complete joint
row selection and original weighting have not been established. Other Table 1
columns quote Fei (2016), Dewaele (2006), and Uchida (2001), not new Miozzi fits.

## Equation and normalization

Equation (3) is standard BM3; Table 1 separately identifies the Vinet alternative.
The 300 K pressure receives a vibrational thermal increment with
`gamma(V)=gamma0*(V/V0)^q` and
`theta(V)=theta0*exp((gamma0-gamma(V))/q)` (Equations 4, 6, 8).
The implementation uses `debye_temperature_law: integrated_gruneisen`, with
zero thermal increment at 300 K. The conventional P63/mmc cell has Z=2, so
molar cm³/mol of Fe is `V_cell*N_A/(2e24)`. The published record now uses
the author-reported **n=2** at that per-Fe volume, based on personal
communication recorded on 2026-10-03. The independent refits use the
physically consistent **n=1** per-Fe basis. A consistent whole-cell basis
would instead use n=2 with doubled molar volumes. These settings are distinct;
see [the normalization explanation](../units.md#atom-count-and-molar-volume-basis).

The supplied PDF has apparent typesetting defects: Equation (5) prints gamma
rather than volume in the denominator and a comma between the energies;
Equation (7) omits the cube on T/theta and uses theta0/T as the upper integration
limit despite the volume-dependent theta in Equation (8). The represented
candidate uses the standard physical MGD expression, `gamma/V_molar` times
the difference in Debye vibrational energies, integrated to theta(V)/T.
Zero-point energy cancels at fixed volume. This interpretation is explicit,
and its failure to reproduce the complete reported fit remains unresolved.
The published gamma is retained; the atom-count setting is attributed to
personal communication rather than introduced as an empirical correction.

## Pressure calibration

Section 2.5 and Supplementary Figure S1 select the **Speziale (2001) thermal MgO
scale** after comparison with Tange (2009). The bundled Speziale catalog record
is only the isothermal BM3 fit; linking that record as a full thermal standard
would misrepresent the source. Paired MgO observations are retained, but pressure
recalculation remains `reference_eos_not_bundled`.

The helium series used the Mao (1986) hydrostatic ruby scale (exponent 7.665)
and a Datchi (1997) SrB4O7:Sm2+ gauge. Their individual readings and combination
are absent. Calibration metadata is therefore `partially_resolved`; per-record
pressure reductions are explicitly `unresolved`. The original reproduction and
code audit use the printed pressure column unchanged. The separate conditional
reconstruction at the end of this note retains both supplied and recalculated
pressures. The discrepancy does not establish that a
particular pressure scale was used incorrectly by the authors.

## Independent reproduction and conditional refits

`scripts/reproduce_miozzi_2020_iron.py` independently implements BM3, Vinet,
and the standard Debye integral using Gauss–Legendre quadrature. This does not
call the Peritheos evaluator. Tests compare the native evaluator to that
independent calculation, verify published coefficients, original data holes,
checksums and inversion, and confirm that strict material loading rejects the
unvalidated records.

| Printed parameterization | Observations | RMS pressure error (GPa) | Residual range (GPa) |
|---|---:|---:|---:|
| BM3 at 298/300 K | 36 | 3.7536 | −7.7591 to +0.5179 |
| Vinet at 298/300 K | 36 | 4.1764 | −8.3823 to +0.2167 |
| Preferred BM3–MGD, author-reported n=2 | 131 | 2.4989 | −6.9970 to +9.4776 |

The thermal residuals still exceed Section 3.3's reported −3 to +3 GPa interval.
The earlier n=1 interpretation gave RMS 6.9826 GPa; the n=2 author setting
substantially reduces that discrepancy without changing the printed coefficients.
The helium series alone lies within approximately 0.52 GPa of the printed BM3
curve, while larger discrepancies occur in the supplied MgO-series rows.
The 0.01 A³ reference-volume rounding difference is too small to explain them.
These are numerical source-reproduction limits, not evaluator disagreement.

Both ordinary pressure residuals and a conditional effective-coordinate-error
objective are fitted. Effective errors combine the squared printed pressure
error, `(dP/dV * dV)^2`, and `(dP/dT * dT)^2`, with derivatives frozen at the
published coefficients. Missing errors contribute no term; one cold observation
with no positive effective error is excluded only from that weighted diagnostic.
Printed error widths are treated as standard deviations solely for this
conditional calculation. They are not claimed to be source-defined sigmas.
Every varied coefficient has conditional standard errors and covariance; thermal
V0, theta0, Tr, and n remain fixed. The original fitting files are not supplied;
the authors' exact weight selection, iteration settings, and treatment of
missing/zero errors have not been established. The thermal diagnostic fits all
four free coefficients simultaneously with V0 fixed rather than replaying the
reported stepwise procedure. It is distinct from the registered five-parameter
Speziale refit below, whose final V0 is free.

EosFit's documented effective-variance method updates weights every
least-squares cycle ([Angel et al., 2014, p. 418](https://www.rossangel.com/Download/2014_Angel_etal_EosFit7.pdf)).
The frozen-weight diagnostic above therefore is not an exact EosFit replay.
A sensitivity check on 2026-10-02 repeatedly recomputed the effective errors
from the fitted EOS and refitted with those weights until the largest change
in scaled parameters was below 1e-8. Keeping the same observations, error
handling, equations and fixed parameters, this gave BM3
(V0=22.47454 A³, K0=163.0332 GPa, K0′=5.50772) and thermal BM3–MGD
(K0=132.94294 GPa, K0′=6.34459, gamma0=2.37984, q=1.25700).
Thus updating the weights alone does not recover the published coefficients.
This check is still conditional and does not reproduce the authors' original
EosFit run. Parameter differences are not evidence that the published fit is
wrong; evaluating the printed coefficients against the supplied rows also
requires confirmation of the input selection, pressure reduction and thermal
normalization.

For the 35-row effective-error cold diagnostic, BM3 gives
V0=22.47155 A³, K0=163.3571 GPa, K0′=5.49869; Vinet gives
V0=22.51847 A³, K0=156.1931 GPa, K0′=5.94292.
The conditional 131-row thermal diagnostic gives K0=132.9967 GPa,
K0′=6.33509, gamma0=2.21755, q=1.11705 at fixed V0=22.81 A³.
These do not recover the published parameter sets. Full conditional fit errors,
covariance, ordinary-fit alternatives and residual metrics are in
[`miozzi-2020-iron-reproduction.json`](../data/miozzi-2020-iron-reproduction.json).
They are diagnostics, not replacement published EOS parameters.

## EosFit code audit (2026-10-02)

The public CrysFML library used by EosFit was inspected at historical revision
`5bb3630d684387230065fa7c5933af2d81d19645` (2020-01-24) and at current
CrysFML2008 revision `f43802268348b0171fee32ce5460bbbc4ca896d7`.
The historical source is
[`Src/CFML_EoS_Mod.f90`](https://code.ill.fr/scientific-software/crysfml/-/blob/5bb3630d684387230065fa7c5933af2d81d19645/Src/CFML_EoS_Mod.f90),
with the Debye implementation in `Src/CFML_Math_Gen.f90`. This historical
revision is not claimed to be the exact library build used by Miozzi et al.

`Pthermal`, `EthDebye`, `Get_DebyeT`, and `Get_Grun_V` agree with the standard
MGD interpretation used here: gamma times the difference of Debye energies,
divided by molar volume; Debye energy is `3*Natom*R*T*D3(theta/T)`.
MGD uses molar cm³/mol and pressure GPa (or kbar). It does not automatically
convert crystallographic cell volumes or supply the number of atoms. For molar
Fe volume, `Natom=1`; the two-atom hcp cell conversion must be performed
separately. Normal MGD uses the integrated Gruneisen Debye-temperature law.
The separate q-compromise model is not the published model with q=0.3.

The unmodified historical thermal routines and their Debye/Chebyshev functions
were compiled with gfortran in a restricted double-precision harness. The
harness supplies the required parameter structure and volume/unit helpers;
the cold pressure uses the upstream BM3 finite-strain expression. Across all
131 observations, its total pressures differ from our independent evaluator
by at most **0.00078148 GPa**. Accounting for the source's rounded R=8.314 and
default-real unit-factor literals reduces the maximum difference to
**3.98e-13 GPa**. The compiled source still gives **6.9830 GPa RMS** against
the supplied pressures, versus our **6.9826 GPa**. Thus the large mismatch
cannot be attributed to the normal MGD equation implementation under these
settings.

This initial audit is a calculation comparison, not an EosFit refinement replay. The public
library provides pressure/parameter derivatives; the original console
least-squares driver and the authors' input/session files were not located.
The downloaded macOS console executable could not run because its XQuartz
libraries are absent. It was extracted in a temporary directory; no system
installation was performed. A later direct console replay is documented below;
the original author input/session and numerical weights remain unavailable.

### Pressure-column evidence

At 300 K, the thermal contribution is identically zero. The paired MgO volumes
in the 21 room-temperature rows permit a pressure check independent of the Fe
fit and its thermal parameters. Using the printed cold parameters from
[Speziale et al. (2001)](https://doi.org/10.1029/2000JB900318),
V0=74.71 A³, K0=160.2 GPa, K0′=3.99, gives **4.0930 GPa RMS** against the
supplied pressure column. The
[Tange et al. (2009), Table 4, Fit3-3BM](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008JB005813)
values (74.698 A³, 160.64 GPa, 4.221) give **0.4492 GPa RMS**.

| MgO-series source row | MgO V (A³) | Supplied P (GPa) | Speziale cold P (GPa) | Tange Fit3-BM3 cold P (GPa) |
|---|---:|---:|---:|---:|
| 1 | 62.625 | 41.5691 | 40.2150 | 41.1468 |
| 43 | 55.861 | 87.1354 | 83.3107 | 86.5626 |
| 72 | 52.706 | 117.8332 | 112.3328 | 117.6668 |

The main paper and supplementary discussion explicitly select Speziale after
comparing both calibrations. The supplied column being much closer to Tange
is evidence of an unresolved pressure-reduction/input discrepancy. It does
not prove which calculation generated that column or establish an author
error: even the Tange calculation does not match it exactly. A different
pressure column would change the inferred cold elasticity and thermal gamma.
This is a stronger explanation to investigate than an optimizer or MGD
implementation difference, but the full original fit has not been recovered.

The reproducible harness is `scripts/audit_miozzi_2020_eosfit.py`; source hashes,
numerical metrics and all 21 pressure comparisons are retained in
[`miozzi-2020-eosfit-code-audit.json`](../data/miozzi-2020-eosfit-code-audit.json).
Run it against a local upstream clone with:

```sh
python scripts/audit_miozzi_2020_eosfit.py --source-checkout /path/to/crysfml
```

## Conditional pressure reconstruction (2026-10-02)

**Recalculating the paired MgO measurements substantially improves the cold Fe
agreement, but does not recover the published thermal fit.** The reconstruction
and its fits remain separate diagnostic artifacts. Neither the original
observations nor the three published EOS records were changed.

The exact implementation meant by the paper's Speziale calibration is not
specified. [Speziale et al. (2001), Sections 2.1 and 5.4](https://duffy.princeton.edu/sites/g/files/toruqf616/files/speziale_et_al-2001-jgrse.pdf)
provides constant-q and volume-dependent-q models.
[Dorogokupets (2010), Equations 1–5 and Table 1](https://doi.org/10.1007/s00269-010-0367-2),
which Miozzi cites when discussing the calibration, distinguishes several
implementations. Its full text was checked in the
[author-uploaded article](https://www.researchgate.net/publication/225361067_P-V-T_equations_of_state_of_MgO_and_thermodynamics).
The Speziale parameter and variable-q pages were also visually inspected.
Three explicit scenarios were calculated for every one of the **116** paired
MgO measurements:

| Reconstruction | Cold parameters | Thermal definition |
|---|---|---|
| Speziale variable-q Debye | V0=74.71 Å³, K0=160.2 GPa, K′=3.99 | gamma0=1.524, q0=1.65, q1=11.8, thetaD0=773 K, Tr=300 K |
| Dorogokupets 2010 Fit #2 | V0=11.248 cm³/mol, K0=160.2 GPa, K′=3.99 | gamma0=1.524, gamma_inf=1.325, beta=11.8, thetaE0=599 K, Tr=298.15 K; Einstein oscillator |
| Speziale constant-q Debye sensitivity | Same as variable-q | gamma0=1.524, q=1.65, thetaD0=773 K, Tr=300 K |

For the first model, `gamma=gamma0*exp(q0/q1*((V/V0)^q1-1))`;
theta is computed by numerically integrating `d ln(theta)/d ln(V)=-gamma`.
The constant-q expression for theta is not applied to a volume-dependent q.
The Dorogokupets scenario uses the published asymptotic-gamma expression and
its integrated theta, with the source gas constant 8.31451 J/(mol K).
The Debye scenarios use 8.31446261815324 J/(mol K). MgO has n=2 atoms per
formula unit and Z=4 formula units per cell. Fe retains n=1 and Z=2.
No fitted calibration coefficient or empirical thermal-energy prefactor was
introduced. The small differences in reference temperature and rounded V0
between scenarios are retained explicitly.

The implementation reproduces three independent pressure checkpoints printed
by Dorogokupets: 17.69 GPa at x=1, T=3000 K; 54.21 GPa at x=0.85,
T=3000 K; and 203.34 GPa at x=0.64, T=3663 K, within their 0.01 GPa
rounding. Tests also check the integrated theta–gamma identity, the Debye
quadrature against adaptive integration, and the q=0 limit.

### What changes before any Fe refit

| Pressure input | Published cold Fe BM3 RMS, 36 rows (GPa) | Published thermal Fe RMS, 131 rows (GPa) | Published thermal Fe RMS, 95 heated rows (GPa) |
|---|---:|---:|---:|
| Supplied pressures | 3.7536 | 6.9826 | 7.9502 |
| Speziale variable-q Debye | 0.8570 | 5.4522 | 6.3898 |
| Dorogokupets Fit #2 | 0.8637 | 5.4431 | 6.3790 |
| Speziale constant-q Debye | 0.8570 | 2.5883 | 3.0125 |

The variable-q reconstruction shifts the supplied MgO pressure by
−5.5004 to +7.9266 GPa (RMS 2.5234 GPa). The cold published curve agrees
much better with the reconstructed pressures. Its thermal residuals still
span −10.7052 to +1.0371 GPa. The constant-q scenario reduces that mismatch
further, but still gives residuals down to −6.0782 GPa. Numerical closeness
does not identify the authors' actual implementation.

### Conditional Fe refits and errors

Both equal-weight pressure fits and effective-coordinate-error fits are saved.
For reconstructed pressures, the weighted residual variance is

```text
sigma_residual² = (dP_Fe/dV_Fe * sigma_V_Fe)²
                + (dP_MgO/dV_MgO * sigma_V_MgO)²
                + ((dP_Fe/dT - dP_MgO/dT) * sigma_T)².
```

This includes the covariance induced by the shared temperature measurement.
The printed pressure error is retained as metadata but is not reused after
recalculation. Other coordinate errors are provisionally independent; no
calibration parameter uncertainties or correlations across observations are
included. Printed error widths are provisionally treated as standard deviations,
because the source does not define their confidence level. Missing errors
exclude an observation only from the weighted fit; they are never replaced by
zero. Original printed zeros remain zero.

The printed-pressure control instead uses independent printed P/V/T errors,
since its pressure–temperature covariance is unavailable. Therefore its
weighted objective and row selection differ from the reconstructed case.
This stricter missing-error handling also differs from the earlier diagnostic
above, which omitted missing variance terms.

Cold fits vary V0, K0 and K′. Thermal fits hold V0=22.81 Å³ and theta0=420 K
fixed. An initial stage varies gamma0 and q with K0=129 GPa and K′=6.2
fixed; a subsequent stage varies K0, K′, gamma0 and q together. Effective
weights are recomputed between solves until changes in scaled parameters are
below 1e-7. The initialization uses frozen starting weights. This follows the
reported sequence at a broad level; the authors' exact cycle settings and
input selection are unavailable. Both all-131-row and heated-95-row thermal
selections are retained.

For the variable-q Debye reconstruction:

| Fit | V0 (Å³) | K0 (GPa) | K′ | gamma0 | q | Fit rows | RMS (GPa) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Cold, equal weights | 22.8063 ± 0.0976 | 128.509 ± 8.117 | 6.2869 ± 0.3064 | — | — | 36 | 0.6511 |
| Cold, effective errors | 22.6421 ± 0.1424 | 146.118 ± 11.534 | 5.6304 ± 0.3154 | — | — | 35 | 0.7353 |
| Thermal, equal weights, all rows | 22.81 fixed | 130.855 ± 1.904 | 6.1424 ± 0.1227 | 1.9626 ± 0.1131 | 0.4660 ± 0.2174 | 131 | 0.9321 |
| Thermal, effective errors, all rows | 22.81 fixed | 134.004 ± 1.411 | 5.9551 ± 0.0893 | 1.8821 ± 0.0928 | 0.3444 ± 0.1949 | 129 | 0.9244 |
| Thermal, equal weights, heated only | 22.81 fixed | 139.459 ± 4.448 | 5.7032 ± 0.2621 | 1.7139 ± 0.2081 | 0.2144 ± 0.4801 | 95 | 0.9705 |
| Thermal, effective errors, heated only | 22.81 fixed | 139.865 ± 4.352 | 5.7131 ± 0.2688 | 1.7795 ± 0.2137 | 0.4391 ± 0.5332 | 94 | 0.9583 |

These are **conditional local standard errors**, not the publication's error
definition. Each varied coefficient, including the two initialization
coefficients, has an error and full covariance. Covariance uses the full-rank
Jacobian at fixed final weights; weighted errors scale by
`max(1, reduced_chi_square)`, and equal-weight errors by `RSS/dof`.
Nonconverged, bounded, rank-deficient or nonpositive-DOF fits are rejected.
The reconstructed all-row weighted thermal fit has reduced chi-square 62.0,
so the printed coordinate errors do not describe its residual scatter under
these assumptions. Tiny nominal errors must not be taken as a complete
uncertainty budget.

The equal-weight cold parameters are close to the publication's
(22.80 Å³, 129 GPa, 6.2). The weighted cold result shows that this agreement
depends on weighting. For the thermal all-row equal-weight fits, the
Dorogokupets model gives gamma0=1.9695, q=0.4839; the constant-q calibration
gives gamma0=2.2505, q=2.0916. Neither recovers gamma0=1.11, q=0.3.
The variable-q result's q is close, while its gamma differs substantially.
The remaining mismatch is therefore not explained by replacing the pressure
column alone or by the audited normal EosFit MGD equation.

A separately labeled normalization sensitivity doubles Fe vibrational energy
at unchanged gamma, theta and molar Fe volume. Even with the variable-q
reconstructed pressures, the published coefficients then give 2.4477 GPa RMS
and a −2.0777 to +6.3418 GPa residual range. This is an intentionally altered
normalization, not the physical n=1 convention, not an adopted correction,
and not evidence of the authors' settings. It does not establish exact parity.

Recovering the original thermal result still requires the authors' actual
EosFit input/session, chosen pressure implementation, thermal normalization,
row selection and weight settings. No reconstruction is promoted into an
executable catalog record; the published records remain `not_reproduced`.

### Reproducible artifacts

- [Row-level reconstructed pressures](../data/miozzi-2020-reconstruction/pressures.csv):
  348 rows, retaining all 116 source rows for each of the three calibrations,
  original pressure/errors, computed pressure/errors, P–T covariance and
  published-Fe residuals. Original missing values remain explicit.
- [Full fit summary](../data/miozzi-2020-reconstruction/summary.json): 24 final
  fits with covariance, initialization results, metrics, fixed parameters,
  model definitions, qualifications and SHA-256 hashes of original CSV inputs.
- [Comparison figure](../data/miozzi-2020-reconstruction/comparison.png): cold
  observations, published thermal residuals, and separate equal-weight fits.
  Shading marks ±3 GPa for comparison with the paper's reported interval.

```sh
python scripts/reconstruct_miozzi_2020_iron.py --plot
python scripts/reconstruct_miozzi_2020_iron.py --check
```

Before upgrading these candidates to validated records, resolve the actual
source pressure reduction and original EosFit-7c fitting inputs/objective.
The tested standard thermal normalization is supported by the code comparison;
the authors' actual Natom/unit settings are still unknown. Published
coefficients, original pressures and validation status remain unchanged. No
author contact has been made in this task.

## Selectable Peritheos refit on the Tange scale (2026-10-02)

At the user's request, a separate executable alternative is now registered:
`iron_miozzi_2020_tange_2009_bm3_mgd_refit`, labeled **Peritheos refit of
Miozzi (2020) hcp-Fe — Tange (2009) MgO pressures, equal weights**.
It is nondefault. The Sakai (2025) default, three published Miozzi coefficient
sets, their `not_reproduced` status, and both original source CSVs are unchanged.

This uses **Miozzi's Fe observations**, not Tange's MgO fitting observations.
Tange provides the pressure calibration. All 116 paired MgO V/T observations
are recalculated using the exact bundled `mgo_b1_tange_2009_vinet` record:
[Tange (2009), Table 4, Fit3-Vinet](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008JB005813).
Its V0=74.698 Å³, K0=160.63 GPa, K′=4.367, theta0=761 K,
gamma0=1.442, a=0.138 and b=5.4 are fixed calibration inputs.
The 15 helium-series pressures remain as supplied. The calibration choice is
explicitly an analyst choice, not a claim about the authors' original reduction.

All **131** central observations receive equal pressure weight, without
excluding rows with missing errors. V0, K0, K′, gamma0 and q are fitted jointly
with physical Fe n=1 and Z=2, Tr=300 K, and theta0=420 K fixed. Unlike the
earlier diagnostic reconstruction, V0 is now free. This change in fit freedom
also affects K0; differences from the fixed-V0 fits must not be attributed
solely to the pressure calibration.

| Fitted parameter | Value | Conditional standard error |
|---|---:|---:|
| V0 (Å³/cell) | 22.578319 | 0.092116 |
| K0 (GPa) | 150.484868 | 8.628764 |
| K′ | 5.801501 | 0.266405 |
| gamma0 | 1.989345 | 0.123305 |
| q | 0.643085 | 0.237983 |

The pressure RMS is **0.984468 GPa**, with residuals −2.740718 to
+3.275262 GPa. There are 126 residual degrees of freedom. The full 5×5
covariance is stored in the record, ordered as `rt_eos.V0`, `rt_eos.K0`,
`rt_eos.K0_prime`, `gamma0`, `q`. It is the SVD inverse of `JᵀJ`, scaled by
`RSS/126`. These conditional local errors assume independent equal-variance
pressure residuals and omit calibration/fixed-parameter systematics,
predictor error, correlations across runs, and model discrepancy. Fixed
parameters have null fit errors. Parameter confidence remains null rather
than assigning the publication's undefined error convention.

All four starting points reach the same optimum, with no active bounds and a
full-rank Jacobian. The independent Tange calculation agrees with the native
evaluator across all paired observations to within 1e-9 GPa, and both Tange
calibration forms reproduce independent Table 5 checkpoints within rounding.
The executable Fe evaluator and its pressure–volume inversion are checked
against the independent refit calculation.

The [refit report](../data/miozzi-2020-tange-refit.json) also retains:

- Tange Fit3-BM3 calibration sensitivity: V0=22.680515 Å³, K0=138.156738 GPa,
  K′=6.350158, gamma0=2.006448, q=0.983225.
- Shared-temperature effective-error sensitivity: V0=22.565880 Å³,
  K0=155.081375 GPa, K′=5.595190, gamma0=1.906665, q=0.520975.
  Only complete-error rows enter this diagnostic; it includes covariance
  induced by shared Fe/MgO temperature and assumes printed widths are sigmas.
- Fixed published V0=22.81 Å³ sensitivity: K0=130.800738 GPa,
  K′=6.419571, gamma0=2.067313, q=0.721134.
- Fixed K0=160 GPa sensitivity: V0=22.483138 Å³, K′=5.525211,
  gamma0=1.933647, q=0.566859; pressure RMS=0.989409 GPa.
- Heated-only sensitivity with V0 fixed at 22.81 Å³, to expose input-selection
  dependence.

The 150.485 GPa free-fit K0 is not a strong determination that hcp Fe is
softer than common reference EOS. For example, the independent Vinet refit
of Dewaele's hcp-Fe observations by
[Morrison et al. (2018), Table 1](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2017JB015343)
gives K0=162.9 ± 4.5 GPa. The present conditional K0 error is 8.629 GPa;
its correlations with V0 and K′ are −0.969 and −0.977, respectively.
Fixing K0 at 160 GPa increases pressure RMS by only 0.004941 GPa (about
0.5%). Thus this dataset supports a family of nearly equivalent EOS around
the optimum, and does not clearly distinguish 150 from 160 GPa. These hcp
reference parameters extrapolate to zero pressure; they are not measured
ambient properties of stable bcc Fe. The fixed-K0 result is retained as a
diagnostic, with K0 omitted from its conditional fitted errors/covariance.

Every varied parameter in each sensitivity has conditional errors and
covariance. Sensitivity results are not an uncertainty distribution and are
not additional catalog selections. The main refit still does not recover
the published thermal gamma0=1.11.

The derived input dataset is `iron_miozzi_2020_tange_2009_vinet_pvt`, bundled in
`peritheos/data/datasets/iron-miozzi-2020-tange-2009-vinet-pvt.csv`.
Each row retains medium, source row/page, original P/error, Fe and MgO
volumes/errors, and temperature/error. Its canonical `pressure_gpa` is the
explicit fit coordinate; `printed_pressure_gpa` preserves the original.
Conditional propagated MgO pressure errors and P–T covariance are retained
but do not weight the preferred fit. He P–T covariance is unknown and blank.
Original missing values and zeros remain explicit. Source and calibration
hashes are pinned in the fit provenance, and the dataset and record link
to each other. EOS recalculation provenance is represented by the dataset's
`pressure_reconstruction` metadata; the typed `pressure_reductions` API
currently supports optical-scale transformations only.

The record's `primary_source_validated` status is scoped to traceability and
numerical reproducibility of **this independent Peritheos refit**. The
paper/record ledgers explicitly retain
`original_publication_reproduction_status: not_reproduced`; parity for this
refit means regeneration of stored coefficients, not recovery of the published
EOS or independent validation as a pressure standard. Actual marginal input
coverage is 10.7–129.307150 GPa and 298–3399 K. Neither the reference-state
extrapolation nor core-condition predictions receive new validation.

```sh
python -m scripts.refit_miozzi_2020_tange --register
python -m scripts.refit_miozzi_2020_tange --check
```

The first command regenerates evidence/derived inputs and registers only the
independent alternative. The second verifies pinned inputs/calibration and
refit parameters. The selection is available from the current library checkout:

```python
from peritheos import get_eos_record

refit = get_eos_record("iron_miozzi_2020_tange_2009_bm3_mgd_refit")
pressure_gpa = refit.pressure(18.0, 2000.0)
```

Installed downstream applications using pinned catalog snapshots require their
normal snapshot refresh to display this new library selection.

## Author-described three-stage EosFit replay (2026-10-03)

On 2026-10-03 the user relayed an in-person conversation with Miozzi:

1. Fit only the room-temperature observations.
2. Fit all observations with the room-temperature coefficients fixed.
3. Use stage two as the starting point and release the fitted coefficients.

The user subsequently clarified that **V0 is free in the final stage and
theta0 remains fixed at 420 K**. This is user-reported author clarification,
not a recovered author input file or a published correction. It supersedes
the earlier uncertainty about this sequence. The paper's sentence describing
fixed V0 remains preserved as a source discrepancy rather than silently edited.

`scripts/replay_miozzi_2020_staged.py` implements these three stages in the
unmodified official EosFit7c 7.60 executable, carrying each fitted result
forward within one process. Stage one fits V0, K0 and K′ to all 36 observations
at 298/300 K. Stage two fits gamma0 and q to all 131 observations with that cold
curve fixed. Stage three releases V0, K0, K′, gamma0 and q. Full MGD uses
theta0=420 K, Tr=300 K, n=1 and cm³/mol Fe throughout; q-compromise is disabled.

All three stages converge in every calibration case below. These are
**equal-pressure-weight diagnostics**; the original numerical weights and
handling of incomplete errors are still unknown. Canonical observations and
source coefficients are unchanged. Recalculated MgO pressures are separate
inputs, and the 15 He pressures stay unchanged in all cases.

| Pressure input | Stage-one K0 (GPa) | Final V0 (Å³) | Final K0 (GPa) | Final K′ | Final gamma0 | Final q | Final RMS (GPa) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Original printed pressures | 136.163 | 22.60259 | 146.810 | 6.05344 | 2.23572 | 1.51423 | 1.27798 |
| Tange Fit3-Vinet recalculation | 136.573 | 22.57832 | 150.486 | 5.80147 | 1.98945 | 0.64308 | 0.98447 |
| Speziale variable-q Debye reconstruction | 128.509 | 22.70312 | 139.419 | 5.86997 | 1.93069 | 0.42959 | 0.92813 |

The Tange result reproduces the registered independent refit, including its
conditional errors and covariance, within console rounding and gas-constant
precision. Thus following the clarified sequence does not materially change
that selection. The Speziale reconstruction recovers the published cold K0
near 129 GPa in stage one, but its final gamma0 still differs substantially
from the published 1.11. The sequence alone does not recover the source thermal
coefficients or establish that our assumed pressure reconstruction and weights
match the authors' actual input.

The [direct-run report](../data/miozzi-2020-eosfit-staged/report.json) retains
each stage's parameters, EosFit conditional errors/covariance, warnings, input
hashes, and executable hash. Each case also retains the input files, complete
macro, stdout/log and three saved EOS files. Stage-two errors condition on the
cold coefficients; they do not propagate stage-one uncertainty. RMS values
here are recomputed from the console's printed coefficients with our gas
constant, and are subject to that output precision. This evidence changes the
documented protocol, while the source-author fit remains not reproduced.

The [full comparison tables](../data/miozzi-2020-eosfit-staged/comparison.md)
include reported/refitted parameter errors, signed absolute and percentage
differences, the cold and thermal-only stages, and residual comparisons on
each pressure input.

## GUI q-compromise diagnostic (2026-10-03)

The user subsequently asked whether the GUI's q-compromise MGD estimation
could explain the discrepancy. The bundled EosFit 7.6 manual
(`EoSFitanddisplay.html`) confirms that the GUI estimates MGD in q-compromise
form for stability and allows switching to the general model after refinement.
This concerns initial estimation, not proof of the authors' final model.
The current manual identifies this model as a v7.6 addition, while the source
paper dates from January 2020 and reports a fitted q. The historical GUI
settings remain unestablished.

Q-compromise holds theta(V)=theta0 and gamma(V)/V=gamma0/V0, so thermal
pressure depends only on temperature. There is no adjustable q. It differs
from full MGD with q fixed to either zero or one. The official console was
rerun with this approximation in both thermal stages, using the same cold
stage, 131 rows, unit pressure weights, theta0=420 K, n=1 and Tr=300 K.
Stage two releases gamma0 only; stage three releases V0,K0,Kprime,gamma0.

| Pressure input | Final V0 (Å³/cell) | Final K0 (GPa) | Final Kprime | Final gamma0 | RMS (GPa) |
|---|---:|---:|---:|---:|---:|
| Supplementary printed pressures | 22.58353 ± 0.11879 | 151.288 ± 10.847 | 5.82850 ± 0.31650 | 1.93913 ± 0.04477 | 1.28915 |
| Tange Fit3-Vinet MgO pressures | 22.60233 ± 0.09375 | 146.697 ± 8.346 | 5.95921 ± 0.25453 | 2.14302 ± 0.03465 | 0.97958 |
| Speziale variable-q Debye MgO reconstruction | 22.74796 ± 0.10415 | 133.059 ± 8.280 | 6.16139 ± 0.27767 | 2.20012 ± 0.03422 | 0.93549 |

Errors are conditional local EosFit standard errors, with complete covariance
retained; they exclude pressure calibration and fixed-parameter uncertainty.
All nine fits converged. The Speziale result has K0 closer to the published
129 GPa but does not recover gamma0=1.11.

A second set of nine converged fits switches from q-compromise to full MGD
before the final stage, retaining the preceding cold curve and gamma0 in
memory and starting q at 0.3. This returns the previous full-MGD optimum
within console convergence and printing precision. For Speziale, the final
values are K0=139.41826 GPa, gamma0=1.93066 and q=0.42951, with RMS
0.928128 GPa. The GUI's starting approximation therefore does not explain
the discrepancy under the tested inputs and equal weights.

The [full q-compromise comparison](../data/miozzi-2020-eosfit-q-compromise/comparison.md)
retains all tables and replay commands. The
[q-compromise report](../data/miozzi-2020-eosfit-q-compromise/report.json) and
[q-compromise-to-full report](../data/miozzi-2020-eosfit-q-compromise-start/report.json)
retain console files, saved-model flags, source hashes and conditional
covariance. Independent pressure calculations match the console within
0.0025 GPa, allowing for gas-constant and output rounding; independent
four-parameter least-squares fits reproduce q-compromise RMS within 1e-6 GPa.
The published EOS and selectable Tange full-MGD refit remain unchanged.

## Author-confirmed full MGD with n=2 (2026-10-03)

The user then relayed two further confirmations from Miozzi: the author model
was **full MGD, not q-compromise**, and EosFit used **n=2 with V0 approximately
6.87 cm³/mol**. This updates the author-reconstruction target. The preceding
n=1 tests used the standard per-mole-Fe normalization and were not faithful
replays of this newly confirmed author setting. The q-compromise tests remain
diagnostic experiments, not an explanation of the authors' final model.
This provenance is a user-relayed in-person account; the original session,
weights and pressure column have not been recovered.

With EosFit's definition of Natom, the consistent per-mole-Fe convention is
n=1 and V0 approximately 6.87 cm³/mol. The two-atom hcp cell has Z=2, so a
per-mole-cell convention requires n=2 and V0 approximately 13.74 cm³/mol.
The Debye energy is proportional to n, while thermal pressure divides by
the molar volume. At fixed gamma and theta, using n=2 with the per-mole-Fe
volume doubles the thermal-pressure amplitude. A reported author setting
can be reconstructed without identifying it as the standard physical
normalization. In full MGD, changing n alone does not give an exact
half-gamma transformation because gamma also controls theta(V).

Both conventions were replayed with the unmodified official console. The
author-setting run holds n=2 and retains the original molar volumes. All
nine stages converged, with the same 36 RT/131 total rows, unit pressure
weights, free final V0,K0,Kprime,gamma0,q, theta0=420 K and Tr=300 K.

| Pressure input | Final V0 (Å³/cell) | Final K0 (GPa) | Final Kprime | Final gamma0 | Final q | RMS (GPa) |
|---|---:|---:|---:|---:|---:|---:|
| Supplementary printed pressures | 22.60890 ± 0.12411 | 146.082 ± 11.387 | 6.07796 ± 0.37080 | 1.12654 ± 0.09609 | 1.58039 ± 0.32995 | 1.27445 |
| Tange Fit3-Vinet MgO pressures | 22.58426 ± 0.09239 | 149.777 ± 8.626 | 5.82394 ± 0.26808 | 1.00371 ± 0.06204 | 0.71483 ± 0.23857 | 0.97879 |
| Speziale variable-q Debye MgO reconstruction | 22.70980 ± 0.09880 | 138.716 ± 8.367 | 5.89337 ± 0.27576 | 0.97443 ± 0.05756 | 0.50272 ± 0.22274 | 0.92278 |

The errors are conditional local standard errors, with covariance retained;
pressure-calibration and fixed-parameter uncertainties are excluded. The
change explains the dominant near-factor-of-two gamma discrepancy, but does
not recover the full published fit. Evaluating the unchanged published
coefficients with n=2 reduces RMS from 6.98261 to 2.49886 GPa on the printed
pressures and from 5.45224 to 2.44770 GPa on the Speziale reconstruction.
Remaining pressure-input, weighting and K0/q differences are unresolved.

A separate nine-stage control sets n=2 and doubles every molar volume,
volume error and initial V0. It returns the earlier n=1 full-MGD optimum
within console convergence/printing precision. Independent calculations
verify the normalization equivalence, all console pressures within 0.0025
GPa, and the author-setting five-parameter optimum and conditional covariance.
The [full n=2 comparison](../data/miozzi-2020-eosfit-n2-author-volume/comparison.md),
[author-setting report](../data/miozzi-2020-eosfit-n2-author-volume/report.json),
and [doubled-volume control](../data/miozzi-2020-eosfit-n2-double-volume-control/report.json)
retain the evidence and replay commands. Existing library records and the
selectable n=1 Tange refit are unchanged by these diagnostic runs.

## Whole-hcp-cell refit with Peritheos (2026-10-03)

At the user's request, fresh fits were performed using Peritheos's own fitting
APIs, with the entire two-atom hcp cell as the molar entity and **n=2**.
The cell volumes remain the measured coordinates in Å³; inside MGD they are
converted to volume per mole of whole cells, rather than per mole of Fe.
For 22.81 Å³/cell the whole-cell starting reference volume is 13.73650 cm³/mol,
or 1.373650 J/bar/mol. Both observed volumes and V0 are twice their per-Fe
molar values. This is distinct from the author-reported n=2 with 6.87 cm³/mol.

Each pressure input was fitted independently on both molar bases, using the
36 RT rows first, all 131 rows with the cold coefficients fixed second, then
all 131 rows with V0,K0,Kprime,gamma0,q free. Theta0=420 K and Tr=300 K remain
fixed. Equal pressure-residual weights and bounds expressed on the same
physical cell basis were used in both runs. All 18 stages converged.

| Pressure input | Whole-cell V0 (Å³/cell) | K0 (GPa) | Kprime | gamma0 | q | RMS (GPa) |
|---|---:|---:|---:|---:|---:|---:|
| Original printed pressures | 22.602632 | 146.806469 | 6.053557 | 2.235644 | 1.514293 | 1.27798169 |
| Tange Fit3-Vinet recalculation | 22.578319 | 150.484853 | 5.801501 | 1.989345 | 0.643085 | 0.98446756 |
| Speziale variable-q Debye reconstruction | 22.703126 | 139.419572 | 5.869960 | 1.930567 | 0.429554 | 0.92812717 |

The independently optimized n=1 and n=2 pressure surfaces differ by at most
6.36e-7 GPa on the 131 observations. The same fitted coefficients under
matching n/volume scaling give identical pressures within numerical precision.
A separate direct whole-cell calculation using n=2, kB and volume in m³/cell
agrees with the Peritheos molar calculation to better than 4.39e-10 GPa.
Thus the physically consistent whole-cell convention reproduces the existing
per-Fe fit; it does not halve gamma0 or resolve the original-publication
discrepancy. Library records and default selections are unchanged.

The [comparison with conditional errors](../data/miozzi-2020-cell-basis/comparison.md)
and [complete fit report](../data/miozzi-2020-cell-basis/report.json) retain
the six fits, all stages, residuals, covariance, source/code fingerprints,
and independent normalization checks. Errors are conditional local standard
errors, excluding calibration, predictor, fixed-theta and inter-row covariance
uncertainties. Reproduce with
`python -m scripts.refit_miozzi_2020_cell_basis`.

## Registered published n=2 record and Speziale n=1 refit (2026-10-03)

The executable `iron_miozzi_2020_bm3_mgd` record now retains the published
coefficients with **n=2** and V0 corresponding to approximately **6.87 cm³/mol**.
This was the setting used in the paper according to **personal communication
with Miozzi et al.**, relayed by the user and recorded on 2026-10-03. The
communication also confirms full MGD, theta0=420 K fixed and final V0 free;
the article's fixed-V0 wording remains documented. No printed coefficient or
reported uncertainty is replaced by a refit. The record remains explicitly
selected and `not_reproduced`: this setting improves the residuals but does
not recover the complete original fit or establish physical normalization.

The separate `iron_miozzi_2020_speziale_2001_bm3_mgd_refit` record uses the
consistent **n=1 per mole Fe** basis. The 116 MgO-series pressures are
recalculated with the Speziale (2001) variable-q Debye calibration described
above; the 15 source He pressures remain unchanged. The calibration is a
documented thermal reconstruction, not the isothermal-only Speziale catalog
record and not a recovered author input file. All 131 central observations
are fitted with equal pressure-residual weights in the three-stage sequence.
V0,K0,Kprime,gamma0,q are free in the last stage; theta0=420 K and Tr=300 K
remain fixed. The resulting coefficients are:

| Parameter | Published record, n=2 | Independent Speziale refit, n=1 |
|---|---:|---:|
| V0 (Å³/cell) | 22.81; printed molar error maps to ±0.06642 | 22.703126 ± 0.098420 |
| K0 (GPa) | 129 ± 1 | 139.419575 ± 8.365237 |
| Kprime | 6.24 ± 0.04 | 5.869960 ± 0.273780 |
| gamma0 | 1.11 ± 0.01 | 1.930567 ± 0.114262 |
| q | 0.3 ± 0.3 | 0.429554 ± 0.221943 |
| theta0 (K) | 420, fixed | 420, fixed |
| Tr (K) | 300 | 300 |
| n | 2, author-reported | 1, normalized per Fe |
| RMS on reconstructed Speziale input (GPa) | 2.447704 | 0.928127 |

Published error confidence is unspecified. Refit errors are conditional local
standard errors; the complete five-parameter covariance includes cold/thermal
cross-correlations but excludes calibration, fixed-theta, predictor and inter-row
uncertainties. The refit is labeled `independent_refit_reproduced`, referring to
its documented numerical procedure, while
`original_publication_reproduction_status` remains `not_reproduced`.

The derived dataset `iron_miozzi_2020_speziale_2001_variable_q_pvt` retains original
pressures/errors alongside recalculated pressures, all original coordinates,
conditional propagated MgO errors and available P-T covariance. Missing errors
remain missing. Equal weights do not use these conditional errors as fit weights.
The [refit report](../data/miozzi-2020-speziale-refit.json) retains each stage,
residuals, covariance, calibration coefficients, bounds and source/code hashes.
Neither independent refit changes the existing iron default.

```python
from peritheos import get_eos_record

published = get_eos_record("iron_miozzi_2020_bm3_mgd")
refit = get_eos_record("iron_miozzi_2020_speziale_2001_bm3_mgd_refit")
# Both record wrappers accept measured hcp cell volumes in Å³.
# Internally, both use volume per mole Fe: published n=2 versus refit n=1.
```

Reproduce and check with
`python -m scripts.refit_miozzi_2020_speziale --check`.
