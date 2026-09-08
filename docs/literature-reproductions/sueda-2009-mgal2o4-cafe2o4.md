# Sueda et al. (2009): CaFe2O4-type MgAl2O4

Audit date: 2026-09-08. **Accepted: two published thermal parameterizations.**

## Source and record scope

Yuichiro Sueda, Tetsuo Irifune, Takeshi Sanehira, Takehiko Yagi, Norimasa
Nishiyama, Takumi Kikegawa, and Ken-ichi Funakoshi (2009), “Thermal equation
of state of CaFe2O4-type MgAl2O4,” *Physics of the Earth and Planetary
Interiors* **174**(1–4), 78–85.
DOI: [10.1016/j.pepi.2008.07.046](https://doi.org/10.1016/j.pepi.2008.07.046).
The DOI contains 2008; the final issue year is 2009.

The published eight-page PDF was read from the user's Zotero library, parent
item `9JHZNNDC`, attachment `FLBPJ9MU`. Its checksum is recorded in each EOS
record's `primary_source_check`; Tables 1–2 and the equations were also checked
as rendered pages. This full-text inspection supersedes the preliminary
abstract-only access hold. No paper or supplement is redistributed.

Both fits extend the existing `mgal2o4_cafe2o4` material:

| Record identifier | Published parameterization | Disposition |
|---|---|---|
| `mgal2o4_cafe2o4_sueda_2009_htbm_2` | Table 2 HTBM-EOS | Accept |
| `mgal2o4_cafe2o4_sueda_2009_mgd_3` | Table 2 MGD-EOS | Accept |

The shared reference isotherm is not duplicated as a third executable record.
The Funamori record remains the equilibrium default. The two Sueda alternatives
are not separate samples or phases. Other Table 2 rows and the non-CF Table 3
parameters are citation comparisons, not new Sueda fits. Table 3's CF density
is an adopted MORB-composition value, not the density of pure MgAl2O4.

## Equation and parameter mapping

Equation (1) is standard BM3, with the fitted 300 K triplet
`V0=240.1(2) A^3`, `K0=205(6) GPa`, and `K0'=4.1(3)`. Table 1's orthorhombic
axes establish the approximately 240 A^3 conventional cell consistent with the
existing Z=4 material. The structure retains its separate Funamori provenance;
Sueda does not supply a new refined atomic-site model. Preserve the printed
volumes rather than replacing them with products of rounded axes.

Both thermal fits fix the static triplet at its published estimates. Table 2's
static errors remain stored as first-stage errors, not refitted thermal-stage
uncertainties. Parenthetical error confidence, covariance, exact residual
objective, and numerical weights are not specified. Missing values remain null.

For HTBM, Equations (2)–(4) define

```text
V0(T) = 240.1 exp[a0 (T-300) + b0/2 (T^2-300^2)]
K0(T) = 205 - 0.030 (T-300)
K0'(T) = 4.1
```

The published coefficients are `a0=1.96(13)e-5 K^-1`,
`b0=1.64(24)e-8 K^-2`, and `dK0/dT=-0.030(2) GPa/K`.
The source explicitly fixes `c0=0` after its fit failed to converge with that
coefficient free. This maps to `thermal_reference_state`,
`thermal_expansion_law=linear_temperature`, and
`reference_volume_law=integrated_expansivity`. In particular, the stored
`alpha0` is the intercept at T=0, not the expansivity at 300 K.

For MGD, Table 2 gives `gamma0=1.73(7)`, `q=2.03(37)`, and
`theta0=1546(104) K`. Equations (5)–(8) give the 300 K referenced thermal
increment, `gamma=gamma0*(V/V0)^q`, and
`theta=theta0*exp[(gamma0-gamma)/q]`. The implementation uses
`debye_temperature_law=integrated_gruneisen` and `n=7` atoms per MgAl2O4
formula unit. The material loader converts cell volume to molar formula-unit
volume using Z=4. The independent reproduction instead uses 28 atoms per
physical cell and the Boltzmann constant, checking this normalization directly.

### Printing errors and explicit interpretation

The rendered page 82 has several errors; none is silently reproduced in code:

- Equation (3) prints `a0` on the left where the surrounding prose and Table 2
  define `alpha0(T)`.
- Equation (4) prints a derivative denominator `partial K`; the prose and
  Table 2 unambiguously specify `partial T`.
- Equation (6) places the energy difference in the denominator. Use
  `gamma/V * [Eth(V,T)-Eth(V,300)]`, the pressure-dimensioned form with zero
  thermal increment at 300 K. The literal expression is singular there.
- The Debye integral's lower limit in the first Equation (7) is printed as
  theta. Use zero: the integration variable is dimensionless, and the literal
  theta limit cannot define the stated thermal energy.
- Two different equations carry number (7); the second is the correctly
  printed integrated-Gruneisen temperature law.

These interpretations follow the stated physical models, dimensional checks,
reference identity, and numerical agreement with Table 1 and Figure 4. They are
explicit audit interpretations, not claimed publisher errata.

## Primary observations and calibration

All 46 Table 1 rows are bundled in
`peritheos/data/datasets/mgal2o4-cafe2o4-sueda-2009-table1-pvt.csv` with a SHA-256
resource link: 21 M114, four M140, five M194, five CF-1, and eleven CF-2 rows.
The table includes 33 Au volumes at compressed states and four ruby comparison
pressures. All printed CF axes, volumes, and parenthetical errors are retained,
including the unusual M194 `221.8(20)` volume error (2.0 A^3) at 2000 K.
Absent errors/calibrants are empty fields, not zero. The repeated 298 K points
and atmospheric `10^-4 GPa` pressures retain their printed values.

Sections 2.1–2.2 select Anderson, Isaak, and Yamamoto's (1989) Au EOS for both
KMA and DAC fitting. The bundled `gold_anderson_1989_bm3_1` is linked explicitly.
Mao et al. (1986) ruby pressures are checks only. Section 2.3 uses Au thermal
expansion from Touloukian et al. (1975) for HTXRD peak-position calibration;
no raw peak positions are printed. KMA temperature fluctuations were kept
within ±5 K, and HTXRD control within ±2 K; these are not invented row-wise
standard deviations or fitting weights.

## Independent reproduction and refit

Run `python scripts/reproduce_sueda_2009_mgal2o4.py`. The script evaluates the
above equations independently of Peritheos and refits with SciPy. It uses
unweighted pressure residuals because the authors do not state their weights.
The thermal stage follows the explicit all-Table-1 selection while fixing the
printed static triplet, rather than substituting this audit's static estimates.

| Fit | Published-curve pressure RMSE | Refit pressure RMSE | Largest thermal coefficient difference / printed error |
|---|---|---|---|
| HTBM, 46 observations | 0.307386 GPa | 0.292911 GPa | 0.73 |
| MGD, 46 observations | 0.288157 GPa | 0.283862 GPa | 0.96 |

HTBM recovers approximately `a0=2.05444e-5`, `b0=1.52832e-8`,
`dK/dT=-0.0301358`. MGD recovers `gamma0=1.77195`, `q=2.38351`,
`theta0=1535.20 K`. These are diagnostics, not replacement EOS records.
The generated refit ledger independently repeats the thermal fits through the
native Peritheos fitter and classifies both as parity.

The 16 explicitly 300 K rows give a diagnostic static triplet
`240.14447 A^3, 204.58272 GPa, 4.09682`, within every printed error.
Section 3.1 states a KMA lower bound of 28.7 GPa, whereas Table 1 and Figure 2
also include 24.8 and 26.6 GPa states. The 16-row reconstruction is documented
as such; the source does not print a row-level static selection mask.

At Table 1's `(218.9 A^3, 1800 K, 32.4 GPa)` state the predicted pressures are
32.004484 GPa (HTBM) and 31.978130 GPa (MGD). At
`(223.9 A^3, 2400 K, 30.1 GPa)` they are 30.917103 and 30.428270 GPa.
The tests allow 0.6 and 1.0 GPa respectively, covering propagated printed
pressure/volume errors conservatively. At ambient pressure and 836 K both
models reproduce `243.8(3) A^3`. Native and independent curves agree across
all 46 rows to better than 2e-10 GPa. Not every point lies inside a single
printed error; the maximum residuals are 0.817103 and 0.651881 GPa.

The experimental marginal range is 0.0001–41.5 GPa and 298–2400 K. KMA,
room-temperature DAC, and ambient HTXRD cover distinct paths; this is not a
rectangular stability field. The abstract's 42 GPa is rounded.

## Candidate and ledger disposition

This DOI has no publication row in the imported LitCurate ledger. The
manual candidate is marked accepted in `docs/material-eos-candidates.md`;
no deposit-derived row or identifier is fabricated. The preliminary
nonproduction hold is removed. The manifest, executable primary-source audit,
refit ledger, paper ledger, source index, and inventory counts now include
the two accepted records and one new primary dataset.
