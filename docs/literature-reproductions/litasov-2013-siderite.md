# Litasov et al. (2013): natural siderite

The 300 K BM3 record **`siderite_fe095mn005_litasov_2013_bm3`** is accepted.
Both published thermal versions remain non-executable because their expansion
intercepts disagree. A separate **`siderite_fe095mn005_litasov_2013_joint_thermal_refit`**
record is executable and explicitly opt-in. The published 300 K record remains
the default. This supersedes the earlier source-access hold.

Primary reference: K. D. Litasov, A. Shatskiy, P. N. Gavryushkin,
I. S. Sharygin, P. I. Dorogokupets, A. M. Dymshits, E. Ohtani, Y. Higo,
and K. Funakoshi (2013), *P–V–T equation of state of siderite to 33 GPa and
1673 K*, Physics of the Earth and Planetary Interiors **224**, 83–87,
[doi:10.1016/j.pepi.2013.07.011](https://doi.org/10.1016/j.pepi.2013.07.011).

The user supplied the complete **Accepted Manuscript**, accepted 18 July 2013,
21 PDF pages. Manuscript page numbers below exclude its cover: Table 1 is
manuscript page 12/PDF page 13; Table 2 spans manuscript pages 13–15/PDF pages
14–16. Its checksum and retrieval history are in the
[source audit](../data/litasov-2013-siderite-source-audit.json).
The final publisher abstract was also checked; the full typeset article and
any official supplementary/correction files were not recovered. No absence of
such files is inferred. The copyrighted PDF is not redistributed.

## Identity and structure

Section 2 identifies natural **Fe0.95Mn0.05CO3**, from the Panasqueira tungsten
mine, Covilhã, Castelo Branco, Portugal. Powder made from the single crystal
was mixed with 5 wt% Au. Pre/post-run microprobe detected no composition change.
This material remains separate from pure FeCO3 and from low-spin siderite.

The conventional hexagonal cell of R-3c (No. 167) has Z = 6 and contents
Fe5.7Mn0.3C6O18. The card's lattice uses the first ambient Table 2 observation,
a = 4.6934 Å and c = 15.3844 Å. Its geometrical volume, approximately
293.485 Å³, agrees with the reported 293.49(3) Å³. The EOS instead uses the
paper's fixed ambient mean V0 = 293.4(1) Å³; these values are not forced equal.

Fractional sites come separately from Effenberger, Mereiter and Zemann (1981),
[Table 3, page 237](https://rruff.info/doclib/zk/vol156/ZK156_233.pdf): cations
on 6b (0,0,0), C on 6a (0,0,¼), and O on 18e (0.27427,0,¼).
Indexed primary-table text supplies these sites and agrees with
[COD 9017361 / AMCSD 0020837](https://www.crystallography.net/cod/9017361.cif).
Direct retrieval of that older primary PDF failed. The deposited structure
has matching Fe0.95Mn0.05 chemistry but is a different crystal from Ivigtut,
Greenland. Fe/Mn share the cation site with occupancies 0.95/0.05. The card
explicitly combines those sites with the Panasqueira lattice; it does not
claim a source-specific atomic refinement or use the database as EOS authority.

## Complete source-owned fit inventory

| Candidate | Source values | Disposition |
|---|---|---|
| 300 K volume BM3 | V0 = 293.4(1) Å³ fixed; K0 = 120(1) GPa; K′ = 3.57(9) | Accepted, published values unchanged |
| Thermal volume BM3 | Same cold coefficients; dK/dT = −0.015(1) GPa/K; a1 = 0.06(14) × 10⁻⁸ K⁻²; a0 = **3.57(9)** × 10⁻⁵ K⁻¹ in the abstract versus **3.77(9)** × 10⁻⁵ in Table 1 | Withheld; both variants evaluated diagnostically |
| a³ axial BM3 | Fixed fictive V0a = 103.4 Å³; Ka = 166(2) GPa, K′a = 14.0(8) | Axial diagnostic, not a physical-volume EOS |
| c³ axial BM3 | Fixed fictive V0c = 3633.8 Å³; Kc = 59(1) GPa, K′c = 2.7 | Axial diagnostic; K′c error is 0.1 in the abstract versus 0.4 in Section 3 |

Table 1 gives K0 = 120(2), whereas the room-temperature abstract gives
120 ± 1 GPa. The RT card retains the latter and records both source widths.
The table summarizes P–V–T properties; neither width is silently treated as
a known confidence interval. Parameter covariance and confidence convention
are not reported. Parenthetical **observation** errors in Table 2 are explicitly
1σ. Comparison columns for Zhang (1998), Lavina (2010), Badaut (2010, DFT), and
Sanchez-Valle (2011, adiabatic moduli) are not additional fits by this paper.

## Equation mapping and fitting protocol

With x = (V0/V)^(1/3), the RT equation is

```text
P = 3*K0/2 * (x^7-x^5) * [1 + 3/4*(Kprime-4)*(x^2-1)].
```

The cited HTBM convention is explicitly given by Litasov et al. (2007),
[Section 3.1, Equations (1)–(4)](https://doi.org/10.1016/j.pepi.2007.06.003):

```text
Tr = 300 K
alpha(T) = a0 + a1*T = (1/V0(T)) * dV0(T)/dT
V0(T) = V0 * exp[a0*(T-Tr) + a1*(T^2-Tr^2)/2]
K0(T) = K0 + dK_dT*(T-Tr)
Kprime(T) = Kprime.
```

Thus a0 is an absolute-temperature intercept and alpha is volumetric.
Peritheos's existing `ThermalReferenceStateEOS`, `linear_temperature`
expansivity and `integrated_expansivity` volume law could express it exactly;
no new model is needed. A source-reported thermal production record is nevertheless withheld
under the [acceptance protocol](../adding-materials-and-eos.md), because the
coefficient conflict has not been resolved by an authoritative correction.

Section 3 explicitly varies K0, K′, dK/dT, a0 and a1 **jointly**, with V0
fixed. The abstract describes RT fitting followed by thermal analysis.
The audit therefore performs the joint five-parameter fit and additionally
tests staged fits holding either the published or independently refitted cold
coefficients. Source weighting, residual definition and any unprinted row
exclusions are not supplied. No rows are discarded to improve agreement.

## Observations and pressure calibration

The checksummed CSV contains **all 111 rows**: runs 1–4 have 66, 12, 6 and 27
rows respectively. It retains T, original pressure, paired Au conventional-cell
volume, siderite a/c/V, their printed 1σ errors, run, row and PDF page.
There are 27 rows at 300 K and 84 heated rows. Across the complete table,
the largest discrepancy between printed V and √3 a²c/2 is only 0.00535 Å³,
consistent with the displayed precision. No synthetic observations are used.

The paper states a pressure-error ceiling of 0.1 GPa, not row-wise pressure
standard deviations. It prints no Au-volume or temperature errors. The
W97Re3/W75Re25 thermocouple junction was at the X-ray position; pressure
correction to its EMF was ignored. The Methods use Au scales from Litasov
et al. (2013a), JAP **113**, 093507, and Sokolova et al. (2013), RGG **54**,
181–199. Dorogokupets and Dewaele (2007) is a comparison, not the adopted scale.

The related catalog `gold_sokolova_2013_holzapfel_4` uses a corrected 2016
realization. Evaluating the printed marker volumes gives 0.166 GPa pressure
RMSE and a maximum 0.716 GPa discrepancy. That compatibility check does not
establish the exact 2013 reduction. Calibration is `partially_resolved`, with
recalculation `not_possible` pending that reduction, despite available marker
observations. No executable exact-scale link or pressure substitution is made.

## Independent reproduction

Run `uv run python scripts/reproduce_litasov_2013_siderite.py` to regenerate
the [numerical audit](../data/litasov-2013-siderite-reproduction.json).
Sample pressure expressions and least-squares fits are independent of the
Peritheos implementation. The catalog is used only for the separately labelled
gold compatibility diagnostic.

The predefined highest-pressure benchmark is Table 2 row 79: V = 241.54(4) Å³,
T = 300 K, P = 33.01 GPa. Published BM3 coefficients predict **32.92407 GPa**,
within the source's 0.1 GPa pressure bound. This is an off-reference primary
observation, not an identity at V0 or a generated expected value.

| Diagnostic | K0 (GPa) | K′ | a0 (10⁻⁵ K⁻¹) | a1 (10⁻⁸ K⁻²) | dK/dT (GPa/K) | Pressure RMSE (GPa) |
|---|---:|---:|---:|---:|---:|---:|
| RT, 27 rows, fixed V0 | 120.2242 | 3.572214 | — | — | — | 0.22731 |
| Joint thermal, 111 rows | 121.9029 | 3.427642 | 3.58784 | 0.101323 | −0.0145979 | 0.22121 |
| Heated 84 rows, published cold fixed | 120 | 3.57 | 3.81618 | −0.003997 | −0.0152890 | 0.22639 |
| Heated 84 rows, refitted cold fixed | 120.2242 | 3.572214 | 3.75572 | 0.083027 | −0.0156258 | 0.22425 |

These fits use equal pressure-residual weights. A separate sensitivity using
only propagated volume errors is included in the JSON; it ignores unknown
row pressure errors and is not identified as the original source objective.
Estimated refit covariances are residual-scaled local estimates, clearly
separate from the unavailable source covariance.

The original RT coefficients give 0.22958 GPa RMSE. The thermal Table 1
and abstract variants give 0.23256 and 0.29937 GPa, respectively, across all
111 observations. None reproduces the manuscript's stated 0.11–0.17 GPa RMS
range. A lower residual alone cannot identify the intended coefficient.
RT coefficient recovery is classified **similar**, not original-fit parity.
The all-cold-row axial diagnostics give Ka = 219.263, K′a = 4.94465 and
Kc = 61.7153, K′c = 2.65811; their disagreement with the printed axial fits
is preserved. No axial material record has been added.

## Executable thermal refit and source evidence

The material's `source.thermal_parameterizations` stores both complete printed
thermal coefficient sets, their errors, source locations and `executable: false`.
They are provenance metadata, not members of `eos_records`, and neither can be
selected with `get_eos_record()`. No typo is declared resolved.

The executable thermal record is a **Peritheos refit**, using the full joint
result rather than mixing a fitted expansion intercept with published cold
coefficients:

| Parameter | Stored refit | Conditional standard error |
|---|---:|---:|
| V0 (Å³) | 293.4, fixed | Not propagated |
| Tr (K) | 300, fixed | Not propagated |
| K0 (GPa) | 121.9028804 | 0.88559 |
| K′ | 3.42764232 | 0.08617 |
| a0 (K⁻¹) | 3.58783926 × 10⁻⁵ | 1.44217 × 10⁻⁶ |
| a1 (K⁻²) | 1.01323442 × 10⁻⁹ | 1.95881 × 10⁻⁹ |
| dK/dT (GPa/K) | −0.0145979110 | 0.00118842 |

Selecting this refit uses its own K0 and K′ at every temperature, including
300 K; it does not switch to the published cold coefficients.

It uses all 111 printed observations with equal pressure-residual weights.
The full five-parameter covariance, including cold/thermal cross-correlations,
is retained in physical units. Errors are local, residual-scaled estimates
conditional on fixed V0 and Tr; they exclude the reported V0 uncertainty,
pressure-scale systematics and model inadequacy. They are not the authors'
reported coefficient uncertainties. The fit's 0.22121 GPa pressure RMSE
measures agreement with its training observations, not holdout predictive skill.

```python
from peritheos import get_eos_record

thermal = get_eos_record("siderite_fe095mn005_litasov_2013_joint_thermal_refit")
volume = thermal.volume(20.0, 1200.0, check_validity=True)
```

The native evaluator matches the independent reference expression on all 111
states and inverts those pressures back to their input volumes. Peritheos's
separate native joint fitter also recovers the independent SciPy fit. At the
high-temperature Table 2 state V = 247.82 Å³, T = 1673 K, the refit predicts
32.35748 GPa against the printed 32.10 GPa; its conditional parameter-only
pressure standard error is 0.30757 GPa. This is an in-sample check, not an
independent validation experiment. Out-of-envelope calls are rejected when
`check_validity=True`.

The ledger's `parity` status for this record means reproduction of the stored
Peritheos refit. It does not mean recovery of either published thermal set.

## Domain and remaining limits

The full observation paths span 0–33.01 GPa and 300–1673 K. The published
record is restricted to **300 K**; the separate refit represents the full
measured temperature envelope. Marginal extrema are not a rectangular
phase-stability guarantee. Section 4 discusses redox-dependent decomposition
and defers the authors' own numerical decomposition boundary to subsequent
papers. It cites a spin transition near 45 GPa with roughly 10% volume collapse,
outside this EOS's measurements. No low-spin extension or invented stability
boundary is represented.

The remaining source work is to resolve thermal a0 from a final table,
correction or author clarification, identify the original fitting objective
and row treatment, and reproduce the exact 2013 Au reduction. The retained
data and diagnostics make those questions reviewable without altering the
published room-temperature record.

