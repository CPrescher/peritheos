# Hirose et al. (2008): gold EOS and Au–MgO cross-calibration

Hirose, K., Sata, N., Komabayashi, T., and Ohishi, Y., *Physics of the Earth
and Planetary Interiors* **167**, 149–154,
[doi:10.1016/j.pepi.2008.03.002](https://doi.org/10.1016/j.pepi.2008.03.002).
The final article was read from Zotero item `RCISM6W3`, PDF attachment
`2DS57FBL`. Its checksum is retained in the dataset provenance.

Three separate gold records preserve Table 2. The preferred thermal result
is `gold_hirose_2008_bm3_fit2`; the alternatives are
`gold_hirose_2008_bm3_fit1` and `gold_hirose_2008_bm3_300k`.
They supplement the existing gold catalog; its default is unchanged.

## Equation and parameters

All three use the explicit third-order Birch–Murnaghan pressure expression
in Section 3.1 with conventional four-atom fcc-cell volume. The fixed
reference values are V0=67.85 Å³ and K0=167 GPa. The room-temperature fit
finds K0′=5.58(2), then holds that value fixed in both thermal fits.

The thermal records use `ThermalReferenceStateEOS` with
`thermal_expansion_law="linear_temperature"`:

\[
V_0(T)=67.85\exp[\alpha_0(T-300)+\tfrac12\alpha_1(T^2-300^2)].
\]

| Parameter | Fit #1 | Fit #2 |
|---|---:|---:|
| α0 (K⁻¹) | 3.179(139) × 10⁻⁵ | 3.824 × 10⁻⁵, fixed |
| α1 (K⁻²) | 1.477(310) × 10⁻⁸ | 1.499(18) × 10⁻⁸, fixed |
| dK/dT (GPa/K) | −0.028(3) | not used |
| β1 (GPa⁻¹ K⁻¹) | — | 1.03(44) × 10⁻⁶ |
| β2 (GPa⁻¹ K⁻²) | — | 3.95(388) × 10⁻¹⁰ |
| β3 (GPa⁻¹ K⁻³) | — | 1.61(97) × 10⁻¹³ |
| a (K⁻¹) | — | 3.61(22) × 10⁻⁴ |

Fit #1 uses K0(T)=167−0.028(T−300) GPa with constant K0′. Fit #2 uses
`bulk_modulus_law="reciprocal_cubic"` and
`kprime_log_coefficient=a`:

\[
K_0(T)^{-1}=167^{-1}+\beta_1(T-300)+\beta_2(T^2-300^2)
+\beta_3(T^3-300^3),
\qquad K_0'(T)=5.58+a(T-300)\ln(T/300).
\]

The paper prints the equivalent unshifted polynomial with
β0=5.639(0)×10⁻³ GPa⁻¹, while explicitly fixing K0(300)=167 GPa. Its rounded
coefficients give 167.003540642 GPa. Peritheos enforces the fixed reference
state, deriving β0=0.0056391269521 GPa⁻¹, which rounds to the printed value.
Both values and the derivation are retained in `source_coefficient_rounding`.
This also makes the thermal increment exactly zero at 300 K.

The α1 error for fit #2 comes from the separately fitted ambient-expansivity
relation on page 153. It was fixed during the EOS regression. Printed
parenthetical errors have no stated confidence convention; no covariance is
invented. The zero printed error on β0 is not treated as an independently
measured zero uncertainty. These mechanical equations do not specify a
caloric free energy.

## Observations and numerical checks

The packaged `gold-hirose-2008-table1.csv` retains all 21 Table 1 rows:
nine 300 K observations, four run-2 heating points, and eight run-3 points.
It preserves both MgO pressure columns, all seven Au pressure columns,
Au/MgO lattice parameters, NaCl B1/B2 lattice parameters where printed,
and all parenthetical errors. Conventional volumes are explicitly derived
as a³, with linearized uncertainty 3a²δa. Missing errors remain missing.

Run 3 uses a single peak per phase and gives no lattice errors. Its printed
pressure errors reflect temperature only; the authors estimate an additional
MgO volume contribution as large as 0.5 GPa. The reported <±10% spatial
variation in temperature is not silently converted into row-wise standard
deviations. The isolated 140 GPa point has reported peak overlap and is retained.

Run `python -m scripts.reproduce_hirose_2008_gold` to reproduce the audit.
The independent Python expression and native implementation agree within
10⁻¹⁰ GPa. Fit #2 matches all 21 printed Au pressures within 0.0563 GPa,
consistent with rounded coefficients and table precision. At aAu=3.7152 Å
and 2070 K, it gives 119.7147 GPa versus the printed 119.7 GPa.

| Diagnostic | Published curve RMSE (GPa) | Published curve MAE (GPa) | Unweighted refit |
|---|---:|---:|---|
| 300 K, all nine rows | 0.670655 | 0.555772 | K0′=5.583345 |
| Fit #1, all twelve hot Hirose rows | 0.459525 | 0.337300 | dK/dT=−0.0283325, α0=3.176726×10⁻⁵, α1≈1.49981×10⁻⁸ |
| Fit #2, all 38 hot observations | 0.447919 | 0.359374 | RMSE=0.445142 GPa, MAE=0.352703 GPa |

Fit #2 uses the complete recovered selection of 38 hot observations:
26 from Fei et al. (2004),
[doi:10.1016/j.pepi.2003.09.018](https://doi.org/10.1016/j.pepi.2003.09.018),
and twelve from Hirose Table 1.
These are all rows of the packaged `gold-fei-2004-table1.csv`, combined with
the twelve hot rows of Hirose Table 1 by `complete_fit2`. The selection matches
the reported count, paired Au–MgO phases and marginal ranges; it includes the
8.64 GPa Fei row because the paper's lower bound of 9 GPa is rounded. No
Fei Table 2 Au–Pt rows are included. The
[Fei source manifest](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/datasets/fei-2004-pressure-scales-source.json)
records the primary-table provenance, selection evidence and source rights.
The source-reported MgO pressures remain the fit targets; the calculated Au
pressure columns are checked separately and are not substituted as observations.

The complete selection gives a published-curve MAE of 0.359374 GPa, close to
the paper's reported 0.35 GPa. The unweighted pressure-residual refit gives
β1≈9.9243×10⁻⁷, β2≈4.3171×10⁻¹⁰, β3≈1.5138×10⁻¹³ and
a≈3.5949×10⁻⁴, with MAE 0.352703 GPa. All four coefficients lie within
the printed error widths. The 300 K and fit #1 diagnostic refits also recover
coefficients within their printed widths. These refits retain the staged
selection and fixed parameters, including fit #2's ambient expansivity.

Complete published PVT recovery does not establish the original fit objective.
The source does not specify residual coordinate, weights or covariance, and
unrounded observations and the underlying Touloukian ambient-expansivity data
are unavailable. The printed-pressure-error sensitivity fit uses those errors
only as relative weights, not established standard deviations; it omits Au
volume and temperature uncertainty and correlated calibration errors. Its MAE
is about 1.847 GPa, illustrating the dependence on weighting. The archived
[reproduction report](../data/hirose-2008-gold-reproduction.json) retains both
diagnostics, with `full_dataset_recovered: true` and
`full_fit_reproduced: false`. Its residual-scaled linearized covariance is
diagnostic, not source covariance. Exact author-fit parity remains unclaimed,
and published coefficients and source uncertainty conventions are unchanged.

## Validity and calibration

Fit #1 is explicitly useful only above 70 GPa and is based on 1340–2330 K,
78.2–119.1 GPa hot observations. Fit #2 incorporates ambient expansivity and
9–26 GPa multi-anvil data; the authors recommend it to 140 GPa and 3000 K,
although the highest measured temperature is 2330 K. Stored marginal ranges
are not rectangular phase-stability guarantees.

The source reports pressures on the Speziale et al. (2001) MgO scale. The
300 K Au record links to executable `mgo_speziale_2001_bm3_2`, and the
explicit `gold_hirose_2008_to_mgo_speziale_2001_300k` edge records that
calibration ancestry in its room-temperature scope. That MgO record now also
bundles the published variable-q Debye thermal EOS. Its gamma law and
Debye-temperature integration prescription are validated and
[author-confirmed](speziale-2001-mgo.md); the original calculator, integration
grid and complete thermal-fit construction remain unrecovered.

The executable published Speziale model does not establish the historical
Hirose pressure reducer. The paired hot MgO states retain discrepancies with
the published model, and the exact later reduction convention remains unresolved
(see the [Speziale pressure-reduction audit](speziale-2001-mgo.md#diagnosis-of-feihirose-pressure-discrepancies)).
The `gold_hirose_2008_to_mgo_speziale_2001_thermal` ancestry edge therefore
remains evidence-only and `executable: false`; the unresolved historical
reduction, rather than absence of a bundled thermal model, is the limitation.
Source-reported pressures and paired lattice observations are preserved.

## Complete pressure-reduction check and Au recalibration sensitivity

The [complete-pair audit](../data/gold-speziale-reduction-audit.json), checked
on 12 October 2026, covers all 26 Fei (2004) Table 1 pairs, all 21 Hirose
(2008) Table 1 pairs, and the three paired Hirose (2006) rows. The other ten
Hirose (2006) rows have no MgO observations and cannot identify its reduction.
The 1340/2330 K heating states are repeated in the 2008 table, so the total
is **48 distinct observations: 38 hot and ten at 300 K**, rather than 50
independent pairs. No calculated Au pressures are used as fit targets. The
[residual table](../data/gold-speziale-reduction-residuals.csv) preserves source
pressure widths and predictions for all 30 tested gamma/theta/reference-temperature
combinations; inputs and evaluator scripts are hashed in the report.

| MgO convention | Fei hot RMS (26 rows, GPa) | Hirose 2008 hot RMS (12 rows, GPa) | All distinct hot RMS (38 rows, GPa) | Maximum hot difference (GPa) |
|---|---:|---:|---:|---:|
| Published Eq. 11 gamma, integrated theta, Tr=300 K | 0.254647 | 1.389455 | 0.808719 | 1.952128 |
| Local-power gamma, integrated theta, Tr=300 K | 0.008997 | 0.027977 | 0.017394 | 0.053917 |
| Local-power gamma, power-law theta, Tr=300 K | 0.007366 | 0.023209 | 0.014396 | 0.048887 |

Here local-power gamma means gamma=gamma0*x^q(x), and power-law theta means
theta=theta0*x^(-gamma(x)), with x=V/V0. Both are explicitly diagnostic
alternatives. The local-power gamma is not Speziale Equation 11, and this
power-law theta is not its thermodynamically consistent integral. Nevertheless,
the power-law combination matches all 38 hot pressures within the intervals
allowed by printed lattice and pressure last digits at fixed temperature and
coefficients. Local-power gamma with integrated theta and Tr=298 K also matches
all 38 within those intervals. These deterministic intervals are not measurement
uncertainties or confidence tests. The expanded comparison strongly supports a
local-power gamma convention, but **does not identify the theta prescription or
reference temperature**. Using the published gamma with alternative theta laws
does not achieve comparable agreement in the tested cases.

At 300 K, every Tr=300 K convention gives the same cold BM3 pressure; those
ten observations cannot distinguish gamma or theta. Their RMS difference is
0.027965 GPa, with maximum 0.045335 GPa. The nine Hirose (2008) differences
are below 0.035 GPa and within their printed pressure widths; the additional
Hirose (2006) cold row has no printed width. Rounded reference coefficients and
historical numerical settings remain possible contributors to these small offsets.

**A recalibration onto the published, author-confirmed Speziale formulation is
numerically consequential and should be identified separately from the published
Au EOS.** Re-evaluating the paired MgO measurements lowers the Fei hot pressure
targets by 0.011–0.362 GPa and the Hirose hot targets by 0.891–1.952 GPa.
Diagnostic unweighted Au refits have therefore been calculated and archived:

| Au diagnostic | Refit to reported pressures | Refit to published-Speziale pressures | Maximum curve shift at the fitted observations |
|---|---|---|---:|
| Hirose fit #2, 38 points | MAE=0.352703 GPa | MAE=0.315781 GPa | 1.862263 GPa lower |
| Hirose fit #1, twelve points | MAE=0.285198 GPa | MAE=0.276620 GPa | 1.901255 GPa lower |

These curve shifts compare two diagnostic refits, not an author fit against a
newly established absolute pressure standard. Hirose thermal refits keep the
published K0-prime=5.58 and the applicable ambient constraints fixed; the
separate 300 K sensitivity is not propagated into them. **Fei EOS work remains
paused pending the author's reply; no Fei-specific Au refits are retained in
this follow-up.** The 26 Fei (2004) measurements remain source inputs to
Hirose's published 38-point selection. The report retains local-power-target
controls for Hirose. Using local-power integrated-theta targets changes its fit #2
curve by at most 0.033075 GPa at these observations.

The published Au records remain historical source results. The separate
registered refit below specifies its calibration, objective and uncertainty
limits. Missing full covariance and the original later reducer remain evidence gaps.
This audit neither replaces source pressures nor establishes that recalculated
values are exact absolute pressures. Other gold EOS records on independent scales
are not implicated by this reduction discrepancy.

Reproduce with `python -m scripts.audit_gold_speziale_reduction`; use `--check`
to verify the saved report and residual table without rewriting them. No shared
material, calibration, ledger or manifest files are regenerated by this script.

## Registered Hirose–Fei refit on the published Speziale calculation

`gold_hirose_2008_fei_2004_speziale_2001_bm3_refit` is an optional **Peritheos refit of
Hirose (2008) and Fei (2004) observations**, on the author-confirmed published
Speziale MgO formulation. It is not a new author result and does not refit Fei's
own EOS. The published Au records and default are preserved.

The [derived dataset](../../peritheos/data/datasets/gold-hirose-2008-fei-2004-speziale-2001-recalculated-pvt.csv)
contains 47 rows with original pressure targets and widths in separate columns.
All new targets are calculated from the paired MgO volumes and source temperatures.
The nine 300 K Hirose rows first determine K0-prime with V0=67.85 Å³ and
K0=167 GPa fixed. That K0-prime is then fixed for the 38-row thermal regression:
all 26 Fei rows plus twelve hot Hirose rows. Ambient alpha0=3.824×10⁻⁵ and
alpha1=1.499×10⁻⁸ remain fixed. Both stages use equal pressure-residual weights.

| Refit coefficient | Value |
|---|---:|
| K0-prime at 300 K | 5.58115032 |
| beta1 (GPa⁻¹ K⁻¹) | 1.09175271×10⁻⁶ |
| beta2 (GPa⁻¹ K⁻²) | 3.96648026×10⁻¹⁰ |
| beta3 (GPa⁻¹ K⁻³) | 1.65378308×10⁻¹³ |
| a (K⁻¹) | 3.58321693×10⁻⁴ |

Thermal RMSE is **0.396463 GPa** and MAE **0.315847 GPa**. This differs
slightly from the fixed-published-Kprime diagnostic above because the cold
anchor is now refitted on recalculated pressures. Cold RMSE is 0.664115 GPa.
The native evaluator, independent equations and Python fallback agree, and the
new calibration ancestry explicitly links to `mgo_speziale_2001_bm3_2`.

The [refit report](../data/hirose-2008-fei-2004-gold-speziale-2001-refit.json) retains
stage-conditional residual-scaled curvature widths and covariance. They exclude
MgO calibration systematics, coordinate uncertainty, source-fixed anchors and
cold-to-thermal covariance propagation. Catalog parameter errors and covariance
are therefore unset. Thermal observations cover 8.629–117.445 GPa and
1273–2330 K; cold observations extend to 140.545 GPa. These are marginal coverage
ranges, not a validated P-T rectangle or a recommendation to extrapolate to 3000 K.
No independent holdout validation or absolute-pressure accuracy is claimed.

Reproduce and check the registered coefficients with
`python -m scripts.refit_hirose_fei_gold_speziale --check`. Registration uses
`--register` and inserts or updates only this record and its derived dataset;
shared global ledgers and manifests require later reconciliation.

The new gold pressure standards can be compared to any other executable gold
EOS through the shared Au volume coordinate, including at high temperature:

```python
from peritheos import get_eos_record, recalculate_xrd_pressure_scale

au = get_eos_record("gold_hirose_2008_bm3_fit2")
pressure = au.pressure(volume=51.28, temperature=2070.0)
comparison = recalculate_xrd_pressure_scale(
    pressure,
    "gold_hirose_2008_bm3_fit2",
    "gold_fei_2007_vinet_2",
    temperature_k=2070.0,
)
print(comparison.target_pressure_gpa)
```

Paired Au/MgO observations are available for future observation-level
recalculation. An executable scale ancestry edge represents the shared pressure
basis; it does not equate Au and MgO volumes or eliminate fit residuals.
