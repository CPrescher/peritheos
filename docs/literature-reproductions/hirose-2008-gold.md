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

| Diagnostic | Published curve RMSE (GPa) | Unweighted refit |
|---|---:|---|
| 300 K, all nine rows | 0.670655 | K0′=5.583345 |
| Fit #1, all twelve hot rows | 0.459525 | dK/dT=−0.0283325, α0=3.176726×10⁻⁵, α1=1.499813×10⁻⁸ |
| Fit #2, twelve local hot rows | 0.635474 | Full refit unavailable |

Both diagnostic refits recover coefficients within the source errors. They
retain the paper's staged selection and fixed parameters, but the source
does not specify residual coordinate, weights, or covariance, so exact solver
parity is not claimed and the published records are not replaced.

Fit #2 used 38 hot observations, including 26 from Fei et al. (2004),
[doi:10.1016/j.pepi.2003.09.018](https://doi.org/10.1016/j.pepi.2003.09.018).
Those 26 rows are not printed in Hirose, were not found in the searched Zotero
library, and the publisher endpoint returned HTTP 403 during this audit.
The Hirose article contains no separate supplementary observation table.
Consequently, the paper's overall 0.35 GPa mean absolute residual is not
claimed reproduced from the twelve available hot points. Table 1's full
calculated-output reproduction supports the executable parameterization.

## Validity and calibration

Fit #1 is explicitly useful only above 70 GPa and is based on 1340–2330 K,
78.2–119.1 GPa hot observations. Fit #2 incorporates ambient expansivity and
9–26 GPa multi-anvil data; the authors recommend it to 140 GPa and 3000 K,
although the highest measured temperature is 2330 K. Stored marginal ranges
are not rectangular phase-stability guarantees.

All pressures were reduced on the Speziale et al. (2001) MgO scale. The
300 K Au record links to executable `mgo_speziale_2001_bm3_2`, and the
explicit `gold_hirose_2008_to_mgo_speziale_2001_300k` edge records that
calibration ancestry. The thermal Speziale scale is not bundled: the
`gold_hirose_2008_to_mgo_speziale_2001_thermal` edge is evidence-only and
`executable: false`. A thermal Au record must not point to the 300 K-only MgO
record as its thermal anchor.

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
