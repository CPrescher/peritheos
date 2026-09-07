# Zhang et al. (2025) hcp-Fe thermal EOS

This reproduction audits the three executable hcp-Fe records taken from Fits 1,
2, and 5 of Zhang et al., *Equation of State Parameters of hcp-Fe Up to
Super-Earth Interior Conditions*, Crystals 15, 221 (2025),
[doi:10.3390/cryst15030221](https://doi.org/10.3390/cryst15030221).

## Primary-source evidence

The version-of-record PDF and official supplementary ZIP were retrieved from
MDPI's static article assets. The article and data are licensed CC BY 4.0. The
derived CSV records these source checksums:

| artifact | SHA-256 |
| --- | --- |
| version-of-record PDF | `f7da2f0446bd76d35b25f363e7f4ee4266db130636b7ef907ddd459af9cc1fe8` |
| official supplementary ZIP | `518ed6f4a17e5dee22b69bfe395fe160e81bf03ea3e9090a8573e251a77ab5a8` |
| `crystals-3484051-supplementary.xlsx` | `a3654047a5cabf03633c9ce280a264fa2c7bfae1c86d4e4c04028e715a0cd4a2` |

The supplement is decisive: Tables S1, S3, and S4 contain every selected P-V-T
row and residuals for all five source fits; S2 contains the self-calibrated
static rows excluded by the authors; S5 contains the variance-covariance
matrices; and S6 gives the precise RMSE values and fit envelope.

## Exact fit selection

The source regression contains 1,313 rows:

| table | rows | observation provenance | temperature provenance |
| --- | ---: | --- | --- |
| S1 | 1,078 | static-compression experiments | reported experimental temperature |
| S3 | 132 | 72 Hugoniot and 60 ramp experimental P-V states | calculated/model-assigned, not measured |
| S4 | 103 | ab-initio calculations | calculated |

All 325 Table S2 rows were excluded because pressure was determined from hcp-Fe
itself. No S2 row is present in the bundled fit dataset. One S1 volume
(`Table S1!D461`) is stored as numeric text; the extractor coerces that exact
cell and records the anomaly instead of silently dropping the row.

The article's text says the static data extend to 3400 K, but S1 contains fitted
rows through 4695 K and supplies fit residuals for them. S6 reports the actual
combined envelope as 10.7-1374 GPa and 298-12000 K. The executable records use
the deposited envelope and document the textual discrepancy.

## Calibration trace

For the 1,078 S1 rows, the fit uses the workbook's “P after correction” value
and corresponding error when populated:

- 499 rows name Fei et al. (2007),
  [doi:10.1073/pnas.0609013104](https://doi.org/10.1073/pnas.0609013104).
- 142 rows name Ye et al. (2018),
  [doi:10.1080/08957959.2018.1493477](https://doi.org/10.1080/08957959.2018.1493477).
- 437 rows retain the compiled pressure. The sheet says these sources either
  already use an internally consistent scale or lack a calibrant developed on
  one.

Methods section 2.1 instead cites Fei et al. (2007) and Dorfman et al. (2012) as
the consistent scales. The row-level workbook attribution is preserved rather
than resolving that conflict by assumption. The final 641 corrected pressures
are therefore source-backed and directly reusable, but their upstream
recalculation is not independently executable: the supplement omits row-wise
calibrant volumes and temperatures. The records consequently use
`partially_resolved` calibration with
`missing_calibrant_observations`.

## Dynamic temperatures are derived inputs

Table S3 calls its temperatures calculated from Zhuang et al. (2021),
[doi:10.1103/PhysRevB.103.144102](https://doi.org/10.1103/PhysRevB.103.144102).
That study obtains the Hugoniot from the Rankine-Hugoniot energy condition and
the ramp path by entropy integration; it does not publish a machine-readable
P-T series or Zhang's interpolation coefficients. In Zhang's workbook, the 62
older Hugoniot temperatures follow a pressure quadratic exactly and the 60
Smith ramp temperatures follow a pressure cubic exactly, but those interpolation
rules are not documented in the paper or supplement. The dataset metadata and
tests record these relations as inferred integrity checks: `T = 0.01491 P^2 +
21.77176 P - 395.42764` for the 62 Hugoniot rows and `T = 4.22e-7 P^3 -
0.00155 P^2 + 2.58948 P + 920.817` for the 60 ramp rows. They reconstruct the
deposited temperatures within `4.4e-12 K` but are not attributed to the authors
as a published method.

The final ten Huang et al. (2022) laser-shock rows are also calculated
temperatures—the Huang article explicitly describes the matching 1073-3919 K
values as calculated—not direct temperature measurements. The CSV labels these
separately as `model_calculated_huang_2022`. Thus experimental dynamic P-V
observations and theory/model-generated temperatures remain distinguishable.

## Model and units

Fits 1 and 2 use a third-order Birch-Murnaghan 300 K reference isotherm; Fit 5
uses Vinet. All three add the source's Mie-Grüneisen-Debye pressure:

`P(V,T) = P_300K(V) + gamma(V) [E(V,T) - E(V,300 K)] / V`,

with `gamma = gamma0 (V/V0)^q` and
`theta = theta0 exp[(gamma0-gamma)/q]`. The independent reproduction uses the
paper's stated `R = 8.314 J mol^-1 K^-1`, `n = 1`, and unweighted pressure
residuals. Source volumes are cm3/mol Fe. They are converted to the conventional
two-atom hcp cell using `2 / 0.602214076`; the same factor is applied to the V0
row and column of each S5 covariance matrix.

The Peritheos runtime uses the modern gas constant rather than the rounded
source value. Across all selected rows this changes published-model pressures
by at most 0.0066 GPa, far below the roughly 4.5-5.1 GPa fit RMSE.

## Independent numerical reproduction

The dedicated script evaluates the equations independently of Peritheos's EOS
classes, reconstructs the deposited residual columns, and runs an unweighted
nonlinear least-squares refit:

| fit | source coefficients `(V0, K0, K0', theta0, gamma0, q)` | independent refit | source RMSE (GPa) | refit RMSE (GPa) |
| --- | --- | --- | ---: | ---: |
| 1 | `(6.756, 174.7, 4.790, 1209, 2.86, 0.84)` | `(6.755761, 174.668655, 4.790198, 1207.3901, 2.857148, 0.837205)` | 4.504302 | 4.504036 |
| 2 | `(6.753 fixed, 175.1, 4.787, 1205, 2.86, 0.84)` | `(6.753 fixed, 175.130777, 4.787028, 1201.5192, 2.856710, 0.840197)` | 4.504408 | 4.504165 |
| 5 | `(6.753 fixed, 151.6, 5.845, 960, 3.51, 1.28)` | `(6.753 fixed, 151.565874, 5.844690, 962.3478, 3.508689, 1.276369)` | 5.084659 | 5.084341 |

The deposited residuals are reconstructed to at most `1.45e-6 GPa`. Covariance
matrices derived from the unweighted Jacobian and `SSE/(N-p)` agree with S5 to
within 0.6% in relative Frobenius norm. The small coefficient/RMSE differences
are consistent with fitting from the rounded published rows and coefficients.

Run the reproducibility checks with:

```text
python scripts/reproduce_zhang_2025_hcp_iron.py --check
pytest -q tests/test_zhang_2025_hcp_iron.py
```

## Reproducibility boundary

The EOS regression is independently reproducible from the deposited final fit
inputs, including source residual and covariance parity. Two upstream reductions
are not fully reproducible from public numerical inputs: the 641 static-pressure
recalculations and the interpolation/generation of most Table S3 temperatures.
Those limitations do not require inventing fit observations, but they prevent a
claim that the entire experimental reduction pipeline is independently
reconstructed. Each record therefore exposes the machine-readable scope status
`final_input_parity_upstream_reduction_partial`; the paper ledger displays that
qualification alongside its coefficient-parity result.
