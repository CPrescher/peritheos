# Grossular: full-inventory fit with estimated errors

**Subsequent controlled study:** the [normalization and Isaak specimen tests](grossular-angel-2022-normalization-specimen.md)
show that normalization alternatives can move θ across the source value, while
K′ remains high in every tested combination. This page retains the earlier
estimated-error sensitivity with normalization held fixed.

User-requested diagnostic and sensitivity study, 12 November 2026. This extends
the [193-row normalization audit](grossular-angel-2022-remaining-channels.md)
by assigning explicitly estimated errors to the remaining **19 Isaak moduli
and seven digitized Thiéblot volumes**. No original error cell or observation
is overwritten. The published EOS and every earlier receipt remain unchanged.

**All 219 non-anchor rows can now be fitted, but the tested error assumptions
do not reproduce Ross's central K′ or θ values.** The predeclared baseline is
K₀ = 167.016 GPa, K′ = 5.23585, θ = 779.754 K and γ = 1.15783. Across all
29 declared cases, θ stays at 773.477–781.417 K and K′ at 5.15975–5.29031.
The published values, 749.825 K and 5.06534, lie outside those ranges.
These ranges describe sensitivity to chosen weights, not confidence intervals
or a demonstration of failure of the published physical parameterization.

## What “full inventory” means here

The earlier 193 rows and all their assumptions are retained exactly. We add
the 26 previously unweighted observations, giving **219 fitted rows** and
**six free parameters**. The nominal S1 total is 220; the remaining Bosenick
datum is a normalization anchor, so it is not added as an independent residual.
The exact author residual count and masks remain unknown.

The added Isaak values are tabulated GGD moduli. The added Thiéblot values are
already processed, digitized Figure 4 points, not newly recovered original
measurements. The existing GR-2B expansion remains a mixed-composition
source-processing replay, not a measured pure-grossular volume series.
“Full” therefore describes coverage of the nominal selected inventory at
these mixed evidence levels; it does not mean author-workflow reproduction.

## Error model declared before solving

The [plan](../data/garnet-angel-2022/grossular/full-sensitivity-plan.json)
was saved before fitting. It defines the baseline, all sensitivity cases,
starts, bounds and stopping rule. No case is selected for closeness to the
published parameters.

**Nineteen Isaak modulus errors.** The baseline linearly interpolates the
three representative source widths: 0.71 GPa at 300 K, 1.02 GPa at 800 K
and 1.20 GPa at 1300 K. The final segment is linearly extrapolated to
1350 K, giving an estimated width of **1.218 GPa** there. The original three
printed-error rows and their weights remain unchanged. Numerical temperature
errors remain missing and omitted. These interpolated widths are newly
estimated weights, not recovered source errors.

**Seven Thiéblot points.** We hypothesize an experimental volume width of
**0.03 cm³/mol** and a temperature width of **1% of T**. The volume scale is
approximately the propagated volume width for a 0.001 Å cell-edge error near
11.85 Å, and ten times the digitization bound. The temperature percentage is
a declared hypothesis for the unknown power-calibration error. Neither is a
reported uncertainty from Thiéblot's paper.

Separately, we assume uniform digitization errors within the recorded bounds
±0.003 cm³/mol and ±2 K, whose variances under this new assumption are the
squared bounds divided by three. We combine those variances with the
hypothetical experimental widths:

\[
s_V^2=(m_V\,0.03)^2+0.003^2/3,\qquad
s_T^2=(m_T\,0.01T)^2+2^2/3.
\]

Thus the digitization bounds are not silently reinterpreted as reported
standard deviations. Shared axis and calibration errors remain unmodeled;
the diagonal approximation is an explicit limitation.

The **27-case factorial grid** independently varies the estimated Isaak Kₛ,
Thiéblot volume and Thiéblot temperature widths by **0.5×, 1× and 2×**.
The uniform digitization terms are held fixed. Two further controls replace
all 19 missing Isaak widths by a constant **0.71 GPa** or **1.20 GPa**, with
the other multipliers at one. These constants test the interpolation shape;
they are not rigorous bounds on the original individual errors.

All 29 cases use the same two declared numerical starts, equations, six free
parameters, normalizations and coordinate-error objective as the 193-row fit.
No additional weights, models, row masks or normalizations were searched.

## Baseline and sensitivity envelope

| Parameter | Published source | Earlier 193 rows | Baseline 219 rows | Range across all 29 cases |
| --- | ---: | ---: | ---: | ---: |
| K₀T, GPa | 167.002 | 166.754 | 167.016 | 166.870–167.147 |
| K′₀T | 5.06534 | 5.33063 | 5.23585 | 5.15975–5.29031 |
| θ₀, K | 749.825 | 780.407 | 779.754 | 773.477–781.417 |
| γ₀ | 1.15602 | 1.15522 | 1.15783 | 1.15540–1.15925 |
| Kono calculated-modulus multiplier | 0.99455 | 0.985270 | 0.987698 | 0.986272–0.990705 |
| Gwanmesia calculated-modulus multiplier | 1.00551 | 1.00127 | 1.00181 | 1.00138–1.00346 |

Baseline relative differences from the source are **+0.00852% in K₀**,
**+3.366% in K′**, **+3.991% in θ** and **+0.1569% in γ**. The K₀ and γ
sensitivity envelopes include their source central values; the K′, θ and
modulus-scale envelopes do not. Parameterwise envelope inclusion does not
establish simultaneous six-parameter agreement.

![Thermal-parameter sensitivity to estimated errors, light theme](../data/garnet-angel-2022/grossular/full-error-sensitivity.svg#only-light)
![Thermal-parameter sensitivity to estimated errors, dark theme](../data/garnet-angel-2022/grossular/full-error-sensitivity-dark.svg#only-dark)

The figure shows only the 27 factorial-grid cases. Columns change the estimated
Isaak modulus widths; each panel independently varies the estimated Thiéblot
volume and temperature widths. The pink outline marks the predeclared unit-
multiplier baseline. Every θ cell remains above the published 749.825 K.
The two constant-Isaak controls are included in the table envelope, not the
figure. Colours indicate conditional point estimates; exact values are in the
case summary CSV.

### Which assumptions matter most?

Holding the other estimated widths at their baseline values gives:

| Varied estimated width | K′ range | θ range, K | Interpretation within this design |
| --- | ---: | ---: | --- |
| Isaak Kₛ, 0.5×–2× | 5.16234–5.28747 | 779.108–779.940 | Strongest effect on K′ and modulus scale factors |
| Thiéblot V, 0.5×–2× | 5.23542–5.23602 | 776.441–781.196 | Strongest single-factor effect on θ |
| Thiéblot T, 0.5×–2× | 5.23585–5.23591 | 779.544–780.373 | Smaller effect at baseline volume width |

The combined grid gives a wider θ span because the temperature and volume
errors enter the same effective pressure variance. The largest θ change over
all declared assumptions is about **7.94 K**, versus a baseline offset of
**29.93 K** from the source. This supports the limited conclusion that these
tested error assumptions do not explain the remaining central-value difference.
It does not rule out different original normalization, covariance, specimen
treatment, selection or errors outside the declared hypotheses.

The two constant-Isaak controls give:

| Estimated missing Kₛ widths | K₀, GPa | K′ | θ, K | γ |
| --- | ---: | ---: | ---: | ---: |
| Constant 0.71 GPa | 167.147 | 5.17920 | 780.233 | 1.15925 |
| Constant 1.20 GPa | 167.078 | 5.22918 | 779.940 | 1.15840 |

## Convergence, residuals and uncertainty scope

**All 58 optimizations converge without hitting bounds.** Starts agree to
better than 3.1 × 10⁻⁸ in the declared parameter scaling. The baseline reduced
weighted objective is **4.42911**, with **213 nominal degrees of freedom**.
Estimated widths and omitted covariance prevent a calibrated χ² probability
interpretation. Reduced objectives from different weight cases are not used
to select a preferred error model.

Conditional local inverse-Jacobian widths are retained for every case,
without reduced-objective rescaling. They quantify local sensitivity under
each assumed diagonal error model; they are not recovered published errors.
The cross-case envelopes are also not statistical confidence intervals.

The baseline forward RMS residuals for the newly included observations are
**0.8704 GPa for the 19 estimated-error Isaak moduli** and
**0.08515 cm³/mol for the seven Thiéblot volumes**. Their unchanged-source
comparisons give 0.6171 GPa and 0.08164 cm³/mol, respectively. These are
unweighted RMS values, not separate fitted objectives. Both curves use the
same fitted modulus multipliers; source scale factors are not independently
reoptimized here. A joint compromise need not improve every dataset's RMS.

## Reproducible artifacts and remaining limitations

- [Baseline 219 inputs and error estimates](../data/garnet-angel-2022/grossular/full-baseline-inputs.csv): explicitly marks the added estimated-error rows and digitized observations.
- [29-case summary](../data/garnet-angel-2022/grossular/full-sensitivity-summary.csv): parameters, factors and convergence.
- [Complete sensitivity receipt](../data/garnet-angel-2022/grossular/full-sensitivity-refits.json): both starts, covariance, traces, plan/input hashes and parameter envelopes.
- [6,351 row-level pressure residuals](../data/garnet-angel-2022/grossular/full-sensitivity-residuals.csv): one 219-row table per case, including effective weights.
- [Baseline source/fit forward comparisons](../data/garnet-angel-2022/grossular/full-baseline-native-comparisons.csv): observed-minus-calculated V or Kₛ, in cm³/mol or GPa, on identical rows. These differ from the minimized pressure-coordinate residuals.

The original Thiéblot table, original assigned errors, exact author masks and
normalization/preliminary-EOS inputs remain unrecovered. The earlier GR-2B
specimen issue, Rodehorst 156 K discrepancy and shared anchor/calibration
covariance qualifications remain. No author was contacted, no catalog record
was written, and no diagnostic is recommended for scientific prediction.

Tests verify selection, the 26 estimated-error assignments, preservation of
all original blank cells and earlier artifact hashes, the complete 29-case
design, both-start convergence and covariance rank. They replay the baseline
and a grid corner, independently check all 219 states for source and fitted
coefficients, and reconcile the exported residuals with every objective.
See the [verification receipt](../data/garnet-angel-2022/grossular/full-sensitivity-verification.json).

Repository-only commands:

```bash
python -m scripts.refit_grossular_angel_2022_full_sensitivity
python -m scripts.plot_grossular_angel_2022_sensitivity
pytest -q tests/test_grossular_full_sensitivity.py
```
