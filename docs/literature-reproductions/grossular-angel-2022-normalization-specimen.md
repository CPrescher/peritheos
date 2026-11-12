# Grossular: normalization and Isaak specimen tests

User-requested controlled diagnostics, 12 November 2026. This follows the
[219-row estimated-error study](grossular-angel-2022-full-sensitivity.md).

**Volume normalization can shift θ by enough to account for its central-value
discrepancy, but this does not identify Ross's original workflow.** The Isaak
specimen treatments alone move θ by a few kelvin under the original plotted
normalization. **No tested combination reaches the published central K′.**
All 40 optimizations, two starts for each of 20 declared cases, converge
without hitting a bound. Estimated errors are held at the preceding baseline;
we do not search weights to improve agreement.

## Separate effects under controlled comparisons

| Change from the 219-row baseline | Rows | K′ | θ, K | θ change, K |
| --- | ---: | ---: | ---: | ---: |
| Published source, for reference | — | 5.06534 | 749.825 | — |
| Baseline plotted processing | 219 | 5.23585 | 779.754 | 0 |
| Isaak GR-2B strain re-referenced to 298 K | 219 | 5.23779 | 780.048 | +0.294 |
| Omit only the 15 GR-2B expansion rows | 204 | 5.24009 | 777.200 | −2.554 |
| Calculated endmember-α curve in place of GR-2B strain | 219 | 5.24191 | 774.317 | −5.437 |
| Tabulated room-condition denominators | 219 | 5.23260 | 779.449 | −0.305 |
| Linearly derived 298 K denominators | 219 | 5.23304 | 763.163 | −16.591 |
| Free Milani/Du volume offsets | 219 | 5.19971 | 740.650 | −39.104 |
| Free all five diffraction-volume offsets | 219 | 5.20154 | 751.073 | −28.681 |

The last two cases add nuisance parameters. They are deliberate relaxation
controls, not recovered source normalization. Under the last control θ differs
from the source by only **+1.248 K**, while K′ still differs by **+0.13620**.
Across all 20 cases, K′ spans **5.19971–5.24191**, entirely above 5.06534.
These are central-value comparisons, not significance tests.

![Normalization and Isaak effects, light theme](../data/garnet-angel-2022/grossular/normalization-specimen-effects.svg#only-light)
![Normalization and Isaak effects, dark theme](../data/garnet-angel-2022/grossular/normalization-specimen-effects-dark.svg#only-dark)

Colour identifies normalization and marker shape identifies the Isaak treatment.
The star and dashed lines show the published central values, not uncertainty
limits. Normalization changes move θ substantially, but all cases remain to
the right of the published K′. Overlapping markers do not indicate identical
inputs; full precision is retained in the case table.

## What was fixed and what was changed

The [plan](../data/garnet-angel-2022/grossular/normalization-specimen-plan.json)
was saved before solving. It declares a **five-by-four factorial design**,
the same two starts for every case, bounds and a fixed compute cap.
All 22 Isaak GGD moduli, all modulus-error estimates, compression calibration,
thermal law, adopted reference volume and other input selections remain fixed.
Only the declared volume processing and the 15 Isaak expansion rows change.

The normalization alternatives are:

| Normalization | Milani denominator, Å³ | Du denominator, Å³ | Bosenick denominator, cm³/mol | Additional free volume scales |
| --- | ---: | ---: | ---: | --- |
| Plotted baseline | ≈1664.05 | ≈1664.26 | 125.293 | None |
| Specific tabulated room measurements | 1663.98 at 293 K | 1664.2 at 296.1 K | 125.293 at 295 K | None |
| Local straight line to 298 K | 1663.99875 | 1664.23504 | 125.30104 | None |
| Free Milani/Du offsets | Baseline processing | Baseline processing | 125.293 | Two |
| Free diffraction offsets | Baseline processing | Baseline processing | Baseline processing | Five: Milani, Du, Bosenick, Rodehorst and Thiéblot |

The linear case uses the nearest bracketing Milani 293/301 K observations,
the first Du 296.1/393.7 K observations, and a three-kelvin extrapolation of
the Bosenick 270/295 K pair. It uses no EOS coefficients. These are calculated
denominators, not new measured anchors or recovered author instructions.
Their displayed calculation precision does not establish experimental precision.
The interpolation also assumes the chosen near-room measurements are mutually
comparable; exact instrument/cycle subgroup processing remains unrecovered.

For fixed-denominator alternatives, both volume and its coordinate width are
multiplied by the same old/new denominator ratio. Temperatures remain unchanged.
Shared anchor errors and uncertainty of the reconstructed denominators remain
unpropagated, as in the reference diagnostic. Rodehorst keeps its actual 298 K
anchor; the already processed Thiéblot series has no new raw denominator.

The free-offset controls retain all baseline processed observations/errors,
and fit calculated-volume multipliers bounded at 0.98–1.02. They do not free
an Isaak expansion scale or any compression scale. The cases have six, eight
or eleven free parameters; omitting Isaak expansion reduces the observations
from 219 to 204. Reduced objectives therefore do not supply a common calibrated
probability or an automatic model-selection score.

## Why a free-scale θ match is not workflow reproduction

With GR-2B kept as plotted, freeing Milani/Du offsets gives multipliers
**1.0000924** and **0.9995208**, respectively. Since a scale multiplies the
calculated volume, its equivalent raw denominator is

\[
D_{\mathrm{equiv}}=D_{\mathrm{plotted}}\,s.
\]

Du's equivalent denominator becomes **1663.4629 Å³**, a change of
**−0.7976 Å³**. Its equivalent observed-minus-calculated offset near the
reference molar volume is about **−0.0601 cm³/mol**. This is about **20 times** the
0.04 Å³ resolution of the denominator recovered from the plotted series.
Milani changes by +0.1538 Å³, about four such resolution units.

Those ratios are not statistical sigma scores. They show that the scale
relaxation is substantially larger than a digitization-reading ambiguity.
It might describe a different calibration/reduction hypothesis, but it does
not reproduce the observed Figure 4 normalization. The five-offset control
also relaxes measured anchors; its near-source θ is not evidence that this
was the author's input processing.

The [offset receipt](../data/garnet-angel-2022/grossular/normalization-specimen-offsets.csv)
retains every fitted multiplier, equivalent denominator and offset separately
from the original observations.

## What the Isaak alternatives mean

[Isaak et al. (1992)](https://doi.org/10.1007/BF00198608), Table 2, gives
measured **GR-2B** strain and an **estimated endmember volumetric expansivity**.
These are different types of evidence. The GGD modulus specimen is unchanged
throughout this study.

The reference-only control linearly interpolates the defining strain ratio
1 at 293 K and the tabulated ratio at 300 K. It gives
**r(298)/r(293) = 1.000092857**, and divides the 15 volume values and their
widths by that ratio. This is a diagnostic reference correction, not a newly
recovered 298 K measurement. Its small θ effect under the plotted normalization
shows that the five-kelvin reference issue alone does not explain the 29.93 K
baseline central-value difference.

The calculated endmember control linearly interpolates the final Table 2 α
column, extrapolates the first 300/350 K segment to 298 K, and calculates

\[
V_{\mathrm{calc}}(T)=125.35\exp\left[\int_{298}^{T}\alpha_{\mathrm{endmember}}(t)\,dt\right].
\]

It is **a calculated curve from a compositional estimate**, not measured
endmember volume. The source estimate assumes Skinner's relative expansivity
difference between andradite and grossular. We introduce no new andradite EOS
or mixing law. Isaak estimates about 2% uncertainty in α, but its temperature
and compositional covariance are unavailable. Keeping the original assigned
0.01 cm³/mol and 2 K comparison weights isolates central-curve changes; it does
not propagate that α uncertainty into measured-volume errors.

Under plotted normalization, this control lowers θ by only 5.44 K and changes
γ from 1.15783 to **1.14407**, farther from the published 1.15602. It does not
resolve the source parameter set. The
[treatment-curve receipt](../data/garnet-angel-2022/grossular/isaak-specimen-treatment-curves.csv)
keeps measured-strain reductions and calculated curves in separate columns.

## Interactions matter

The following table gives θ for every declared combination, in kelvin:

| Normalization | GR-2B plotted | GR-2B referenced to 298 K | Omit GR-2B expansion | Calculated endmember α |
| --- | ---: | ---: | ---: | ---: |
| Plotted | 779.754 | 780.048 | 777.200 | 774.317 |
| Tabulated room anchors | 779.449 | 779.735 | 778.953 | 774.047 |
| Linear 298 K anchors | 763.163 | 763.329 | 760.220 | 757.492 |
| Free Milani/Du offsets | 740.650 | 740.955 | 746.549 | 735.875 |
| Free five diffraction offsets | 751.073 | 759.252 | 766.854 | 765.235 |

For example, omitting GR-2B expansion changes θ by **−2.55 K** under plotted
normalization, but **+15.78 K** with all diffraction offsets free. The latter
near-source θ therefore depends on retaining the mixed-composition expansion
series. Treating normalization and the specimen issue as independent additive
corrections would be misleading. Paired main effects and factorial interactions
are recorded for all six shared EOS/modulus-scale coefficients.

## Artifacts, verification and limits

- [20-case parameter table](../data/garnet-angel-2022/grossular/normalization-specimen-summary.csv).
- [Full fit receipt](../data/garnet-angel-2022/grossular/normalization-specimen-refits.json), with both starts, nuisance parameters, covariance, normalization constants and paired effects.
- [4,305 pressure residuals and effective errors](../data/garnet-angel-2022/grossular/normalization-specimen-residuals.csv), one row table per primary-start case.
- [Source/fit native V and Kₛ comparisons](../data/garnet-angel-2022/grossular/normalization-specimen-native-comparisons.csv), with calculated endmember curves explicitly marked. Residual sign is observed/diagnostic-input minus calculated, in cm³/mol or GPa; both curves share the case's fitted nuisance factors.

Tests verify the controlled changes, all 22 GGD moduli and unchanged error
background, denominator calculations, exact integration of the initial linear
α segment, baseline replay, both-start agreement and covariance rank. They
replay the eleven-parameter calculated-curve control, compare its full set of
states with independent forward Debye quadrature, and reconcile all exported
objectives and paired effects. All historical source and audit hashes remain
unchanged. See the [verification receipt](../data/garnet-angel-2022/grossular/normalization-specimen-verification.json).

These results establish **strong θ sensitivity to volume processing** and
**persistence of the K′ central-value discrepancy under the tested hypotheses**.
They do not recover the original normalization, resolve all covariance or
composition assumptions, or establish a replacement EOS. No catalog record
or published coefficient was changed, and no author was contacted.

Repository-only commands:

```bash
python -m scripts.refit_grossular_angel_2022_normalization_specimen
python -m scripts.plot_grossular_angel_2022_normalization_specimen
pytest -q tests/test_grossular_normalization_specimen.py
```
