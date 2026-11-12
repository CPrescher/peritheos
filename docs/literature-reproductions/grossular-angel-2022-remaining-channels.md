# Grossular: the remaining 73 channels and normalization

Investigation: 12 November 2026. This follows the
[147-channel diagnostic](grossular-angel-2022-expanded-fit.md).
**All 73 excluded channels now have an explicit disposition.** We recovered
additional error information and reconstructed plotted normalization, enabling
a **193-row source-processing diagnostic**. A faithful full 220-channel fit
remains unavailable: 26 observations still lack usable individual errors,
and one nominal channel is a normalization anchor.

The larger diagnostic **does not converge toward all of Ross's coefficients**.
K₀ and γ become closer, but K′ and θ move farther away. Selection and
normalization both changed, so this is not an isolated test of adding rows.
No catalog record, published coefficient or scientific classification changed.
No author was contacted; the pyrope request remains waiting for a reply.

## Disposition of the 73 channels

| Dataset | Nominal channels | Newly established evidence | Treatment in the diagnostic |
| --- | ---: | --- | --- |
| Isaak expansion | 15 | S1 assigns 2 K and 0.01 cm³/mol or better; Figure 4 matches scaled GR-2B strain with a 293 K reference | Include 15 solely as a replay of the plotted reduction; mixed composition explicitly retained |
| Bosenick–Geiger expansion | 20 | Primary PDF recovered; all table values/errors visually checked; temperature accuracy ±3 K down to 100 K and ±5 K below | Include 19 non-anchor rows; 295 K measurement is the normalization anchor |
| Thiéblot expansion | 7 | S1 specifies preliminary-EOS 300/298 K scaling; seven processed points already digitized | Exclude: original measurements, experimental errors and preliminary coefficients remain missing |
| Rodehorst expansion | 9 | Table 4A gives σ(a) = 0.001 Å and a(298 K) = 11.850 Å | Include nine non-anchor rows with propagated cell errors; temperature errors remain missing |
| Isaak moduli without printed errors | 19 | Primary text explicitly says errors are printed at representative temperatures only | Exclude: cannot reconstruct individual errors from the three printed examples |
| Conrad moduli | 3 | Cᵢⱼ errors bound σ(Kₛ) even with unknown covariance; S1 assigns σ(P) = 0.1 GPa | Include three with a conservative covariance-compatible upper error bound |
| **Total** | **73** | **46 additional diagnostic rows, 26 still unweighted, one anchor** | **147 + 46 = 193 rows** |

The [channel audit](../data/garnet-angel-2022/grossular/remaining-channel-audit.json)
records evidence, missing inputs and assumptions separately. These counts do
not establish the author's exact residual count: the S1 ambient-anchor footnote
still conflicts with several original/used counts.

## Bosenick: primary errors and anchor recovered

The final [Bosenick–Geiger (1997) publication](https://doi.org/10.1029/97JB01612)
was downloaded through the publisher's ePDF reader. Table 2 on p.22651
confirms all 20 grossular volume/error tokens in the existing transcription.
The parenthetic volume errors derive from one-standard-deviation cell errors.
The source uses Nₐ = 6.02214 × 10²³ mol⁻¹; its tabulated molar volumes remain
unaltered, and that normalization constant cancels in the volume ratio.

Page 22650 gives temperature accuracy of **±3 K down to 100 K** and **±5 K
below 100 K**, with collection drift below 0.1 K. Accuracy is not declared to
be one standard deviation. Using it as a diagonal temperature weight is a
conditional choice in this diagnostic.

All **19 blue marker paths** in Angel Figure 4 match original rows 2–20 when
normalized by **Vₘ(295 K) = 125.293 cm³/mol**. With the fixed tabulated anchor,
differences remain below the 0.003 cm³/mol digitization bound. The 295 K
anchor is absent from those marker paths. This supports anchor omission;
it does not prove the exact fit-project mask or uniquely resolve S1's 20-used
statement.

The 1997 authors excluded cooling points below 150 K and their problematic
20 K datum from their own thermal-expansion fits, and averaged other repeats.
Angel's plotted series contains the separate measurements, including those
low-temperature points. This replay follows the plotted rows and preserves
the disagreement between these selection conventions.

## Volume normalization: recovered method and plotted constants

Angel (2022), p.3, explicitly gives the diffraction reduction:

\[
Y_i = 125.35\,\frac{V_i}{V_a}\quad\mathrm{cm^3/mol},
\]

where Vₐ is that dataset's own room-condition measurement. Refining a volume
scale is the stated alternative for datasets lacking such a reference. Our
earlier fits instead freely refined two volume scales, which differed from
this source processing.

The [new digitization receipt](../data/garnet-angel-2022/grossular/normalization-recovery.json)
and [77 matched marker coordinates](../data/garnet-angel-2022/grossular/figure4-normalization-points.csv)
identify the plotted normalization independently of the EOS coefficients:

| Series | Matched markers | Denominator inferred from the plotted points | Plot residual RMS, cm³/mol |
| --- | ---: | ---: | ---: |
| Bosenick | 19 | 125.293 cm³/mol, agreeing with measured 295 K anchor | 0.000805 |
| Milani | 35 | approximately 1664.05 Å³ | 0.000836 |
| Du | 13 | approximately 1664.26 Å³ | 0.000647 |

The Milani and Du denominator resolution is approximately **0.04 Å³**, from
the declared 0.003 cm³/mol digitization bound. This is a shared resolution
estimate, not an experimental standard deviation and not reduced by √N.
These constants reproduce plotted processing; the original anchor interpolation
or preliminary EOS calculation remains unidentified. The new diagnostic fixes
these reconstructed constants, rather than optimizing them to match Ross.

Rodehorst has a measured 298 K anchor, so we use
**125.35 × [a(T)/11.850]³**. Nine of the ten plotted points agree with the
literal table within the digitization bound. The **156 K point differs by
0.0331 cm³/mol**, equivalent to about 0.001 Å in cell edge. We preserve the
primary Table 4A value **11.842 Å**; it is not repaired to the plotted value.
The source input may therefore contain a transcription or plotting difference.

Shared anchor errors would induce covariance even if raw cell errors were
independent. For a non-anchor pair i,j, first-order propagation gives

\[
\operatorname{Cov}(Y_i,Y_j)=
\left(\frac{125.35}{V_a}\right)^2\delta_{ij}\sigma_i^2
+\frac{Y_iY_j}{V_a^2}\sigma_a^2.
\]

Original calibration correlations and the author's propagation convention are
unknown. The current diagonal diagnostic omits this shared covariance and
uncertainty in the plotted denominators. Its parameter widths are consequently
conditional, not reconstructed author uncertainty intervals.

## What can and cannot be recovered for the moduli

[Isaak (1992)](https://doi.org/10.1007/BF00198608), p.110, says errors in
Tables 3a/b are listed at representative temperatures and include density
uncertainty and modal-frequency standard deviations. GGD Table 3b prints Kₛ
widths of 0.71, 1.02 and 1.20 GPa at 300, 800 and 1300 K. It does not supply
the underlying mode-level observations/covariance needed for the remaining
19 individual errors. Polynomial fit coefficient errors are not independent
observation errors. Interpolation would be an assumption, not data recovery,
so those 19 rows remain outside the weighted fit.

[Conrad (1999)](https://doi.org/10.2138/am-1999-0321) supplies one-standard-
deviation C₁₁/C₁₂ errors, but not their covariance. Since
Kₛ = (C₁₁ + 2C₁₂)/3, any admissible correlation gives

\[
\frac{|\sigma_{11}-2\sigma_{12}|}{3}
\leq\sigma_K\leq
\frac{\sigma_{11}+2\sigma_{12}}{3}.
\]

The three selected rows have admissible intervals **0.010–0.030**,
**0.0567–0.0900** and **0.0200–0.0467 GPa**. We use the upper bound as a
conservative diagnostic width and S1's assigned 0.1 GPa pressure error.
That bound is not the recovered original Kₛ error. The two highest-pressure
outliers remain excluded; ambient temperature is explicitly adopted as 298 K.

## Declared 193-row fit and outcome

The [plan](../data/garnet-angel-2022/grossular/normalized-refit-plan.json) was
saved before solving. It keeps the same BM3/q-compromise equations and fixed
reference conventions. Four EOS coefficients and the two modulus scales are
free, matching the six nonzero covariance indices in the official EOS file.
The two volume scales from the earlier diagnostic are now fixed by the
normalization reconstruction. This does not establish identical author inputs.

We retain the earlier 147 coordinate-error conventions and add only the 46
rows described above. Isaak expansion uses **S1's assigned 2 K and 0.01 cm³/mol**
and the plotted GR-2B strain reduction. It is explicitly a mixed-composition
source-processing replay, not a physical fit to 15 measured pure-grossular
volumes. No compositional correction is invented.

Both declared numerical starts converge to the same solution without hitting
a bound. The conditional reduced weighted objective is **4.76095**, with
**187 nominal degrees of freedom**. It is not directly comparable as a
goodness-of-fit statistic with the 147-row result: both input rows and volume
normalization changed, and the weights lack complete covariance semantics.

| Parameter | Published source | Earlier 147 rows | Normalized 193 rows | 193-row relative difference from source |
| --- | ---: | ---: | ---: | ---: |
| K₀T, GPa | 167.002 | 167.353 | 166.754 | −0.149% |
| K′₀T | 5.06534 | 4.98267 | 5.33063 | +5.237% |
| θ₀, K | 749.825 | 740.183 | 780.407 | +4.079% |
| γ₀ | 1.15602 | 1.16467 | 1.15522 | −0.0688% |
| Kono modulus multiplier | 0.99455 | 0.996303 | 0.985270 | −0.933% |
| Gwanmesia modulus multiplier | 1.00551 | 1.00547 | 1.00127 | −0.422% |

Full precision, local covariance, optimizer traces, effective errors and all
193 fitted residuals are retained in the
[fit receipt](../data/garnet-angel-2022/grossular/normalized-refits.json).
The [forward comparison](../data/garnet-angel-2022/grossular/normalized-native-comparisons.csv)
contains original row identities, processed observations and both curves on
identical P/T coordinates. Residual sign is **observed minus calculated**,
in cm³/mol for volumes or GPa for moduli. Both curves use the same fitted
modulus multipliers to isolate coefficient effects; these source predictions
do not use the source's original modulus multipliers. The minimized pressure
residuals are retained separately in the fit receipt.

Gwanmesia's fit RMS becomes **0.911 GPa**, versus 0.678 GPa in the earlier
147-row diagnostic. Du's volume RMS remains **0.0585 cm³/mol** under the
fixed reconstructed normalization. The larger constrained input set exposes
tensions; more rows do not establish closer agreement. No alternative weights,
normalizations or models were searched after obtaining this result.

## Access, verification and next evidence needed

The [access receipt](../data/garnet-angel-2022/grossular/remaining-source-access.json)
records the recovered Bosenick PDF hash and failed download attempts. RRUFF's
human-confirmation gate had disappeared when revisited after user approval,
but the browser showed a blank PDF surface and no Rodehorst PDF was captured.
The publisher explicitly reports no access. Its indexed primary Table 4A
remains usable, but numerical temperature errors remain unrecovered.

A full author-workflow replay still requires **19 Isaak modulus errors**, the
**original seven Thiéblot values and experimental errors**, the original
normalization/preliminary-EOS inputs, resolution of the Rodehorst 156 K value,
and the exact fit masks and controls. S1's nominal 220 count is not forced
into 220 independent residuals by retaining an anchor or inventing errors.

Tests preserve every historical artifact hash, check the normalization and
Conrad covariance bounds, replay both earlier diagnostics and the new fit,
and compare all 46 added states with independent forward Debye quadrature.
See the [verification receipt](../data/garnet-angel-2022/grossular/normalized-verification.json).

These are repository-only maintenance tools:

```bash
# Optional PyMuPDF and a privately obtained final Angel PDF required:
python -m scripts.recover_grossular_angel_2022_normalization --pdf /path/to/article.pdf
python -m scripts.refit_grossular_angel_2022_normalized
pytest -q tests/test_grossular_normalization_recovery.py
```
