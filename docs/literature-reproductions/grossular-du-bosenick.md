# Grossular: Du processing and Bosenick anchor/covariance

User-requested investigation, 13 November 2026. This follows the
[θ influence study](grossular-angel-2022-theta-influence.md). The additional
independent **2% Conrad modulus-error allowance remains diagnostic and fixed**.
All 219 fit rows, including all three Conrad points and all 19 non-anchor
Bosenick rows, remain in every case. No published coefficient or source error
was changed.

**The tested preprocessing and covariance assumptions do not close the
remaining central θ discrepancy.** The two prior diagonal baselines replay.
Separating the Bosenick anchor moves θ by about 17.3 K, but propagating its
shared volume uncertainty moves θ upward toward 850 K. The tested Du changes
also leave θ high. All 32 cases, two starts each, converge without bounds.
Across the whole declared study θ is **763.879–895.616 K**, above the
published **749.825 K**. These are central-value comparisons, not significance
tests or proof that all possible original covariance assumptions are excluded.

## Main comparisons

The table uses the earlier **Milani/Du linear 298 K background**. The
Bosenick denominator is independently selected as its measured 295 K value
or calculated 298 K value; it no longer changes implicitly with the background.

| Bosenick anchor and error treatment | θ, K | K′ | Conditional unscaled θ width, K |
| --- | ---: | ---: | ---: |
| Measured 295 K, original diagonal weights | 781.850 | 5.06619 | — |
| Calculated 298 K, original diagonal weights | **764.605** | **5.06502** | 9.89 |
| Calculated 298 K, propagated volume-anchor covariance | **849.605** | **5.07157** | 14.97 |
| Calculated 298 K, volume plus independent raw-T propagation | 851.127 | 5.07172 | 15.04 |
| Calculated 298 K, volume plus fully shared standardized raw-T errors | 858.408 | 5.07222 | 15.53 |
| Volume + independent raw T, with Du a³/σₐ reconstruction | **821.662** | **5.06490** | 15.23 |

The last row changes Du's representation while using Bosenick's
volume-plus-independent-T covariance. The widths are conditional inverse-
information widths, without residual rescaling, under the adopted diagonal/raw
error assumptions. They are not established experimental confidence intervals.
The [complete table](../data/garnet-angel-2022/grossular/du-bosenick-summary.csv)
retains every parameter and both backgrounds. Objectives across different
covariance assumptions are not a calibrated model-selection ranking.

![Du and Bosenick controls, light theme](../data/garnet-angel-2022/grossular/du-bosenick.svg#only-light)
![Du and Bosenick controls, dark theme](../data/garnet-angel-2022/grossular/du-bosenick-dark.svg#only-dark)

Top panels isolate Bosenick's anchor and covariance with Du's listed
temperatures. Bottom panels vary Du processing with Bosenick's propagated
volume and independent raw-T covariance. Lines connect categorical controls;
they do not interpolate physical models. Starred calibration-stage cases
contradict the preprint's stated processing and remain error controls only.

## Du temperatures were already corrected

The available [Du author preprint](https://academiccommons.columbia.edu/doi/10.7916/D8V415Z5/download),
methods p. 7, explicitly uses NaCl-derived temperatures throughout. Its
Eq. 1 maps thermocouple temperature T₁ to NaCl temperature T₂, **both in °C**:

\[
T_2=-3.6\times10^{-5}T_1^2+0.982T_1-0.565.
\]

The source describes a difference of about 25 °C near 900 K. Table 2's listed
temperatures therefore must not receive that correction again as a source
replay. The experiment's thermocouple control/stability statements and its
0.2 °C table error do not supply raw NaCl observations or calibration covariance.

Two explicitly unsupported stage-error controls reconstruct T₁ by inverting
Eq. 1 or apply Eq. 1 a second time to the listed T₂. They use Celsius for the
calculation, then return to kelvin, and propagate the existing 0.2 K width
through the mapping. The calculated Du 298 K denominator is recomputed when
the background requires it. Neither stage control closes θ:

| Du treatment, Bosenick calculated 298 K anchor and independent raw-T covariance | θ, plotted background, K | θ, linear 298 K background, K |
| --- | ---: | ---: |
| Listed NaCl temperatures | 851.495 | 851.127 |
| Reconstructed thermocouple temperatures | 895.616 | 892.096 |
| Accidental second correction | 804.091 | 806.484 |
| Listed temperatures, fully shared 0.2 K error | 848.446 | 848.162 |
| Du a³ and propagated σₐ | 821.213 | 821.662 |

The final publisher table remains unverified: the current full-text request
was forbidden, other public retrievals were inaccessible/gated, and the cached
institutional preprint remains the evidence basis. Inaccessible is not absent.
No fitted temperature correction or additional temperature-error width was
introduced to force agreement.

## Du's lattice and volume precision are not identical

The alternate column-based control computes V=a³ and σ(V)=3a²σ(a). For the
plotted background it preserves the first observation's normalized ratio;
for the linear background it computes the 298 K denominator from the first
two lattice-derived volumes and listed temperatures. This is a calculated
precision/error-propagation control, not recovery of an author input file.

Under the linear background, the central normalized volumes change by at
most **0.00514 cm³/mol**, while coordinate widths become **0.843–1.696 times**
those from the rounded volume-error column. The combined representation change
lowers θ by about 29.5 K with the full covariance, but leaves it near 821.7 K.
Since both central values and widths change, that fit shift cannot be uniquely
assigned to one of them. It does not repair the much larger previously observed
middle-temperature volume departures or recover the published θ.

## What Bosenick actually measures

[Bosenick and Geiger (1997)](https://doi.org/10.1029/97JB01612), Table 2,
gives the room anchor **125.293(16) cm³/mol at 295 K**, not at 298 K. The
calculated alternative is

\[
D_{298}=X_{295}+\frac{298-295}{295-270}(X_{295}-X_{270})
       =1.12X_{295}-0.12X_{270}=125.30104\;\mathrm{cm^3/mol}.
\]

With diagonal errors, changing only this anchor lowers θ by **17.31 K**
under the plotted Milani/Du background and **17.25 K** under the linear
background. This isolates the anchor contribution more directly than the
preceding simultaneous normalization comparison, but does not account for
the remaining approximately 14–15 K gap.

The final primary paper, p. 22650, reports temperature accuracy of 3 K down
to 100 K and 5 K below, alongside much smaller within-scan variation. It
also describes cooling/heating discrepancies below 150 K and measurements
after holding the sample at 20 K for eight hours. These statements do not
establish Gaussian errors or a cross-point covariance. Neither repeat set
was discarded. The source's counts and rows are preserved; the original
room measurement remains the normalization anchor, not a new fit channel.

## Propagating the anchor without double counting

Let X contain all 20 raw molar volumes, D=wᵀX, and let each fitted non-anchor
volume be Yᵢ=125.35 Xᵢ/D. For the measured anchor w=(1,0,…); for the
extrapolated anchor w=(1.12,−0.12,0,…). Its volume Jacobian is

\[
J_{ij}=\frac{125.35}{D}\delta_{ij}-\frac{Y_i}{D}w_j.
\]

The **270 K measurement is both a numerator and part of D** in the second
case. Adding an independent anchor error to each row would miss that
dependence and its cross terms. We instead propagate the entire Jacobian
against the declared independent printed raw-volume variances. A common
multiplicative raw-volume calibration cancels in the ratios, as checked by
JX=0; it must not be added again as independent anchor noise.

For the calculated anchor, the full temperature control also differentiates
D with respect to its 295 K and 270 K temperatures. The pressure residual
Jacobian includes both direct ∂P/∂T and the denominator-induced Y change.
For the measured-295 room-reference approximation, D has no T dependence
in the declared processing. These are different hypotheses, not one implicit
interpretation of a room-temperature measurement.
For the fixed measured-295 denominator, independent raw-T propagation gives
the same temperature term as the volume-only control, so those two fits coincide.

The source's 3/5 K accuracy bounds are used as **comparison widths**. One
control treats their raw-temperature effects as independent; the other fully
correlates their standardized errors. The latter is an extreme calibration
control with those marginals, not a recovered constant temperature offset.
All pressure covariance blocks are recomputed at each IRLS iteration and
whitened by Cholesky factors. No jitter or eigenvalue clipping was used.

Propagating the shared volume anchor weakens the cold series' common-offset
constraint. The upward shift toward 850 K is consistent with allowing Du
more relative influence, as the preceding Bosenick-omission diagnostic also
gave high θ. It is not evidence that the chosen covariance model is the
original author's or a preferred physical fit.

## Scope of the result

These controls are more explicit about normalization covariance, but they
do not recover raw Si/calibration correlations, thermal-equilibration/run
errors, Du's denominator covariance or the original EosFit error columns.
Raw-volume independence and the temperature-correlation limits remain
assumptions. Consequently **the tested choices do not explain the central
θ discrepancy**, while other unrecovered author inputs remain unresolved.
The close K′ values in several controls do not validate θ or the full EOS.

The additional independent 2% Conrad error allowance remains a separate,
retained diagnostic; it has not become measured uncertainty or a catalog
default. No source observations, coefficients or classifications changed.
No further author message was drafted or sent. These fits are not recommended
for predictions or calibration.

## Reproducible artifacts and checks

The [predeclared plan](../data/garnet-angel-2022/grossular/du-bosenick-plan.json)
contains 16 anchor/covariance cases and 16 single Du-treatment controls,
two starts each. [All fits and covariance blocks](../data/garnet-angel-2022/grossular/du-bosenick-refits.json),
[pressure residuals](../data/garnet-angel-2022/grossular/du-bosenick-residuals.csv),
[transformed inputs](../data/garnet-angel-2022/grossular/du-bosenick-inputs.csv)
and [anchor Jacobians](../data/garnet-angel-2022/grossular/du-bosenick-anchor-jacobians.json)
are retained separately from originals. There are **7,008 pressure residuals**
and 1,120 affected-series/Conrad input rows.

Tests check the complete ratio and pressure Jacobians by finite differences,
Celsius calibration inversion and derivatives, common-scale cancellation,
dense GLS versus Cholesky whitening, stationary scores, inverse-information
covariance, baseline parity, case replays, all exported objectives and all
nine previous artifact manifests. See the [verification receipt](../data/garnet-angel-2022/grossular/du-bosenick-verification.json)
and [primary source evidence](../data/garnet-angel-2022/grossular/du-bosenick-source-evidence.json).

Repository-only commands:

```bash
python -m scripts.refit_grossular_du_bosenick
python -m scripts.plot_grossular_du_bosenick
pytest -q tests/test_grossular_du_bosenick.py
```
