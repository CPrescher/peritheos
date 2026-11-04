# Bykova (2018): published joint coesite-IV/V experimental EOS

Record: `coesite_iv_v_bykova_2018_experimental_bm2_reference`.
Audit: 24 October 2026. Source: Bykova et al., *Nature Communications* 9, 4789,
[doi:10.1038/s41467-018-07265-z](https://doi.org/10.1038/s41467-018-07265-z).

The **published experimental curve is reproduced**, with all source coefficients
and all 16 tabulated observations preserved. Its large between-run scatter and
unavailable regression protocol remain explicit. This is an empirical fit to
**both coesite-IV and coesite-V**, not an independently established EOS for each
phase, a phase mixture model or a recommended pressure calibrant. It is stored
with coesite-IV for discovery, with both phases identified in `phase_scope`;
the existing theoretical coesite-IV refit remains separate and unchanged.

## Published parameters and reference convention

Figure 1's caption and Supplementary Figure 9 report a second-order
Birch–Murnaghan EOS referenced to **32.7 GPa**:

| Source quantity | Printed value | Executable parameter |
| --- | ---: | --- |
| Reference pressure | 32.7 GPa | `P0 = 32.7`, fixed reference coordinate |
| V32.7/Z | 23.44(3) Å³ per SiO2 | `V0 = 375.04`, error 0.48 Å³ per Z=16 cell |
| K32.7 | 254(9) GPa | `K0 = 254`, error 9 GPa |

Both phases have Z=16. The parentheses are preserved as printed error widths;
their confidence and covariance are unspecified. The model parameter names
V0 and K0 follow the library's reference-state convention: for `BM2Reference`
they mean values **at P0**, not ambient values. The source's printed V/Z and
uncertainty are also retained without conversion in `published_reference_state`.
Room temperature is represented by nominal 298 K with an explicit assumption;
it is not a measured temperature or a thermal model.

With eta = (V0/V)^(1/3), the represented law is

\[
P(V)=P_0+\frac32 K_0(\eta^7-\eta^5).
\]

This convention is identified by a quantitative comparison with the source's
drawn red curve. P(V0)=P0, K(V0)=K0 and dK/dP at V0 is 4. It is not generally
equivalent to taking an arbitrary ambient BM2 and simply relabelling its
parameters at another pressure. The new explicit `BM2Reference` model preserves
that distinction. Ordinary `BM2` retains its existing zero-pressure semantics.
Readers without this model must reject its identifier; a v3 compatibility
export does not make it executable in an older reader.

## Independent published-curve check

The original article PDF has SHA-256
`776ca4ca61872b7540c1ebaa70c29046f84a41c2d0700725d009eba203bd3cb5`.
Figure 1 is on PDF page 2. Its red fitted curve is a vector path: a short line,
a cubic Bezier segment and another short line. The exact path commands,
control points, 0.475-point stroke width, and axis ticks are archived in
`peritheos/data/datasets/coesite-iv-v-bykova-2018-figure1-vector.json`.
These are **published fitted-curve checkpoints**, not measured observations.

Affine coordinate maps use all eight labelled pressure ticks (0–70 GPa) and
eight volume ticks (20–34 Å³ per formula unit). Sampling the archived path
provides **23 checkpoints within the observed 29.9–70.4 GPa envelope**.
The unchanged published model differs from that drawn curve by pressure RMS
**0.0247 GPa**, maximum **0.0427 GPa**. Every checkpoint's residual interval
contains zero using half the source stroke width plus the maximum axis-map
residual: ±0.1055 GPa horizontally and ±0.3943 Å³ per conventional cell
vertically. These graphical bounds are propagated through the exact equation;
they are not standard deviations or confidence intervals. The Bezier path is
itself a drawn approximation to the authors' analytic curve.

This validates the equation mapping and printed coefficients at the figure's
resolution. It does not reproduce the regression used to obtain them.

## All experimental observations and their residuals

Supplementary Table 2 (PDF page 13) supplies **12 coesite-IV and four coesite-V
rows**, spanning 29.9–70.4 GPa. The new CSV retains source row numbers,
phase, sample/run identifiers, original pressure/volume/error tokens, decoded
error widths, and the separately printed V/Z column. The latter has its own
rounding; it is not substituted for cell V/16. No outlier is removed.
The supplement SHA-256 is
`9d4bcb36ae6a89b19f5bbcc3fc671f3c1bfc6530a76915382c6132d3d907e459`.

At each printed **cell volume**, pressure residual is P(published model) minus
P(table). The pressure RMS over all 16 rows is **4.0466 GPa**, maximum
**11.8546 GPa** at source row 13 (coesite-V, 44.5(5) GPa,
V=346.4(8) Å³, DAC#1 c1). The source's Figure 1 caption explicitly notes
substantial coesite-III/IV volume scatter, possibly from different deviatoric
stresses between runs. That statement does not justify deleting any row or
changing its printed value; the coesite-V residual is retained independently.

An **analyst diagnostic**, fitting all 16 rows with equal pressure weights and
P0=32.7 GPa fixed, yields Vref≈378.272 Å³ (23.6420 Å³ per formula unit),
Kref≈212.381 GPa and pressure RMS≈3.6587 GPa. It does not recover the printed
reference values. Original weights, objective, row selection and whether the
reference volume was fixed to the datum or refined are unspecified. This
diagnostic is not a replay of a known author protocol, a replacement EOS, or an
additional catalog parameterization. Its coefficients are retained in the
numerical audit only. No statistical combined-error or author-fit parity claim
is made from the unspecified printed error convention.

| Parameter | Published | Printed error width | Diagnostic refit |
| --- | ---: | ---: | ---: |
| Vref, Å³ per cell | 375.04 | 0.48 | 378.2718 |
| Kref, GPa | 254 | 9 | 212.3808 |
| Pref, GPa | 32.7 | Not a fitted-error estimate | 32.7 (fixed) |

Diagnostic uncertainty is not inferred here. The dashed curve and open grey
residual markers show this diagnostic on exactly the same observations as
the published curve; it remains unadopted.

![Published experimental curve, PDF curve checkpoints, and all observations with residuals](assets/coesite-iv-v-bykova-2018-experimental.svg#only-light)
![Published experimental curve, PDF curve checkpoints, and all observations with residuals](assets/coesite-iv-v-bykova-2018-experimental-dark.svg#only-dark)

The source equation is validated against its published curve; agreement with
the observations is much poorer and is reported separately. Use is restricted
to reproducing the reported empirical curve by explicit selection. It is
**not recommended for quantitative pure-phase predictions or pressure
calibration**, and is not made an automatic default.

## Pressure calibration and remaining evidence

The Methods identify Ne as both pressure medium and diffraction pressure
standard. The exact reference EOS version and paired Ne measurements were not
recovered; no calibration is guessed from the software or publication date.
Source-observation pressure recalculation is therefore unavailable.

The original regression inputs, row mask, objective, weights, reference-volume
constraint and covariance are needed for author-workflow reproduction. The
graphical check does not supply them. No outreach was sent.

## Reproduce or use

Repository-only audit, leaving all published parameters unchanged:

```bash
uv run python -m scripts.reproduce_bykova_2018_experimental --check
uv run --group docs python -m scripts.reproduce_bykova_2018_experimental --check --plot
```

Installed-package evaluation of the published curve:

```python
from peritheos import get_eos_record

record = get_eos_record("coesite_iv_v_bykova_2018_experimental_bm2_reference")
pressure = record.pressure(375.04)  # 32.7 GPa, not zero
volume = record.volume(50.0)  # conventional-cell Å³
```

[Full numerical comparison](../data/coesite-iv-v-bykova-2018-experimental.json) ·
[Separate theoretical-refit audit](coesite-iv-bykova-2018.md) ·
[Material parameters and source data](../material-summaries/coesite_iv.md)
