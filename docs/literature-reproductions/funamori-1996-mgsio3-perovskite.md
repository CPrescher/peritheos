# Funamori et al. (1996): MgSiO3 perovskite thermoelastic audit

## Audit outcome and primary source

This audit covers all four LitCurate rows grouped under Nobumasa Funamori et
al., *Thermoelastic properties of MgSiO3 perovskite determined by in situ X ray
observations up to 30 GPa and 2000 K*, **Journal of Geophysical Research 101**,
8257--8269 (1996),
[doi:10.1029/95JB03732](https://doi.org/10.1029/95JB03732).

The scientific authority was the final 13-page publisher article, locally
archived as `95jb03732.pdf` with SHA-256
`0f1a72b523d8ba242bf4baadbd0c31cc2eb7a488d312e52fd05442d9e1a8eb13`.
LitCurate was used only to find the candidate rows.

The paper contributes **zero production EOS records** in this audit. The two
apparently source-reported LitCurate rows are not two author-generated EOS
fits. They split one normalized, combined-data thermoelastic model into two
pseudo-records by substituting the separate ambient volume calibration used
for each experimental run. Encoding both would duplicate one physical model;
encoding either one would silently select a volume normalization that the
authors did not designate as the model's unique `V0`.

The genuine Funamori et al. thermal EOS remains a documented future candidate.
It is not currently representable by Peritheos's supported thermal terms, and
an exact refit would also require the numerical lower-pressure observations
from Wang et al. (1994) and Utsumi et al. (1995) that the authors combined with
their own data.

## Material and experiment

The authors crushed flux-grown enstatite single crystals, compressed the
powder to about 25 GPa, and heated it to 1200 K to transform it through the
ilmenite phase to pure MgSiO3 perovskite. The phase is the conventional
orthorhombic perovskite now called bridgmanite. Two in situ energy-dispersive
X-ray diffraction runs used the MAX90 MA8 apparatus with sintered-diamond
second-stage anvils. Run 2 replaced the cobalt binder on the X-ray-side anvil
with silicon carbide to reduce exposure time.

The observed conditions span 21--29 GPa and 293--2000 K. Run 1 used NaCl plus
Au pressure markers and seven alternative pressure scales. Run 2 used MgO and
two alternative scales. The paper ultimately chooses the Decker (1971) NaCl
coordinate, `ND1`, for the global combined-data fit so that all included data
share one pressure scale. The `AA1` Au scale is the authors' representative
coordinate for their separate thermal-expansion discussion, while `MJ2` is a
successful independent run-2 check. The exact calibrant EOS implementations
and raw diffraction energies are not all available here, so full pressure
recalculation is not possible.

## What the source actually parameterizes

The paper explicitly says that uncertainty in the diffraction angle or MCA
energy-channel relation creates a small ambient-volume offset between runs.
It therefore normalizes *every observation by the ambient volume of its own
run* before all analyses:

| calibration role | source construction | volume (A3) |
|---|---|---:|
| Run 1 normalization | one post-run ambient observation in Table 4a | 162.32 |
| Run 2 normalization | mean of the two post-run ambient values 162.33 and 162.05 A3 in Table 4b | 162.19 |

These are measurement calibrations, not two zero-pressure volumes obtained by
separate EOS regressions. The paper never reports a run-1 `K0`, a run-2 `K0`,
or two independent BM parameter tables.

The actual source model is a single third-order Birch--Murnaghan thermal EOS:

\[
P(V,T)=\frac{3K_{T,0}}{2}
\left[\left(\frac{V_{T,0}}{V}\right)^{7/3}
-\left(\frac{V_{T,0}}{V}\right)^{5/3}\right]
\left\{1+\frac{3}{4}(K'_{T,0}-4)
\left[\left(\frac{V_{T,0}}{V}\right)^{2/3}-1\right]\right\}.
\]

It fixes `K300,0=261 GPa` and `K'300,0=4` to the Mao et al. (1991) result,
fixes the 300 K ambient expansivity to the Wang et al. (1994) result, and fits
the remaining thermal coefficients to the authors' selected observations plus
the lower-pressure Wang and Utsumi data. Table 6 gives the preferred global
model:

| coefficient | value | source role |
|---|---:|---|
| `K300,0` | 261 GPa | fixed from Mao et al. (1991) |
| `K'300,0` | 4 | fixed from Mao et al. (1991) |
| `(dKT,0/dT)P` | -0.0280(32) GPa/K | fitted; authors expand the realistic uncertainty to +/-0.017 GPa/K |
| `a0` | 1.982 x 10^-5 K^-1 | constrained through the 300 K Wang expansivity |
| `b0` | 0.818(257) x 10^-8 K^-2 | fitted |
| `c0` | 0.474(118) K | fitted |

The accompanying definitions are
`KT,0=K300,0+(dKT,0/dT)P(T-300)`, `K'T,0=K'300,0`,
`alphaT,0=a0+b0*T-c0*T^-2`, and
`VT,0=V0*exp(3*integral(alphaT,0 dT))`. This is not a static BM3 fit with a
newly measured `K0`. Nor is it equivalent to two static BM3 records merely
because its run-normalized measurements require two calibration volumes.

## Numerical convention check

The standard BM3 reference curve can still be checked against the near-room-
temperature Table 4 observations. Using each run's own normalization, the
published fixed `K0=261 GPa`, `K0'=4`, and the pressure coordinate closest to
the paper's preferred analysis gives:

| run | T (K) | V (A3) | pressure scale | observed P (GPa) | calculated P (GPa) | residual (GPa) |
|---|---:|---:|---|---:|---:|---:|
| 1 | 293 | 150.51 | ND1 | 23.10 | 22.933889 | -0.166111 |
| 2 | 293 | 150.55 | MJ2 | 23.28 | 22.561706 | -0.718294 |
| 2 | 293 | 149.93 | MJ2 | 23.85 | 24.009529 | +0.159529 |
| 2 | 293 | 149.88 | MJ2 | 23.70 | 24.127528 | +0.427528 |

With `K0'=4` fixed, a conditional unweighted fit to those four rounded points
gives `K0=261.698443 GPa` and pressure RMSE `0.428970 GPa`. This supports the
standard BM3 algebra and the plausibility of the adopted modulus. It does not
turn the two run normalizations into two EOS fits, and it cannot reproduce the
global thermal regression, which deliberately used a larger combined dataset
and fitted thermal rather than static coefficients.

## Disposition of every same-DOI LitCurate row

| LitCurate identifier | reported candidate | disposition | reason |
|---|---|---|---|
| `litcurate_b925041100f09082` | Run 1: `V0=162.32 A3`, `K0=261 GPa`, `K0'=4` | **rejected as a duplicate pseudo-record / run normalization** | `162.32 A3` calibrates run 1. The modulus and derivative were fixed from Mao et al. (1991), not fit to run 1, and Table 6 reports no run-specific EOS. |
| `litcurate_3548aee7dcd347b9` | Run 2: `V0=162.19 A3`, `K0=261 GPa`, `K0'=4` | **rejected as a duplicate pseudo-record / run normalization** | `162.19 A3` is the mean of two ambient run-2 measurements. It normalizes the run-2 observations but does not define a second physical EOS. |
| `litcurate_292dc6dfe6cd7bf0` | Mao et al. (1991), `(Mg0.9Fe0.1)SiO3`: `K0=261 GPa`, `K0'=4` | **citation trace only** | The source explicitly imports these fixed static coefficients. The incomplete comparison belongs under the Mao primary publication and is a different composition. |
| `litcurate_d71096ad980e2a61` | Jeanloz and Hemley (1994) consensus: `K0=260 GPa`, `K0'=4` | **citation trace only** | This is a literature consensus cited to justify the adopted constraints, not a Funamori measurement or complete independent volumetric EOS. |

## Future-work boundary

No existing production record uses DOI `10.1029/95JB03732`; therefore this is
not a catalog-level duplicate of another source. A future implementation
should represent the *one* combined-data thermal model and preserve the two
run volume normalizations inside its fit dataset. It should not promote the
two LitCurate rows independently. That work requires either a thermal model
capable of the paper's empirical `alpha(T)` and linear `K0(T)` laws or an
explicitly documented equivalent mapping, plus the complete Wang and Utsumi
input tables and the source's weighting/exclusion protocol.
