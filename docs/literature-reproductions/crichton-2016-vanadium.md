# Crichton et al. (2016): bcc vanadium

Crichton, W. A., Guignard, J., Bailey, E., Dobson, D. P., Hunt, S. A., and
Thomson, A. R., *High Pressure Research* **36**, 16–22 (2016),
[doi:10.1080/08957959.2015.1123256](https://doi.org/10.1080/08957959.2015.1123256).
Audit date: 2026-09-19.

## Source access and acceptance

The full author-accepted manuscript is available from
[UCL Discovery](https://discovery.ucl.ac.uk/id/eprint/1478098/).
Its Results and Figure 2, together with the independently published EosFit
methodology, establish the main experimental parameterization. The publisher's
publication-of-record abstract was inspected in a browser and agrees with all
five printed modulus/thermal coefficients and their errors. The final full text
and figures require subscription access; direct full-text, PDF and supplement
requests returned HTTP 403. This audit does **not** claim inspection of the
complete typeset publication. UCL identifies its deposited version as accepted
and the publication date as 21 January 2016.

Crossref metadata has no correction/update relation. The publisher article
navigation, UCL deposit, Manchester institutional record, and DOI/title searches
yielded no accessible official supplement or correction. The CrossMark browser
widget did not return usable content. Thus no correction was found, not a proof
that none exists. Zotero's local API was available: exact-title and Crichton
searches, and a vanadium search in Methods → EOS Library (`JT8V6LUL`), found no
matching paper/attachment. No library write or author contact was performed.
The [source inventory](../data/crichton-2016-vanadium-source-audit.json) records
URLs, checksums, candidate dispositions and retrieval limitations.
The [author-request entry](../dataset-requests.md#crichton-et-al-2016-bcc-vanadium)
lists the missing inputs, their purpose, and the current outreach status.

| Source-owned candidate | Decision |
|---|---|
| Unconstrained 62-point thermal BM3 | Accepted as `vanadium_bcc_crichton_2016_bm3_thermal`; independent plotted-curve reproduction passes. |
| BM3 with fixed K′ = 3.5; K = 160.5(3.5) GPa | Withheld: the alternative fitted V0 and thermal coefficients are not supplied. |
| BM3 with fixed K′ = 4; K = 157.9(3.5) GPa | Withheld for the same missing values; no coefficients borrowed from the unconstrained fit. |

The other EOS and elastic numbers in Results are explicitly attributed to
earlier experimental or computational papers. They are comparisons, not fits
owned by this paper. Likewise h-BN, NaCl and Au are supporting materials, not
new source-owned EOSs. No computational fit or aggregate EOS is added.

## Material and volume

The EOS sample is elemental Advent vanadium foil, initially 0.020 mm thick;
purity is not quantified. The bcc lattice is the measured V phase, separate
from the co-loaded NaCl/Au powder and h-BN capsule. The first compression/heating
cycle up to about 600 K showed anomalous strain, attributed to foil fabrication
and removed by annealing. The text does not give a row-level exclusion list.
The final recovered ambient point was explicitly **not** fitted.

Results (manuscript p.5) gives a fitted relative reference volume **1.005(1)**,
equal at printed precision to the recovered value **1.005 = 27.995 Å³**;
Figure 2's caption confirms this equality. The executable V0 is therefore
**27.995 Å³ per conventional bcc cell**, containing two V atoms (Z = 2).
The plotted denominator is instead 27.995/1.005 = 27.855721393 Å³. Replacing
V0 by that denominator would introduce a spurious 0.5% shift. The relative
uncertainty 0.001 is preserved explicitly; the source does not print a separate
absolute fitted-volume error, so the record's absolute V0 error remains null.

Crystallography is separately sourced from James and Straumanis (1960),
*Journal of the Electrochemical Society* **107**, 69, through the attributed
[COD 9012770 CIF](https://qiserver.ugr.es/cod/9012770.cif)
(AMCSD 0014111). The original CIF is bundled and checksummed. At 298.15 K,
a = b = c = 3.0241 Å, all angles 90°, space group Im-3m (229), with a fully
occupied V 2a site at (0,0,0). Body centering generates the second atom at
(½,½,½), giving Z = 2 and a³ = 27.6559412875 Å³. This structural cell is
not forced to equal the relaxed foil's EOS reference volume. No inherited
database EOS coefficient supplies the thermal record.

## Exact equation and provenance

The publication names third-order Birch–Murnaghan, EOS-FIT5.2 (Angel, 2002),
and a linear alpha(T). The primary software-methodology paper
[Angel, Gonzalez-Platas and Alvaro (2014)](https://www.rossangel.com/Download/2014_Angel_etal_EosFit7.pdf),
pp.410 and 412, explicitly states that earlier EosFit versions used the
exponential integral of alpha0 + alpha1 T and documents linear K0(T).
Thus the exact existing `BM3` + `ThermalReferenceStateEOS` mapping is:

\[
\begin{aligned}
\alpha(T)&=\alpha_0+\alpha_1T,\\
V_0(T)&=V_0\exp\left[\alpha_0(T-300)
             +\tfrac12\alpha_1(T^2-300^2)\right],\\
K_0(T)&=K_0+(\partial K_0/\partial T)(T-300),\\
P(V,T)&=\tfrac32 K_0(T)(x^7-x^5)
 \left[1+\tfrac34(K'_0-4)(x^2-1)\right],\qquad
x=[V_0(T)/V]^{1/3}.
\end{aligned}
\]

| Parameter | Value | Source / status |
|---|---:|---|
| V0 | 27.995 Å³ | Results p.5 and Figure 2, derived absolute normalization of fitted 1.005(1) |
| K0,300 | 150.4 ± 6.2 GPa | Results p.5 and publisher abstract; fitted |
| K′0 | 5.5 ± 1.0 | Same; fitted, no temperature variation specified |
| alpha0 | 4.8(6) × 10⁻⁵ K⁻¹ | Same; fitted absolute-temperature intercept |
| alpha1 | −2.4(9) × 10⁻⁸ K⁻² | Same; fitted |
| dK0/dT | −0.0446(7) GPa K⁻¹ | Same; fitted zero-pressure modulus slope |
| Tr, Pr | 300 K, 0 GPa | Experimental reference and Results K0,300 convention |

All six independent coefficients were free in one unweighted fit of 62 P-V-T
states. No covariance or confidence level for the printed errors is reported.
`alpha0` is not alpha at 300 K and is not linear lattice expansivity. The
reference-temperature equivalent intercept would be 4.08e-5 K⁻¹; its uncertainty
cannot be propagated independently without the missing covariance. The stored
absolute-temperature convention preserves both published coefficients/errors.
Neither the Berman polynomial nor a linear volume law is used.

The printed coefficients imply alpha(300)K0 = 0.00613632 GPa/K, consistent
with the approximate 0.00616 in the abstract within parameter rounding/errors,
and K0(1000 K) = 119.18 GPa. Alpha's extrapolated zero at 2000 K supports the
absolute-temperature intercept; it is not a valid extrapolation limit.

## Observations, pressure calibration and domain

The accepted manuscript has no numerical data table. Figure 2 displays only
points within ±20 K of the 300, 450, 700, 900 and 1000 K isotherms. The bundled
CSV preserves every distinguishable symbol: **29 filled source symbols**, ten
open Ming–Manghnani (1978) comparison symbols, and the excluded recovered star.
This is an approximate partial projection of the source's 62-state dataset,
not the original table. Overlap can hide additional points. Exact measured
temperatures, measurement uncertainties and paired marker readings are absent.

Digitization uses the unchanged 950×1017 JPEG extracted from manuscript p.11:
`pdfimages -f 11 -l 11 -j accepted.pdf fig2`. Pixel centres are recorded in each
CSV. Calibrations are x = 159 at 0 GPa, x = 849 at 12 GPa; y = 201 at V/Vref
= 1.02, y = 894 at 0.94. Both axes use a conservative three-pixel reading
bound (0.05218 GPa and 0.00034633 in relative volume). These are digitization
bounds, **not** source standard deviations. The ±20 K interval is a selection
bin, not temperature measurement uncertainty. Two partly overlapping 900 K
symbols near 6 GPa are separately read, with the same conservative bound.

The independent curve CSV has 24 readings across all five isotherms, from
0.50 to 11.36 GPa. It records calculated line centres separately from symbols;
the merged 900/1000 K lines at the highest-pressure reading cannot be read
separately, so only red is used there. Curve data never enter the proxy fit.

Pressures **and** temperatures were determined by PTX-Cal cross-calibration of
simultaneous NaCl/Au diffraction lattices, accepting |PNaCl−PAu| < 1e-6 GPa
while scanning temperature and averaging returned solutions. The quoted
approximately 0.5 K solution width is numerical solver resolution, not a
physical error bar. The method cites Crichton and Mezouar (2002), but the
accessible paper does not identify the NaCl/Au coefficient versions. No
specific bundled calibrant is substituted. Re-reduction is unavailable because
both paired readings and exact calibration versions are missing.

Validity bounds are the reported marginal 0–11.5 GPa, 300–1000 K experimental
range, not a sampled rectangular domain or a bcc phase-stability guarantee.
The cited bcc/rhombohedral transitions near 60–69 GPa and reentrant bcc phase
are outside the represented data; they do not justify extending the bounds.
Oxidation, chemical reaction, and the initial stressed foil are outside scope.

## Reproduction and qualified refit

Run:

```bash
uv run python scripts/reproduce_crichton_2016_vanadium.py --check
uv run pytest -q tests/test_crichton_2016_vanadium.py
uv run python scripts/validate_primary_eos_refits.py --check
```

The standalone NumPy/scipy equation and bracketing inverse are independent of
Peritheos. Across the 24 source curve readings, maximum relative-volume error
is **0.00056695**, below the 0.001 tolerance set by the reported reference-ratio
precision, rounded coefficients and three-pixel digitization bound. Native
pressure agrees with the independent expression to about 1.1e-13 GPa.

The Conclusion's 300–1000 K expansions (2.7% near ambient, 1.2% at 10 GPa)
are **not exactly reproduced**: the printed coefficients give 2.2939% and
1.0255%. Figure 2 supports those coefficients. The ambient plotted coordinate
at 1000 K is 1.02805 relative to the pre-relaxation denominator, suggesting a
possible normalization contribution to the prose discrepancy; this is an
inference, not a published correction. We do not change alpha0, alpha1 or V0
to match the prose percentages or falsely list them as passed benchmarks.

An independent simultaneous six-parameter, unweighted pressure-residual fit
uses all 29 distinguishable filled symbols at their **nominal** temperatures.
It excludes open comparison symbols, the recovered star and all curve samples.
The approximate optimum is:

| Parameter | Published | Nominal-temperature proxy |
|---|---:|---:|
| V0 (Å³) | 27.995 | 28.0460 |
| K0 (GPa) | 150.4 | 150.443 |
| K′0 | 5.5 | 4.75936 |
| alpha0 (K⁻¹) | 4.8e-5 | 3.54924e-5 |
| alpha1 (K⁻²) | −2.4e-8 | −6.83787e-9 |
| dK0/dT (GPa/K) | −0.0446 | −0.0400178 |

Pressure RMS falls from 0.24264 to 0.18248 GPa. This does not recover all the
printed thermal coefficients. Missing rows, nominal temperatures and overlapping
symbols prevent calling this original-fit parity or estimating original-fit
parameter errors. The ledger therefore says **direct refit unavailable**, with
the proxy diagnostic attached. No opt-in replacement refit is justified by
these incomplete observations. The original published coefficients remain the
only executable record. The full numerical report is
[bundled here](../data/crichton-2016-vanadium-reproduction.json).
