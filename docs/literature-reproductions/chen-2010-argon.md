# Chen et al. (2010): fcc argon supporting study

**Outcome: supporting study; published EOS not reproduced.** The paper and
its useful density observations are included, but no executable Chen EOS record
is added. Substituting its 2 GPa elastic constants for zero-pressure BM3
parameters would be incorrect. The standard BM3 reconstructed from those
local constants also fails the source curve. This is an unresolved source
parameterization problem, not a claim that the experiment is invalid.

## Source and physical interpretation

Bin Chen, A. E. Gleason, J. Y. Yan, K. J. Koski, Simon Clark and Raymond
Jeanloz, *Elasticity, strength, and refractive index of argon at high
pressures*, Physical Review B **81**, 144110 (2010),
[DOI 10.1103/PhysRevB.81.144110](https://doi.org/10.1103/PhysRevB.81.144110).
Published 14 April 2010. The verified five-page publisher PDF was retrieved
from [APS's full-text service](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.81.144110/fulltext).
SHA256: `17a0bcd30604bd81f125894a18a45a91255b235fb80643cddd4bcbc6e5bd01f1`.
The PDF is external to git at
`/Users/clemens/Documents/Peritheos/papers/chen-2010/chen-2010.pdf`.
Title, all authors, DOI, volume/article, and all five page numbers were checked.

The authors measured polycrystalline **fcc** Ar using two Brillouin scattering
angles. XRD confirmed fcc at 4 GPa. The headline range is **1.3–30 GPa** at
room temperature. Table I specifically states **2 GPa and 290 K**, whereas
the strength discussion on page 4 says 300 K. The datasets retain 290 K as
the table reference without pretending the experiment was exactly isothermal
to a known temperature precision. No thermal model is supplied.

Pressure comes from ruby fluorescence, Experimental section reference **14**:
Mao, Bell, Shaner and Steinberg (1978), J. Appl. Phys. 49, 3276. Reference 13
is Mao et al. (1986), but is not the calibration cited for these measurements.
Raw ruby wavelengths are absent, so pressure-scale recalculation is unavailable.
The Figure 4 caption gives typical pressure uncertainty 0.2–0.7 GPa and
velocity uncertainty about 2%; no confidence convention or covariance is stated.

Density is in g/cm³. Converted volumes are conventional **four-atom fcc cells**:

\[
V_{\rm cell}[\mathrm{\AA^3}]
=\frac{4(39.948)}{0.602214076\,\rho[\mathrm{g/cm^3}]}.
\]

The mass conversion uses 39.948 g/mol and the exact Avogadro constant. It does
not add experimental precision. Density remains the source quantity.

## What is actually published

The source prints the acoustic relationships
\(K_S=\rho(V_P^2-4V_S^2/3)\), \(G=\rho V_S^2\), and
\(K_S=(C_P/C_V)K_T\). Its Eq. (4) is

\[
\frac{d\rho}{dP}=\frac{C_P/C_V}{V_P^2-4V_S^2/3},\qquad
\frac{C_P}{C_V}=1+0.25(\pm0.11)\exp[-0.197(\pm0.041)P].
\]

Here P is GPa and velocities are km/s. Integration starts at
**ρ = 2.02 g/cm³, P = 1.3 GPa**. The authors mention numerical integration
and iterations, but supply no full velocity table, interpolation specification,
integration covariance or fit weights. The script implements the printed
pointwise equation without constructing an invented velocity interpolant.

The paper names a third-order Eulerian finite-strain/Birch–Murnaghan EOS,
citing Birch (1977, 1978), but does **not** print its pressure equation or
unrounded coefficient set. It reports:

| Quantity | Value | Location |
|---|---:|---|
| ρ at 2 GPa | 2.18 ± 0.06 g/cm³ | Page 3 text |
| KT at 2 GPa | 15.1 ± 1.1 GPa | Table I |
| dKT/dP at 2 GPa | 5.4 ± 0.3 | Table I |
| KT d²KT/dP² at 2 GPa | −7.3 ± 1.2 | Page 4 text |
| Extrapolated ρ at zero pressure | 1.52 ± 0.05 g/cm³ | Page 3 text/Figure 5 |

All 14 own-study values recovered from Table I and adjacent prose/captions are
in `argon-fcc-chen-2010-reported-values.csv`, including the adiabatic and shear
constants. Comparison rows attributed to other studies are not Chen observations.
The uncertainty symbols are retained; they are not assumed to be one sigma.

## Recovered Figure 5 data

The black squares are **densities derived by integrating Brillouin velocities**.
They are not XRD measurements. The gray crosses are Ross et al. (1986) XRD
comparison data; gray triangles are Shimizu et al. (2001). Neither gray series
is included in the Chen dataset.

`scripts/digitize_chen_2010_argon.py` extracts the vector geometry in the
checksummed PDF, page index 2, using PyMuPDF. Main-axis calibration in PDF
points (origin at page top left):

- P = 0 and 25 GPa: x = 374.885986328125 and 513.6329956054688.
- ρ = 1.5 and 4.0 g/cm³: y = 230.6690673828125 and 81.11602783203125.
- Paths 297–444: black filled squares, side approximately 4.646 points.
- Paths 291–296: matching error-bar segments.
- Path 445: the continuous published fitted line, 261 segments/262 vertices.

The square centers yield **80 distinct PDF positions**, spanning approximately
1.232–26.062 GPa. Repeated graphics commands within 0.005 PDF points are
deduplicated, with rendering multiplicity retained in a separate column.
Closely spaced distinct centers are preserved. Their count is **not** an
original sample count or a claim of statistically independent observations.
The lower plot-derived endpoint differs from the rounded headline 1.3 GPa;
neither was silently altered.

All 80 centers have matched horizontal and vertical graphical error bars.
Halfwidths are retained as plot uncertainties, with confidence convention
unspecified. Conservative additional coordinate-reading bounds are 0.05 GPa
and 0.01 g/cm³; these are not fitted statistical errors. The source errors,
integration correlations, shared anchor and common calibration must be
considered before using an independent-error regression.

The **262 curve vertices are a separate `published_fit_curve` dataset**. They
must render as a line labeled as a published fit, never as measured points or
as a first-principles theory dataset. PDF coordinates are included for both
assets so the extraction is auditable. Figure-derived values are not the
original source tables. Source article/figure rights remain APS's; CC0 applies
only to contributor-held factual transcription and arrangement rights.

## Why an executable published EOS is withheld

We tested an actual nonzero-pressure local reference, not a shifted shortcut.
Every standard BM3 pressure expression can be written, with
\(x=\rho/\rho_2\), as

\[
P(x)=c_3x^3+c_{7/3}x^{7/3}+c_{5/3}x^{5/3}.
\]

Let \(D=x\,d/dx\). Then \(K=DP\) and \(KK'=D^2P\).
Solving these three constraints at x = 1 for the printed P = 2, K = 15.1,
K′ = 5.4 gives **c = (32.5325, −47.415, 16.8825) GPa**. This is a
mathematical reconstruction of rounded local constraints, not a published
coefficient list. The existing native BM3 evaluator reproduces the polynomial
within 8×10⁻¹⁴ GPa over the audit grid.

The stable zero-pressure branch gives ρ₀ = **1.674647 g/cm³**, K₀ =
**2.579673 GPa**, K₀′ = **9.081605**, and cell V₀ = **158.445809 Å³**.
Its difference from the printed ρ₀ is 3.09 times the printed 0.05 uncertainty;
this ratio is **not a significance test**, because parameter correlations and
confidence conventions are unknown. It predicts ρ(1.3 GPa) = 2.066389 rather
than the integration anchor 2.02, and KK″(2 GPa) = −5.703635 rather than −7.3.
More decisively, against the 80 drawn densities it gives pressure RMS error
**6.880 GPa**, reaching **19.068 GPa**.

The published line itself gives ρ(2 GPa) ≈ **2.18017 g/cm³**, agreeing with
the text. A local cubic interpolation over 1–3 GPa gives K(2 GPa) ≈
**11.157 GPa**, which differs substantially from the table's 15.1 GPa.
This derivative is only a diagnostic of the drawn line, not an independently
measured modulus.

We also explicitly tested the alternative convention
P = 2 GPa + BM3(V; Vref = Mcell/2.18, Kref = 15.1, K′ref = 5.4).
This pressure-increment formulation gives KK″ = **−7.248889** at 2 GPa,
consistent with the rounded printed −7.3. Thus the second-derivative result
suggests that this alternative reference convention may have been used.
However, it still predicts ρ₀ = **1.633737 g/cm³** and misses the plotted
densities by **6.251 GPa RMS**, with maximum pressure error **16.972 GPa**.
The paper does not print the pressure expression needed to resolve its
convention. Neither interpretation recovers the source density curve; the
alternative is recorded as a diagnostic, not silently adopted as a solution.

We therefore do not export the diagnostic coefficients as a Chen EOS or
silently choose a different equation. The study status is **not_reproduced**.
No full original-data refit or validated parameter covariance is claimed.

## Reproduction and missing inputs

Run `python -m scripts.reproduce_chen_2010_argon` to regenerate
`docs/data/chen-2010-argon-reproduction.json` from the bundled assets.
Run `python scripts/digitize_chen_2010_argon.py /absolute/path/to/chen-2010.pdf`
to reproduce source extraction (optional PyMuPDF dependency).
Tests are in `tests/test_chen_2010_argon.py`.

Checked: complete publisher paper and reference list, APS article metadata,
author-uploaded complete text on ResearchGate, and the local Zotero library.
No supplement is cited in the paper or linked in the accessed APS metadata;
no original numerical observation table was found. A source equation with
its exact reference convention and unrounded coefficients, or original
velocity/density data and fitting procedure, would resolve the ambiguity.
No authors were contacted and no purchase was made.

Packaged supporting-study metadata: `peritheos/data/studies/argon-fcc-chen-2010.json`.
All datasets have empty `used_by_eos_records`; they are retained for comparison
and future audit, without attaching them to unrelated Dewaele or Ono fits.
