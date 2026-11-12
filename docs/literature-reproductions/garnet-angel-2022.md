# Pyrope and grossular: Angel et al. (2022)

Audit and integration: 8 November 2026. Source: Angel, Gilio, Mazzucchelli and
Alvaro, *Garnet EoS: a critical review and synthesis*, Contributions to
Mineralogy and Petrology **177**, 54,
[doi:10.1007/s00410-022-01918-5](https://doi.org/10.1007/s00410-022-01918-5).

**Correspondence status:** Clemens confirmed on 12 November 2026 that he sent
the pyrope input/analysis request to Ross Angel. **Waiting for Ross's reply.**
See the [request tracker](../dataset-requests.md#angel-et-al-2022-pyrope-and-grossular)
and [user-confirmed receipt](../data/garnet-angel-2022/pyrope/correspondence.json).
The recovery and verification receipts below describe work performed before
that email; their historical contact statements are retained in that context.

**12 November follow-up:** 99 of the 102 previously missing pyrope channels
have now been recovered, including 26 digitized channels. Two expanded
conditional fits use 182 and 191 of the 194 channels listed in S2. The
[expanded recovery and fit](#expanded-pyrope-recovery-and-fit-12-november-2026)
still do not reproduce the published thermal coefficients. Three temperature
channels, the increased Hartwig weights and exact author input reductions
remain unresolved. No author was contacted and no catalog coefficients changed.

The subsequent [user-requested Hartwig weight sensitivity](#user-requested-hartwig-weight-sensitivity)
shows that increased weighting explains much of the thermal-parameter shift,
but no tested multiplier reproduces both published thermal central values.
Direct same-coefficient property checks at 27 states agree with official
EosFit; the original fitting workflow remains unrecovered.

## Grossular follow-up: 12 November 2026

The [grossular observation audit](grossular-angel-2022.md) and its
[recovery follow-up](grossular-angel-2022-recovery.md) now cover 213 tabulated
or table-derived channels plus seven digitized Thiéblot plotted points under
the nominal 220-channel S1 inventory. All 52 Gwanmesia modulus rows are recovered.
The Isaak plotted expansion series matches GR-2B strain despite S1's GGD
composition label. Original Thiéblot values/errors, exact normalization and fit
inputs remain incomplete. Nine same-coefficient official EosFit property checks
pass. The separately declared 95-channel partial refit remains diagnostic:
K0=167.373 GPa, K0′=4.96395, θ0=744.384 K and γ0=1.16605. No source coefficients
changed and the full author workflow is not reproduced.

The subsequent [147-channel fit](grossular-angel-2022-expanded-fit.md) adds
the 52 Gwanmesia moduli with explicitly declared relative-error weights.
It converges to K0=167.353 GPa, K0′=4.98267, θ0=740.183 K and γ0=1.16467.
Three coefficients move slightly closer to the source; θ moves farther away.
The Gwanmesia calculated-modulus multiplier is 1.00547 versus published 1.00551.

The [remaining-channel investigation](grossular-angel-2022-remaining-channels.md)
subsequently recovers additional source errors and normalization. A declared
193-row diagnostic converges to K0=166.754 GPa, K0′=5.33063, θ0=780.407 K
and γ0=1.15522. K′ and θ move farther from the source; the original fit remains
unreproduced and no weights or normalization are tuned after this result.

The subsequent user-requested [full-inventory sensitivity study](grossular-angel-2022-full-sensitivity.md)
fits all 219 non-anchor rows with explicit estimates for the 26 missing-error
observations. Its baseline is K0=167.016 GPa, K0′=5.23585, θ0=779.754 K
and γ0=1.15783. All 29 declared cases converge, with θ0=773.477–781.417 K
and K0′=5.15975–5.29031. These tested weights do not reproduce the source
central θ0 or K0′; all earlier scientific qualifications remain.

The [normalization/specimen factorial test](grossular-angel-2022-normalization-specimen.md)
subsequently isolates these input assumptions at fixed error weights.
Normalization changes move θ by tens of kelvin; the 15 Isaak expansion rows
alone have smaller effects under plotted processing. No tested combination
reproduces the central K′, which remains 5.19971–5.24191. The close θ match
with relaxed offsets is a diagnostic control, not recovered author processing.

## Decision and comparison scope

The two published **q-compromise pressure parameterizations** are accepted:
`pyrope_angel_2022_bm3_q_compromise_debye` and
`grossular_angel_2022_bm3_q_compromise_debye`. The existing Milani (2015)
pyrope record remains unchanged and retains its automatic default position.
Grossular is a new material with separately sourced diffraction structure.

Source-equation checks pass against independent SciPy quadrature, the paper's
Table 1 thermal expansivities, and the actual official EosFit7c calculator.
These establish the represented pressure equations and numerical conventions;
they do not independently validate the experiments or reproduce the original
joint regression. The reconstruction ledger therefore uses `not_refittable`,
with `source_equation` scope in each record, rather than coefficient parity.

The exact mixed volume/adiabatic-modulus input project is not supplied by the
four official EOS supplements or the supplementary PDF. The latter is a
selection inventory, not numerical observations. Some upstream measurements
are public, including the separately retained Milani (2015) pyrope PV table;
that table is not the pressure-rescaled, volume-normalized 2022 fit input.
No reconstructed observations or replacement fit are added to the catalog.
The explicitly conditional diagnostic recoveries and refits below remain
separate from published material records.

## Source versions and coefficients

The authoritative numerical inputs are official supplementary files **1**
(grossular) and **2** (pyrope). The unmodified files, including **all 59
parameter tokens and the full 59 × 59 covariance**, are retained under
[data/garnet-angel-2022/sources](../data/garnet-angel-2022/sources/pyrope-official.eos).
The [source manifest](../data/garnet-angel-2022/source-manifest.json) identifies
all downloaded versions and checksums.

| Quantity | Pyrope | Grossular | Role |
|---|---:|---:|---|
| Molar V0, cm³/mol | 113.13 | 125.35 | Fixed adoption from Holland–Powell (2011); no error supplied |
| K0T, GPa | 169.3 ± 0.3 | 167.00 ± 0.20 | Fitted |
| K0′ | 4.55 ± 0.05 | 5.07 ± 0.08 | Fitted |
| θ0, K | 771 ± 28 | 750 ± 13 | Fitted |
| γ0 | 1.185 ± 0.012 | 1.156 ± 0.006 | Fitted; printed Table 1 error widths |
| Reference temperature, K | 298 | 298 | Fixed |
| Atoms per formula unit | 20 | 20 | Formula normalization |
| Conventional-cell formula units Z | 8 | 8 | Crystallographic normalization |

The source EOS files supply finer central values, retained literally in
EOSMAT: pyrope K0=169.339, K0′=4.55314, θ0=771.096, γ0=1.18470;
grossular K0=167.002, K0′=5.06534, θ0=749.825, γ0=1.15602. These are not
independent Peritheos refits.

The author's current downloads have March2023 filenames. Their represented
pressure coefficients are identical to the official supplements. Grossular's
saved **shear modulus** changes from 107 to 106.42 GPa; this is outside the
integrated pressure model. We retain the official files and do not import a
shear model. Pyrope's 59 parameter values are unchanged. The paper's full-MGD
alternatives in Table 2 remain distinct published comparisons; they are not
silently substituted for the selected q-compromise records.

## Equation, units and capability boundary

For molar volume Vm in cm³/mol and Debye thermal energy ED in J/mol,

\[
P(V,T)=P_{\rm BM3}(V;V_0,K_0,K'_0)
 +\frac{\gamma_0}{1000 V_{m0}}
 [E_D(\theta_0,T)-E_D(\theta_0,298)],
\]

\[
E_D(\theta,T)=9nRT(T/\theta)^3
\int_0^{\theta/T}\frac{x^3}{e^x-1}\,dx.
\]

Both θ0 and γ/V are constant. **There is no q parameter.** A full MGD model
with q fixed to zero or one does not reproduce both conditions. This is why
Peritheos supplies the explicit `QCompromiseDebye` class in Python, Rust and
EOSMAT instead of selecting an approximate existing MGD setting.

EOSMAT uses conventional-cell Å³, with
`Vcell = Vm × 8 × 10^24 / NA` and exact `NA = 6.02214076e23 mol⁻¹`.
The internal molar volume is `Vm/10` J/bar/mol. Twenty atoms belong to one
Mg3Al2Si3O12 or Ca3Al2Si3O12 formula unit; Z=8 is not the oscillator count.
The adopted EOS reference volumes are separate from measured structural cells.

This approximation provides P–V–T, derivatives, inversion, fitting and
parameter covariance propagation. It does not supply a full thermodynamic
potential, heat capacities, entropy, shear elasticity or compositional mixing.
Its thermal pressure is volume-independent, so KT at fixed V is independent
of T and the reference increment is exactly zero at 298 K. Caloric operations
are explicitly unsupported. This boundary avoids presenting the simultaneous
constant-θ and constant-γ/V approximations as an integrated quasiharmonic
free-energy model.

## Numerical checks

The repository-only command `python -m scripts.audit_garnet_angel_2022` uses
independent SciPy quadrature in SI/molar units and compares native evaluation
at 20 mathematical states for each material: V/V0=0.85, 0.95, 1.00, 1.03 and
T=90, 298, 700, 1000, 1700 K. The maximum differences are 1.28 × 10⁻¹⁰ GPa
(pyrope) and 1.13 × 10⁻¹⁰ GPa (grossular). These are calculated checkpoints,
including deliberate extrapolation tests, **not observations or measured
phase-stability limits**. Full-precision values are in the
[numerical audit](../data/garnet-angel-2022/audit.json).

Independent differentiation gives α298=2.2564254 × 10⁻⁵ K⁻¹ for pyrope and
2.0474093 × 10⁻⁵ K⁻¹ for grossular, agreeing with Table 1 values
2.256(5) and 2.047(3) × 10⁻⁵ K⁻¹.

The unchanged official Mac EosFit7c **7.60, 13 May 2021** executable was run
on both unchanged source files in a temporary authenticated display. X11
libraries/server were extracted and relocated locally; the EOS calculator
executable was not modified and no system installation was performed.
The [calculator receipt](../data/garnet-angel-2022/eosfit/checks.json) records
its checksum, all 40 requested states, returned pressures and command files.
`python -m scripts.run_garnet_angel_2022_eosfit --executable /path/to/eosfit7c`
replays the checks when the external console and its display are available.

EosFit uses the rounded gas constant (the binary contains the float32 8.314
literal, consistent with the independently audited CrysFML source). Peritheos
retains its documented R=8.31446261815324. The largest raw pressure offsets
are 0.000432 GPa and 0.000387 GPa. A separately labelled comparison using the
source gas constant reduces the maxima to 0.0000703 and 0.0000601 GPa.
The declared 0.0001 GPa calculator tolerance covers 0.00005 GPa output
rounding plus 0.00005 GPa for its finite-precision evaluation. No material
coefficient or atom count is adjusted to improve agreement, and bitwise
identity is not claimed.

![Calculated EOS comparisons; pressure on horizontal axis](../data/garnet-angel-2022/calculated-comparison.svg#only-light)
![Calculated EOS comparisons; pressure on horizontal axis](../data/garnet-angel-2022/calculated-comparison-dark.svg#only-dark)

Lines are independent source-law calculations at 298 and 1000 K. Points are
EosFit calculated checkpoints, not diffraction observations or refit residuals.
The bottom panels show EosFit minus the source-constant calculation, in GPa.
No observation residual RMS is available for the unrecovered joint fit.
Recreate the figure with the repository-only command
`python -m scripts.plot_garnet_angel_2022`.

## Covariance and uncertainty qualifications

The four-parameter covariance order is KT0, K0′, θ0, γ0. Its positive-definite
submatrix is copied directly from source indices 2, 3, 11, 18. It retains the
marginal covariance of these parameters from the joint fit, including their
correlations; excluded nuisance scale-factor parameters are not treated as new
independent uncertainties. V0 was fixed, so its zero covariance row does not
establish that its adopted value is exact or experimentally error-free.

Printed Table 1 errors remain unchanged alongside the finer covariance.
Grossular's printed γ0 width is 0.006, while the source-file covariance gives
sqrt(0.00002598)=0.0050971. That difference is recorded rather than silently
resolved. The source confidence convention and absolute/residual scaling are
not inferred. Propagated errors are conditional on the source model, fixed V0,
selections, weighting and calibration, and do not include total experimental
or inter-study systematic uncertainty.

## Data selection, calibration and validity

The source simultaneously fitted V and adiabatic KS, with KS/KT conversion
inside EosFit, and assessed consistency using heat-capacity data. It normalized
volume series to fixed endmember V0 and refined some group scale factors.

For grossular, supplementary S1 retains selected room-pressure expansion,
seven pressure-volume points, Bass/Isaak/Conrad elasticity and the Kono and
Gwanmesia high-P/T moduli; the latter have scale factors near 0.995 and 1.006.
Several simultaneous high-P/T volume series were rejected as too scattered;
they must not be represented as observations directly fitted by this EOS.
Supplementary S1 labels the seven-point Milani compression row as 2015,
whereas the main article and grossular EOS-file comments cite 2017. That
source-label discrepancy is retained; no pyrope measurements from 2015 are
reassigned to grossular.

For pyrope, S2 includes selected Zhang helium/neon compression, Du and Milani
PV, Milani and Hartwig–Galkin expansion, and selected elastic series. The
source increased Hartwig's weights; the exact factors are not recovered.
The raw Milani (2015) table already in Peritheos contains 12 PV rows; the
2022 inventory reports 11 used after normalization. No mask or pressure
reduction is invented to make those rows into the original fit.

Milani's quartz pressures were recalculated from Angel (1997) to Scheidl
(2016). Other source series use ruby/Au/NaCl/MgO standards; complete paired
readings and joint reductions remain unavailable. Pressure provenance is
therefore partially resolved, with no replacement universal scale.

The marginal selected-constraint envelopes are 0–33.4 GPa/100–1700 K for
pyrope and 0–16.9 GPa/20–1650 K for grossular. They combine separate series,
not measured rectangles or stability maps. No transition, decomposition,
nonhydrostatic-stress or composition correction is included.

## Diffraction identity and preservation

Pyrope's existing Armbruster-based structure and all existing source records,
errors and observations are unchanged. Grossular's nominal Fe-free structure
comes from Wang et al. (2019), *Minerals* 9, 691,
[doi:10.3390/min9110691](https://doi.org/10.3390/min9110691), independently of
its EOS provenance. Tables 3 and 5 were visually verified. The quenched
single-crystal cell is a=11.8418(3) Å; source V=1660.55(11) Å³ agrees with
a³=1660.5546 Å³. This is not an in-situ 3 GPa cell.

In Ia-3d (230), Ca24c, Al16a, Si24d and O96h give Ca24Al16Si24O96, hence
Z=8. The source fixes Ca=(1/8,0,1/4), Al=(0,0,0), Si=(3/8,0,1/4);
O=(0.0384,0.5453,0.3487), each coordinate esd 0.0001. Full nominal
occupancies follow the Fe-free endmember model. Table 4's measured small
stoichiometric deviations do not establish a refined defect model.

The [integration preservation receipt](../data/garnet-angel-2022/integration-preservation.json)
binds the new scientific projections. Removing only those exact additions
recovers the preceding full scientific catalog fingerprint. Historical
scientific receipts remain unchanged.

## Remaining inputs

The complete original joint project is still needed for author-workflow and
coefficient/uncertainty reproduction: all numerical V/KS inputs, included-row
identities, normalization anchors, revised pressures, refined scale groups,
assigned errors, Hartwig weights, fitting controls and confidence convention.
The official article, all five official supplements, current author parameter
page and both author EOS downloads were checked. The Milani (2017) university
manuscript URL was checked but not retrieved (access failed); no assertion
that its numerical tables are absent is made. No author has been contacted.
See the [dataset request](../dataset-requests.md#angel-et-al-2022-pyrope-and-grossular).

Both source publications are CC BY 4.0. Exact source terms, numerical
transcription scope, attribution and exclusions are recorded in the
[scoped notice](../licensing.md#garnet-angel-2022).

## Pyrope follow-up: recovered observations

The follow-up on 8 November 2026 recovered **33 Milani expansion rows and
40 Zou high-pressure/high-temperature rows**. This is a partial observation
audit of the unchanged source EOS, not a reconstruction of its joint fit.
The pressure-equation acceptance above does not establish agreement within
every upstream measurement's reported precision. No coefficients, defaults,
catalog classifications or grossular data were changed in this follow-up.

### Original tables and recovery limits

The official [Milani (2015) supplementary DOC](https://ars.els-cdn.com/content/image/1-s2.0-S0024493715001097-mmc1.doc)
contains Appendix A, Table A1: all 33 temperature, lattice-parameter and
cell-volume rows from 291 to 1073 K. Parentheses explicitly denote standard
deviations; temperature errors are not tabulated. The table has measurements
at 291 and 303 K, with no 298 K measurement. The final primary article was
read from an existing local Zotero attachment. The DOC was extracted with
macOS textutil and independently rendered with LibreOffice; all rows were
checked visually. The [transcription](../data/garnet-angel-2022/pyrope/milani-2015-table-a1.csv)
retains original tokens, values, errors and row order. It complements the
previously retained compression Table A2; neither table is silently replaced
by normalized 2022 inputs.

An [author-hosted copy of Zou (2012b)](https://www.mattwhitaker.net/docs/zou2012b.pdf),
[doi:10.1063/1.4736407](https://doi.org/10.1063/1.4736407), provides all
40 rows in Table I on article page 014910-4 (PDF page 5). The
[transcription](../data/garnet-angel-2022/pyrope/zou-2012b-table-i.csv)
preserves all nine columns and parenthetic errors. Its printed `V0` column
contains the volume at each experimental state, not a constant ambient
volume. The velocity-derived moduli, density and velocities are not separate
independent observations. Row 37's density token `3.8086(1)` is preserved
literally. The confidence convention of these errors is not newly assigned.

Download hashes, extraction checks, attribution, source terms and failed
retrievals are in the [source manifest](../data/garnet-angel-2022/pyrope/sources.json).
Publisher DOC/PDF files stay in private audit output; these factual numerical
transcriptions do not grant rights in the publisher documents. These are
repository audit datasets, not newly registered installed-package datasets.

The [Hartwig–Galkin publisher page](https://link.springer.com/article/10.1007/s10973-020-09396-2)
and its PDF endpoint initially returned subscription-preview HTML. A later
check of the publisher's **separate Table 3 URL recovered all 12 measured
rows**, with the qualifications below. The complete Sinogeikin–Bass (2002) primary article,
[doi:10.1016/S0012-821X(02)00851-8](https://doi.org/10.1016/S0012-821X(02)00851-8),
was recovered. Its tables give fitted elastic constants and temperature
slopes; the individual pyrope KS(T) points appear in Figure 2 rather than a
numerical 12-row table. They have not been digitized or synthesized from the
reported slope. The Zhang (1998) publisher PDF endpoint returned preview
HTML; this follow-up did not recover its primary numerical table.

The included-row counts below are the authors' Table S2 inventory, not a
newly inferred joint input file. Counts for volume and modulus channels can
refer to the same experimental state.

| Source/channel | S2 rows used | Recovered input status |
|---|---:|---|
| Zhang (1998), helium/neon PV | 19 / 14 | Primary numerical table and exact normalization not recovered |
| Sinogeikin–Bass (2000), KS(P) | 17 | Scale recovered; numerical input/scale-group mask not recovered |
| Sinogeikin–Bass (2002), KS(T) | 12 | Primary article recovered; points graphical, not transcribed |
| Zou (2012b), volumes | 7 | All 40 source rows recovered; only seven 300 K volumes selected in S2 |
| Zou (2012b), KS(P,T) | 40 | All source rows recovered; no modulus scale fitted here |
| Du (2015), PV | 9 | Original rows and excluded-row identity not recovered here |
| Milani (2015), PV | 11 | All 12 original rows already retained; exact revised pressures unavailable |
| Milani (2015), VT | 33 | All 33 original rows recovered; exact 2022 normalization/errors unavailable |
| Chantel (2016), KS / volumes | 10 / 10 | Room-temperature selection described; numerical inputs not recovered here |
| Hartwig–Galkin (2021), VT | 12 | All 12 Table 3 rows recovered; source-used column/corrections and increased-weight multiplier unavailable |

### Hartwig–Galkin Table 3 recovered

The publisher's [Table 3](https://link.springer.com/article/10.1007/s10973-020-09396-2/tables/3)
is publicly accessible even though the article/PDF endpoints return a
subscription preview. The [literal transcription](../data/garnet-angel-2022/pyrope/hartwig-galkin-2021-table-3.csv)
retains all seven columns. Temperature, lattice parameter `a` and `Vmeas`
are kept separate from the source's interpolated/calculated `Vcalc`, `ΔV`,
expansivity and heat capacity. No calculated column is counted as additional
measured data. No measurement errors are tabulated. The measured row
temperatures are **103–650 K**, resolving the tabulated coverage rather than
substituting the abstract's 100–630 K or S2's rounded 100–650 K.
The [source receipt](../data/garnet-angel-2022/pyrope/hartwig-galkin-2021-source.json)
records the original HTML checksum, units and source terms.

The table has apparent internal inconsistencies. For example:

- At 297.5 K, `Vmeas=113.1639` and `Vcalc=113.2946` cm³/mol imply
  ΔV=−0.1307, whereas the printed ΔV is +0.0005 cm³/mol. The Vcalc
  token is repeated in the next row.
- At 600 K, `a=11.4827 Å` gives Vm=113.97064 cm³/mol for Z=8 and
  exact Avogadro normalization, whereas the printed `Vmeas=114.121`
  cm³/mol. Their difference, 0.15036 cm³/mol, exceeds printed rounding.
  The printed `Vcalc=113.9767` and ΔV=−0.0055 are also inconsistent
  with that Vmeas.

Every printed token is retained; none is repaired using a smooth curve,
another column or Table 4. The numerical audit reports all rounding and
arithmetic checks. It does not establish whether an error is confined to the
publisher HTML, appears in the full PDF, or was corrected in the 2022 input.
That distinction still requires an authoritative original/corrected version.
Tables 1/2 are literature comparisons, and Tables 4/5 contain calculated
extrapolations/moduli; they do not replace the measured Table 3 values.

A bounded comparison uses the literal Vmeas column, the unchanged Angel EOS
at 1 bar and one profiled constant molar-volume scale. The scale on the
source's measurement basis is 113.17435 cm³/mol at 298 K; it is not a
replacement for the adopted physical V0. The 12-row RMS residual is
**0.01884 cm³/mol**, with maximum absolute residual **0.04105 cm³/mol**.
This is 0.628 times S2's **initial assigned** 0.03 cm³/mol V error on a
volume-only RMS basis. S2's initial 1% T errors are recorded separately;
this comparison does not reproduce an effective-variance fit or the later
increased weights. It cannot identify the actual weight multiplier, original
normalization or which table values entered the 2022 project.

The exact 2022 EosFit joint input remains unrecovered. The author parameter
page provides a final `.eos` file, and the public Pavia repository record
lists the article PDF; neither exposed the fitted-data file or weight
multiplier in the checked listings. The local official EosFit example files
contain no pyrope/Hartwig project. Crossref's fulltext link still returned a
preview; OpenAlex listed no open full-text location, which is a discovery
result rather than proof that none exists. Browser access to the publisher
failed with `ERR_INTERNET_DISCONNECTED`, so institutional-session access
could not be established. No author was contacted.

### Source adiabatic conversion and scale direction

The source operational conversion is
KS = KT(1 + αγT), with **γ = γ0 V/V0** for q-compromise, while θ remains
constant. This follows the CrysFML `Get_Grun_V` implementation, and was
checked against the unchanged official EosFit7c 7.60 calculator at
(0 GPa, 298 K), (10 GPa, 1000 K), and (19.86 GPa, 1700 K).
Independent V, KT, KS, α and γ agree within the declared output/finite-
precision tolerances, with the calculator's rounded gas constant harmonized
for that comparison only. The
[receipt, commands and full output](../data/garnet-angel-2022/pyrope/eosfit-properties/checks.json)
are retained. This audit-only conversion does not extend the library's
pressure model to a caloric potential or establish quasiharmonic consistency.

Source parameter **53 = 0.97828**, with covariance marginal width
0.00587707, belongs to the Sinogeikin–Bass (2000) modulus scale. The main
article rounds it to 0.978(6), while S2 gives 0.977(6). The saved source file
resolves the numerical input, while the printed difference remains recorded.
The methods specify multiplying **calculated moduli** by this factor to
compare with the original observations. Displaying observations on the EOS
scale therefore requires division by 0.97828 (multiplication by 1.02220223).
Applying 0.97828 directly to the already low observed moduli would reverse
the correction. The exact high-P group membership still requires the input
project; the ambient datum must not automatically be assigned this scale.

### Quantitative comparisons and what they establish

Run the repository-only command
`python -m scripts.audit_pyrope_angel_2022_observations` to regenerate the
[full residual audit](../data/garnet-angel-2022/pyrope/audit.json) and figure.
All four EOS coefficients and adopted V0 are fixed to the source values.
The observation calculations use Peritheos's documented modern gas constant.

| Comparison | RMS residual | Largest absolute residual | Scope |
|---|---:|---:|---|
| Milani VT, 33 rows | 0.43056 Å³/cell | 1.24876 Å³/cell | One constant instrumental volume scale profiled; reported volume std only |
| Zou KS, 40 rows | 0.81543 GPa | 1.89291 GPa | Each measured P,T; no fitted modulus scale |
| Zou volumes, seven 300 K rows | 0.11771 GPa | 0.15941 GPa | P residual using published V/V0 × adopted molar V0 |

For Milani, the sole fitted quantity is an instrumental cell-volume scale,
1504.41429 Å³ at 298 K, minimizing diagonal squared volume residuals with
the reported volume standard deviations. It is **not** a new physical V0,
a refitted EOS, or a recovered author normalization. No temperature errors
are assigned, and normalization-anchor/cross-row covariance is unavailable.
The standardized volume-residual RMS is 3.316 on this limited error basis.
The residuals have a systematic temperature trend that a constant scale
cannot remove: the directly measured V1073/V303 is **1.02236415**, while
the unchanged EOS predicts **1.02323217**. A constant instrumental scale
cancels from this comparison. Its diagonal, volume-only propagated error
is 0.00012998; missing temperature errors, calibration covariance and EOS
parameter uncertainty prevent treating this as a joint-fit significance.
It establishes tension at the reported volume precision, rather than exact
agreement of the source EOS with every retained expansion datum.

For Zou, KT and KS are calculated at the measured P,T using the EOS volume,
not the high-temperature diffraction volumes excluded by S2. The KS residual
RMS is smaller than the typical reported modulus errors. A separately
labelled diagonal propagation of the reported P/KS errors and S2's assigned
2% T errors gives a standardized RMS of 0.536. This is a diagnostic in KS,
not EosFit's pressure-residual objective or a reproduced χ². The seven
room-temperature volume checks use each published V/V0 directly; rounding
prevents recovery of a unique exact original ambient V0 and its uncertainty.
No extra scale or uncertainty is invented to remove their pressure offset.

![Pyrope observation residuals](../data/garnet-angel-2022/pyrope/observation-residuals.svg#only-light)
![Pyrope observation residuals](../data/garnet-angel-2022/pyrope/observation-residuals-dark.svg#only-dark)

The bars show only the reported volume or modulus errors. The EOS is fixed;
only the one stated Milani volume scale is profiled. Pressure and assigned
temperature errors are included in the separate numerical diagnostic, not
these plotted bars.

This completed partial audit supports the Zou modulus subset and exposes
expansion-data tension. It does **not** reproduce the source's joint
χ² = 2.24, coefficient uncertainties or workflow. A compromise joint fit
need not minimize each individual dataset. The source-used Hartwig values and
weight multiplier, normalization anchors, revised quartz pressures, other
selected inputs and refinement controls remain required for that replay.
There is no demonstrated failure of the implemented source pressure law;
there is also no basis here to promote the record to full observation/workflow
reproduction or make a stronger recommendation. Further work should recover
those specific inputs, rather than tune new weight/model hypotheses to
eliminate these residuals.

## Conditional fits of the recovered pyrope subset

At the user's request, four bounded source-model refits were performed before
author contact. They are numerical reproduction investigations, **not new
catalog EOSs or replacements**. The [fit plan](../data/garnet-angel-2022/pyrope/refit-plan.json)
was saved before refinement; [full results and row-level residuals](../data/garnet-angel-2022/pyrope/subset-refits.json)
retain the inputs, source hashes, convergence traces, covariance assumptions
and fitted instrumental scales.

The model is unchanged BM3 plus q-compromise Debye, with fixed
Vm0=113.13 cm³/mol, n=20 and Tr=298 K. Source-model fits use the EosFit
rounded gas constant, float32(8.314), rather than altering any atom count.
The separate expansion fits fix the published KT0=169.339 GPa and
K0′=4.55314; the joint recovered-subset fit refines both along with θ0 and
γ0. The published thermal coefficients are θ0=771.096 K and γ0=1.18470.

| Constructed fit | Data channels | θ0 (K) | γ0 | Conditional reduced χ² |
|---|---:|---:|---:|---:|
| Milani VT alone | 33 | 302.976 | 1.075144 | 0.80484 |
| Hartwig VT alone | 12 | 815.661 | 1.258032 | 0.05811 |
| Both expansion tables | 45 | 495.476 | 1.097269 | 1.22307 |
| Expansion + Zou PV/KS | 92 | 483.318 | 1.075930 | 0.64631 |

The joint fit gives **KT0=166.73850 GPa** and **K0′=4.790252**.
Its local, conditional unscaled widths are 0.71485 GPa, 0.07050,
35.551 K and 0.008060 for KT0, K0′, θ0 and γ0, respectively. These
are conditional on the assumptions below, not the source covariance or
total systematic errors. All four fits converged within the declared
iteration/evaluation limits, without reaching a parameter bound.

![Conditional thermal coefficient fits](../data/garnet-angel-2022/pyrope/subset-coefficients.svg#only-light)
![Conditional thermal coefficient fits](../data/garnet-angel-2022/pyrope/subset-coefficients-dark.svg#only-dark)

### Explicit subset and objective

The 92 fitted data channels comprise 33 Milani volumes, 12 literal Hartwig
`Vmeas` values, seven Zou room-temperature volumes and 40 Zou adiabatic
moduli. A volume and modulus at the same state are separate fitted channels,
not independent experimental runs. Supplementary S2 tabulates 194 selected
volume/modulus channels in total. The missing 102 in this replay are the
Zhang PV, Sinogeikin–Bass P–KS and T–KS, Du PV, recalibrated Milani PV,
and Chantel room-temperature PV/KS inputs. The raw Milani compression table
is available, but it is not silently substituted for its revised 2022
pressures.

The objective is squared residuals in **pressure**, using iteratively updated
effective variance:

σeff² = σP² + (∂P/∂Y σY)² + (∂P/∂T σT)²,

where Y is observed volume or KS. For a KS datum, the source's operational
KS(V,T) law is inverted at the measured KS and T to calculate P. It never
uses the excluded high-temperature diffraction volume. The effective weights
are frozen within each nonlinear solve and recomputed between solves.
No parameter derivative of the weight is introduced into that inner objective.

The finite assumptions are:

- Milani uses the reported volume standard deviations. Its unreported
  temperature/pressure errors are **omitted variance terms**, not asserted
  physical zero errors.
- Hartwig uses S2's initial 1% T and 0.03 cm³/mol V errors. **No increased
  Hartwig multiplier is supplied or fitted.** No lattice-derived or smoothed
  values replace the literal measured-volume column.
- Zou uses the tabulated P, relative-volume and KS errors, with S2's assigned
  2% T errors. The seven volumes and their errors are V/V0 and its error
  multiplied by the adopted molar V0. No extra modulus scale is fitted.
- Each expansion series has one fitted constant instrumental volume scale.
  This is a conditional alternative to the unrecovered original normalization
  anchors, not their reconstruction. The source physical V0 stays fixed.
- All coordinate-error and inter-row covariances are unavailable and omitted.
  No unreported uncertainty floor, alternative mask, weight search or new
  thermal model is introduced.

The four reported reduced χ² values therefore belong only to these stated
subsets and error protocols. They cannot be compared as reproductions of
the article's joint χ²=2.24. Residual-scaled widths, where applicable,
use max(1, χ²/dof), following the documented EosFit convention; the full
unscaled covariance is retained separately. No source confidence level is
assigned to either.

![Recovered-subset pressure residuals](../data/garnet-angel-2022/pyrope/subset-residuals.svg#only-light)
![Recovered-subset pressure residuals](../data/garnet-angel-2022/pyrope/subset-residuals-dark.svg#only-dark)

Gray crosses replay the published EOS with the same two instrumental scales
profiled under the stated protocol. Blue points use the joint-subset refit;
bars are the corresponding effective pressure errors. The four data groups
retain their individual coordinates; mixed-temperature data are not joined
as an isotherm.

### Official EosFit replay

The unchanged official **EosFit7c 7.60** executable was also given this exact
constructed 92-channel subset, its coordinate errors and volume-scale groups.
The [receipt](../data/garnet-angel-2022/pyrope/subset-eosfit/checks.json)
links the saved inputs, row map, commands, complete console output and
diagnostic `.eos` output. It reports numerical convergence with:

| Coefficient | Independent IRLS | EosFit saved output |
|---|---:|---:|
| KT0 (GPa) | 166.738503 | 166.743 |
| K0′ | 4.790252 | 4.78982 |
| θ0 (K) | 483.318298 | 483.401 |
| γ0 | 1.075930 | 1.07597 |
| Conditional reduced χ² | 0.646308 | 0.6467 |

The nearby numerical results agree within 1% of the conditional marginal
widths plus half the saved parameter's last digit; reduced χ² agrees within
0.001. These tolerances account for default shift/esd termination, single-
precision evaluation and saved-file precision. In particular, the instrumental
scales are saved to only five decimal places. The raw differences and rounding
allowances are retained; bitwise agreement or equal covariance is not claimed.

The console warns about some zero coordinate errors: these are solely its
input encoding for the omitted P/T variance terms described above. Every
total effective pressure variance is positive in the independent check.
No source measurement error is rewritten to suppress the warning.

Repository-only reproduction commands:

```console
python -m scripts.refit_pyrope_angel_2022_subset
python -m scripts.plot_pyrope_subset_refits
python -m scripts.run_pyrope_subset_eosfit --executable /path/to/eosfit7c
```

The last command requires an available external console/runtime/display.
It does not install software or obtain missing original inputs.

**Finding:** the conditional recovered-subset fit is numerically confirmed
by the official console, but it does not reproduce the published coefficients.
The expansion tables constrain different thermal parameters under their
reported/initial error assumptions. The missing input channels, original
normalization and increased Hartwig weighting prevent attribution of the
coefficient difference to any one cause. The literal Hartwig table conflicts
remain unresolved. No weight multiplier was tuned to recover the publication,
and the resulting coefficients are retained solely as diagnostics, without
catalog integration or a recommendation for use. The original source record
and its parameterization-only acceptance remain unchanged; this partial
refit does not establish failure of the complete author fit.

## Integration checks and unrelated baseline drift

Both new records pass source checks, schema/loading, unit conversion,
inversion, covariance preservation, native fitting and Python/Rust checks.
All 685 prior ledger rows are preserved. A fresh read-only replay of those
prior rows exposes existing Wang (2004) akimotoite and Sata (2010) cementite
ledger drift, together with policy-text drift; the exact identifiers and
fields are in the [baseline receipt](../data/garnet-angel-2022/baseline-ledger-drift.json).
These findings predate the garnet entries and are not silently reconciled or
reclassified by this curation.

The [verification receipt](../data/garnet-angel-2022/verification.json) records
the full Python run and exact passing rechecks, the 17 garnet-specific checks,
Rust tests, lint/type checks, strict documentation build and clean-wheel
example. The common-ledger metadata drift remains separately identified.

## Expanded pyrope recovery and fit, 12 November 2026

The requested recovery obtained **99 of the 102 previously missing channels**.
There are now **191 of the 194 channels listed in supplementary S2**, including
26 digitized channels. This is a near-complete conditional reconstruction,
not the original author project or a full 194-channel fit. The original
published coefficients and catalog records remain unchanged. Ross Angel
has not been contacted.

### Newly recovered inputs

| Source group | Additional fitted channels | Recovery and limitations |
|---|---:|---|
| Zhang (1998), helium / neon PV | 19 / 14 | All 39 Table 1 rows retained, including ambient and five explicitly excluded neon rows; 33 selected. Common ambient anchor 1502.9 Å³. |
| Du (2015), PV | 9 | Author preprint Table 4, pages 26–27, from Columbia. All ten pyrope rows retained; ambient 1503.2 Å³ used for normalization. Nine selected points corroborated by the 2022 vector plot. Final publisher table not directly recovered. |
| Milani (2015), PV | 11 | Existing original compression/decompression table retained; nonambient pressures recalculated with the author-distributed old/new quartz EOS files. |
| Chantel (2016), PV / KS | 10 / 10 | Published Table 1 recovered through its indexed primary-paper text. All 16 rows retained, but only the ten nonambient 298 K states selected. KS calculated from the reported density and velocities. |
| Sinogeikin–Bass (2000), KS(P) | 17 | Original open-square series and error bars extracted from Angel (2022) Fig. 8 vectors; **digitized plotted observations**, not the original numerical table. Includes two ambient points. |
| Sinogeikin–Bass (2002), KS(T) | 9 | **Digitized observations**: nine clearly readable C11/C12 pairs from primary Fig. 2; KS=(C11+2C12)/3. Three further channels listed by S2 remain unresolved. |
| **Total additional** | **99** | **73 from tables/recalibration, 26 digitized** |

The [recovery manifest](../data/garnet-angel-2022/pyrope/expanded-input-recovery.json)
records source URLs, versions, checksums, units and error qualifications.
Raw transcriptions preserve the printed tokens; neither plotted fitted lines
nor invented duplicate points fill the three-channel gap. The overlapping
sample symbols and absence of an original numerical table prevent reliable
recovery of every Sinogeikin–Bass (2002) C11/C12 pair.

The Zhang table was recovered from the primary article's public full text on
[ResearchGate](https://www.researchgate.net/publication/225719091_Hydrostatic_compression_and_crystal_structure_of_pyrope_to_33_GPa).
The Du preprint is available through
[Columbia Academic Commons](https://academiccommons.columbia.edu/doi/10.7916/D87P8XSM).
The primary Chantel table is indexed at the
[RRUFF paper endpoint](https://www.rruff.net/odr/view/downloadfile/82036);
direct PDF retrieval returned HTTP 403. The complete Sinogeikin–Bass (2002)
article is retained privately from the
[geoscience repository](https://repository.geologyscience.ru/bitstream/handle/123456789/30263/Sino_02.pdf?isAllowed=y&sequence=1).
Publisher documents and source-figure crops are not redistributed with these
factual numerical audit transcriptions.

### Reduction checks and remaining choices

The [declared protocol](../data/garnet-angel-2022/pyrope/expanded-refit-plan.json)
uses the same source BM3 plus q-compromise Debye model, fixed Vm0=113.13
cm³/mol, Tr=298 K, n=20 and adiabatic-modulus conversion as the earlier subset.
The four EOS coefficients, two expansion-series scales and one Sinogeikin–Bass
(2000) modulus scale are free. The source model is not replaced, and the
increased Hartwig weight multiplier is neither supplied nor fitted.

For PV tables, volumes and numerator errors are normalized as
Vm0·Vcell/Vambient. Zhang neon uses the common ambient volume tabulated with
the helium data. Chantel volumes use Vm0·ρambient/ρ, with ρambient=3.565 g/cm³.
Reference-anchor errors and induced correlations are unavailable and omitted.
Room-temperature-only sources without numerical T are conditionally encoded
at 298 K; this is not a recovered temperature measurement. Missing coordinate
errors remain missing and are omitted from the effective-variance sum.

Chantel provides **two pressure columns**. Extracted vector coordinates in
both 2022 figures match the **Au-calibrated pressures**, 2.3–23.8 GPa, rather
than the internally derived 2.0–25.2 GPa column. S2's quoted 25.2 GPa range
alone would select the wrong column. Both columns are preserved in the
transcription, and the final diagnostic fits use the calibrated column.
The internally derived pressures depend on Chantel's fitted elasticity and
must not be described as independent pressure measurements.

The [figure cross-check](../data/garnet-angel-2022/pyrope/expanded-figure-crosschecks.json)
finds Du's nine preprint points agree with plotted P to 0.02 GPa and normalized
V to 0.004 cm³/mol. It also retains unresolved differences: the conditional
Chantel normalization differs from Fig. 7 by about 0.045 cm³/mol; Milani
volumes differ by up to 0.008 cm³/mol and transformed pressures by up to
0.042 GPa. No anchor or measurement is adjusted to force the plotted values.
These are additional reasons not to claim recovery of the original reduction.

For Milani, the pressure conversion inverts the 1997 quartz BM3 for an
**equivalent calculated quartz volume**, then evaluates the 2016 quartz BM4
at that absolute volume. Both evaluations at all eleven nonambient states
were checked against the unchanged official EosFit calculator:
[22 native PCAL checks](../data/garnet-angel-2022/pyrope/expanded-eosfit/quartz-calibration/checks.json).
All agree within 0.0006 GPa, justified by three-decimal output rounding plus
single-precision allowance. The original measured quartz lattice parameters
were not recovered. Calibration coefficients are held fixed; original
pressure errors are transformed by the new/old derivative ratio. Calibration
covariance is not added as independent measurement noise. EosFit's PSHIFT
command opened an inaccessible dialog in the headless runtime; **PSHIFT's
interactive workflow itself was not replayed**.

The digitized Sinogeikin–Bass (2000) inputs preserve vertical half-error-bar
spans and bounds of approximately 0.01 GPa in P and 0.05 GPa in KS for
coordinate resolution; extra saved digits retain the deterministic extraction,
not experimental precision. Fig. 8 displays **all 17 points scaled**, including
both ambient points. That displayed mask is used conditionally, while S2's
wording refers to high-pressure scaling. The plotted original/scaled ratios
are about 0.9763–0.9766, which also differ slightly from the saved source EOS
scale 0.97828; they are not used to overwrite that source coefficient.
Original fit groups remain unavailable.

For Sinogeikin–Bass (2002), pixel centers and error-bar endpoints are retained
separately. Approximate digitization bounds are 2 K in T and 0.3 GPa per elastic
constant. The experimental T error is encoded as 5 K, as stated for the same
experiment in Fig. 3; the C11/C12 error bars have no independently recovered
confidence level. Their covariance is missing, so the conditional propagated
KS error is sqrt(σ11²+4σ12²)/3. That assumption is not an author covariance.
The nine resolved pairs are not sufficient to recover the original twelve
measurements or their complete statistical dependence.

### Fit results

The [complete numerical results](../data/garnet-angel-2022/pyrope/expanded-refits.json)
retain every prepared observation, effective pressure error, residual,
conditional covariance, convergence trace and input hash. The pressure
objective and iterated effective-variance procedure are unchanged from the
92-channel fit. Both expanded fits converge without reaching parameter bounds.

| Coefficient | Published source | 182-channel independent fit | 191-channel independent fit | 191-channel official EosFit |
|---|---:|---:|---:|---:|
| K0T (GPa) | 169.339 | 169.420 | 169.510 | 169.511 |
| K0′ | 4.55314 | 4.54676 | 4.53504 | 4.53482 |
| θ0 (K) | 771.096 | 492.954 | 496.101 | 498.981 |
| γ0 | 1.18470 | 1.097532 | 1.098700 | 1.09915 |
| Sinogeikin–Bass (2000) KS scale | 0.97828 | 0.981177 | 0.981252 | 0.98135 |
| Conditional reduced χ² | 2.24 in article, different unrecovered protocol | 2.03608 | 1.94316 | 1.9454 |

The 191-channel residual-scaled conditional widths are approximately 0.286 GPa
for K0, 0.0410 for K0′, 48.2 K for θ0 and 0.00825 for γ0. They omit unknown
weight choices, normalization/calibration covariance and other systematic
errors, and do not define statistical significance against the source fit.
The profile calculation with the four published EOS coefficients fixed and
only the three scales free gives reduced χ²≈3.593 on these same 191 inputs.
Neither reduced χ² is a replay of the article's reported 2.24.

![Expanded conditional coefficient fits](../data/garnet-angel-2022/pyrope/expanded-coefficients.svg#only-light)
![Expanded conditional coefficient fits](../data/garnet-angel-2022/pyrope/expanded-coefficients-dark.svg#only-dark)

Adding the recovered compression data brings the cold coefficients close to
the published values. Adding the nine resolved temperature-dependent modulus
channels changes θ0 only from about 493 K to 496 K under this protocol. The
published thermal coefficients are still not obtained. This does not identify
which missing author settings, corrections or observations explain the gap,
and does not establish that the publication is wrong.

![Expanded conditional pressure residuals](../data/garnet-angel-2022/pyrope/expanded-residuals.svg#only-light)
![Expanded conditional pressure residuals](../data/garnet-angel-2022/pyrope/expanded-residuals-dark.svg#only-dark)

### Official solver comparison and reproducibility

The unchanged EosFit7c 7.60 executable was also run on both constructed expanded
inputs. It reports convergence, with fresh saved outputs and full input row
maps retained under [expanded-eosfit](../data/garnet-angel-2022/pyrope/expanded-eosfit/191/checks.json).
The same predeclared comparison threshold as the earlier subset—1% of the
conditional marginal width plus saved half-last-digit, and 0.001 in reduced
χ²—**is not passed** for these expanded cases. The threshold was not relaxed.
For 191 channels the θ difference is about 2.88 K, or 0.083 of its unscaled
conditional width; the reduced χ² difference is about 0.00224.

Inspection of the available CFML EoS source shows that `Get_Pressure_K` uses
a modulus residual stopping tolerance of 0.001 times its target modulus.
Small native/independent pressure and weight differences are present in the
replay. This coarser inversion is a plausible contributor to the solver gap,
not an isolated demonstration of its entire cause. Both implementations give
the same broad thermal discrepancy, but exact solver, coefficient, covariance
or author-workflow parity is not claimed. The tests preserve the failed strict
comparison explicitly; passing them does not convert it into numerical parity.

Repository-only commands:

```bash
python -m scripts.refit_pyrope_angel_2022_expanded
python -m scripts.plot_pyrope_expanded_refits
python -m scripts.digitize_pyrope_angel_2022 --pdf /path/to/original-angel-2022.pdf
python -m scripts.run_pyrope_expanded_eosfit --executable /path/to/eosfit7c \
    --quartz-source-dir /path/to/recovered-quartz-eos-files
```

Digitization additionally requires PyMuPDF; official replay requires the
external executable and its functioning display/runtime. These maintenance
scripts are not installed-package user APIs. Source PDFs, quartz EOS downloads,
failed-access receipts, raster inspection crops and the rejected initial
pressure-column trial remain in the private audit output. No author message
was sent, no replacement EOS was installed, and grossular was not modified by
this follow-up.

## User-requested Hartwig weight sensitivity

On 12 November 2026 the user explicitly requested testing whether the increased
Hartwig weights could explain the thermal-coefficient discrepancy before any
author contact. This separate diagnostic extends the earlier fixed-weight
protocol. It does **not** recover the multiplier actually used by Angel et al.,
change any measured value or replace a catalog EOS.

### Matched test

The [predeclared sensitivity plan](../data/garnet-angel-2022/pyrope/hartwig-weight-plan.json)
uses the same 182- and 191-channel inputs, source equations, calibration,
normalization, seven fitted parameters and baseline coordinate errors. Only
the twelve Hartwig rows receive an inverse-variance multiplier w. The fixed
grid is 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000 and 10000. The four orders
of magnitude test sensitivity, not a presumed plausible range of author
weights. No multiplier is continuously optimized to match the publication.

For a Hartwig row, its effective pressure variance is divided by w by dividing
**both** the assigned volume error and temperature error by sqrt(w) in copied
diagnostic rows. This preserves their relative variance contributions. Thus
500× weight is equivalent to 22.36× smaller assigned errors, not 500× smaller
errors. At w=500 the equivalent volume error is 0.001342 cm³/mol and T error
is 0.04472%; at w=1000 they are 0.000949 cm³/mol and 0.03162%. These are
weighting devices, not newly established experimental uncertainties. If the
authors changed only one coordinate error or individual row weights, that
would be a different protocol, still unrecovered.

Each grid point was fitted from the published coefficients and independently
from the preceding grid solution. All **44 fits** converge without reaching
parameter bounds. The two starts give θ differences below 0.001 K and γ
differences below 0.000001; this checks these two starts, not global uniqueness.
The [complete results](../data/garnet-angel-2022/pyrope/hartwig-weight-sensitivity.json)
retain both fits, traces, conditional covariances and parameter differences.
Baseline files and their hashes are unchanged.

### Results

| Hartwig weight multiplier | K0 (GPa) | K0′ | θ0 (K) | γ0 | Conditional weighted reduced χ² |
|---:|---:|---:|---:|---:|---:|
| 1 | 169.510 | 4.53504 | 496.101 | 1.098700 | 1.94316 |
| 10 | 169.505 | 4.52844 | 619.704 | 1.120361 | 2.28498 |
| 100 | 169.502 | 4.52014 | 678.552 | 1.148357 | 4.20754 |
| 500 | 169.464 | 4.50985 | 738.539 | 1.194131 | 8.66512 |
| 1000 | 169.406 | 4.50952 | 764.679 | 1.215148 | 11.65205 |
| 10000 | 168.262 | 4.62719 | 807.571 | 1.242835 | 40.29597 |
| Published coefficients | 169.339 | 4.55314 | 771.096 | 1.184700 | 2.24 for the original, unrecovered protocol |

These are the 191-channel independent results. The 182-channel path is very
similar: w=500 gives θ0=738.504 K, γ0=1.193990; w=1000 gives θ0=764.645 K,
γ0=1.214975. The nine resolved Sinogeikin–Bass temperature channels therefore
do not materially change this particular weight sensitivity. This does not
establish the influence of the three unresolved channels.

**Increased weighting explains much of the direction and magnitude of the
thermal-parameter shift**, but no tested grid point reproduces both published
thermal central values. At 500× weight γ is near the published value while
θ remains about 33 K lower. At 1000× weight θ is only about 6 K lower while
γ is about 0.030 higher. This is a conditional consistency test, not evidence
that the authors used one of these factors or that weighting alone accounts
for the full publication. Untested intermediate multipliers are not ruled out
by this discrete grid; no continuous matching search was performed.

Changing w changes the statistical objective, so weighted χ² values across w
must **not** be used as a fit-quality ranking under one error model. For a
matched descriptive comparison, every solution was also evaluated with the
original, unchanged coordinate errors. At w=1, Hartwig's pressure residual RMS
is 0.05265 GPa and Milani's is 0.01490 GPa. At w=1000 they become 0.01337 GPa
and 0.07078 GPa, respectively: fitting Hartwig more closely increases the
Milani RMS by about 4.75×. The largest Milani residual increases from 2.49
to 12.94 times its conditional effective pressure error. Unknown Milani
coordinate errors and systematic covariance remain outside that calculation.

![Hartwig weight sensitivity](../data/garnet-angel-2022/pyrope/hartwig-weight-sensitivity.svg#only-light)
![Hartwig weight sensitivity](../data/garnet-angel-2022/pyrope/hartwig-weight-sensitivity-dark.svg#only-dark)

Gray bands show marginal widths from the original source covariance, with no
assigned confidence level. They are descriptive comparison scales, not
acceptance criteria. The coefficient path shows why agreement in θ alone
cannot establish reproduction of the thermal EOS.

### Official EosFit checks

The same w=500 and w=1000 inputs were also fitted inside the unchanged official
**EosFit7c 7.60** executable, with fresh saved outputs:

| Multiplier | Independent θ0 (K) / γ0 | Official EosFit θ0 (K) / γ0 |
|---:|---:|---:|
| 500 | 738.539 / 1.194131 | 738.265 / 1.19420 |
| 1000 | 764.679 / 1.215148 | 764.826 / 1.21534 |

The [official receipts](../data/garnet-angel-2022/pyrope/hartwig-weight-eosfit/500/checks.json)
retain the exact constructed data, commands and saved files. Both native fits
converge and confirm the direction and magnitude of the coefficient shifts.
The earlier strict numerical-comparison criterion remains **not passed**;
for example the w=1000 weighted reduced χ² is 11.65205 independently and
11.7396 in the official console. These differences are retained, not hidden
by rounding the coefficients. They do not change the broad sensitivity result
but preclude exact fitting-solver or covariance parity claims.

The actual author multiplier, any nonuniform error changes, source-used
Hartwig corrections, original normalization and three missing channels remain
unresolved. No author was contacted and no scientific catalog record was
changed by this test.

Repository maintenance commands:

```bash
python -m scripts.audit_pyrope_hartwig_weights
python -m scripts.plot_pyrope_hartwig_weights
python -m scripts.run_pyrope_hartwig_weights_eosfit --executable /path/to/eosfit7c
```

### Are these the same EOS equations as Ross used?

There is strong evidence for the **published EOS formulation**. The model is
read from the unchanged official pyrope `.eos` file: BM3 reference compression
at 298 K, plus q-compromise Debye thermal pressure with constant θ and γ/V,
Vm0=113.13 cm³/mol and n=20 atoms per formula unit. In molar units,

\[
P(V,T)=P_{\rm BM3}(V)+\frac{\gamma_0}{1000V_{m0}}
 \left[E_D(\theta_0,T)-E_D(\theta_0,298\,{\rm K})\right],
\]

with P in GPa, Vm0 in cm³/mol and Debye energy in J/mol. For the source's
operational modulus conversion, γ(V)=γ0·V/V0 and KS=KT(1+αγT). A full MGD
with volume-dependent θ must not be substituted for this q-compromise model.
The source gas constant is also matched in the diagnostic fits.

The earlier 20-state pressure checks match official EosFit to 0.0001 GPa.
To examine the low-T range and changed coefficients directly, a further
[same-coefficient calculator check](../data/garnet-angel-2022/pyrope/hartwig-weight-eosfit/equation-checks/checks.json)
uses three EOS coefficient sets: the published set and the saved native
w=500 and w=1000 diagnostic sets. Each is evaluated at P=0 with T=103,
150, 297.5, 600, 650 and 1073 K; at 10 GPa/300 K and 10 GPa/1000 K; and
at 19.86 GPa/1700 K. **All 27 states pass** the previously used calculator
bounds: 0.0001 cm³/mol in V, 0.001 GPa in KT and KS, 2×10⁻⁹ K⁻¹ in α,
0.0001 in γ, and 0.001 K in constant θ. The θ bound permits single-precision
output at hundreds of kelvin. Maximum KT and KS differences are about
0.000287 and 0.000261 GPa, respectively. Independent SciPy quadrature and
root finding are compared to native evaluations, not to generated observations.

This verifies the represented pressure/property equations over these checks.
It cannot establish that Ross used this exact executable build or recover his
private data reductions, weights, original fit groups, uncertainty handling or
solver settings. The separate weighted fitting comparisons still fail the
strict numerical-parity threshold. Matching the EOS formulation and reproducing
the author's regression are distinct claims; only the former has passed these
calculator checks. The coefficient shifts also occur when the sensitivity
fits are performed inside official EosFit itself, so they are unlikely to be
caused solely by a different independently implemented EOS equation.
