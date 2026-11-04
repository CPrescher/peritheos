# Fei 2007 gold author fit evidence follow up

Evidence checked on 12 October 2026. Source-exact reproduction remains
unresolved. The existing staged equal-pressure-weight replay retains
**K0-prime = 5.995326** from 37 Dewaele cold rows and **q = 0.619149** from
28 reported-MgO hot rows. These are conditional numerical parity within the
printed widths of 6.00(2) and 0.6(3), with V0/K0/gamma0/theta0 adopted and
fixed. No new fit was performed for this follow-up.

The additional evidence strengthens the identification of the two Hirose
observations and documents why later reanalyses cannot establish the original
MgO reducer. It supplies neither the historical reducer nor Fei's regression
objective, weights, selection file or covariance. The
[existing Au reproduction](fei-2007-gold.md) and its
[numerical report](../data/fei-2007-gold-reproduction.json) remain the authority
for the conditional replay.

## Author evidence and source selection

[Fei et al. (2007)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890468/),
pages 9183-9184, explicitly describes the revised-ruby cold anchor and a thermal
refit of Fei (2004) and Hirose (2006), with gamma0 = 2.97 and theta0 = 170 K.
It does not specify the Au residual coordinate, relative weights or uncertainty
calculation. The plotted isotherms are calculated outputs, not additional
observations or evidence of the optimizer.

**Hirose's later paper independently identifies the high-pressure inputs.**
[Hirose et al. (2008)](https://doi.org/10.1016/j.pepi.2008.03.002), section 2.2
on page 151, says the first heating cycle of run 3 measured 1340 and 2330 K
and was previously reported in Hirose (2006). Page 153 states that Fei (2007)'s
megabar thermal-pressure prediction was based on that first cycle. Table 1
repeats the same Au/MgO lattice parameters and reported MgO pressures:

| T (K) | a(Au) (angstrom) | a(MgO) (angstrom) | Reported MgO P (GPa) |
| --- | --- | --- | --- |
| 1340 | 3.7237 | 3.7827 | 106.5(12) |
| 2330 | 3.7254 | 3.7767 | 117.9(21) |

These are repeated observations, not two new independent data points. This
author statement supports the existing two-row Hirose selection, but does not
confirm all 26 Fei rows or their weights. Hirose (2008)'s separate 38-point
Au fit is a different regression; its row count and fit quality cannot define
the 2007 objective.

The [Speziale author confirmation](../data/speziale-2001-author-confirmation.json)
records a user-supplied reply on 9 October 2026, with no archived original
email headers. It confirms the variable-q gamma law in Speziale Eq. 11 and
numerical integration of d ln(theta) = -gamma d ln(V), generally using
trapezoidal integration. The original calculator is unavailable. This is direct
author recollection about the 2001 prescription; its scope does not extend to
Fei/Hirose's later implementations or Fei's Au fit.

## Evidence about later MgO implementations

[Dorogokupets (2010)](https://doi.org/10.1007/s00269-010-0367-2),
*P-V-T equations of state of MgO and thermodynamics*, explicitly distinguishes
at least four implementations called Speziale MgO. The introduction separates
the Fei (2004)/Hirose (2008)/Komabayashi (2008)/Tateno (2009) tabulations from
the Wu, Zha and Gao versions. Figures 3 and 4 and Table 1 fit 4 analyze the
former group. That fit is a reconstruction from published pressures using
another parameterization, not recovered Fei/Hirose code. Its agreement cannot
identify their gamma or theta algorithm. The footnote discussing 51.041 and
24.619 J mol^-1 K^-1 heat-capacity normalizations is evidence of another
implementation detail to ask about, not confirmation of the 2007 Au prefactor.
The existing Speziale audit already documents this paper's general convention
warning; the useful follow-up is the explicit separation of calculator families.

**A later published recalculation must be distinguished from the original
reported pressures.** In the
[Tange et al. (2010) correction](https://doi.org/10.1029/2010JB007959),
paragraph 2, the authors say that their original Figure 11 recalculated Fei
and Hirose data and boundaries with Speziale's MgO EOS. The correction replaces
them with the original boundaries for comparison. Paragraph 1 also corrects
a plotted theoretical boundary. Thus the original 2009 Figure 11 is evidence
of Tange's recalculation, not an archived output of Hirose's reducer. The
correction does not identify either implementation's numerical code or Fei's
Au objective.

[Ye et al. (2017)](https://doi.org/10.1002/2016JB013811), section 3.3,
page 3457, Eq. 3 prints the exponential variable-q law:

\[
\gamma(V)=\gamma_0\exp\left[\frac{q_0}{q_1}
\left((V/V_0)^{q_1}-1\right)\right].
\]

The equation was verified against the PDF page image, not only extracted text.
It agrees with Speziale Eq. 11 and the author-confirmed prescription. It does
not support identifying the alternative local-power law as Speziale's printed
law. The paper's later scale comparisons and refinements do not recover the
2004-2007 reducer, Debye-temperature integration or original Au fitting weights.
Its supplement inventory consists of supporting information and four data
sets for its own experiments; no historical Fei optimizer is identified there.
The recovered [publisher-version PDF](https://www.osti.gov/servlets/purl/1368221)
is a reference source, never an Au fit input.

## Diagnostic inference retained without promotion

The existing [Speziale equation audit](speziale-2001-mgo.md) finds that
substituting q(V) into gamma0*(V/V0)^q(V), a different gamma law, closely
matches later reported MgO pressures. With integrated theta and Tr = 298 K,
its Fei 26-row RMS is 0.003264 GPa, versus 0.254647 GPa for printed Eq. 11
with integrated theta and Tr = 300 K. These are inherited diagnostics, not
new fits or newly recovered source evidence. Alternative theta conventions
also agree at the printed precision, so the tabulated pressures do not uniquely
identify the entire historical algorithm.

The direct Speziale reply and the later printed exponential equation constrain
what may be called the source equation. The pressure matches constrain possible
later implementations. Neither establishes who used which implementation,
whether any code contained an error, or whether Fei refitted against reported
pressures rather than recalculated unrounded MgO states. Reported targets,
published Au coefficients and conditional parity therefore remain intact.

## Supplement and calculator recovery limits

The original
[Hirose (2006) publisher page](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005GL024468)
lists `grl20669-sup-0001-t01.txt` as a 1012-byte tab-delimited Table 1.
Both its download action and direct supplement route returned HTTP 403 on
12 October 2026. The existing displayed-table transcription remains usable;
the author-file bytes, unrounded values and any additional columns are
unverified. A table supplement, even if recovered, would not by itself supply
the reducer or regression objective.

Gao's calculator is a specific archival lead. Gao's
[2010 thesis](https://www.ideals.illinois.edu/items/18462), Appendix B,
documents a Java Pressure Scale program with modified MGD support, EOS-file
inputs and uncertainty outputs. Only indexed appendix documentation and the
repository record were accessible in this search; the full thesis download
returned HTTP 403. The original `netfiles.uiuc.edu/liligao2/` host no longer
resolved. A bounded Wayback prefix query recovered only a 2006 personal
homepage, with no calculator binary, EOS file or source code. This does not
prove that other archives lack it. Even a recovered Gao version would require
evidence linking that version to the Fei/Hirose reductions before it could
establish historical parity.

The Fei lab publication page and searched publisher/source routes supplied no
Au regression file or pressure-reduction code. Current PMC retrieval returned
a challenge page and Europe PMC full-text XML returned HTTP 500; the existing
cached original-paper text was used for the Au passages. These access results
are bounded recovery outcomes, not proof that no supplement or private file
exists. The IOP conference-paper route and Springer PDF route returned HTML
instead of PDF and were not archived as publications. No additional original
Au numerical observations, covariance or weights were recovered.

## Precise author questions

The highest-value request is a copy of the original Au fitting input and
calculator/workbook, with its version and a saved output. If those are
unavailable, the following answers would narrow independent gaps:

1. **Selection and stages:** Which exact cold and hot rows entered the 2007
   Au regressions? Were the cold inputs the 37 Dewaele rows only, or did they
   include the new 2007 points and the Hirose cold datum? Did all 26 Fei Table 1
   rows and both 1340/2330 K Hirose rows enter the thermal fit? Were any shock
   states included, and if so which original states and reductions?
2. **Objective and weights:** Was the minimized residual in pressure, volume,
   thermal pressure, or another coordinate? Were weights equal per row,
   inverse variance, per run/isotherm, or adjusted manually? Were weights
   parameter-dependent, and how were Au volume and shared-temperature errors
   treated? A formula or original per-row weight column would resolve this.
3. **Fixed and fitted anchors:** Were V0 = 67.850, K0 = 167, gamma0 = 2.97 and
   theta0 = 170 fixed exactly as printed? Was the hot stage anchored to an
   unrounded cold K0-prime, or to 6.00? Please provide unrounded values and
   the temperature assigned to Dewaele's original 298 K observations.
4. **MgO reducer identity:** Were the reported Table 1 MgO pressures used
   directly, or were they recalculated for the 2007 fit? Which calculator,
   spreadsheet or executable and version produced them? What gamma(V),
   theta(V), reference-temperature subtraction, gas/heat-capacity constant,
   molar normalization and MgO anchors did it use? Was Fei's reducer shared
   with Hirose's? In particular, how should the MgO mass-specific heat-capacity
   entry in Fei (2004) Table 3 be interpreted?
5. **Reducer checkpoints:** At the printed Hirose inputs
   (aMgO = 3.7827 angstrom, T = 1340 K) and
   (aMgO = 3.7767 angstrom, T = 2330 K), what pressures did that historical
   program return before rounding? Can the corresponding unrounded lattice
   parameters and temperatures be supplied? These would distinguish code
   output from reproduction of rounded table values.
6. **Au thermal convention:** Did the Au fitting code use the printed
   theta = theta0*(V/V0)^(-gamma(V)) law or an integrated constant-q theta
   law? Did it use modern molar 3R or the upstream mass-specific Au prefactor
   0.125 J/(g K)? Please supply the Debye energy definition and reference
   subtraction actually used, rather than only an EOS name.
7. **Published error widths:** What confidence/coverage convention produced
   K0-prime = 6.00(2) and q = 0.6(3)? Were these curvature errors,
   residual-scaled errors, a covariance propagation, run scatter or an
   estimated range? Did they include uncertainty in the cold anchor, gamma0,
   thermometry and MgO calibration? Is the covariance or normal matrix retained?

## Scope and verification

The numerical values above were checked against the current saved reports and
source transcriptions; the repeated Hirose states were cross-checked against
Hirose (2008) Table 1. Source links, the recovered PDF signature and SHA-256,
and local report links were checked. No fitting, material registration,
calibration/coverage edit, shared ledger/manifest regeneration, Pt work,
commit, push or author contact was performed.

The two existing tests for source-table columns/missing errors and reported
MgO hot-target selection passed (2 passed). The saved numerical report's
SHA-256 is unchanged; no diagnostic fits were rerun. Local links, archived
source hashes, repeated Hirose table values and new-file whitespace passed
the follow-up checks; `git diff --check` also passed.

The follow-up adds this report and a Fei-specific source directory with the
recovered Ye PDF, bounded archive-query response and
[source manifest](../data/fei-2007-gold-follow-up-sources/source-manifest.json).
Third-party source rights remain with the respective rightsholders; no CC0
claim is made for their articles, source tables or code. There are no shared
generated files from this follow-up requiring reconciliation. Any concurrent
Au changes in the shared checkout belong to their owning chats and are
preserved.
