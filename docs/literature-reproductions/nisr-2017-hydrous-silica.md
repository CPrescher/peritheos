# Nisr et al. (2017): hydrous silica and the dry reference

Audit date: 2026-09-19. Primary authority: C. Nisr, K. Leinenweber,
V. Prakapenka, C. Prescher, S. Tkachev and S.-H. Dan Shim, *Journal of
Geophysical Research: Solid Earth* **122**, 6972–6983,
[10.1002/2017JB014055](https://doi.org/10.1002/2017JB014055).
Use the **26 October 2017 corrected publication of record**, not the accepted
manuscript. The terminal Erratum corrects the dry sample's `c` and `V0` in
Table 2; it does not revise the EOS tables or claim to update the supplement.

## Disposition of every candidate

| Source-owned experimental fit | Decision | V0 (Å³/cell) | K0 (GPa) | K0′ | Fit rows |
|---|---|---:|---:|---:|---|
| Dry stishovite, Table 3 | Add to existing `sio2_stv_andr` | 46.569 fixed | 312(2) | 4.59 fixed | S1, 20 |
| Hydrous stishovite, Table 3 | New composition/phase card | 47.191 fixed | 257(9) | 4.59 fixed | S2, 7 |
| Hydrous CaCl2 type, Table 4 | Separate new phase card | 47.23(36) fitted | 286(18) | 4 fixed | S4, 8 |

Identifiers are `stishovite_nisr_2017_dry_bm3`,
`hydrous_stishovite_nisr_2017_bm3`, and
`hydrous_silica_cacl2_nisr_2017_bm2`. The dry fit uses a separately synthesized
specimen and new measurements co-compressed with the hydrous sample on an Au
scale. It is not the already catalogued Andrault (2003) fit: only its fixed
pressure derivative is adopted from that paper. The existing default dry
record stays unchanged. All three new records preserve the published values;
the diagnostic refits below are not additional executable records.

The remaining five rows in Tables 3–4 cite Andrault (2003),
Bolfan-Casanova (2009), or Lakshtanov (2005). They are comparisons, not fits
owned by this paper. In particular, the dry CaCl2-type row `334(7), 4, 46.31`
belongs to Andrault (2003). No computational EOS, aggregate EOS, transition
EOS, or thermal EOS is added. None of the three source-owned fits is withheld.

## Source recovery and version control

Zotero item `FKB5RJG9` is outside Methods → EOS Library (`JT8V6LUL`). Its
attachment `XMRTZ8UU` is the corrected journal PDF; `7ILRMC37` is an accepted
manuscript and was not used as parameter authority. The publisher HTML was
also checked for its correction notice and Tables 1–4. The journal's official
[Supporting Information S1](https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2017JB014055&file=jgrb52249-sup-0001-Supplementary.pdf)
was downloaded through the browser after direct HTTP retrieval returned 403.
The publisher still serves a five-page file whose footer says DRAFT,
June 8, 2017. That version is recorded without silently repairing it.

The source manifest `peritheos/data/datasets/nisr-2017-source-manifest.json`
records PDF SHA-256 hashes, byte counts, version and retrieval provenance.
The primary numeric transcription is
`peritheos/data/datasets/nisr-2017-silica-tables-s1-s4.csv`.
It retains every row, decimal precision, lattice parameter, printed sigma and
available background-subtracted Rwp. Counts are S1=20, S2=7, S3=8, S4=8.
The eight S3 entries are four pressures with two separately refined phases,
not eight aggregate volumes. No run identifiers occur in these tables.
Tetragonal `b` and absent Rwp values stay empty rather than being fabricated.
Corrected Table 2 is separately bundled as an ambient structural dataset.
No separate data license is stated; the CSV is an attributed factual numeric
transcription, and publisher PDFs are checksummed rather than redistributed.

## Equation, reference state and uncertainty

Section 3.3 specifies Birch–Murnaghan fits, cites Birch (1978), and fixes the
pressure derivatives. With `eta=(V0/V)^(1/3)`, the exact mapping is

`P = 3 K0 / 2 * (eta^7 - eta^5) * [1 + 3/4 (K0′ - 4) (eta² - 1)]`.

The stishovite records use native BM3. Fixing `K0′=4` makes the high-pressure
record exactly BM2; the derivative is implicit model configuration, not a
fitted coefficient with zero error. There are no temperature-dependent terms.
The reference is 300 K and ambient pressure, which the source treats as zero
(the difference from 1 bar is 0.0001 GPa). Volumes are the conventional
primitive tetragonal or orthorhombic cell, with two nominal formula units.
They are not volumes per atom, formula unit or mole.

Tables 3–4 define parentheses as estimated one-sigma errors. Fixed-parameter
errors are stored as missing, with fixed status separately recorded. Measured
ambient errors in Table 2 are retained as observation errors. No covariance,
fit statistic, exact dependent variable, weighting rule, or fitting software
is reported. The CaCl2 phase is not quenchable; its fitted V0 is an
extrapolation and must not be substituted for a measured ambient cell.

Three independent ambient numbers must remain distinct:

- Corrected Table 2: dry `a=4.1805(9)`, `c=2.6647(9)`, `V=46.569(11)`;
  this volume agrees with the dry EOS reference in Table 3.
- Official S1: dry `a=4.18614`, `c=2.666`, `V=46.612` with
  `sigma(V)=0.0178`. This ambient observation
  was not overwritten by the correction. It has zero sensitivity to K0 in
  a volume-residual fit when V0 is fixed, but a small effect in a
  pressure-residual fit because its V differs from V0.
- Hydrous Table 2: `V=47.198(21)`, versus S2 and the fixed Table 3
  EOS `V0=47.191`. The S2 lattice product is about 47.198696 Å³.
  Preserve the tabulated volume as an observation and the Table 3 value
  as the EOS parameter; do not replace either by the product.

## Composition, crystallography and phase limits

The hydrous sample was synthesized at 9 GPa and 723 K for 49 h with excess
water; the dry sample at 10 GPa and 1473 K for 1 h. These are synthesis
conditions, not EOS temperature coverage. Both hydrous phase cards preserve
Table 1's nominal `Si0.954O2H0.184`. Charge balance is consistent with four
protons per silicon vacancy. Standard atomic masses give a water-equivalent
mass fraction **2.8102 wt%**, not exactly the reported central **3.2 wt%**;
the difference is within the quoted ±0.5 wt%. No formula is reverse-engineered
to force agreement.

The water estimate is indirect. Section 3.1 and Text S1 use a cell-volume
calibration updated in the companion Raman paper,
[10.2138/am-2017-5944](https://doi.org/10.2138/am-2017-5944), from Spektor's
thermal-analysis traces, volumes, phase fractions and private communication.
The original Spektor calibration would give 2.7 wt%. Text S1 does not publish
numerical calibration coefficients or the underlying communication. Consequently,
the 3.2±0.5 wt% calibration cannot be independently rebuilt from this paper's
files. The corrected Table 2 volumes do reproduce the stated approximately
1.3% expansion (1.350684%). Raman OH bands establish structural OH
qualitatively; they are not a quantitative assay. Retention of this water
content throughout compression is an interpretation rather than a measured
high-pressure composition. The cards state these limitations.

Structural provenance is separate from EOS parameter provenance: source
Sections 2 and 3.2, Figures 2a/3b, and S2/S4 supply the refinement symmetry,
indexed reflections and cell metrics. Andrault et al. (2003),
[10.2138/am-2003-2-307](https://doi.org/10.2138/am-2003-2-307), was also read
from Zotero attachment `JPM8IJSR` for the distinct rutile/CaCl2 structural
context; it does not establish the hydrous coordinates. The conventional
rutile topology has two cation positions and four oxygen positions, hence
Z=2 for the nominal formula, retained on distortion to CaCl2 type.
The introduction's `P42/mmn` spelling is treated as a transposition of
the canonical rutile group `P42/mnm` (136), consistent with the source's
rutile assignment and cited crystallography.

The authors explicitly did not refine atomic coordinates. Nor does their
XRD distinguish a Pnnm Si–O framework from a hydrogen-ordered
δ-AlOOH-related subgroup. The hydrous cards therefore use the protocol's
**documented peak-position fallback**, with no fictitious H sites or
fully occupied Si model:

- Rutile: P42/mnm (136); S2 ambient metric. Five hkl labels from Figure 2a
  (101, 200, 210, 112, 301) give d spacings calculated from that metric.
- CaCl2 type: Pnnm (58), explicitly a framework assignment; S4's
  47.57 GPa cell, matching Figure 3b. Its eleven explicitly indexed hkl
  labels give d spacings calculated from the measured metric.

Peak intensities are uniformly 1 as **display markers only**. They are not
claimed as measured or simulated intensities. The fallback supports peak
positions, not intensity analysis or hydrogen ordering. The measured structural
reference cell intentionally differs from the EOS's extrapolated V0, and
its reference pressure is explicit. Isotropic EOS rescaling is not a model
of the pressure-dependent axial ratios reported in Figure 4.

The hydrous low-pressure fit uses S2, 0–22.747 GPa (text: approximately
23 GPa); the high-pressure fit uses S4, 42.225–62.934 GPa. Neither is valid
in the 28–42 GPa coexistence interval. Broad peaks, no diffraction measurements
between about 22 and 28 GPa, and different media (He for XRD, Ar for Raman)
limit comparison with the companion Raman transition at 24–28 GPa. Do not
encode 28 GPa as an exact equilibrium transition or the finite interval as
a demonstrated thermodynamic two-phase boundary. Kinetic/first-order
interpretations remain the authors' hypotheses.

S1's dry volumes are refined as tetragonal throughout and support the
published dry reference fit. The paper discusses dry distortion near 60 GPa;
the full fitted envelope is therefore not a guarantee of stable rutile.
There is no independently tabulated dry two-phase aggregate to promote.

## Pressure calibration and unavailable inputs

Section 2 uses gold foil, individual Au XRD peaks, Fei et al. (2007), and
helium medium in the shared DAC chamber. The specimen volumes come from GSAS
Rietveld refinements of Dioptas-integrated patterns. The article does not
choose between Fei's published BM3 and Vinet Au parameterizations. The scale
is therefore partially resolved and no exact executable standard link is
invented. Neither the paper nor S1–S4 supplies paired gold volumes, gold
lattice parameters, pressure errors, or raw diffraction frames. Observation
level re-reduction remains unavailable. No pressure errors are inferred from
the stated pressure adjustment precision of approximately 0.5–1 GPa.

## Deterministic reproduction and independent refits

Run `uv run python scripts/reproduce_nisr_2017_silica.py`.
[Machine-readable results](../data/nisr-2017-silica-reproduction.json) include
three objectives, residuals, residual-scaled covariance, standard errors and
high-pressure observation benchmarks. The implementation evaluates the
Eulerian strain equation independently of Peritheos; SciPy least squares
starts from K0=300 GPa and, when free, V0=47 Å³. It uses all rows in the
specified source table, excludes S3 entirely, fixes the published derivatives
and stishovite V0, and fits the two CaCl2 parameters simultaneously. No
staged source fitting sequence is stated.

| Fit | Published | Independent unweighted P residual fit | Published P RMSE |
|---|---|---|---:|
| Dry stishovite | K0=312(2) | 312.135508 | 0.647669 GPa |
| Hydrous stishovite | K0=257(9) | 256.916727 | 0.853457 GPa |
| Hydrous CaCl2 | V0=47.23(36), K0=286(18) | 47.267767, 284.217610 | 0.190944 GPa |

Every free coefficient is within one reported source sigma. This establishes
coefficient parity under an explicit diagnostic objective, not knowledge of
the unpublished weights. The principal refit P RMSEs are respectively
0.647464, 0.853450 and 0.180529 GPa. High-pressure observations independently
check the curves: at V=43.9177 Å³ hydrous BM3 gives 21.796653 GPa versus
22.747 observed; at V=40.2441 Å³ hydrous BM2 gives 63.081538 versus
62.934 observed; dry V=40.3876 gives 61.701073 versus 62.934 observed.
Their volume residuals are below twice the corresponding source-table volume
RMS scatter (0.137035, 0.016530 and 0.059486 Å³ respectively).
These empirical tolerances are not the much smaller Rietveld volume sigmas;
the paper's printed EOS does not fit all individual points within those sigmas.

Unweighted volume residuals give dry K0=310.878937, hydrous rutile
K0=254.857999, and CaCl2 V0=47.257632/K0=284.779339. Using positive printed
volume sigmas gives 311.036019, 257.883636, and 47.170658/289.410543.
The S4 44.857 GPa row reports sigma(V)=0: it remains in the principal and
unweighted-volume fits and is excluded only from the sigma-weighted
sensitivity calculation. The S3 zero sigma is irrelevant to these selections.
No finite sigma is invented, and no zero is assigned infinite weight.
Diagnostic errors and covariance are computed from the Jacobian and residual
variance, not presented as the authors' uncertainties. No refit record is
needed because the source values are reproducible and the original weights
are unknown.

## Validation

The focused tests check all row groups, zero errors, checksums, corrected
reference choices, exact equation mapping, phase-gap rejection, independently
refitted coefficients, source observation residuals, native pressure/volume
round trips, and metadata-preserving material reconstruction. The global
primary EOS refit validator calls this independent reproduction and registers
three parity outcomes. Generic Python/Rust bundle tests cover schema,
serialization and native construction. Required contributor-check outcomes
are reported with the completed change.
