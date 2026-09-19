# Sakai et al. (2011): experimental NaCl-B2 on four Pt scales

## Disposition and primary authority

Sakai, Ohtani, Hirao and Ohishi, *Journal of Applied Physics* **109**, 084912,
[doi:10.1063/1.3573393](https://doi.org/10.1063/1.3573393), publishes eight
experimental parameterizations in Table II. All eight are accepted: BM3 and
Vinet separately calibrated against Matsui (2009), Fei (2007), Dorogokupets–Oganov
(2007), and Holmes (1989) Pt. These are fits to measured NaCl compression even
where a Pt calibrant incorporates theoretical information. No computational
NaCl EOS or borrowed comparison curve is added. The existing Sakai 2025,
Heinz 1984, Shen 2026 and Sakai 2014/Yokoo-Pt records remain unchanged.

The publication-of-record PDF was recovered from Zotero item `WZC4B8LJ`,
attachment `CEQD3YIJ`, and checked visually against Tables I–II. The publisher
HTML snapshot and live article were also inspected. No visible supplementary
file was found (the page contains a hidden generic Supplementary Material
control). Crossref has no correction/update relation and a focused correction
search found none. This is a retrieval result as of 2026-09-19, not a guarantee
that no later correction can exist. The upstream [Sata et al. (2002) publication
of record](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.65.104114/fulltext)
was recovered directly from APS. Source hashes, URLs and retrieval details are
in `peritheos/data/datasets/nacl-b2-sakai-2011-source.json`.

## Identity, equations and uncertainties

The existing `nacl_b2` material is reused: pure NaCl, CsCl-type B2, cubic
Pm-3m (221), Z=1, occupied Na 1a and Cl 1b sites. The existing structural
reference cell is retained; the EOS V0 values are separate extrapolated
references, not measured ambient B2 lattice volumes. Table I's NaCl volumes
are cubic conventional-cell volumes in Å³; Pt is a four-atom conventional
fcc cell. No factor of four is applied to the NaCl-B2 volume.

Section III.C prints the standard third-order Birch–Murnaghan and Vinet
expressions for P300K; the existing exact implementations are reused.
With f=((V0/V)^(2/3)-1)/2 and x=(V/V0)^(1/3), these are
P=3K0 f(1+2f)^(5/2)[1+3(K0′-4)f/2] and
P=3K0(1-x)x^(-2) exp[3(K0′-1)(1-x)/2].
The reference is P=0 at 300 K. No thermal EOS is inferred from laser annealing.

| Pt scale | V0 (Å³), fixed after extrapolation | BM3 K0 (GPa), K0′ | Vinet K0 (GPa), K0′ |
|---|---|---|---|
| Matsui | 37.73 ± 4.05 | 47.00 ± 0.46, 4.10 ± 0.02 | 40.40 ± 0.54, 5.04 ± 0.04 |
| Fei | 37.73 ± 3.88 | 47.80 ± 0.45, 4.06 ± 0.02 | 41.36 ± 0.53, 4.96 ± 0.04 |
| Dorogokupets–Oganov | 37.76 ± 4.21 | 47.24 ± 0.47, 4.14 ± 0.02 | 40.37 ± 0.55, 5.11 ± 0.04 |
| Holmes | 38.01 ± 5.26 | 43.42 ± 0.50, 4.36 ± 0.02 | 36.07 ± 0.54, 5.49 ± 0.04 |

V0 errors come from Section III.B, not Table II's `fix` entries. The large
printed errors are retained literally, including the explicitly printed
37.73(4.05). They combine the g–G fit error and B1 reference-volume error.
K0 and K0′ are fitted with V0 fixed; these error sources are not silently
combined. The confidence level, regression weights and covariance are not
reported; missing metadata is null, not zero.

Section III.B uses V01=179.42(16) Å³ for the B1 conventional cell as the
numerical strain reference. Figure 4's g coordinates confirm that reference
choice. Its printed G denominator has exponent −5/2, which is inconsistent
with its plotted normalized stresses and the stated quadratic/BM3 equivalence.
The diagnostic uses G=P/[3(1+2g)^(+5/2)], equivalently
G=(P/3)(1+2g)^(-5/2), and explicitly records this interpretation. It does not
change the unambiguous Section III.C production equations. Using a consistently
rescaled B1 reference shifts the g coordinate but does not change the physical
zero-pressure volume obtained from an exact quadratic transformation.

## Observations, calibration and domain

Two checksummed CSVs retain every relevant primary row and its reported errors:

- `nacl-b2-sakai-2011-table1.csv`: all 27 rows, both lattice parameters and
  volumes with errors, and all four original pressure columns.
- `nacl-b2-sata-2002-table1-sakai-input.csv`: all 29 rows with Pt and MgO
  diffraction readings, original marker pressures and errors, NaCl spacings
  and volumes with errors, annealing flags, and Sakai's explicit exclusions
  of rows 8 and 11. No discarded source observation is erased.

A source-printed lattice error of 0.000 is preserved; it is not interpreted
as infinite statistical weight. The confidence convention is unspecified.
Sata marks 16 rows as measured after annealing; its own EOS uses those rows,
whereas Sakai explicitly names only exclusions 8 and 11. The default diagnostic
uses all remaining 27 Sata rows, and tests annealed-only selection separately.

The [Matsui primary paper](https://doi.org/10.1063/1.3054331), Equation (6) and
Table II, gives the 300 K Vinet Pt scale V0=60.38 Å³, K0=273 GPa, K0′=5.20.
The [Fei primary Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC1890468/)
gives Pt Vinet V0=60.38 Å³, K0=277 GPa, K0′=5.08. The
[Dorogokupets–Oganov primary Table I](https://uspex-team.org/static/file/PressureScales-PRB-2007.pdf)
gives V0=9.091 cm³/mol, K0=276.07 GPa and K0′=5.30 with its full 298.15 K
reference model; the existing exact record is evaluated at 300 K. The new
NaCl records link that exact bundled calibration and the separately audited
300 K Matsui and Fei reference-isotherm records described below. The original
refit reproduction keeps independent Matsui/Fei expressions as equation checks.

Re-reducing the rounded Sakai Pt volumes differs from its Matsui, Fei and
D&O pressure columns by at most 0.1582, 0.1454 and 0.1638 GPa respectively,
consistent with 0.01 Å³ and 0.1 GPa table rounding. The exact bundled Holmes
Equation (11) implementation differs by up to 0.701 GPa. A commonly rounded
Vinet candidate 60.38/266/5.81 reduces that maximum to 0.159 GPa, but the
paper does not specify that implementation. This is retained as a calibration
ambiguity; no exact Holmes catalog link is asserted. The diagnostic uses the
primary Holmes Equation (11) model and original Sakai pressure column, not
an adjusted calibration selected to improve agreement.

The 27 new observations span 53.7–304.4 GPa on Matsui, 54.0–301.8 on Fei,
54.7–313.8 on D&O and 54.5–333.1 on Holmes. These scale-specific bounds are
the conservative executable validity envelopes. Lower-pressure Sata points
also enter the source fit; their reconstructed envelope is reported separately
in the reproduction JSON. No claim of measurement at 0 or 364 GPa is made.
The paper's 50–364 GPa curve comparison and Earth-center discussion extend
beyond its new measurements. B2 is unquenchable and V0 is extrapolated.
The sample was annealed around 1500–2000 K; diffraction/EOS values here are
at 300 K after annealing.

The stress assessment is for Pt, using its gamma plot and published elastic
compliances. Assuming alpha=1 gives a lower-bound |t|=2.4 GPa at 304 GPa.
Only Pt (111) and (200) were available there, rather than all three peaks.
The reported lattice deviations are within ±0.16% over the series, ±0.08%
at 304 GPa, with reflection-derived pressures varying within 2.5%. Comparison
with helium-medium studies is evidence of relatively low deviatoric stress,
not proof of a hydrostatic sample or a universal 2.4 GPa pressure uncertainty.

## Deterministic reproduction and independent refits

Run `uv run --python 3.9 python scripts/reproduce_sakai_2011_nacl_b2.py --check`.
The checked output uses Python 3.9.6, NumPy 2.0.2 and SciPy 1.13.1;
optimizer last digits can vary with the numerical environment.
The result is `docs/data/sakai-2011-nacl-b2-reproduction.json`.

For each scale the independent refit retains the 27 original Sakai pressures
and reconstructs the 27 selected Sata pressures from the equal-weight mean
of sqrt(3)d111 and 2d200 for Pt, cubed. That averaging rule is explicit and
is not claimed to be the authors' unpublished rule. Stage one fits an
unweighted quadratic in g–G and selects the physically relevant zero-pressure
root. Stage two fixes V0 and fits K0/K0′ in pressure space. A separate
conditional fit fixes the published V0, so upstream extrapolation differences
can be distinguished from final-regression differences. Reported uncertainties
are retained but not used as weights in the default diagnostic because the
source objective is unspecified. Volume-error weighting and annealed-only
selection are alternative diagnostics, not opt-in catalog records.

| Scale | Conditional BM3 K0, K0′ | Conditional Vinet K0, K0′ | Staged V0 (Å³) |
|---|---|---|---|
| Matsui | 47.14270, 4.09757 | 40.56908, 5.03092 | 37.85144 |
| Fei | 47.95634, 4.05604 | 41.54273, 4.95128 | 37.84272 |
| D&O | 47.41983, 4.13437 | 40.57182, 5.10033 | 37.90054 |
| Holmes | 43.61800, 4.35058 | 36.27065, 5.47382 | 38.22321 |

All conditional K0/K0′ estimates lie within their respective printed error
widths; all staged V0 estimates lie within the large printed errors. This is
conditional numerical parity, not exact recovery of the undisclosed source
weights, uncertainty calculation or marker averaging. The 54-point unweighted
pressure RMS residuals range from 1.67 to 2.14 GPa. Weighting by the propagated
NaCl volume errors moves fitted K0 appreciably, so no diagnostic refit is
promoted to a new scientific recommendation.

The production curves agree with the independently coded Section III.C
expressions to better than 1e-10 GPa. All 27 independent Table I observations
for each curve fall inside a tolerance constructed from three printed volume
error widths (without assigning a statistical confidence to those widths).
Both the low-pressure and highest-pressure observations are checked.
At the Matsui BM3 volume for 364 GPa, Vinet gives 358.99456 GPa, reproducing
the separately reported approximately 5 GPa difference. Other BM3 scale
shifts are −3.47076 (Fei), +12.35965 (D&O) and +39.13535 GPa (Holmes), versus
reported −4, +12 and +37. The Holmes central value is visibly different;
the source result is compatible with independently propagated half-last-digit
rounding of both sets of Table II coefficients (36.179–42.122 GPa). These
rounding bounds are numerical compatibility intervals, not confidence bounds.

All eight candidates pass the acceptance gate through complete equation and
parameter mapping, primary observation coverage and independent numerical
benchmarks. No source-owned candidate is withheld. Remaining reproducibility
limits are the unspecified weighting/averaging and the Holmes calibration
variant, all retained in the audit and refit ledger.

## Executable recalibration and comparison with published variants

The calibration links now resolve to `platinum_matsui_2009_vinet_300k` and
`platinum_fei_2007_vinet_300k`, in addition to the existing complete D&O model.
These two Pt records implement only the exact 300 K reference isotherms needed
here. Their cell convention is fcc Pt, Z=4, V0=60.38 Å³. They do not implement
either paper's thermal model. No additional NaCl fit is created by conversion.

Matsui's Equation (6) and Table II give K0=273 GPa and K0′=5.20; no coefficient
errors or covariance are supplied. Table III independently gives 235.96 GPa at
V/V0=0.70 and 300 K, recovered within its 0.005 GPa rounding half-width.
The source's shock constraints extend to 290 GPa, and its proposed standard
extends to 300 GPa. The record's 0–300 GPa validity is a reference-model scope,
not a claim of static 300 K measurements across that range. Seven existing
Holmes shock rows remain linked as upstream constraints. A static fit to those
hot Hugoniot observations would be invalid; reconstructing Matsui's complete
thermal reduction is outside this reference-isotherm addition. Its refit-ledger
classification is therefore `not_refittable`.

Fei's Table 1 gives V0=60.38(1) Å³, K0=277 GPa and K0′=5.08(2) for Vinet;
the source does not identify an uncertainty confidence level or covariance.
The source's explicit reference temperature is 300 K. Although its rounded
cold coefficients equal Dewaele (2004), that existing record retains its
298 K reference. The Fei record preserves the 300 K convention and the
Dewaele revised ruby calibration for the cold subset. Its conservative
0–93.6 GPa validity follows the existing 36 Dewaele observations.
Those rows and all 42 upstream Fei (2004) P-V-T rows remain bundled and linked.
The latter are context for the joint thermal fit, not a new 300 K dataset.
A conditional 36-row cold refit, with V0 and K0 fixed, gives K0′=5.08448 and
0.221 GPa pressure RMSE. Its classification is capped at `similar`: the complete
joint thermal optimization and Au re-reduction have not been reconstructed.
Both new isotherms also reproduce all 27 independently printed Sakai Pt
pressure reductions within 0.17 GPa, consistent with the rounded input table.

Run:

```bash
uv run python scripts/compare_sakai_2011_pressure_scales.py --check
```

The output is [`sakai-2011-pressure-scale-comparison.json`](../data/sakai-2011-pressure-scale-comparison.json).
It evaluates every directed pair among Matsui, Fei and D&O, for both NaCl
equation forms, with identity controls. Two additional Matsui-to-Holmes
comparisons are diagnostic only. All eight published NaCl variants are covered.

For each NaCl volume, the public API evaluates the source NaCl EOS, inverts
its linked Pt EOS to obtain a **virtual** Pt volume, then evaluates the target
Pt EOS there. It compares that pressure with the published target NaCl fit at
the same NaCl volume. The virtual Pt volume is not a recovered measurement.
No new regression or replacement coefficient set is involved.

Each comparison samples 1,001 equally spaced NaCl volumes across the
intersection of the two sample-fit envelopes. These are sampled maxima,
not analytic bounds. A second result retains only states inside both sample
and both Pt validity envelopes and verifies them with `check_validity=True`.
Thus sample-range results can explicitly include Pt-standard extrapolation.
The percentage denominator is the published target pressure at each state.

Using Matsui as a display baseline, without changing any material default:

| Published target variant | Maximum absolute difference (GPa), sample envelopes | Maximum relative difference | Maximum absolute difference (GPa), all linked envelopes |
|---|---:|---:|---:|
| Matsui BM3 / Vinet, identity controls | 0 | 0% | 0 |
| Fei BM3 | 0.2063 | 0.0683% | 0.0164 |
| Fei Vinet | 0.1298 | 0.0741% | 0.0502 |
| D&O BM3 | 0.3158 | 0.1006% | 0.2001 |
| D&O Vinet | 0.2609 | 0.1246% | 0.1486 |
| Holmes BM3, diagnostic only | 0.0472 | 0.0761% | 0.0472 |
| Holmes Vinet, diagnostic only | 0.4121 | 0.3417% | 0.3805 |

Across **all directed resolved pairs**, the largest sampled absolute difference
is 0.3158 GPa and the largest relative difference is 0.1475% (rounded upward).
The largest absolute difference with every linked envelope enforced is
0.2001 GPa. Fei's strict subsets extend only to approximately 93.5 GPa;
Matsui-to-D&O strict subsets extend to approximately 242.5 GPa, where the D&O
300 K volume-ratio bound becomes limiting. Full per-pair ranges, extrema,
paths, RMS differences and excluded-state counts are retained in the JSON.

These discrepancies support using recalibration for a workflow that tolerates
about 0.5 GPa of difference from the published alternatives over the sampled
resolved ranges. They do not support a uniform 0.1 GPa equivalence claim.
They are differences between rounded fitted curves, not an accuracy estimate
for either pressure standard, an experimental uncertainty, or a confidence
interval. Converting a fitted curve is not mathematically equivalent to
refitting the recalibrated measurements, and it need not preserve BM3/Vinet form.
The eight source records remain useful for reproducing cited published results;
normal recalibration does not require selecting another stored NaCl fit.

For example, this calculation stays inside all linked envelopes:

```python
from peritheos import recalculate_eos_pressure_scale

result = recalculate_eos_pressure_scale(
    source_eos_record="nacl_b2_sakai_2011_matsui_pt_bm3",
    sample_volume=[23.0, 22.5],  # NaCl-B2 conventional cell, Å³
    target_standard_eos_record="platinum_fei_2007_vinet_300k",
    sample_temperature_k=300.0,
    check_validity=True,
)
print(result.target_pressure_gpa)
```

Holmes remains explicitly unlinked on the Sakai records. Choosing the audited
Holmes Equation (11) as a target is an executable comparison with that model;
it is not a resolution of the paper's unspecified Holmes implementation.
