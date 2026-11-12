# Grossular: observation recovery and fit audit

Audit: 12 November 2026. This extends the [garnet source-equation audit](garnet-angel-2022.md)
for `grossular_angel_2022_bm3_q_compromise_debye`. The source is Angel et al.
(2022), [*Garnet EoS: a critical review and synthesis*](https://doi.org/10.1007/s00410-022-01918-5).
The published coefficients remain unchanged.

## Outcome and comparison scope

The latest [remaining-channel investigation](grossular-angel-2022-remaining-channels.md)
recovers Bosenick temperature errors and plotted normalization, and bounds
Conrad modulus errors. Its 193-row source-processing diagnostic gives
K₀ = 166.754 GPa, K′ = 5.33063, θ₀ = 780.407 K and γ₀ = 1.15522.
Nineteen Isaak modulus errors and seven Thiéblot experimental errors remain
unrecovered; one Bosenick row is an anchor. It does not reproduce the full fit.

The [subsequent 147-channel fit](grossular-angel-2022-expanded-fit.md) adds all
52 recovered Gwanmesia moduli. It gives K₀ = 167.353 GPa, K′ = 4.98267,
θ₀ = 740.183 K and γ₀ = 1.16467. K₀, K′ and γ move slightly closer to the
source than the historical 95-channel result below; θ moves farther away.
This remains a conditional diagnostic, not the full author fit.

**The published equation passes numerical checks; the full original regression
cannot yet be reproduced.** The [subsequent recovery](grossular-angel-2022-recovery.md)
adds all 52 Gwanmesia table rows and seven digitized Thiéblot plotted points,
bringing coverage to **213 tabulated/table-derived plus seven digitized
channels** under the nominal 220-channel S1 inventory. It also traces the
Isaak plotted expansion series to GR-2B, rather than the GGD specimen named
in S1. Original Thiéblot values/errors, exact normalization, some assigned
errors and original fit-project selections remain unavailable. The earlier
161-channel receipts below are preserved as historical results.

A separately declared **95-channel partial refit** converges from two starts,
giving K₀ = 167.373 GPa, K′ = 4.96395, θ₀ = 744.384 K and γ₀ = 1.16605.
These are diagnostic coefficients, not a replacement EOS or a reproduction of
the author's complete fit. The source record remains a source-parameterization
assessment with its existing qualifications.

The pyrope request remains **waiting for Ross Angel's reply**. That email asked
about pyrope; this audit has not sent a grossular request.

## Published equation and retained parameters

The represented model is BM3 plus the q-compromise Debye thermal pressure:

\[
P(V,T)=P_{\mathrm{BM3}}(V)+
\frac{\gamma_0}{1000V_{m0}}
\left[E_D(\theta_0,T)-E_D(\theta_0,298)\right].
\]

Here V is molar volume in cm³/mol, Debye energy is J/mol, and pressure is GPa.
Both θ and γ/V are constant; **there is no fitted q**. The fixed reference
volume is 125.35 cm³/mol, adopted from Holland–Powell (2011), with 20 atoms per
formula unit and Z = 8 for the conventional cell. The adopted EOS volume is
distinct from a measured diffraction cell. The [earlier equation audit](garnet-angel-2022.md#equation-units-and-capability-boundary)
explains the capability and thermodynamic-potential boundary.

| Parameter | Official EOS file | Printed Table 1 width |
| --- | ---: | ---: |
| K₀T, GPa | 167.002 | 0.20 |
| K′₀T | 5.06534 | 0.08 |
| θ₀, K | 749.825 | 13 |
| γ₀ | 1.15602 | 0.006 |

The complete original 59-parameter file and 59 × 59 covariance remain intact.
Only indices **2, 3, 11, 18, 53, 54** have nonzero covariance diagonals: the
four EOS coefficients and the two modulus scale factors. The unlabeled value
at parameter 56 has zero covariance and is not evidence for another fitted
scale. The source-file γ₀ width is 0.00509706, retained separately from the
printed 0.006. Neither confidence convention is invented.

EosFit multiplies **calculated** moduli by a scale to compare them with the
observed values. The source scales are **0.99455 for Kono** and **1.00551 for
Gwanmesia**, not their reciprocals. The latter's source covariance gives a
width of 0.00249546; S1 prints 1.006(2), while the main text prints 1.006(3).
Those statements remain distinct.

## Recovered observations and selection inventory

The [recovery manifest](../data/garnet-angel-2022/grossular/input-recovery.json)
records the earlier source locations, access URLs, original-evidence hashes where downloaded,
transcription hashes, missing errors and qualifications. The
[access log](../data/garnet-angel-2022/grossular/source-access-checks.json)
distinguishes inaccessible files, metadata-only responses and unavailable tables.
An unsuccessful download is not a finding that the data do not exist.

| Source and channel | Nominal S1 used | Recovered | Evidence and qualification |
| --- | ---: | ---: | --- |
| Isaak 1992 expansion | 15 | 15 | [Table 2 transcription](../data/garnet-angel-2022/grossular/isaak-1992-table-2.csv); observed GR-2B strain, specimen discrepancy below |
| Bosenick–Geiger 1997 expansion | 20 | 20 | [Table 2](../data/garnet-angel-2022/grossular/bosenick-geiger-1997-table-2.csv); indexed primary author-uploaded text, including repeats; PDF not pixel-checked |
| Thiéblot 1998 expansion | 7 | 7 digitized | [Follow-up extraction](grossular-angel-2022-recovery.md#thieblot-seven-plotted-points-original-table-still-missing); original table/errors unavailable |
| Rodehorst 2002 expansion | 9 | 9 | [Table 4A](../data/garnet-angel-2022/grossular/rodehorst-2002-table-4a.csv); indexed primary text; 298 K anchor retained separately; PDF not pixel-checked |
| Du 2015 expansion | 13 | 13 | [Preprint Table 2](../data/garnet-angel-2022/grossular/du-2015-preprint-table-2.csv); scan visually checked, final publisher table not directly recovered |
| Milani 2017 expansion | 35 | 35 | [Official deposit Table 2](../data/garnet-angel-2022/grossular/milani-2017-table-2.csv); XLSX/PDF agree |
| Milani compression, labeled 2015 in S1 | 7 | 7 | [Official 2017 Table 3](../data/garnet-angel-2022/grossular/milani-2017-table-3.csv); eighth row is ambient anchor |
| Bass 1989 ambient Kₛ | 1 | 1 | [Author abstract values](../data/garnet-angel-2022/grossular/bass-1989-ambient.csv); primary full PDF not recovered |
| Isaak 1992 Kₛ | 22 | 22 | [Table 3b](../data/garnet-angel-2022/grossular/isaak-1992-table-3b.csv); GGD specimen; only three explicit errors |
| Conrad 1999 Kₛ | 3 | 3 | [Table 6 Cᵢⱼ](../data/garnet-angel-2022/grossular/conrad-1999-table-6.csv); Kₛ derived, original two excluded high-P rows also retained |
| Kono 2010 Kₛ | 36 | 36 | [Table 1](../data/garnet-angel-2022/grossular/kono-2010-table-1.csv); all eight original columns retained |
| Gwanmesia 2014 Kₛ | 52 | 52 | [Final publisher Table 1](../data/garnet-angel-2022/grossular/gwanmesia-2014-table-1.csv); full follow-up source receipt |
| **Nominal total** | **220** | **213 tabulated/derived + 7 digitized** | **Exact author input values, errors and residual count remain uncertain** |

The original numerical tables were obtained from the
[MSA 2017 deposit](http://www.minsocam.org/MSA/AmMin/TOC/2017/Apr2017_data/AM-17-45855.zip),
[Isaak primary paper](https://doi.org/10.1007/BF00198608),
[Kono primary paper](https://doi.org/10.1007/s12583-010-0112-2),
[Conrad primary paper](https://doi.org/10.2138/am-1999-0321),
[Du author preprint](https://academiccommons.columbia.edu/doi/10.7916/D87P8XSM),
[Bosenick–Geiger author upload](https://www.researchgate.net/publication/235764549_Powder_X_ray_diffraction_study_of_synthetic_pyrope-grossular_garnets_between_20_and_295_K),
[Rodehorst primary table](https://rruff.info/doclib/am/vol87/AM87_542.pdf), and
[Bass author abstract](https://experts.illinois.edu/en/publications/elasticity-of-grossular-and-spessartite-garnets-by-brillouin-spec/).
Publisher PDFs stay in private audit output; attributed factual observations
are retained here. Accessibility is not treated as a source license.

### Inventory and reduction problems

**Isaak expansion uses a different specimen.** The
[primary paper](https://doi.org/10.1007/BF00198608), pp.107–109 and Table 2,
identifies the observed dilatometry specimen as **GR-2B: 76.1 mol% grossular,
22.1 mol% andradite**. Table S1 instead labels the 15 expansion channels with
the **GGD** composition, about 96.53% grossular, which belongs to the elasticity
specimen. Isaak's last Table 2 column estimates end-member **expansivity** using
an additional compositional assumption; it is not measured end-member volume.
We preserve observed strain, calculated strain and estimated end-member α in
separate columns. The observed strain values were interpolated by about 2 K
onto the published temperature grid. Values extrapolated above 1000 K are
not added as observations. The [follow-up vector comparison](grossular-angel-2022-recovery.md#isaak-specimen-identity-traced-to-the-plotted-series)
identifies the plotted series as consistent with directly scaled GR-2B strain.
Isaak p.109 explicitly confirms a 293 K strain reference. The exact fit-project
reduction remains unconfirmed; these 15 channels stay outside the diagnostic
end-member refit.

**The compression citation/range is inconsistent.** S1 says “Milani 2015,
0–8.4 GPa”; the article and EOS comments cite Milani 2017. The official 2017
Table 3 contains eight rows, including ambient, ending at **7.674 GPa**.
The conditional replay uses that grossular table and retains the discrepancy.
It does not borrow the Milani 2015 pyrope observations. Original paired quartz
volumes are missing: the recalibration infers equivalent quartz volumes from
rounded published pressures, then evaluates the Scheidl (2016) BM4 scale.

**220 is an inventory sum, not a recovered residual count.** S1's footnote
states that normalization removes an ambient datum from a diffraction set.
Yet some original/used counts are identical (35/35, 20/20, 13/13). For
compression and Rodehorst the counts explicitly decrease by one. Anchors,
sub-instrument groupings and their induced covariance cannot be reconstructed
from these counts alone. Coverage identifies upstream or plotted channels
under the S1 inventory; it does not prove they entered the author fit unchanged.

**Selection and error details remain incomplete.** Bosenick–Geiger's own
discussion excludes low-temperature cooling measurements and a 20 K point,
whereas S1 lists all 20 channels as used. Rodehorst appears in S1 but not in
the final EOS comments. Isaak Table 3b prints Kₛ errors only at 300, 800 and
1300 K. Other errors remain blank, rather than being interpolated as factual
source inputs. Conrad supplies Cᵢⱼ errors but no C₁₁/C₁₂ covariance, so we
compute Kₛ = (C₁₁ + 2C₁₂)/3 without inventing an independent Kₛ error.

S1 excludes the 36 Kono and 52 Gwanmesia **volume** constraints while retaining
their modulus constraints. Those volume columns are not added to the refit.
Dachs heat capacities helped select a consistent preliminary Cp-EOS; S1 does
not identify them as additional Debye residual channels in the final joint
V/Kₛ fit. No Cp observation count is invented.

## Fixed-source residual comparisons

![Recovered grossular residuals against unchanged published coefficients](../data/garnet-angel-2022/grossular/source-residuals.svg#only-light)
![Recovered grossular residuals against unchanged published coefficients](../data/garnet-angel-2022/grossular/source-residuals-dark.svg#only-dark)

The expansion panels profile one constant volume scale per dataset with
reported volume errors only. They test curve shape without conflating
instrument offsets with expansion. They do not replay the original
normalization, temperature weighting or shared anchor covariance. Asterisks
mark indexed-only transcriptions. The modulus panel uses the published Kono
scale. The compression panel fixes the measured ambient normalization anchor,
recalibrates pressures as described above and explicitly adopts 298 K.

| Comparison | Rows | RMS observed minus source |
| --- | ---: | ---: |
| Milani expansion, profiled scale | 35 | 0.02473 cm³/mol |
| Du expansion, profiled scale | 13 | 0.04294 cm³/mol |
| Bosenick–Geiger expansion, all tabulated rows | 20 | 0.02531 cm³/mol |
| Rodehorst expansion, anchor omitted | 9 | 0.05641 cm³/mol |
| Milani compression, recalibrated/normalized | 7 | 0.01318 cm³/mol |
| Kono moduli, calculated source scale 0.99455 | 36 | 1.59581 GPa |
| Isaak GGD moduli | 22 | 0.60450 GPa |
| Bass ambient modulus | 1 | 0.21966 GPa |
| Conrad selected derived moduli | 3 | 2.07113 GPa |

Du has a volume-error-only standardized RMS of **4.92**, versus **1.40** for
Milani. These conditional residuals document scatter/shape tension under the
specified comparison; they do not identify an error in the published EOS or
establish its full-fit χ². Source-table version, weighting and normalization
remain relevant limitations.

Independent quadrature and official **EosFit7c 7.60** agree at nine identical-
coefficient states for V, Kₜ, Kₛ, α, γ and θ, including 20 K and 16.91 GPa /
1650 K. The [receipt and raw console output](../data/garnet-angel-2022/grossular/eosfit-properties/checks.json)
retain the previously declared printed-precision tolerances. This extends the
earlier pressure-only checks; it validates equations and conventions, not
the original regression or experimental accuracy.

## Declared 95-channel partial refit

The [plan](../data/garnet-angel-2022/grossular/refit-plan.json) was written
before solving. It includes 35 Milani expansion, seven Milani compression,
13 Du expansion, 36 Kono modulus, three Isaak modulus values with explicit
errors and one Bass ambient modulus. No missing observations or errors are
fabricated. Indexed-only expansion tables remain fixed-source comparisons.

The objective is an iteratively reweighted pressure residual, propagating the
available P/T/V or Kₛ errors into an effective pressure variance. Missing error
terms are omitted, not known zero. Two freely fitted instrumental volume
scales and a Kono modulus scale are nuisance parameters. This differs from
the author's fixed pre-normalization and incomplete source selection.
Shared anchor/calibration covariance is omitted. No weights or selections
were tuned to obtain publication agreement.

| Quantity | Published source | Conditional partial refit |
| --- | ---: | ---: |
| K₀T, GPa | 167.002 | 167.372989 |
| K′₀T | 5.06534 | 4.963951 |
| θ₀, K | 749.825 | 744.384257 |
| γ₀ | 1.15602 | 1.166055 |
| Kono calculated-modulus scale | 0.99455 | 0.996923 |

The fitted Milani and Du volume multipliers are 0.9993928 and 0.9989457.
Both numerical starts converge to the same result, without hitting bounds.
The conditional reduced χ² is **4.36512**, with **88 degrees of freedom**.
It belongs to this 95-channel, seven-parameter objective and is not the
unrecovered full source fit's statistic. Local diagonal-error covariance is
retained in the [fit receipt](../data/garnet-angel-2022/grossular/subset-refits.json);
its widths are conditional, not reconstructed published uncertainty intervals.

This subset supports a thermal temperature near the published 750 K and a
similar pressure curve, with modest coefficient shifts. It does **not**
establish full coefficient or author-workflow reproduction. The audit stops
here without searching alternative EOS laws or unknown weights.

## Integration, verification and resumption

The audit retains attributed source observations and diagnostics separately
from material records. The grossular record's copied pyrope/Hartwig wording
is corrected to name the grossular source gaps. No EOS coefficient, source
covariance, structure, default or scientific classification changes. No
diagnostic refit is registered as a catalog EOS.

The checks cover source hashes/tokens, missing-error preservation, inventory
counts, the official property receipt, independent pressure/quadrature and
modulus inversion, the declared subset selection, convergence and public
material loading. See the [verification receipt](../data/garnet-angel-2022/grossular/verification.json)
for commands and actual results.

These are **repository-only maintenance tools**, run from a checkout with its
development dependencies; they are not installed-package API examples:

```bash
python -m scripts.audit_grossular_angel_2022_observations
python -m scripts.refit_grossular_angel_2022_subset
python -m scripts.plot_grossular_angel_2022_audit
pytest -q tests/test_grossular_observation_audit.py
# Requires a separately installed working EosFit executable/runtime:
python -m scripts.run_grossular_angel_2022_eosfit --executable /path/to/eosfit7c
```

The [follow-up recovery](grossular-angel-2022-recovery.md) supplies the 52
Gwanmesia rows and seven digitized Thiéblot points. To resume a full reproduction,
recover the original Thiéblot table/errors; confirm the Isaak fit input; recover
exact normalization anchors, point masks, coordinate errors, scale groups and
original fit controls. Paired quartz readings would remove the rounded-pressure
reconstruction assumption. The current pyrope email is not an already-sent
request for these grossular files.
