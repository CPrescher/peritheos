# Walker (2002): verified NaCl ancestry and remaining reproduction limits

Audit date: 2026-10-08. Both `kcl_walker_2002_bm3_2` (B2) and
`kcl_b1_walker_2002_bm3_linear_thermal` (B1) use the NaCl-B1 BE2 thermal
equation of Birch (1986). This attribution is now verified directly in
Walker's final article. Complete independent NaCl pressure replay remains
blocked by missing verified reference-isotherm parameters and reduction
conventions. Published KCl coefficients, errors, reference states, defaults,
and both raw CSV resources are unchanged.

Run `python scripts/audit_walker_2002.py --check`. The
[numerical report](../data/walker-2002-reproduction.json) records the independent
fits, all 69 bundled NaCl input rows, their actual Kelvin temperatures,
conventional cell volumes, known thermal increments, and explicit nulls for
unavailable independently replayed pressures. No executable Birch scale or
ready observation-pressure reduction is registered.

## Primary evidence actually inspected

The final [Walker article](https://doi.org/10.2138/am-2002-0701) was read from
the locally held PDF at `/Users/clemens/Zotero/storage/8FEA25CP/am-2002-0701.pdf.pdf`.
Its SHA256 is `01585fbd045709f3c08516aeedb30dcb553aec1d3073542c5450a72a6c3ff4b6`,
matching the variants source manifest. Full pages 806-810, including tables,
footnotes and the displayed strain definition, were inspected; pages 806-809
were also rendered for visual verification.

* Page 806, Experimental method: both pellets share pressure assigned by
  the NaCl BE2 thermal EOS of Birch (1986), and temperature measured by the
  type-S thermocouple. The B1 geometry measures NaCl and KCl together;
  the B2 geometry alternates illumination. Six B1 NaCl-only spot checks test
  gradients; they do not supply KCl fit volumes.
* Page 807: ESDs are spectrum/unit-cell refinement errors. No additional
  NaCl-EOS or temperature errors are assigned. The source says temperature
  was usually stable within one degree, which is not a published statistical
  temperature error to propagate.
* Page 808: the stated objective is the sum of squared pressure residuals
  between experimental pressure and pressure calculated from observed volume.
  Excel Solver performs the downhill minimization. The B1 volume is measured
  rather than fitted, K0' is fixed at 5 because the compression range does not
  resolve it, and all Table 1 sample data are used.
* Page 809, Table 3: the staging footnote applies to the **B2** rows. The bold
  row fits cold elastic coefficients first and then thermal expansion; the
  italic row is the separate simultaneous four-parameter fit. It does not
  prescribe a B1 cold-first fit.
* Page 810 and the abstract: the thermal products and their printed errors
  are better identified than their component parameters. Individual elastic
  errors and author covariance are unavailable.

The [original Birch (1986) publisher abstract](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JB091iB05p04949)
was inspected on the web and in the browser. It verifies

```text
P(V,t) = P25(V) + 0.0286 * (t - 25) kbar, with t in Celsius.
Equivalent thermal term = 0.00286 * (T - 298.15) GPa.
```

Its stated construction covers 25-500 Celsius and 0-300 kbar. It explicitly
says the earlier reference isotherm has small adjustments. The ePDF page
displayed an access-denial/institutional-login/purchase screen. No full text
or adjusted numerical coefficient table was recovered. Searches of the
local source holdings found the 1978 precursor, not the original 1986 paper.
Publisher article/PDF attempts and focused DOI/title searches supplied no
accessible original coefficient table. No purchase or author contact occurred.

The original [Birch (1978) precursor](https://doi.org/10.1029/JB083iB03p01257)
was inspected locally, including the rendered page 1258. Its PDF SHA256 is
`4461ec624cbd39f74af821d0302d3486a7ca82b8eede33921ec73b52c8f67770`.
Equations 8-9 and Table 1 distinguish BE1's linear reduced-strain polynomial
from BE2's quadratic polynomial. With compression-positive strain,

```text
f = ((V0/V)^(2/3) - 1)/2
P_BE2 = 3*K0*f*(1+2*f)^(5/2)*(1 + a*f + b*f^2)
a = (3/2)*(K0' - 4)
b = (1/6)*(9*K0*K0'' + 9*K0'^2 - 63*K0' + 143)
```

Thus BE2 is the fourth-order energy/three-elastic-parameter finite-strain
form, **not BM2**. This precursor verifies the convention; it cannot supply
the unverified adjusted 1986 coefficients. Neither a BM3 replacement nor
1978 coefficients are treated as Walker's exact calibrant.

## NaCl pressure replay: available inputs and missing steps

The ledger preserves all 30 Table 1 entries, including the derived reference
anchor and six NaCl-only checks, and all 39 bundled paired Table 2 sample
states. It converts raw Celsius by +273.15 and lattice a to conventional B1
cell volume a^3 (Z=4). It calculates only the verified thermal increment,
using **Birch's 298.15 K reference**, which differs from **Walker's KCl
296.15 K reference**. A thermal increment is not a complete pressure replay.
All 600 Celsius rows are flagged as above the original abstract's stated
500 Celsius limit; 23/24 Celsius rows are flagged as below its 25 Celsius
lower limit. These flags record source coverage, not a recovered
author extrapolation protocol.

Table 2 also prints separate calibrant ambient anchors r34439,
a=5.6414(14) A at 23 Celsius, and r35101, a=5.6468(4) A at 24 Celsius.
They are retained in the report, without adding fake KCl observations or
pressures. Table 1 gives a=5.6444(4) A at 23 Celsius as its reference anchor;
the two 36 Celsius pressures are bracketed imposed zeroes. These distinct
anchors make reference-lattice/run normalization a material question.

The remaining exact inputs/questions are:

1. What are the final Birch (1986) adjusted P25 BE2 K0, K0', K0'' (or a,b),
   reference volume/density, units, and rounding? Recover the original
   equation and numerical table, not a later unqualified catalog entry.
2. Did Walker use one absolute Birch reference cell, normalize each loading
   to its printed ambient NaCl anchor, or correct an energy/angle offset?
   How were the 23/24/36 Celsius anchors mapped to Birch's 25 Celsius state?
3. Was the 0.0286 kbar/K term extrapolated unchanged to 600 Celsius, or did
   the spreadsheet use temperature-specific tabulated isotherms/interpolation?
4. What unrounded NaCl a/T observations, pressure-reduction spreadsheet,
   and exact propagation prescription generated the printed pressures/ESDs?

Until these are verified, fit-to-printed-NaCl-pressure coefficients would be
an inverse diagnostic, not independently recovered Birch calibration.
No such inverse fit is used to manufacture an executable ancestry link.
Missing NaCl measurements are not the blocker.

## B1 fit: objective recovered, exact coefficients still differ

The 23 measured sample rows are fit with fixed conventional-cell
V0=249.53 A3, K0'=5 and Tr=296.15 K. The derived reference anchor is a fixed
input, not an extra measured residual. With beta=alpha0*K0, the source
pressure objective is linear in K0 and beta. A direct least-squares solve
therefore finds its global minimum for these inputs. A separate nonlinear
K0/alpha0 solve verifies the same minimum without changing parameterization
or adopting undocumented weights.

| Quantity | Published | Joint rounded-table fit |
|---|---:|---:|
| K0 (GPa) | 17.7 | 17.18012994 |
| beta (GPa/K) | 0.00195 +/- 0.00005 | 0.001872169315 |
| alpha0 (1/K) | 0.00011, rounded | 0.000108972943 |
| Pressure RMSE (GPa) | 0.07387843 | 0.06063608 |
| Sum of squared pressure residuals (GPa2) | 0.12553452 | 0.08456488 |

K0 differs by 2.94% and beta by 3.99%. The beta difference is about 1.56
times the printed error width; that width has no published confidence
interpretation, so it is not labeled a sigma or used to prove parity.
The fitted alpha0 rounds to the printed 0.00011, but this does not reproduce
the separately reported thermal product or elastic coefficient.

**A source-backed alternative materially improves agreement:** Figure 1 on
page 808 labels B1 V0 as **37.50 cm3/mol**, whereas Tables 1 and 3 print
**249.53 A3 and 37.57 cm3/mol**. Converting the figure's rounded molar
volume to a Z=4 conventional cell using the exact Avogadro constant gives
249.0808601 A3. Repeating the same 23-row joint objective at this reference
gives **K0=17.6834218 GPa**, **beta=0.00192714587 GPa/K**, and
RMSE=0.05382736 GPa. K0 rounds to 17.7 GPa, alpha0 rounds to 0.00011/K,
and beta differs by 1.17%, within the printed +/-0.00005 GPa/K width.
This identifies a primary-source reference-volume inconsistency capable of
accounting for much of the coefficient drift. It does not prove which volume
the Excel workbook used or establish exact beta/error reproduction. The
figure label is rounded; its conversion is a diagnostic, not a more precise
measurement. The catalog preserves the tabulated fixed V0=249.53 A3.

The printed BE1 pressure expression has minus signs inconsistent with its
positive-compression f. The standard BM3 signs are retained explicitly.
Fitting the literal signs gives K0=22.13552 GPa and beta=0.001777568 GPa/K,
RMSE=0.07983052 GPa, so it does not fix the discrepancy. Fixing published
K0 and fitting only beta gives 0.001836526 GPa/K and likewise does not recover
the published product; this is a conditional diagnostic, not a recovered
B1 staging rule. The B2 footnote cannot authorize undocumented B1 staging.

The Excel input workbook, initial values, stopping tolerance, any additional
constraints, and unrounded data are unavailable. These are possible missing
inputs, not established causes. The figure/table volume discrepancy above is
observed, while its role in the author's solver remains an inference.
No exact author solver or uncertainty replay
is claimed. The published coefficients are preserved.

## B2 staging and Ma boundary

An independent algebraic cold fit with V0=53.53 A3 and all eight 23-24 Celsius
rows gives K0=23.77539804 GPa and K0'=4.41587367. The thermal fit to all 39
rows then gives beta=0.002766245675 GPa/K and RMSE=0.10347680 GPa. This
confirms the existing staged reproduction independently of the native EOS
evaluator. The cold stage treats 24 Celsius as the reference isotherm;
the thermal stage retains its actual 297.15 K. Excluding that row as a
sensitivity check gives K0=23.56826768, K0'=4.50395014 and
beta=0.002787997669. The source does not deposit the exact cold-stage
temperature treatment.

Ma's downstream recalibration investigation is owned separately. This
follow-up verifies its reported Birch ancestry from Walker directly and
records the seven-23/one-24 Celsius distinction needed for any future 300 K
correction. It does not change Ma coefficients, fit selection or pressure
coordinates, and it does not extend a room-temperature recalibration to the
heated Walker observations.

## Integration and verification

This change is based on variants commit `30c5416e71ab3309572647deab444c52a83195a5`.
When combining with calibration commit `3ffb836a234eac78402be12aa99b7236430eda79`,
retain the latter's independent Dewaele/Campbell work and use this follow-up's
Walker attribution, unresolved replay notes and dedicated report. Recompute
manifest counts from the final merged records; do not copy either branch's
whole counts over the other. The older calibration report's claim that Walker
methods were unavailable is superseded by the primary evidence above.

Focused validation covers independent B1/B2 objectives, generated-report
reproducibility, actual Celsius offsets, paired source identities and hashes,
unchanged published errors, schema/API round trips, and rejection of both
unavailable Walker dataset reductions. Article PDFs are not redistributed.
Any contributor CC0 dedication applies only to contributor transcription,
normalization, metadata and arrangement, not third-party article rights.

Validation on this worktree: 110 tests passed across `test_walker_2002_audit`,
`test_kcl_variants`, `test_dataset_pressure` and `test_pressure_calibrations`;
all 143 `test_eosmat` tests passed. The final six Walker tests were rerun after
completing the source-domain flags. Both dedicated numerical reports pass
`--check`, Ruff passes for the new script/tests, and `git diff --check` is clean.
A direct comparison with the task base verifies that every KCl EOS/thermal
coefficient, published error, default and raw dataset resource is unchanged.
