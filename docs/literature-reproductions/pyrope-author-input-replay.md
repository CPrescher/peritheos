# Pyrope: Ross Angel's supplied inputs and fit replay

14 November 2026. Clemens supplied Ross's reply, two processed data files,
an EOS file, the final fitting log and a six-page dataset description.
**All 194 fitted observations, their assigned coordinate errors and their
scale groups are now available.** Replaying these inputs removes the large
thermal-parameter discrepancy of the earlier conditional reconstruction.
Original attachments remain in the private archive; source coefficients and
earlier recovery/refit artifacts remain unchanged.

**Subsequent numerical investigation:** the
[residual/weight/Jacobian audit](pyrope-numerical-differences.md) isolates the
remaining small differences. The audited console uses finite one-error shifts
for coordinate weights and a coarser modulus inversion. Its captured numerical
Jacobian and weights reproduce its errors in Python to single-precision
rounding. An analytic-derivative replay prototype preserves Python's original
objective and results while accelerating this local fit.

## Direct comparison, including errors

The reference is the **25 October 2021 final log**, not the slightly earlier
coefficient snapshot inside the supplied EOS file. Our official EosFit console
and Ross's log both identify version **7.60**. We loaded the unchanged DAT
files, selected their volume/modulus columns, assigned group 3 to rows
118–132, fixed V₀ and groups 1–2, and refined the four EOS coefficients plus
the group-3 scale using all P, T and value errors.

All ± values below are **local residual-rescaled parameter standard errors**,
conditional on the supplied inputs and the effective-variance fitting protocol.
They are not bounds on total experimental or calibration error.

| Parameter | Ross's final log | Our EosFit 7.60 replay | Our independent Python replay |
| --- | ---: | ---: | ---: |
| K₀, GPa | 169.33849 ± 0.31290 | 169.33817 ± 0.31334 | 169.33983 ± 0.31220 |
| K′ | 4.55313 ± 0.04503 | 4.55321 ± 0.04510 | 4.55290 ± 0.04488 |
| θ, K | 771.10760 ± 27.93594 | 771.05811 ± 27.94607 | 770.86593 ± 27.94708 |
| γ | 1.18471 ± 0.01160 | 1.18468 ± 0.01160 | 1.184645 ± 0.011598 |
| Sinogeikin–Bass high-P scale | 0.97828 ± 0.00588 | 0.97828 ± 0.00588 | 0.978234 ± 0.005855 |
| Reduced χ² | 2.2335 | 2.2337 | 2.234252 |

**EosFit does not reproduce every printed digit of Ross's final log.** Its
largest differences remain below 0.003 times the corresponding reported
standard error. The independent fit differs by less than 0.009 times each
reported error. These ratios describe numerical closeness; they are not
independent-sample significance tests. Both Python starting points converge
to the same solution without parameter bounds being reached.

Ross explains that single precision and a shift/error stopping rule prevent
strict rerun identity. The native replay also terminates after one small-shift
cycle. Its final 194 observed-row/error/group token sequences match his log
exactly at the log's printed precision. Residual differences are small and
compatible with numerical/stopping effects; we have not isolated every
remaining contribution or demonstrated identical executable builds.

The supplied March2023-named EOS file contains all **59 identical parameter
values and the identical full 59×59 covariance matrix** already in the
archived official EOS. Its title adds the citation; the bytes differ. It was
saved on 8 October 2021, so its covariance and coefficients must not silently
replace the later final-log snapshot in this comparison.

## What the paper and supplement did not determine uniquely

The published inventory identified the included studies and counts. It did
not uniquely specify every processed value, assigned error or group membership.
The following compares Ross's files with our latest **191-channel conditional
reconstruction**, preserving that reconstruction as a historical diagnostic.

| Item | Earlier reconstruction | Actual supplied fit input |
| --- | --- | --- |
| Dataset completeness | 191 of 194 channels; 9 readable Sinogeikin–Bass 2002 modulus points | All 194: 115 volume and 79 modulus observations; all 12 Sinogeikin–Bass 2002 points |
| Hartwig weighting | σ(V)=0.03 cm³/mol and σ(T)=1% T; uniform whole-row multiplier tested separately | σ(V)=0.01 cm³/mol; temperature errors remain approximately 1%, rounded to the listed values |
| Milani thermal data | Original converted molar volumes; unknown temperature errors omitted; free normalization scale | Already normalized volumes about 0.097% smaller; σ(T)=7 K for all 33 rows; no free thermal-data scale |
| Thermal-series scales | Free Milani and Hartwig scales as conditional controls | All volume rows share fixed scale 1; only the high-P modulus scale is free |
| Sinogeikin–Bass 2000 scaling | One scale applied to all 17 digitized points | Only 15 pressure-dependent points scaled; two ambient points remain unscaled |
| Sinogeikin–Bass 2002 values/errors | 9 points digitized through C₁₁/C₁₂; σ(T)=5 K comparison width | 12 explicit modulus values; σ(T)=3 K; actual modulus errors 1.59–1.71 GPa |
| Zou high-T modulus temperatures | 2% T from the supplementary description | 5 K at 300 K and 20 K at higher T |
| Chantel errors | Missing P/T errors omitted; modulus errors propagated from velocity/density errors assuming independence | Explicit approximately 5% P, 2% T and 2.5% modulus errors; separately assigned volume errors |
| Room-temperature datasets | 298 K where room T was not uniquely specified | Zhang and Milani PV: 295±2 K; Du PV: 296±2 K |
| Exact normalization/calibration | Reconstructed anchors and recalculated quartz pressures | Final preprocessed values supplied; Milani PV pressures differ by up to 0.0469 GPa and Chantel volumes by up to 0.0532 cm³/mol from our reconstruction |
| Covariance convention | Unscaled and residual-scaled widths reported separately; original convention unresolved | Log explicitly reports residual-rescaled covariance and errors |

Reducing Hartwig's volume error from 0.03 to 0.01 increases the **volume-error
contribution's inverse variance by nine**, not necessarily the complete row
weight by nine because temperature errors also enter. Our earlier diagnostic
scaled both volume and temperature errors together. That was not the actual
author treatment. The supplied dataset description still lists the initial
0.03 cm³/mol estimate, making the final DAT file essential.

Several central values and errors were only approximately recoverable from
figures or rounded primary tables. For example, Zou's supplied PV volumes
differ by up to 0.00757 cm³/mol from our reconstruction from rounded V/V₀.
These are processed author inputs, not new measurements or corrections to
the original publication. The combined replay establishes their overall
effect; it does not uniquely attribute the fit shift to any single change.

## Independent check of our errors

There are **194−5=189 degrees of freedom**. We compute local inverse-
information covariance with the final effective coordinate weights frozen,
then multiply covariance by χ²ᵣ=2.234252 and widths by **√χ²ᵣ=1.494741**.
For θ, the unscaled width is **18.69693 K** and the reported rescaled width
is **27.94708 K**, versus **27.94607 K** in the EosFit replay.

An independent five-point pressure Jacobian at three step sizes and a QR
inverse-information calculation reproduce our three-point/SVD widths to
better than 2×10⁻¹⁰ relative. The covariance is positive definite; the
inverse-information identity is satisfied to about 10⁻¹⁴. The θ–γ correlation
is **+0.86108**, consistent with Ross's printed +86.1%. Our rescaled widths
differ from the native replay by at most **0.493%** across all five parameters.
The saved native covariance independently yields its printed errors once
the covariance-file and log rounding are included.

This checks our implementation of the matched local error convention. It
does not establish raw-coordinate independence, total systematic uncertainty
or a Gaussian confidence interpretation. The approximate q-compromise KS
conversion remains the source's operational convention.

## Fixed-parameter calculation and remaining scope

Following Ross's proposed separation of calculation and fitting, we also
checked an exactly representable **proposed** coefficient set: V₀=112 cm³/mol,
K₀=170 GPa, K′=4.5, θ=800 K, γ=1.25, n=20 and Tᵣ=298 K. Ross has not yet
agreed to this particular set. At 12 states spanning 0/10/30 GPa and
100/298/1000/1700 K, independent volume, Kₜ and source-convention Kₛ agree
with EosFit within the previously declared native printed-precision
tolerances: 0.0001 cm³/mol for volume and 0.001 GPa for moduli. This uses
EosFit's gas-constant convention; no library default was changed.

We can now replay the **supplied processed regression** closely, including
its uncertainty scale. Ross explicitly retains some upstream preprocessing
and exclusion decisions outside these attachments. Exact original anchors,
raw calibration covariance and all reasons for exclusions are therefore
not reconstructed. Full experimental validation and an end-to-end replay
from raw studies are separate claims. No replacement EOS was installed,
no grossular input changed and no response was sent to Ross.

## Artifacts and repeatability

- [Author-input audit and both independent fits](../data/garnet-angel-2022/pyrope/author-2026-11-14-audit.json)
- [All 194 input rows with original tokens](../data/garnet-angel-2022/pyrope/author-2026-11-14-inputs.csv)
- [Author log's calculated values and weights](../data/garnet-angel-2022/pyrope/author-2026-11-14-log-calculations.csv), calculated checks, not additional observations
- [EosFit comparison, saved covariance and fixed-parameter checks](../data/garnet-angel-2022/pyrope/author-2026-11-14-eosfit.json)
- [Verification receipt](../data/garnet-angel-2022/pyrope/author-2026-11-14-verification.json)

Repository-only commands (the second requires the external official console
and an already configured runtime/display):

```bash
python -m scripts.audit_pyrope_author_inputs --source-dir <private-attachments>
python -m scripts.run_pyrope_author_eosfit --source-dir <private-attachments> --executable <eosfit7c> --output <private-replay-directory>
pytest -q -W error tests/test_pyrope_author_inputs.py
```
