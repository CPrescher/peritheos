# KCl pressure-calibration capability audit

Audit date: 2026-10-08. Scope: existing KCl records, observation pressure
coordinates and calibration provenance. Published EOS coefficients, original
CSV rows, parameter errors and thermal terms are unchanged.

## Dewaele (2012): reported-pressure normalization is ready

The [primary paper](https://doi.org/10.1103/PhysRevB.85.214105), experimental
methods and Tables I-II, assigns its ruby pressures to
Dorogokupets-Oganov (2007). That scale is already executable as
`ruby_dorogokupets_oganov_2007`. Both complete sample tables are bundled:

| Record | Dataset | Rows | Volume convention | Source-fixed parameters |
|---|---|---:|---|---|
| `kcl_b2_dewaele_2012_vinet_3` | `kcl_dewaele_2012_table1_compression` | 123 | Formula-unit volume equals B2 cell volume | V0=54.5 A^3 |
| `kcl_b1_dewaele_2012_vinet_1` | `kcl_b1_dewaele_2012_table2_compression` | 30 | Formula-unit volume multiplied by four for B1 cell | K0'=5.5 |

The old `missing_calibrant_observations` note incorrectly claimed that sample
P-V rows were needed. Both records now have `recalculation.status: ready`,
qualified specifically for normalization of **reported source-scale pressures**
and conditional reference-isotherm refits. Each table's `pressure_gpa` column
declares its scale and an `as_reported` all-row reduction for its Dewaele EOS.
`resolve_dataset_pressure` verifies the CSV checksum and returns this coordinate
without changing printed pressures or assigning independent uncertainties.

For the DOR equation, write x=lambda/lambda0-1:

```text
P_DOR = 1884*x*(1+5.5*x)
x = 2*(P_DOR/1884)/(1+sqrt(1+22*P_DOR/1884))
P_Mao = (1904/7.665)*((1+x)^7.665-1)
```

This inverts the corrected ruby signal at the precision of the printed
pressures. It does not recover original R1 peak fits, raw wavelengths,
reference wavelengths, unrounded values, or experimental corrections. The
paper reports general pressure/volume uncertainty limits, not row-wise errors
or covariance; these limits are not converted to standard deviations. The
unchanged published EOS ancestry already supplies a path from either Dewaele
record through DOR to other ruby scales and existing ruby/XRD bridges.

## Executable normalization and refit diagnostic

```python
from peritheos import (
    get_material_document,
    resolve_dataset_pressure,
    recalculate_ruby_pressure,
)

document = get_material_document("kcl")
rows = resolve_dataset_pressure(
    document,
    "kcl_b2_dewaele_2012_vinet_3",
    "kcl_dewaele_2012_table1_compression",
)
mao_pressure = recalculate_ruby_pressure(
    rows.pressure_gpa,
    rows.source_calibration["identifier"],
    "ruby_mao_1986",
)
```

Run `python -m scripts.audit_kcl_pressure_normalization` to regenerate
[`kcl-pressure-normalization.json`](../data/kcl-pressure-normalization.json),
or add `--check` to detect stale results. It retains every original row index
and resource checksum, converts all pressures, and fits the reference Vinet
isotherm with the source-fixed parameters. Fits minimize unweighted pressure
residuals because the original weighting and row-wise uncertainties are not
available. Derived coefficients are diagnostic output, not bundled EOS variants.

| Phase / pressure coordinate | V0 (conventional cell, A^3) | K0 (GPa) | K0' | RMSE (GPa) |
|---|---:|---:|---:|---:|
| B2 / reported DOR | 54.5 fixed | 17.24990 | 5.87310 | 0.80025 |
| B2 / normalized Mao | 54.5 fixed | 18.25133 | 5.53566 | 0.70690 |
| B1 / reported DOR | 249.37910 | 17.32529 | 5.5 fixed | 0.05348 |
| B1 / normalized Mao | 249.39910 | 17.44862 | 5.5 fixed | 0.05390 |

Twenty-three normalized B2 rows exceed Mao's declared 80 GPa pressure range;
the output marks their indices and retains them explicitly as extrapolations.
The diagnostic fits use the 298 K observations. Dewaele's published B2 equation
represents these with a 300 K reference term, and its alpha_KT comes from MD.
Neither a new thermal slope on Mao nor uncertainty parity can be inferred from
these static rows. A high-temperature EOS curve transformed through the graph
is a pressure-coordinate normalization, not a refit of new thermal observations.

## Walker (2002): Birch (1986) ancestry verified; exact replay unavailable

The final [Walker article](https://doi.org/10.2138/am-2002-0701), page 806,
explicitly identifies the NaCl-B1 BE2 thermal EOS of Birch (1986) for both
phases. This supersedes the initial audit's unidentified-scale assessment.
Both Walker records now have `partially_resolved` calibration provenance and
`reference_eos_not_bundled` replay status. BE2 includes a quadratic strain
polynomial; substituting BM2 or unverified precursor coefficients is unsupported.

All 30 B1 entries and 39 paired B2 sample states preserve NaCl lattice a,
pressure/ESDs, sample volumes and actual Celsius temperatures. The eight cold
B2 rows comprise seven at 23 degC and one at 24 degC. The
[dedicated Walker report](walker-2002-kcl.md) verifies the source objectives,
reference-volume discrepancy and available Birch thermal increment, while
leaving missing adjusted 1986 elastic coefficients, reference/run normalization
and high-temperature treatment explicit. No complete independently replayed
pressure or executable calibration edge is manufactured.

The derived reference anchor and six B1 calibrant-only checks are retained as
source evidence, not additional measured KCl residuals. Ma's room-temperature
recalibration must not be applied to all heated Walker rows. Source Table 2
ESDs omit NaCl EOS and additional temperature errors.

## Campbell-Heinz (1991): Mao (1978) ancestry verified

The recovered final methods, pages 495-496 and reference 5 on page 499,
identify Mao et al. (1978) ruby fluorescence. All 14 reported Table 1 mean
pressures now declare an executable `as_reported` coordinate on
`ruby_mao_1978`. Each mean has five spatial readings; the reported errors are
spatial standard deviations rather than standard errors of the mean.

The [source-gap audit](kcl-tateno-campbell-source-gaps.md) documents conditional
recovery of the weighted normalized-stress regression and the remaining
missing individual ruby readings, numerical fit weights and covariance.
Nonlinear conversion of a reported mean cannot recover the mean of five
separately transformed readings. The published coefficients and explicit
Campbell-Dewaele composite volume construction remain unchanged.

## Verification

Tests cover schema/export preservation, all 153 Dewaele row identities and
checksums, independent DOR-to-Mao algebra through the existing ancestry graph,
source-fixed refit constraints, generated audit reproducibility, manifest
recalculation counts, and rejection of unsupported Walker reductions. Campbell's verified reported
pressure contract is covered by the source-gap regression tests.
