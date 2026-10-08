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

## Walker (2002): observations exist; exact scale remains unverified

The bundled [Walker primary table](https://doi.org/10.2138/am-2002-0701)
contains 39 B2-KCl observations, including NaCl lattice a and its ESD,
NaCl-derived pressure and its ESD, KCl volume, and temperature. Eight rows are
at 23 degC. There is no need to recover missing NaCl observations before a
future verified reference-EOS re-reduction.

[Ma (2024), section 3](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024JB028819)
identifies Walker's original NaCl-B1 scale as Birch (1986) and recalculates
room-temperature rows using Matsui (2012), correcting the small temperature
difference to 300 K. This is a strong discovery lead, not verification of the
original methods. The MSA PDF and publisher article could not be retrieved in
this recheck; the author-upload page exposes figures and a reference list but
not the original methods text. No local original copy was found. Neither the
Birch (1986) nor Matsui (2012) EOS is bundled in the current NaCl-B1 card.

Consequently the Walker calibration remains `unresolved`, with an explicit
unresolved dataset reduction and a revised note naming the actual evidence and
blockers. No guessed `reference_eos_record` or executable graph edge is added.
Before enabling it, verify the original methods and exact NaCl equation,
parameter convention, cell volume a^3 (four formula units), reference
temperature, validity, and thermal correction. A Ma room-temperature subset
reduction must not be applied to all 39 heated rows. Original Table 2 ESDs do
not include NaCl scale uncertainty or additional temperature error.

## Campbell-Heinz (1991): primary methods still needed

All 14 sample observations are bundled. The accessible
[primary abstract](https://www.sciencedirect.com/science/article/abs/pii/002236979190181X)
does not identify the pressure gauge or calibration equation; no accessible
primary methods evidence was recovered in this recheck. The scale remains
unresolved and its dataset reduction now carries that explicit reason.
No ruby scale is inferred from the DAC apparatus, publication date, or later
papers. The composite B2 reference-volume construction remains unchanged.

## Verification

Tests cover schema/export preservation, all 153 Dewaele row identities and
checksums, independent DOR-to-Mao algebra through the existing ancestry graph,
source-fixed refit constraints, generated audit reproducibility, manifest
recalculation counts, and rejection of unresolved Walker/Campbell reductions.
