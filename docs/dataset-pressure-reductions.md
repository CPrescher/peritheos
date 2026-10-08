# Dataset pressure provenance and EOS reductions

A pressure coordinate belongs to a reported observation table and a declared
calibration. An EOS can use that coordinate unchanged or after a documented
transformation. A dataset link alone does not establish either choice.

These optional, additive format-3 fields make the distinction executable:

| Location | Meaning |
|---|---|
| `datasets[].columns[].pressure_scale` | `calibration_record` identifies the scale of this pressure value column; `reference` and `source_location` establish the attribution. |
| `datasets[].pressure_reductions[]` | One declaration per linked `eos_record`, within this dataset. |
| `status: as_reported` | Use `pressure_column` unchanged, apart from unit conversion. Its scale must equal `target_calibration_record`. |
| `status: transformed` | Use `pressure_column` as the source coordinate and apply the named `convention` to obtain `target_calibration_record`. |
| `status: unresolved` | `notes` explains why a coordinate cannot be established. No pressure column, target, row scope, or convention is permitted. |
| `row_scope: all_rows` | Every existing dataset row, in its original order. Required for resolved declarations. Subset and mixed-gauge reductions are not currently executable in this representation. |
| `reference`, `source_location`, `notes` | Each reduction requires scientific provenance; notes are optional except for unresolved declarations. |

There is no implicit default pressure column or calibration. Missing fields in
older documents remain valid and mean unresolved. `used_by_eos_records` keeps
its broader meaning: a table may inform a self-consistent reference model
without supplying a simple experimental pressure reduction. Multiple EOS
reductions of one dataset share its observation identities; they are alternative
analyses, not independent experiments. Group by dataset identifier and preserve
row identity as `(dataset identifier, zero-based row index)`.

The JSON Schema validates field shapes. Python and Rust validators additionally
reject unknown calibration, column, or EOS references, unlinked EOS records,
duplicate reductions for one EOS, non-pressure source columns, incompatible
as-reported scales, unsupported conventions, and conflicts with the EOS
`pressure_calibration`. A resolved reduction currently requires the EOS to have
one resolved calibration method matching its target. Ambiguous multi-method
cases must remain unresolved. No catalog-wide mapping is inferred.

## Transformation convention

The supported convention is `same_corrected_ruby_r1_ratio`. It means:

1. Convert the original reported pressure units to GPa.
2. Invert the source ruby calibration to the corrected R1 ratio `lambda/lambda0`.
3. Evaluate the target ruby calibration at that same ratio.

Use `recalculate_ruby_pressure` (Python) or the existing
`pressure_calibration::recalculate_ruby_pressure` (Rust). Do not use the EOS
sample volume to construct the observation pressure, interpolate between
scales, or change the ratio to accommodate the EOS reference temperature.
The convention cannot recover unrounded measurements or undo temperature,
stress, composition, or peak-fitting effects.

The data migration distinguishes three calibrations:

| Calibration ID | Equation (GPa), with `r = lambda/lambda0` |
|---|---|
| `ruby_mao_1986` | `(1904 / 7.665) * (r**7.665 - 1)` |
| `ruby_dewaele_2004` | `(1904 / 9.5) * (r**9.5 - 1)` |
| `ruby_dorogokupets_oganov_2007` | `1884 * (r - 1) * (1 + 5.5 * (r - 1))` |

The printed 2004 revised pressure is **not** DOR. The 2019 DOR reduction inverts
the printed classical Mao pressure, then evaluates DOR. The as-reported 2004
revised reduction preserves the separately rounded printed values, including
the Cu pair `37.0 / 37.1 GPa`; it does not recompute or repair that pair.

## Python consumption

```python
from peritheos import get_material_document, resolve_dataset_pressure

document = get_material_document("aluminum")
eos_id = "aluminum_dewaele_2019_dor_vinet"

# Discover links by metadata; do not derive dataset names from EOS suffixes.
linked = [d for d in document.get("datasets", []) if eos_id in d["used_by_eos_records"]]
for metadata in linked:
    result = resolve_dataset_pressure(document, eos_id, metadata["identifier"])
    pressure_gpa = result.pressure_gpa
    original_pressure = result.reported_pressure
    original_unit = result.reported_pressure_unit
    row_indices = result.row_indices
    source_id = result.source_calibration["identifier"]
    target_id = result.target_calibration["identifier"]
    extrapolated = result.source_pressure_extrapolated

    # Select volume by quantity/role; preserve its declared basis.
    volumes = result.dataset.find_columns(quantity="atomic_volume", role="value")
    assert len(volumes) == 1
    atomic_volume = result.dataset[volumes[0].name]
    cell_volume = atomic_volume * document["formula_units_per_cell"]
    assert len(pressure_gpa) == 40
```

`resolve_dataset_pressure(document, eos_record, dataset_identifier, *,
resource_root=None, check_validity=False) -> DatasetPressure` validates the
material and verifies CSV checksums through the existing loader. It returns
read-only pressure arrays, original reported pressures and units, row indices,
the original loaded dataset, reduction metadata, both full calibration
documents, and `eos_reference_temperature_k`. Unknown datasets and absent or
unresolved reductions raise `DatasetError`; malformed declarations raise
`EosmatError`. A UI should display the unresolved reason instead of guessing.

`source_pressure_extrapolated` and `target_pressure_extrapolated` are masks for
the respective calibration pressure ranges; `None` means no range is declared.
The default keeps all observations. `check_validity=True` rejects any pressure
outside either declared range. This is a pressure-only check, not certification
of temperature, hydrostaticity, or all other calibration conditions. The full
calibration validity metadata remains available for display. Aluminum has 12
rows above Mao's published 80 GPa range; its DOR coordinate reaches about
155.29074838 GPa. They remain valid *transcribed observations*, with the source
scale extrapolation explicitly qualified.

Calibration reference conditions, EOS reference temperature, and observation
temperatures have separate meanings. The current calibration library gives
Mao and DOR `validity.temperature_k = [298, 298]`; the 2019 EOS reference is
300 K. The source table does not report row temperatures. Neither 298 K nor
300 K is assigned to those rows, and this conversion applies no temperature
correction. `DatasetPressure` preserves calibration documents and the EOS
reference separately; observation temperatures, when present, belong to the
original dataset columns.

`uncertainty_treatment = "not_propagated"` explicitly records that this API
performs coordinate conversion only. Original uncertainty columns, roles,
notes, and EOS coefficient uncertainties remain accessible and unchanged.
Shared scale errors must not become independent per-row errors, and the
printed volume bound must not become a standard deviation. Unreported
pressure errors and covariance remain unreported.

## Migration and Studio integration

The Dewaele (2012) KCl tables also declare DOR source-scale columns and
as-reported reductions: all 123 B2 and 30 B1 rows. Their reported-pressure
normalization and conditional reference-isotherm refits are ready. See the
[KCl calibration capability audit](literature-reproductions/kcl-pressure-calibration.md)
for executable diagnostics, volume conventions, raw-wavelength and uncertainty
gaps, and explicitly unresolved Walker/Campbell reductions. The B2 table's
separate Chidester link has no declared reduction here; this change does not
infer the pressure treatment of that combined fit.

The six Dewaele (2004) Table I datasets now carry both column attributions:
Al (40 rows), Cu (42), Au (37), Pt (36), Ta (36), and W (42). The migration
adds 12 as-reported 2004 EOS reductions, plus 10 2019 Mao/DOR reductions for
Al, Cu, Pt, Ta, and W. Gold's 2019 EOS uses the separate Takemura (2008) dataset
and is outside this migration. The linked Pt Dorogokupets–Oganov (2007) model
and Fei (2007) optimization remain explicitly unresolved for simple all-row
pressure reduction. Other catalog metadata is untouched.

For exact identifiers, read each dataset's `pressure_reductions`; tests verify
all 22 resolved mappings, both aluminum reductions over all 40 rows, the
printed revised distinction, existing audit fits, and Python/Rust round trips.
Scientific provenance is the
[2004 six-metal reproduction](literature-reproductions/dewaele-2004-six-metals.md)
and the [2019 audit](literature-reproductions/dewaele-2019-static-dac-metals.md),
including its independent numerical regressions. Raw CSV resources, measured
values, coefficients, uncertainties, and existing EOS records are unchanged.

Studio currently pins Peritheos `fa322960c7fd1e2e43b94c833b07f81dad35e7cd`.
This upstream change does not update that pin. After review, repin to a commit
containing this change and regenerate any metadata exports while preserving
these optional dataset/column fields. A separately reviewed metadata supplement
is an alternative for a consumer retaining the old pin; it must include the
explicit mappings and matching schema/validation support. The old pin alone
cannot resolve them.

Studio can then replace its interim source adapter with the metadata-driven
lookup above. For non-Python consumers: validate the document, resolve the
named column and calibration IDs, dispatch only the declared convention through
the existing calibration API, and retain the same original pressures, units,
row identities, provenance, and validity flags. Rust eosmat loading and export
preserve and validate these fields; this change does not add a Rust CSV loader.
Python exports of selected EOS records retain only those records' pressure
reductions, alongside their filtered dataset links.
Keep absent or unresolved mappings visibly unresolved. Do not merge alternative
EOS records or count their shared rows as independent experiments.
