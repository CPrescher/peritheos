# Source observation views

An EOS record, a source table and the observations to plot are different things.
One table can contain several samples, phases, pressure gauges, volume bases,
fit masks, and calculated comparison values. Fitting eligibility must not decide
whether its numerical observations are available for display.

EOSMAT v3 datasets accept optional `observation_views`. A view explicitly names
its EOS records, evidence kind, row selection, pressure and sample-volume
coordinates, source temperature column, uncertainty relationships and optional
fit selection. Evidence distinguishes experimental, digitized, shock,
experimental reduction, simulation, model output and source comparison.

The JSON definition is `peritheos/data/observation-view-v1.schema.json`, also
embedded in the EOSMAT v3 schema. Python validation checks references, units,
uncertainty roles and ambiguous record/table views. Existing files remain valid
without this extension; other consumers can preserve it without resolving it.

## Source-preserving migration

`peritheos/data/observation-views.json` contains the same contracts for existing
tables, without rewriting their declarations, CSVs or EOS coefficients. SHA-256
binds each complete dataset declaration (including the resource checksum) and
material identity/cell metadata. A stale binding is an error. Migration of an
existing display interpretation is not new publication validation or reproduction
of the authors' fitting workflow.

Prefer embedded views for new contributions. Do not declare both an embedded
and registry view for the same record/table pair. A view can explicitly restore
an omitted legacy record/table link.

```python
from peritheos import get_material_document, resolve_dataset_observations

data = resolve_dataset_observations(
    get_material_document("nacl_b2"),
    "nacl_b2_shen_2026_vinet_2",
    "shen_smith_2026_table_s1_simultaneous_volumes",
)
print(len(data.row_indices))  # 131 measurements
print(data.pressure_gpa, data.volume_a3)
```

`plan_dataset_observations` returns a validated contract or a specific unresolved
reason. Automatic plans are restricted to unambiguous typed, single-phase P–V
columns with an explicit conventional-cell basis. Complex tables need views.

The resolver supports reported pressure columns, source-declared constant
pressure, explicit ruby-scale conversion and pressure calculated from paired
calibrant volumes. Volumes support explicit multiplicative or reciprocal-density
conversions and a complete lattice metric with declared constraints. Calibrants
use the existing full `EOSRecord.pressure` evaluator, including nonstandard
cold/thermal terms, rather than a second EOS formula.

Results retain original zero-based row indices, calibrant row pairs, omission
reasons, source hashes and the full calibration record. Arrays are read-only.
Missing sample temperatures remain NaN; a calibrant calculation temperature does
not become a measured sample temperature. Missing T never erases valid P–V rows.
Standard errors remain standard errors. Calibration/lattice errors are not
propagated without a specified model. Reciprocal-density errors use first-order
coordinate propagation with nominal mass and Z treated exact. An absent fit
selection returns `selected_for_fit=None`, rather than declaring every row
excluded.

For Shen–Smith (2026), pairing requires experiment **and** run number. Missing or
duplicate calibrants produce explicit omissions without averaging. Cu uses the
full Fratanduono (2020) EOS at 298 K; no 298-to-300 K correction is applied.
NaCl-B1/B2 retain all 67/131 observations. Separately labelled audited fit masks
contain 55/125 rows. The B2 above-33-GPa table statement versus approximate
30–35-GPa coexistence interval remains documented; this is not an exact recovered
author mask or a new fit.

## Coverage and consumer exports

```sh
python scripts/audit_observation_coverage.py --output /tmp/observation-coverage.json
python scripts/export_observation_views.py --revision <immutable-commit> \
  --output /tmp/observation-views.json
```

The census inspects every stored EOS, including experimental/hybrid records
without datasets. It distinguishes plottable experimental evidence, unresolved
coordinate interpretations and missing numerical sources. Model grids do not
count as experimental observations.

Native observables remain available through `Dataset`: V/V0 without a published
normalization, V–T without pressure, acoustic data and paired marker/sample
volumes must not be converted into absolute P–V measurements using a selected
EOS's V0, invented pressure or another phase's basis. Consumers can show these
original coordinates separately and explain why the EOS overlay is unresolved.

The exporter checks every binding against an immutable Git revision and reports
incompatible declarations. Scientific revision and registry content hash are
separate: a local metadata export does not upgrade a consumer's engine, tables
or scientific release pin. Tests validate all registered views and resolve every
distinct interpretation against bundled observations. Catalog changes must
rerun this census and review affected contracts.
