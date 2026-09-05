# LitCurate candidate workflow

LitCurate is used as a discovery and structured-extraction tool upstream of
Peritheos. Its results never become executable equations automatically.

The public lower-mantle deposit contains 1,334 reported rows from 205 papers.
The checked-in candidate ledger preserves 1,309 rows from 204 papers and classifies
their review state without copying LitCurate's evidence excerpts. It is
attributed to Shakya et al. and the deposited CC-BY-4.0 dataset:
<https://doi.org/10.5281/zenodo.22118629>.

The 25 rows for DOI `10.2138/am-2024-9562` are intentionally filtered during
import. A primary-source audit found that its temperature-indexed BM3 rows are
slices of a first-principles thermal calculation and cannot be represented as
independent static records by the current model registry.

## Rebuild the candidate ledger

Download `data.xlsx` from the dataset record and run:

```bash
python scripts/import_litcurate_eos_candidates.py \
  /path/to/data.xlsx \
  docs/data/litcurate-eos-candidates.json
```

The same command accepts LitCurate's native merged `database.json` or a flat
CSV export. XLSX input is read with the Python standard library and does not
add a spreadsheet package to Peritheos.

The ledger records the input filename and SHA-256 checksum, dataset DOI and
license, original worksheet row, verbatim reported EOS values, current
Peritheos DOI overlap, current documented-backlog overlap, and deterministic
review classifications.

## Review buckets

| Bucket | Meaning |
|---|---|
| `review_first` | Source-reported, DOI not bundled, explicit existing equation family, named phase, and complete V0/K0/Kp values |
| `primary_source_audit` | Existing family candidate with incomplete core parameters or unresolved phase identity |
| `equation_audit` | Equation family or order is ambiguous and must be recovered from the primary paper |
| `model_work` | The reported formulation may need a new family, adapter, or complete thermal treatment |
| `citation_trace` | The extracting paper copied the value; locate the underlying primary publication |
| `already_bundled_source` | The extracting paper's DOI already occurs in an executable Peritheos record |
| `manual_triage` | Missing DOI or source-origin identity |
| `not_an_eos_fit` | LitCurate explicitly classified the row as not being an EOS fit |

`review_first` means suitable for expert review, not suitable for automatic
execution. Even an explicit BM3 label must be checked against the printed
equation, volume basis, reference state, fixed parameters, uncertainties,
pressure calibration, and primary observations.

## Run a new LitCurate campaign

Keep LitCurate and its `runs/` directory outside this repository. The
Peritheos-specific extraction schema and prompt live under
`curation/litcurate/`. They request the additional evidence required by the
[material and EOS contribution policy](adding-materials-and-eos.md), including
complete equations, higher-order and thermal coefficients, crystallographic
identity, pressure calibration, primary-data locations, conflicts, and explicit
blockers.

After export, rebuild this ledger and assign one primary publication to one
isolated review task. Only the task that completes the primary-source and
numerical audit may create or modify an `.eosmat` record.
