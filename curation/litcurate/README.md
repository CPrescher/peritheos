# LitCurate configuration for Peritheos

This directory contains the Peritheos-specific extraction contract used after
LitCurate has discovered and screened a paper. LitCurate remains an external
curation tool; it is not a Peritheos runtime dependency.

Start from LitCurate's `configs/examples/minerals_eos.yaml`, then replace the
extraction schema and prompt with:

```yaml
extraction:
  schemas:
    - name: peritheos_eos_candidate
      path: /absolute/path/to/peritheos/curation/litcurate/peritheos_eos_candidate.schema.json
      prompt: /absolute/path/to/peritheos/curation/litcurate/peritheos_eos_candidate.md
      version: "1"
      empty_list_field: eos_entries
```

Keep the LitCurate `runs/` directory outside this repository. Import its final
`database.json` with:

```bash
python scripts/import_litcurate_eos_candidates.py \
  /path/to/database.json \
  docs/data/litcurate-eos-candidates.json
```

The importer also accepts the public `data.xlsx` deposit from
<https://doi.org/10.5281/zenodo.22118629>. Its output is a discovery ledger,
not an executable material library. Every promoted record must independently
satisfy `docs/adding-materials-and-eos.md`.
