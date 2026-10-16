# Corrected Zhao diopside P V T dataset

Use `diopside-zhao-1998-table1-corrected-pvt.csv` for continued pressure–volume–
temperature analysis of Zhao et al. (1998), DOI 10.4131/jshpreview.7.25. It
contains all 32 printed rows in their original order, with the supported
pressure correction applied to row 15. This is a Peritheos analysis dataset,
not an author-issued erratum or a replacement for the original transcription.

## Correction and provenance

Table 1 prints row 15 as 6.95 GPa, 1072 K and 431.66(6) Å³. Its matching
Figure 2 marker reads approximately 4.5803 GPa; the corrected dataset stores
**4.58 GPa**, avoiding unjustified graphical precision. The independent
graphical reading bound is **±0.05 GPa**. This reading is consistent with
4.5630 GPa from the published EOS, but the replacement was digitized from the
plot rather than calculated from the EOS. The original exact experimental
pressure remains unavailable.

The CSV's `pressure_gpa` is the value to use for analysis. It changes only
on row 15. `source_pressure_gpa` preserves every original pressure, including
6.95 GPa on that row. `pressure_correction_applied`, `pressure_source` and
`correction_note` make the replacement explicit. All temperatures, cell
volumes and their printed errors remain unchanged.

`pressure_graphical_reading_bound_gpa` is a digitization reading bound, **not
an experimental standard uncertainty**. The original source does not print
row-wise pressure or temperature errors, so `pressure_error_gpa` and
`temperature_error_k` remain blank. A weighted fit must not silently treat
the graphical bound as an experimental sigma. The original Decker (1971)
NaCl pressure scale is retained; paired NaCl readings are unavailable.

The inconsistent source lattice dimensions are not included in this P–V–T
view because no defensible corrected a, b, c or β values are available. They
remain in `diopside-zhao-1998-table1.csv` and the original numerical notation
in `diopside-zhao-1998-table1-original.txt`. `quality_notes` marks rows with
known source lattice/volume discrepancies; use the printed cell-volume column
for the EOS analysis.

Row 32 repeats the printed ambient values of row 1 and is marked with
`duplicate_of_source_row=1`. No row is automatically excluded. Choose and
report how to handle this repeated entry when refitting, especially if V0
is varied.

## Reproduction and access

From the repository root, regenerate the plot analysis and corrected dataset:

```bash
python -m scripts.diopside_check_plots
python -m scripts.diopside_correct_dataset
```

`python -m scripts.diopside_correct_dataset --check` verifies the export.
The digitization coordinates and source-PDF hash are in
`docs/data/diopside-plot-check.json`. Its independent refit uses precisely the
rounded 4.58 GPa pressure exported in this CSV. The material document
`peritheos/data/materials/diopside.eosmat` registers this resource as
`diopside_zhao_1998_table1_corrected_pvt`, with its checksum and column units,
and recommends it for refitting.

Published EOS coefficients remain unchanged. Source rights remain unspecified;
`diopside.LICENSE.md` limits the CC0 dedication to contributor transcription,
normalization and metadata rights, including this documented analysis correction.
