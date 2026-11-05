# Optional Peritheos datasets

`peritheos-datasets` adds selected numerical tables to the Peritheos library.
Install it through `python -m pip install "peritheos[datasets]"`. The installer
selects matching versions of both distributions.

The normal `peritheos` package already contains the calculation library,
structures, EOS parameters and dataset descriptions. Use it on its own for
applications such as Dioptas and for fitting your own observations. Installing
this optional package adds values; it does not install a second EOS implementation
or download additional source files at runtime.

## What is included in 0.18.0?

| Content | Amount and location |
| --- | --- |
| Material dataset entries | 451, with declarations in `peritheos_datasets/declarations/` |
| Supporting-study entries | 2 Barker tables, in `peritheos_datasets/study-declarations/` |
| Numerical CSV resources | 392, under `peritheos_datasets/data/datasets/` |
| Small inline tables | 14, within the declaration JSON files |
| Source notices and inclusion decisions | `peritheos_datasets/notices/`, `manifest.json` and `redistribution-review.json` |

The **453 entries** can share resource files. They include measurements, numerical
figure readouts, explicitly calculated tables and a partial source-lineage graph;
these are not 453 independent experiments. Original identifiers, units, errors,
attribution and scientific qualifications are retained.

This package contains no original papers, publisher artwork, source workbooks,
author software or detailed research archives. It has no separate crystal
structure or EOS parameter library: those remain in core.

## What is not included?

The original Sakai Table S8 lineage resource and the Campbell Tables S4/S5
integration resource are excluded. They account for five material entries,
because the Campbell file is referenced by four materials. Separately named
partial lineage, polynomial and relative-pressure calculations are included;
they do not reproduce the complete original tables or author workflows. Their
approximation and missing-baseline limits remain in their declarations/notices.

Barker's two resources are supporting-study data, accessed with
`get_study_dataset()`. They do not define a crystal phase or an executable
Barker EOS, and the original unreproduced-model status remains explicit.

## Licenses and reuse

Some datasets have an explicit source or author license, such as CC BY 4.0.
Others are numerical transcriptions included after reviewing the particular
content. **Including a transcription is not a new license for its source.**
The original source notices remain, including NC/ND notices for 21 entries.
Contributor CC0 statements cover only the contributor's own work.

For an individual entry, read its declaration and source notice. `manifest.json`
lists exactly what is installed; `redistribution-review.json` explains the
reviewed inclusion/exclusion decisions. `admit` means selected for this package,
not “open license granted”. The numerical-transcription notice and the
source-specific notices retain the detailed reasoning and limits.

These decisions concern general redistribution of this optional package. Debian
or another distributor may assess them separately. Packaging the normal
`peritheos` library does not require packaging this collection.

Code remains MIT; data retain their individual recorded terms. Keep source
attribution, applicable license texts and descriptions of changes. The
[package comparison and current inventory](https://github.com/CPrescher/peritheos/blob/main/docs/research/license-inventory.md)
explain the two distributions and identify the preserved detailed reports.
