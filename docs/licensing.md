# Packages and licensing

Since **0.18.0**, Peritheos has two Python distributions. For Dioptas, start with
**`peritheos` alone**. The optional package supplies observation tables and other
numerical datasets for researchers who want those inputs.

## What each package contains

| Content | `peritheos` | `peritheos-datasets` |
| --- | --- | --- |
| Python/Rust library, EOS calculations and fitting your own data | Yes | Uses the core library |
| Material descriptions, crystal structures and EOS parameters | 234 material documents and 685 EOS records | No separate EOS or structure library |
| Structural CIF files and calculated reference diffraction patterns | 10 COD CIFs and 10 calculated patterns | No CIFs or separate pattern library |
| Dataset descriptions, citations and checksums | Yes; descriptions do not include observation values | Complete declarations for selected datasets |
| Experimental tables, numerical figure readouts and selected calculated datasets | No | 451 material entries and 2 supporting-study entries |
| CSV resources | None | 392 files, shared by some entries; 14 further tables are embedded in JSON |
| Original papers, publisher artwork, source workbooks and detailed research archives | No | No |
| License texts, source notices and provenance needed for the packaged content | Yes | Yes |

These are the contents of the **released 0.18.0 wheels and source distributions**.
The Git repository contains additional research evidence. A GitHub-generated
“Source code” ZIP/tarball is a repository snapshot, not the smaller PyPI source
distribution. The [package inventory](research/license-inventory.md) identifies
the release archives and their installed paths.

The core package still contains scientific data: structures, coefficients,
uncertainties and calculated diffraction lines. Here, “without datasets” means
**without the observation tables and other numerical datasets supplied by the
optional package**. The legacy directory name `peritheos/data/datasets/` also
holds the core CIF files and their notices; that name does not mean experimental
tables are installed.

A core dataset description can still record the absent table's source license,
including an NC/ND label. That metadata is not evidence that its numerical table
is shipped, or a license for the whole library. Use the package contents and
manifest to determine which resources are present.

## Installation

For calculations, structures, EOS parameters and fitting your own observations:

```bash
python -m pip install peritheos
```

To add the selected bundled datasets:

```bash
python -m pip install "peritheos[datasets]"
```

The `datasets` extra asks the installer for the matching version of
`peritheos-datasets`. It does not select another implementation of Peritheos.
Data access does not download missing files at runtime. See
[loading datasets](datasets.md) for the API and missing-data messages.

## Which licenses apply

| Component | What to retain when redistributing it |
| --- | --- |
| Peritheos code and contributor-written scripts | MIT license in `LICENSE`. A script's paper citation identifies its scientific source; it does not apply the paper's license to the script. |
| Rust dependencies and standard library | The supplied third-party license and attribution bundle, starting at `peritheos/data/third_party/rust/NOTICE.md`. |
| The ten core COD CIF files | Their CC0 source notice, original attribution and recorded changes. See `peritheos/data/datasets/cod-cif-sources.LICENSE.md`. |
| Other structures, EOS coefficients and calculated patterns | Their references and the supplied structure/parameter/pattern notices. These document the numerical content and its sources; they do not grant rights in an entire paper or database. |
| Optional datasets with an explicit source or author license | The license and attribution recorded for that dataset, including descriptions of changes. |
| Optional numerical transcriptions without an applicable source-data license | The original source notices and our documented reason for including the numerical content. Inclusion is not a new CC BY or CC0 license. |

A paper's license, permission for an author file and a contributor's dedication
can cover different things. A contributor CC0 statement covers only that
contributor's work. Source terms and our decision to include a numerical
transcription are therefore recorded separately. The
[package inventory](research/license-inventory.md) explains these cases in plain
language and links to the individual records.

The package metadata uses `LicenseRef-Peritheos-Distribution` because an archive
contains components with different terms. This is a pointer to the supplied
notices, not an additional license. Peritheos code remains MIT.

## What should a Debian maintainer check?

For **core Peritheos**, the intended packaging route is:

1. Use the official `peritheos` source distribution for 0.18.0, rather than the
   full Git repository snapshot. The observation tables and research archives
   have already been left out of that distribution.
2. Build and test Dioptas with the core package. Do not request its `datasets`
   extra or add `peritheos-datasets` as a dependency unless the application needs
   those bundled observations.
3. Retain the MIT, Rust, COD and other source notices for the components actually
   packaged, using Debian's normal copyright-file conventions. If a particular
   EOS or structure source needs further assessment, the
   [component index](data/license-inventory/eos-structure-review.csv) identifies
   its material, record and source.

We are not asking Debian to review the optional observation collection merely
to package the core library. We also do not have a confirmed list of core files
that Debian must remove. The detailed review records questions about source
collections, including the Holland parameter extracts; it does not establish
that individual coefficients are prohibited from reuse.

Packaging **`peritheos-datasets`** would be a separate task. It contains numerical
transcriptions admitted on a documented assessment as well as explicitly
licensed data. In particular, 21 entries retain NC/ND notices from their sources;
we have not relabeled those sources as freely licensed. A distributor considering
that package needs to decide whether it accepts the documented basis for those
specific numerical extractions. This is not needed for the core-only route.

## More detail

- [Current package inventory](research/license-inventory.md): exact counts,
  locations, exclusions and how to read the supplied records.
- [Preserved 0.17.3 audit and later source-tree notes](research/license-inventory-0.17.3.md):
  the earlier archive investigation, before the split. Its `release` and `local`
  rows do not describe the contents of the 0.18.0 packages.
- [EOS and structure review](research/eos-structure-compilation-rights.md) and
  [numerical dataset review](research/numerical-dataset-redistribution.md):
  the detailed source investigations, retained as evidence.
