# Peritheos: core package contents and licenses

This is the core `peritheos` distribution, introduced with the package split in
0.18.0. Install it with `python -m pip install peritheos`. It supplies the library,
structures and EOS parameters needed for calculations and applications such as
Dioptas. The optional observation package is not a dependency of this install.

## What is included?

- The Python/Rust calculation and fitting library.
- 234 material documents with 685 EOS records, including coefficients, units,
  uncertainties, reference states, source citations and scientific qualifications.
- Crystal structures, 10 COD CIF files and 10 calculated reference diffraction
  patterns. Structural and proxy limitations remain in their material records.
- Dataset descriptions, references and checksums, so users can discover source
  observations without having their values installed.
- The applicable license texts, source notices and review records.

**There are no experimental observation tables or CSV resources in this
package.** Original author workbooks, publisher artwork and detailed research
archives are also excluded. The wheel and the official core source distribution
use the same scientific-content boundary. A full Git repository snapshot contains
additional research evidence and is not the same archive.

The directory `peritheos/data/datasets/` in core holds the CIFs and their notices;
it does not contain the optional observation tables. The material documents also
contain scientific numbers, such as EOS coefficients and calculated peak lines.
“Without datasets” means without the optional numerical tables, not without any
scientific data.

## Which licenses apply?

| Content | License or source information |
| --- | --- |
| Peritheos code | MIT, in the supplied `LICENSE`. |
| Rust dependencies and standard library | Full license texts and attributions in [the Rust notice bundle](third_party/rust/NOTICE.md). |
| Ten COD CIFs | CC0 source declarations, with original attribution and recorded changes in [the COD notice](datasets/cod-cif-sources.LICENSE.md). |
| EOS parameters and other structural values | Scientific references and numerical-content review in [the compilation notice](datasets/eos-structure-compilations.LICENSE.md) and `eos-structure-rights-review.json`. No blanket license for the cited papers or source collections is asserted. |
| Calculated reference patterns | Inputs, contributor rights and model limits in [the pattern notice](datasets/calculated-patterns.LICENSE.md) and `calculated-pattern-inputs.json`. Six inputs are COD CIFs, one model is independently authored and three use attributed numerical structural parameters. |
| Material definitions migrated from Dioptas | Original MIT attribution in [the Dioptas notice](datasets/dioptas-materials.LICENSE.md). |

If rebuilding with another Rust version, check its standard-library notices too.
The supplied bundle records the compiler used for the official wheels.

Source license statements and our assessment of numerical reuse are different
records. The detailed review retains questions about the Holland parameter
extracts; it does not establish that those numerical coefficients must be
removed. An absent source-data license is not by itself a finding that a
numerical parameter is prohibited from reuse.

## Adding observation data

`python -m pip install "peritheos[datasets]"` adds the matching
`peritheos-datasets` distribution. In 0.18.0 it supplies 451 material entries and
two supporting-study entries, using 392 CSV files and 14 inline tables. It has
its own source notices and selection review. It does not include every dataset
described by the core catalog, and access does not download missing data at
runtime.

For the package comparison, source archive names and packaging guidance, see
[packages and licensing](https://github.com/CPrescher/peritheos/blob/main/docs/licensing.md).
The [current inventory](https://github.com/CPrescher/peritheos/blob/main/docs/research/license-inventory.md)
links the detailed evidence without treating the old research checkout as the
contents of either current package.
