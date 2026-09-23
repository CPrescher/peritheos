# Tetragonal CaSiO3: source-specific space-group assignments

Review date: 2026-09-22. The source-specific distinction is between a published
indexing/refinement model and a uniquely established structure. The generic
label "tetragonal CaSiO3" does not resolve that distinction.

## Primary evidence

- **Shim, Jeanloz, and Duffy (2002)**, *Geophysical Research Letters* 29,
  2166, [DOI 10.1029/2002GL016148](https://doi.org/10.1029/2002GL016148).
  Figures 1-2 use P4/mmm for indexing and the one-formula-unit volume basis.
  Paragraph 15 reports almost equal Rietveld residuals: 6.9% for P4/mmm and
  7.1% for I4/mmm. Paragraph 19 explicitly requires further work to confirm
  the assignment and investigate lower symmetry. The shortened pseudocubic
  c-axis differs from the later octahedral-rotation interpretation.
  [Author-hosted article](https://duffy.princeton.edu/sites/g/files/toruqf616/files/shim_et_al-2002-grl.pdf).
- **Chen et al. (2018)**, *American Mineralogist* 103, 462-468,
  [DOI 10.2138/am-2018-6087](https://doi.org/10.2138/am-2018-6087).
  Figure 3 and p. 465 favor the octahedral-rotation model in Ne: the
  background-subtracted residual is 2.2% for I4/mcm versus 3.0% for P4/mmm.
  However, the discussion on p. 465 explicitly leaves the exact space group
  uncertain among competing rotation structures. I4/mcm is the chosen
  refinement and EOS model, not a claimed unique solution. The authors
  discuss stress conditions as a possible explanation for earlier differences.
- **Sun et al. (2022)**, *American Mineralogist* 107, 110-115,
  [DOI 10.2138/am-2021-7913](https://doi.org/10.2138/am-2021-7913).
  Page 113 identifies an additional 211 reflection at 148-199 GPa that is
  explained by I4/mcm among the four tested tetragonal candidates. Reinspection
  of the authors' 2016 patterns also finds it at 24-124 GPa. The abstract
  retains a qualification that other possibilities may exist. Figure 3
  distinguishes the older P4/mmm-indexed data from the I4/mcm data.
  [Author-hosted article](https://www.jsg.utexas.edu/lin/files/SunCaPvAM2022.pdf).

## Catalog treatment

`casio3_perovskite_tetragonal` contains the Shim (2002) EOS. Its
`source.indexing_model` now explicitly records **P4/mmm**, Z=1, and the
non-unique nature of that model. The top-level space-group field remains
unassigned. Copying I4/mcm onto this card would imply an unsupported
reassignment of the original samples and obscure the source's cell convention.

`ca_perovskite_tetragonal` already carries **I4/mcm (No. 140)**, Z=4,
for the later adopted structure model. This is the supported modern model
in these reviewed sources, with the qualifications above. Its conventional
cell volume is four times its formula-unit volume. A future consolidation
would need explicit volume normalization and preservation of the different
structural interpretations, not merely a shared space-group label.

No EOS coefficients, volume bases, or existing space-group assignments were
changed by this review.

## Retrieval

Zotero items: Shim `UJJ6ED5H`, Chen `RS252GDC` (PDF attachment `6ZWSNJHV`),
and Sun `9ZCZ3PSB` (linked article `VYYH5SJ6`). Chen's original PDF was read
from the local library. Shim's item has no attachment; the original paper
was retrieved from Duffy's Princeton site. Sun's linked publisher PDF was
retrieved from Lin's University of Texas site. The relevant printed pages
were visually checked. No Zotero records were modified.
