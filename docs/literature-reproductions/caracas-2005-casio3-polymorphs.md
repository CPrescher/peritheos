# Caracas et al. (2005): static CaSiO3 perovskite tilt polymorphs

## Outcome

The primary article reports nine phase-specific BM3/BM4 pairs. All 18 remain
preserved in the source-audit fixture, but production retains only the four
fits for the ideal `Pm-3m` reference and the lowest-energy `I4/mcm` branch.
The other seven exploratory tilt branches are not exposed as standalone
materials because the source does not provide diffraction-ready structures.

R. Caracas, R. Wentzcovitch, G. David Price, and J. Brodholt, “CaSiO3
perovskite at lower mantle pressures,” *Geophysical Research Letters* **32**,
L06306 (2005),
[doi:10.1029/2004GL022144](https://doi.org/10.1029/2004GL022144).
The [UCL author copy](https://discovery.ucl.ac.uk/id/eprint/99227/1/2004GL022144.pdf)
and publisher HTML were checked directly.

## Source interpretation

Sections 2-4 describe static plane-wave pseudopotential DFT calculations up to
about 160 GPa.  Table 1 defines nine distinct octahedral-rotation structures;
Table 2 gives a complete third- and fourth-order Birch-Murnaghan fit for each.
The printed volume is explicitly `V0/Z`, in A3 per CaSiO3 molecule. Peritheos
converts it to the established conventional-cell basis on the existing
`Pm-3m` and `I4/mcm` material cards.

The source does not publish the E(V) grid, weights, coefficient covariance, or
coefficient uncertainties.  Its approximately 20 meV/molecule uncertainty is
an energy-accuracy estimate and is not copied into parameter-error fields.
No crystallographic lattice constants or optimized fractional coordinates are
invented for the seven nonproduction branches.

## Candidate-by-candidate disposition

LitCurate's generic structure labels are corrected from Table 1.  In
particular, I4/mcm and I4/mmm are CaSiO3 perovskite tilt variants, not
hollandite phases.

| LitCurate ID | Source structure/model | Disposition |
|---|---|---|
| `litcurate_192489e3e54304ee` | Pm-3m BM3 | ACCEPT; already present, retained. |
| `litcurate_d3cc1be1f5314d12` | I4/mcm BM3 | ACCEPT; structure corrected. |
| `litcurate_9ff9925ba424b41d` | Imma BM3 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_d173c5f08b283a85` | R-3c BM3 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_a86085849dfddfac` | P4/mbm BM3 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_30b66389486a3c14` | I4/mmm BM3 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_7812ee572cfb90b8` | Im-3 BM3 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_6a8edbf747ae3e36` | P42/nmc BM3 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_25ceab08a1621ad0` | Pnma BM3 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_4a8451fa3faa3bb0` | Pm-3m BM4 | ACCEPT. |
| `litcurate_61ad352d02ba2bc8` | I4/mcm BM4 | ACCEPT; structure corrected. |
| `litcurate_3b31d38ef03927d1` | Imma BM4 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_eb82dd40a95765b8` | R-3c BM4 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_15a3f8137aceab19` | P4/mbm BM4 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_cf5ba9b21a0f6468` | I4/mmm BM4 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_26e8348ba9113699` | Im-3 BM4 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_3cbd9436503036c1` | P42/nmc BM4 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_0ae9750e49c24da9` | Pnma BM4 | AUDIT ONLY; exploratory branch without a complete structure. |
| `litcurate_4bd92877f00bf9ff` | averaged experimental range | REJECT under this DOI; paragraph 4 summarizes cited experiments rather than fitting a new curve. |
| `litcurate_9ef919ec49605b7a` | averaged theoretical range | REJECT under this DOI; paragraph 4 summarizes cited theory and is not a Caracas fit. |

## Numerical reproduction

Table 2's BM3 coefficients independently reproduce the densities stated in
paragraph 16.  At zero and 130 GPa the Pm-3m curve gives 4.326 and 5.769
g/cm3, matching 4.32 and 5.77; I4/mcm gives 4.331 and 5.780 g/cm3, matching
4.33 and 5.78.  This also explains why BM3 is listed before the same-data BM4
alternative on each material card.

The deterministic check is
`scripts/reproduce_caracas_2005_casio3_polymorphs.py`; focused tests preserve
the complete 18-fit Table 2 coefficient matrix while confirming that only the
four selected `Pm-3m` and `I4/mcm` records instantiate through the production
catalog.
