# Caracas et al. (2005): static CaSiO3 perovskite tilt polymorphs

## Outcome

The primary article reports nine phase-specific BM3/BM4 pairs.  All 18 are
complete source-generated static EOS parameterizations.  One cubic BM3 record
was already in the catalog; this audit adds the cubic BM4 and the remaining
eight pairs, for **17 net-new production records**.

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
The printed volume is explicitly `V0/Z`, in A3 per CaSiO3 molecule.  Peritheos
therefore preserves a one-formula-unit volume basis on the indexing-only cards,
or converts it to the established conventional-cell basis where an existing
Pm-3m or I4/mcm card supplies that basis.

The source does not publish the E(V) grid, weights, coefficient covariance, or
coefficient uncertainties.  Its approximately 20 meV/molecule uncertainty is
an energy-accuracy estimate and is not copied into parameter-error fields.
No crystallographic lattice constants or optimized fractional coordinates are
invented for the new phase cards.

## Candidate-by-candidate disposition

LitCurate's generic structure labels are corrected from Table 1.  In
particular, I4/mcm and I4/mmm are CaSiO3 perovskite tilt variants, not
hollandite phases.

| LitCurate ID | Source structure/model | Disposition |
|---|---|---|
| `litcurate_192489e3e54304ee` | Pm-3m BM3 | ACCEPT; already present, retained. |
| `litcurate_d3cc1be1f5314d12` | I4/mcm BM3 | ACCEPT; structure corrected. |
| `litcurate_9ff9925ba424b41d` | Imma BM3 | ACCEPT. |
| `litcurate_d173c5f08b283a85` | R-3c BM3 | ACCEPT. |
| `litcurate_a86085849dfddfac` | P4/mbm BM3 | ACCEPT. |
| `litcurate_30b66389486a3c14` | I4/mmm BM3 | ACCEPT; structure corrected. |
| `litcurate_7812ee572cfb90b8` | Im-3 BM3 | ACCEPT. |
| `litcurate_6a8edbf747ae3e36` | P42/nmc BM3 | ACCEPT. |
| `litcurate_25ceab08a1621ad0` | Pnma BM3 | ACCEPT. |
| `litcurate_4a8451fa3faa3bb0` | Pm-3m BM4 | ACCEPT. |
| `litcurate_61ad352d02ba2bc8` | I4/mcm BM4 | ACCEPT; structure corrected. |
| `litcurate_3b31d38ef03927d1` | Imma BM4 | ACCEPT. |
| `litcurate_eb82dd40a95765b8` | R-3c BM4 | ACCEPT. |
| `litcurate_15a3f8137aceab19` | P4/mbm BM4 | ACCEPT. |
| `litcurate_cf5ba9b21a0f6468` | I4/mmm BM4 | ACCEPT; structure corrected. |
| `litcurate_26e8348ba9113699` | Im-3 BM4 | ACCEPT. |
| `litcurate_3cbd9436503036c1` | P42/nmc BM4 | ACCEPT. |
| `litcurate_0ae9750e49c24da9` | Pnma BM4 | ACCEPT. |
| `litcurate_4bd92877f00bf9ff` | averaged experimental range | REJECT under this DOI; paragraph 4 summarizes cited experiments rather than fitting a new curve. |
| `litcurate_9ef919ec49605b7a` | averaged theoretical range | REJECT under this DOI; paragraph 4 summarizes cited theory and is not a Caracas fit. |

## Numerical reproduction

Table 2's BM3 coefficients independently reproduce the densities stated in
paragraph 16.  At zero and 130 GPa the Pm-3m curve gives 4.326 and 5.769
g/cm3, matching 4.32 and 5.77; I4/mcm gives 4.331 and 5.780 g/cm3, matching
4.33 and 5.78.  This also explains why BM3 is listed before the same-data BM4
alternative on each material card.

The deterministic check is
`scripts/reproduce_caracas_2005_casio3_polymorphs.py`; focused tests confirm
all 18 records instantiate, give `P(V0)=0` and `K(V0)=K0`, preserve the exact
Table 2 coefficient matrix, and keep the new cards free of inferred
coordinates.
