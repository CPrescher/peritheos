# Source-specific metal EOS batch (47 records)

This batch adds 47 metal equation-of-state records selected to favor widely
used pressure standards and experimentally constrained fits that were absent
from the catalog. The exact machine-readable
inventory is
`peritheos/data/datasets/experimental-metal-eos-batch-47.csv`, and the
idempotent import is implemented by
`scripts/build_experimental_metal_eos_batch_47.py`.

## Composition

| Source family | Records | Character |
|---|---:|---|
| Dewaele (2019), Table 1 | 30 | Static DAC X-ray fits under two published pressure scales |
| Dewaele et al. (2004), Table II | 8 | Static DAC X-ray fits under classical/revised ruby scales |
| Dorfman et al. (2012), Table 2 | 6 | Co-compression pressure-standard fits against MgO |
| Zhang et al. (2025), Table 2 | 3 | High-P/T hcp-Fe thermal fits to mixed experimental/theoretical data |
| **Total** | **47** | **44 experimental; 3 mixed-data thermal** |

## Primary sources

- Dewaele, *Equations of State of Simple Solids ... in the Mbar Range*,
  [doi:10.3390/min9110684](https://doi.org/10.3390/min9110684).
  Both Mao-1986 and Dorogokupets-2007 pressure-scale reductions are retained as
  correlated alternatives. Table 1 states that the unprinted 95% errors for
  the latter equal those printed for the former.
- Dewaele, Loubeyre, and Mezouar, *Equations of state of six metals above
  94 GPa*,
  [doi:10.1103/PhysRevB.70.094112](https://doi.org/10.1103/PhysRevB.70.094112).
- Dorfman et al., *Intercomparison of pressure standards ... to 2.5 Mbar*,
  [doi:10.1029/2012JB009292](https://doi.org/10.1029/2012JB009292).
- Zhang et al., *Equation of State Parameters of hcp-Fe Up to Super-Earth
  Interior Conditions*,
  [doi:10.3390/cryst15030221](https://doi.org/10.3390/cryst15030221).

## Scope and caveats

The 2019 hcp-Pb rows receive a dedicated phase card rather than being mixed
into ambient fcc Pb.

The three hcp-Fe thermal records use BM3 or Vinet reference curves coupled to
the Mie-Gruneisen-Debye model. They span the source's compiled 0--1374 GPa and
300--12000 K envelope. Because the source fit combines static and dynamic
experiments with theoretical points, these records are labeled `thermal_mixed`
rather than purely experimental.

## Verification

The generator asserts an exact 47-record selection before writing either the
CSV or material cards. The dedicated tests verify the source-family counts,
record identifiers, schema loading, reference-state pressure, Baonza scalar
and vector behavior, and nonzero thermal response of the hcp-Fe models:

```text
uv run --frozen python scripts/build_experimental_metal_eos_batch_47.py
uv run --frozen pytest -q tests/test_experimental_metal_eos_batch_47.py
```

After import and removal of benchmark-only or structurally incomplete records,
the complete bundled catalog contains 208 material documents and 488 EOS
records. Every record is included in the generated primary-source
audit and refit-feasibility ledger.
