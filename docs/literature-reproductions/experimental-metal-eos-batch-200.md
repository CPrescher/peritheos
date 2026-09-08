# Experimental metal EOS batch (200 records)

This batch adds 200 metal equation-of-state records selected to favor widely
used pressure standards and experimentally constrained fits that were absent
from the catalog. It contains no DFT-only records. The exact machine-readable
inventory is
`peritheos/data/datasets/experimental-metal-eos-batch-200.csv`, and the
idempotent import is implemented by
`scripts/build_experimental_metal_eos_batch_200.py`.

## Composition

| Source family | Records | Character |
|---|---:|---|
| Sun et al. (2010), Table 3 | 111 | Fits to compiled experimental compression data |
| Delta-project experimental reference table | 42 | Experimental reference coefficients corrected for zero-point and thermal effects |
| Dewaele (2019), Table 1 | 30 | Static DAC X-ray fits under two published pressure scales |
| Dewaele et al. (2004), Table II | 8 | Static DAC X-ray fits under classical/revised ruby scales |
| Dorfman et al. (2012), Table 2 | 6 | Co-compression pressure-standard fits against MgO |
| Zhang et al. (2025), Table 2 | 3 | High-P/T hcp-Fe thermal fits to mixed experimental/theoretical data |
| **Total** | **200** | **197 experimental or experimental-reference; 3 mixed-data thermal** |

The 111 Sun records comprise Vinet, MRS3, SMS3, and SMS4 fits for 26 metals,
plus Baonza pseudospinodal fits for seven high-use pressure-standard metals:
Cu, Mo, W, Ag, Pt, Ta, and Au. Supporting Baonza EOS evaluation was added to
both the Python and Rust runtimes.

## Primary sources

- Sun, Wu, Guo, and Cai, *Two Universal Equations of State for Solids*,
  [doi:10.1515/zna-2010-1-202](https://doi.org/10.1515/zna-2010-1-202).
  Table 3 supplies each low-pressure envelope and model coefficients; Table 2
  supplies the experimental molar reference volume.
- Lejaeghere et al., *Error estimates for solid-state density-functional
  theory predictions*,
  [doi:10.1080/10408436.2013.772503](https://doi.org/10.1080/10408436.2013.772503),
  and the CC BY 4.0
  [Delta-project archive](https://archive.materialscloud.org/record/2023.133).
  The archive's `history/history/exp.txt` is explicitly the experimental table
  with zero-point and thermal corrections.
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

The 42 Delta reference rows publish coefficients but not the experimental
compression range. Their validity interval is therefore deliberately stored
as the reference state only; the catalog does not invent an extrapolation
range. The Sun compilation does not resolve every original phase or pressure
calibration, so those fits remain on the existing phase-unresolved legacy
cards. The 2019 hcp-Pb rows receive a dedicated phase card rather than being
mixed into ambient fcc Pb.

The three hcp-Fe thermal records use BM3 or Vinet reference curves coupled to
the Mie-Gruneisen-Debye model. They span the source's compiled 0--1374 GPa and
300--12000 K envelope. Because the source fit combines static and dynamic
experiments with theoretical points, these records are labeled `thermal_mixed`
rather than purely experimental.

## Verification

The generator asserts an exact 200-record selection before writing either the
CSV or material cards. The dedicated tests verify the source-family counts,
record identifiers, schema loading, reference-state pressure, Baonza scalar
and vector behavior, and nonzero thermal response of the hcp-Fe models:

```text
uv run --frozen python scripts/build_experimental_metal_eos_batch_200.py
uv run --frozen pytest -q tests/test_experimental_metal_eos_batch_200.py
```

After import, the complete bundled catalog contains 286 material documents and
814 EOS records. Every record is included in the generated primary-source
audit and refit-feasibility ledger.
