# Cohen and Lin (2014): static FeSiO3 Pv, PPv, and PPv-II

## Outcome

All three source-reported LitCurate candidates are accepted as distinct static
Vinet records. They describe three different FeSiO3 structures, not repeated
fits or experimental-run splits. The records are stored in phase-specific
material cards:

- `fesio3_bridgmanite_cohen_lin_2014_vinet_1`
- `fesio3_post_perovskite_cohen_lin_2014_vinet_1`
- `fesio3_post_perovskite_ii_cohen_lin_2014_vinet_1`

## Primary source

R. E. Cohen and Y. Lin, “Prediction of a potential high-pressure structure of
FeSiO3,” *Physical Review B* **90**, 140102(R) (2014),
doi:[10.1103/PhysRevB.90.140102](https://doi.org/10.1103/PhysRevB.90.140102).
The audit used the [authoritative UCL repository
copy](https://discovery.ucl.ac.uk/id/eprint/1461221/1/PhysRevB.90.140102.pdf),
SHA-256 `c4bd0265b83f2859c3bb2158f91dd20be064771a03cabbcb6bb164ad4ea5a711`.

The calculations are static PAW-PBE GGA+U (`U=6.0 eV`) Quantum ESPRESSO
calculations in 20-atom cells, using a 4x4x4 k-point mesh and 80 Ry cutoff. The
paper states that all three phases use high-spin antiferromagnetic states.

## EOS basis and reproduction

Table III explicitly labels the equation Vinet and prints `V0`, `K0`, and
`K0'` for Pv, PPv, and PPv-II. Its volumes are per FeSiO3. The material cards
use conventional 20-atom cells (`Z=4`), so only volume is multiplied by four;
bulk modulus and its derivatives are unchanged.

The underlying eight static energies correspond to -10, 0, 25, 50, 75, 100,
125, and 150 GPa and appear only in Figure 3. Digitizing that small energy plot
would add false precision. More usefully, Table III provides three independent
100 GPa checks for every fit: `V100`, `K100`, and `K100'`.

| Phase | V0 (A3/formula) | K0 (GPa) | K0' | Source V100 | Reproduced V100 | Source K100 | Reproduced K100 | Source K100' | Reproduced K100' |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Pv | 44.31 | 225 | 4.42 | 34.27 | 34.28075 | 597 | 596.8896 | 3.34 | 3.33886 |
| PPv | 44.90 | 189 | 4.73 | 33.98 | 33.94965 | 579 | 577.8706 | 3.47 | 3.46808 |
| PPv-II | 45.45 | 195 | 4.67 | 34.49 | 34.50403 | 580 | 580.5804 | 3.44 | 3.44258 |

The PPv `K100` difference is 1.13 GPa (0.20%) and is consistent with the
integer-rounded `K0` and `K100` printed in the source. The other checks round
directly to Table III. `scripts/reproduce_cohen_lin_2014_fesio3.py` evaluates
the Vinet equation independently and performs these checks deterministically.

## Structure and scope

Table I supplies complete 100 GPa Cmcm PPv and Cmmm PPv-II lattices and
fractional coordinates, which are stored without modification. It does not
tabulate the Pv structure, so no Pv lattice, space-group setting, or coordinates
are invented. Table I site multiplicities establish `Z=4` for both tabulated
cells.

The mathematical fit spans the eight static calculation pressures, but this is
not a phase-stability interval. In particular, the authors find that Cmmm
PPv-II distorts to C2/m at low pressure. The Cmmm structure card is therefore
anchored to the explicitly tabulated 100 GPa state.

## Candidate dispositions

| LitCurate identifier | Disposition | Reason |
|---|---|---|
| `litcurate_9b4780a09d66d440` | ACCEPT | Source Table III Pv Vinet; independent V100/K100/K100' parity. |
| `litcurate_ceaea2f0eec4299a` | ACCEPT | Distinct source Table III PPv Vinet; Cmcm structure tabulated. |
| `litcurate_7596b5eb1cb53baf` | ACCEPT | Distinct source Table III PPv-II Vinet; Cmmm structure tabulated. |

There are no citation-reported comparison rows under this DOI in the current
LitCurate ledger.
