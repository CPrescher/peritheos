# Hama and Suito (1996): static-lattice MgO Vinet comparison

Primary source: J. Hama and K. Suito, “The search for a universal equation of state correct up to very high pressures,” *Journal of Physics: Condensed Matter* **8**, 67–81 (1996), <https://doi.org/10.1088/0953-8984/8/1/008>.

Sections 1–2, Table 1, Figure 2b, and the abstract were checked directly. Table 1 reports the static-lattice QSM MgO reference values `V0 = 123.747 bohr3` per formula unit, `K0 = 157 GPa`, and `K0' = 4.37`; Section 2 uses those values to compare the Vinet curve with the theoretical pressure–volume relation. The stored `V0 = 73.34965396218291 A3` is the exact conventional-cell (`Z=4`) conversion using the CODATA bohr. The paper reports agreement for diatomic solids to approximately 1 TPa or `V/V0 = 0.35`.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 154 | `litcurate_c6a0c7a005a7ff43` | ACCEPT | Source-owned, phase-specific static MgO Vinet application with a complete reference triple. The ledger's generic “APW” method label is corrected to the source's Table 1 QSM attribution. |
| 155 | `litcurate_5d29aec553e70970` | REJECT | MgSiO3-perovskite values are explicitly taken from reference [29], not calculated or fitted in this paper. |

The reproduction script converts the printed volume basis, evaluates five compression anchors through the exact Vinet equation, and the focused test checks the bundled implementation and inverse pressure–volume round trip.

Result: **1 production record and 1 rejected citation row**.
