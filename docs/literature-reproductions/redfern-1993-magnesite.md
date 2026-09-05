# Redfern et al. (1993): natural magnesite

Primary source: S. A. T. Redfern, B. J. Wood, and C. M. B. Henderson,
“Static compressibility of magnesite to 20 GPa: Implications for MgCO3 in the
lower mantle,” *Geophysical Research Letters* **20**, 2099–2102 (1993),
<https://doi.org/10.1029/93GL02507>.

The room-temperature powder-diffraction study used natural
`Mg0.991Fe0.008Mn0.001CO3` (Harwood collection specimen 2212). It reports two
Birch–Murnaghan sensitivity fits on the same conventional hexagonal cell
(`Z=6`, `V0=279.4(2) Å3`): `K0=142(9) GPa` with `K0'=4` fixed, and a full fit
with `K0=151(7) GPa`, `K0'=2.5`. Both are executable records, not separate
phases. The accessible primary record did not expose numerical P–V rows, so the
reproduction script uses analytical pressure checkpoints and does not invent
digitized points.

## LitCurate disposition

| Row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 52 | `litcurate_71122768f103ee77` | ACCEPT | Complete fixed-derivative Birch–Murnaghan sensitivity fit. |
| 53 | `litcurate_bd03302a0d214f4a` | ACCEPT | Complete unconstrained Birch–Murnaghan fit. |
| 54 | `litcurate_383afd63c513fd64` | REJECT | The source identifies `148/3.2` as a graphical finite-strain estimate, not an EOS fit. |
| 55 | `litcurate_b94b146ddc4a305e` | REJECT | MgO comparison value is citation-reported and incomplete in this source. |

Result: **2 production records, 2 rejected rows**.
