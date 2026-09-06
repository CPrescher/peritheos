# Zhu et al. (2025): internally consistent Pt, Au, and MgO pressure standards

## Primary-source audit

- Citation: X. Zhu, Y. Ye, W. Sun, T. Katsura, and R. Caracas, “Internally consistent pressure-volume-temperature equations of state among platinum, gold, and MgO: Implications for phase transition boundaries in the mantle,” ESSOAr preprint version 1 (2025), DOI [`10.22541/essoar.176236186.65259830/v1`](https://doi.org/10.22541/essoar.176236186.65259830/v1).
- Primary preprint PDF: <https://d197for5662m48.cloudfront.net/documents/publicationstatus/288995/preprint_pdf/85b400a506980803732aafbd5816f7fb.pdf>
- Author data and scripts: Mendeley Data, version 3, DOI [`10.17632/6kxnhc2g73.3`](https://doi.org/10.17632/6kxnhc2g73.3), CC BY 4.0.
- Audited 2026-09-05. This is explicitly a non-peer-reviewed preprint; records preserve that status in their citation source.

Primary Equation 1 and Table 1 give three exact 300 K Vinet branches:

| Material | V0 (Å³ conventional cell) | K0 (GPa) | K0' |
|---|---:|---:|---:|
| Pt | 60.38 | 277.3 | 5.230(33) |
| Au | 67.85 | 167 | 5.897(22) |
| MgO | 74.71 | 160.3 | 4.182(19) |

The fixed/reference `V0` and `K0` values are marked in Table 1 as adopted from Ye et al. (2017), while the derivative and thermal terms belong to this internally consistent optimization. The complete high-temperature EOS adds a generalized Grüneisen law `gamma(V)=gamma0{1+a[(V/V0)^b-1]}` and a volume-dependent `T²` excess free-energy/pressure term. No existing single Peritheos thermal component represents both terms together. Consequently, only the exact, independently executable 300 K Vinet components are production records; the wider P–V–T branches are held rather than approximated.

## Exhaustive LitCurate disposition

| Candidate | Origin | Decision | Reason |
|---|---|---|---|
| `litcurate_b1a3d7a39cb0fae8` | source | **accepted after primary correction** | LitCurate correctly transcribed MgO `V0=74.71` and `K0=160.3`, but column-shifted Pt’s `K0'=5.230` into MgO. Primary Table 1 unambiguously gives MgO `K0'=4.182(19)`. The corrected MgO branch is stored. |
| `litcurate_fdde71245d253100` | citation | **rejected here** | `K0'=4.367` is a cited comparison value without a complete same-source parameterization; it belongs to the cited primary work. |

Primary-source expansion adds the co-optimized Pt and Au 300 K branches, which LitCurate omitted. Thus two ledger rows yield one corrected production record, and the complete same-paper audit yields three production records total. This is not padding: the three curves represent different materials within the paper’s explicitly coupled pressure-scale family.

## Reproduction

Run `uv run python scripts/reproduce_zhu_2025_pressure_standards.py`. The script verifies the exact Table 1 coefficient triples, evaluates each Vinet curve at three compressions, and performs inverse volume recovery. The Mendeley deposit contains the authors’ broader optimization/calculation scripts, but numerical parity for unsupported high-temperature terms is not claimed by these deliberately isothermal records.
