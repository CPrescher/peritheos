# Funamori et al. (1998) MgAl2O4 high-pressure phase audit

## Scope and primary source

This audit covers every LitCurate candidate attached to Funamori et al.,
“High-pressure transformations in MgAl2O4,” *Journal of Geophysical Research:
Solid Earth* **103**, 20813–20818 (1998),
[doi:10.1029/98JB01575](https://doi.org/10.1029/98JB01575).

The primary evidence is the free-access publisher PDF.  The authors compressed
99% pure polycrystalline spinel mixed with 5 wt% Pt, laser heated each charge
at 2000–3000 K, and measured diffraction at room temperature.  Pressures after
heating were obtained independently from Mao et al. (1978) ruby fluorescence
and the Holmes et al. (1989) Pt EOS.  The source uses their average as the
pressure attached to each in-situ cell.

A catalog scan found no pre-existing EOS record with this DOI.

## Phase identity and fit independence

The two source parameterizations are not arbitrary run-specific fits:

- Run 2 contains a CaFe2O4-type orthorhombic MgAl2O4 phase.  The source gives
  its recovered ambient cell and a separate in-situ cell at 33.4 GPa.
- Runs 3 and 4 establish a higher-pressure CaTi2O4-type MgAl2O4 phase.  TEM on
  run 3 confirms the MgAl2O4 composition and narrows its symmetry; the source
  assigns `Cmcm`.  Run 4 supplies recovered ambient and 64.3 GPa cells.

Section 3.3 applies a second-order Birch finite-strain EOS separately to the
two structures.  Thus `K0_prime=4` is assumed by equation order.  Each fit has
only an ambient/high-pressure pair, and the authors explicitly say more points
would be needed to estimate `K0_prime`.  These are sparse but independently
reported phase EOS parameterizations.

No target atomic coordinates were refined.  The authors calculated patterns
using CaFe2O4 and CaTi2O4 prototype coordinates.  The material cards therefore
store measured target lattices and structure assignments but do not invent
Mg/Al/O atom sites.

## Data reconstruction and numerical check

The paper does not print a separate P-V table.  All four source states are
recoverable without figure digitization: Sections 3.1–3.2 and Tables 1–2 print
the ambient and high-pressure axes, ambient volumes, and both high-pressure
calibrations.  High-pressure volumes in the bundled CSVs are direct products
of the printed axes; their uncertainties are propagated from the printed axis
standard deviations.

| Phase | Source BM2 `(V0 A3, K0 GPa)` | High-P state `(P GPa, V A3)` | Fixed-V0 BM2 refit K0 | Published-curve residual |
|---|---|---|---:|---:|
| CaFe2O4 type | `(240.3(2), 211(6))` | `(33.4, 212.406277266)` | 211.436871 GPa | −0.069011 GPa |
| CaTi2O4 type | `(240.3(4), 206(3))` | `(64.3, 195.540421098)` | 206.401375 GPa | −0.125040 GPa |

Both refitted moduli are comfortably within the published uncertainties.  The
small residuals arise from multiplying rounded lattice axes rather than from a
phase or EOS mismatch.

## Disposition of every same-DOI candidate

| LitCurate identifier | Candidate | Disposition |
|---|---|---|
| `litcurate_ca9a93f28c82d5bc` | Run-2 `V0=240.3`, `K0=211`, `K0'=4` | **ACCEPT** as the CaFe2O4-type MgAl2O4 BM2; independent source phase and two source states verified. |
| `litcurate_45e8795108756c17` | Runs-3/4 `V0=240.3`, `K0=206`, `K0'=4` | **ACCEPT** as the CaTi2O4-type `Cmcm` MgAl2O4 BM2; independent source phase and two source states verified. |
| `litcurate_cc1efc12c6e1c6e0` | Yutani et al. (1997) CaFe2O4-type `K0=241`, `K0'=4` | **REJECT as a Funamori-source record**; it is a citation-reported comparison and lacks `V0` in LitCurate.  It remains eligible for audit under the Yutani primary paper. |

Final and net-new production count for DOI `10.1029/98JB01575`: **2 records**.
No source-reported candidate is held.

## Reproduction assets

- `peritheos/data/materials/mgal2o4_cafe2o4.eosmat`
- `peritheos/data/materials/mgal2o4_cati2o4.eosmat`
- `peritheos/data/datasets/mgal2o4-cafe2o4-funamori-1998-text-pv.csv`
- `peritheos/data/datasets/mgal2o4-cati2o4-funamori-1998-text-pv.csv`
- `scripts/reproduce_funamori_1998_mgal2o4.py`
- `tests/test_funamori_1998_mgal2o4.py`

Run `python scripts/reproduce_funamori_1998_mgal2o4.py` for the independent
lattice-product, curve, and fixed-V0 diagnostics.
