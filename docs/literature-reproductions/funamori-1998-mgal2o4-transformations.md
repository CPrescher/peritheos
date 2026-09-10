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

A catalog scan found no pre-existing EOS record with this DOI before these two
records were added.

## Source construction and generic-refit classification

The two source parameterizations are phase-specific, but neither is a
regression fit:

- Run 2 contains a CaFe2O4-type orthorhombic MgAl2O4 phase.  The source gives
  its recovered ambient cell and a separate in-situ cell at 33.4 GPa.
- Runs 3 and 4 establish a higher-pressure CaTi2O4-type MgAl2O4 phase.  TEM on
  run 3 confirms the MgAl2O4 composition and narrows its symmetry; the source
  assigns `Cmcm`.  Run 4 supplies recovered ambient and 64.3 GPa cells.  A
  CaTi2O4-type phase was also identified in run 2, but its peaks were too weak
  and overlapped to determine a reliable unit-cell volume.

Section 3.3 states the exact one-point construction.  For each polymorph the
authors hold the separately measured recovered `V0` fixed, assume
`K0_prime=4`, and calculate `K0` from the only usable compressed state as the
normalized pressure

`F = P / [3 f (1 + 2 f)^(5/2)]`, where
`f = [(V/V0)^(-2/3) - 1] / 2`.

For a second-order Birch EOS, `F=K0`.  There is therefore one informative
finite-pressure observation for one remaining unknown and zero independent
residual degrees of freedom.  The authors themselves say that more points
would provide a more robust `K0` and allow `K0_prime`, which they assumed to be
4, to be estimated.  The generic ledger consequently remains
`not_refittable`: there is no statistically independent fit to rerun.

No target atomic coordinates were refined.  The authors calculated patterns
using CaFe2O4 and CaTi2O4 prototype coordinates.  The material cards therefore
store measured target lattices and structure assignments but do not invent
Mg/Al/O atom sites.

## Search for additional observations

The publisher article explicitly reports that there are no supplementary
materials.  DOI metadata has no related-data object, OpenAlex exposes no
repository full text, and targeted searches of author/repository records,
general research-data repositories, and theses found no source-associated
compression series.  The only author-uploaded copy found is the article
itself.  Later papers do publish multi-point measurements, but those are new
experiments and cannot be substituted as Funamori et al. fit inputs.

The article does not print a separate P-V table because the authors did not use
a longer compression series.  The complete source construction is nevertheless
recoverable without figure digitization: Sections 3.1–3.2 and Tables 1–2 print
both recovered
ambient cells, both usable compressed cells, and ruby- and Pt-derived
post-heating pressures.  High-pressure volumes in the bundled CSVs are direct
products of the printed axes; their uncertainties are propagated from the
printed axis standard deviations.

## Analytical reproduction

| Phase | Source BM2 `(V0 A3, K0 GPa)` | Compressed state `(P GPa, V A3)` | Fixed-V0 endpoint `K0` | Ruby/Pt `K0` half-range | Published-curve residual |
|---|---|---|---:|---:|---:|
| CaFe2O4 type | `(240.3(2), 211(6))` | `(33.4, 212.406277266)` | 211.436871 GPa | 6.013923 GPa | −0.069011 GPa |
| CaTi2O4 type | `(240.3(4), 206(3))` | `(64.3, 195.540421098)` | 206.401375 GPa | 2.888977 GPa | −0.125040 GPa |

The endpoint calculation reproduces the published moduli.  Applying it to the
two separately printed calibrant pressures reproduces the reported `K0`
uncertainties after rounding: about 6 and 3 GPa.  The small curve residuals
arise from multiplying rounded lattice axes.  These are transcription and
equation-parity checks, not fit validation.

## Production decision and external context

Both records remain `record_kind: published`, are defaults for their distinct
polymorph material cards, and remain executable production records.  This is
justified because the source equation, coefficients, fixed-parameter protocol,
complete construction observations, cell basis, and provenance are all known.
Their scientific-validation blocks now label them `not_refittable` and warn
that they are historical endpoint curves, not multi-point EOS fits.

Later independent experiments support the scale of the CaFe2O4-type modulus
but also show why users should prefer multi-point results: Irifune et al.
(2002), [doi:10.1007/s00269-002-0275-1](https://doi.org/10.1007/s00269-002-0275-1),
obtained `K0=213(3) GPa` from three compressed observations, while Sueda et al.
(2009), [doi:10.1016/j.pepi.2008.07.046](https://doi.org/10.1016/j.pepi.2008.07.046),
explicitly described the Funamori value as based on one 33.4 GPa cell and an
ambient cell and obtained a room-temperature multi-point result near
`K0=205(6) GPa`.  For the CaTi2O4 type, Ono et al. (2006),
[doi:10.1007/s00269-006-0068-z](https://doi.org/10.1007/s00269-006-0068-z),
published a separate 0–91.2 GPa table and a BM2 fit of `V0=238.9(9) A3`,
`K0=219(6) GPa`.  None of those later rows belongs in the Funamori datasets.

## Disposition of every same-DOI candidate

| LitCurate identifier | Candidate | Disposition |
|---|---|---|
| `litcurate_ca9a93f28c82d5bc` | Run-2 `V0=240.3`, `K0=211`, `K0'=4` | **ACCEPT** as the source-reported CaFe2O4-type MgAl2O4 historical endpoint BM2; complete construction evidence verified, but not independently refittable. |
| `litcurate_45e8795108756c17` | Runs-3/4 `V0=240.3`, `K0=206`, `K0'=4` | **ACCEPT** as the source-reported CaTi2O4-type `Cmcm` historical endpoint historical endpoint BM2; complete construction evidence verified, but not independently refittable. |
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
lattice-product, fixed-V0 endpoint, pressure-bracket, and curve diagnostics.
