# Kawai and Tsuchiya (2012) NaMg2Al5SiO12 NAL and CF audit

## Scope and source integrity

This audit covers every LitCurate candidate attached to Kawai and Tsuchiya,
“Phase stability and elastic properties of the NAL and CF phases in the
NaMg2Al5SiO12 system from first principles,” *American Mineralogist* **97**,
305–314 (2012), [doi:10.2138/am.2012.3915](https://doi.org/10.2138/am.2012.3915).

The primary source was the final typeset article linked by the official MSA
2012 table of contents.  The current publisher endpoint is protected by a web
application firewall, so the byte-identical MSA URL was retrieved from the
Internet Archive's 2012-07-21 capture.  Its SHA-256 is
`8273dc07c0ba3dca62724d490e035e945ee6c4123d7c6319f140afc44c7ec089`.
The official table of contents has no active data-deposit link for this paper;
its deposit markup is commented out.  No official supplement or numerical
DFT table was found.

A pre-implementation scan of all canonical `.eosmat` records found no existing
record with this DOI.  The two records introduced by this audit are therefore
not duplicates.

## Source-reported phases and calculations

The paper studies one exact bulk composition, NaMg2Al5SiO12, in two structures:

- NAL is hexagonal `P63/m`, with the general formula `XY2Z6O12`.  Eight
  distinct 1×1×3, 63-atom ordered supercells model Na/vacancy and Al/Si
  configurations.  The volume curve spans 0–50 GPa.
- CF is calcium-ferrite type, orthorhombic `Pbnm`, with the general formula
  `XY2O4`.  It is not an orthorhombic perovskite.  Twelve 1×1×3, 84-atom
  ordered supercells model the cation configurations.  The volume curve spans
  0–150 GPa.

The calculations used LDA in PWSCF/Quantum ESPRESSO, a 50 Ry plane-wave
cutoff, a 2×2×2 supercell k-point grid, and fully relaxed variable cells at
constant pressure until residual forces were below `5.0e-5 Ry/a.u.`.  These
are static 0 K calculations.  There is no experimental pressure gauge or
pressure calibration.

Table 1 reports standard third-order Birch–Murnaghan coefficients.  It gives
no equation, weighting, covariance, confidence level, or parameter errors.
The standard Eulerian BM3 convention is independently confirmed by the plotted
curves.  No coefficient is marked fixed for either current-study row.

## Volume basis

Table 1 labels `V0` only as cm³/mol, even though its current-study row labels
use the 12-oxygen composition.  Figure 2c–d volumes are exactly three times the
Table 1 scale, showing that the table silently uses the conventional
four-oxygen mineral formula normalization.

Using `1 cm3/mol = 0.602214076 A3` per formula unit:

| Phase | Printed O4 molar V0 | Conventional-cell contents | Peritheos V0 |
|---|---:|---|---:|
| NAL | 35.8 cm³/mol | one NaMg2Al5SiO12 (three O4 units) | 178.34189581447112 Å³ |
| CF | 35.2 cm³/mol | four O4 units = 4/3 NaMg2Al5SiO12 | 233.80390065807762 Å³ |

The non-integer `formula_units_per_cell=4/3` in the CF material is intentional:
the human-readable material formula is the exact 12-O study composition, while
the primitive `Pbnm` conventional cell contains four O4 units.  The 84-atom
1×1×3 supercell independently verifies three 28-atom conventional cells.

## Vector digitization and numerical reproduction

The official PDF page was converted to SVG with Poppler.  Marker centers were
read from the vector paths, not estimated from a raster.  The visible
current-study points are six red circles for NAL and eight dark-blue triangles
for CF.  All cited experimental comparison markers were excluded.

The CSV files retain the PDF point coordinates, plotted 12-O molar volumes,
Table 1 O4-normalized volumes, conventional-cell volumes, explicit curation
uncertainties, and fit-selection flags.  The 0.05 GPa and 0.06 cm³/mol
digitization bounds are conservative curation estimates covering frame/stroke
placement; they are not author-reported uncertainties.

The rounded published curves and an unweighted volume-residual diagnostic
refit give:

| Phase | Points | Published-curve volume RMSE | Maximum | Diagnostic refit `(V0 A3, K0 GPa, K0')` | Refit RMSE |
|---|---:|---:|---:|---|---:|
| NAL | 6 | 0.100424430486 Å³ | 0.112884772937 Å³ | `(178.233853876982, 217.562095743266, 4.088521901863)` | 0.003452426412 Å³ |
| CF | 8 | 0.151787368860 Å³ | 0.198344317864 Å³ | `(233.944765268382, 213.205584668798, 4.132190010213)` | 0.030775628572 Å³ |

These refits verify equation order, cell-volume basis, and the rounded Table 1
coefficients.  They are not replacements for the author fits: the source does
not disclose whether the published fit used one configuration, all
configurations, or configuration-averaged volumes, and the plotted coordinates
are rounded presentation data rather than the original DFT grid.

## Crystallographic proxies

The paper illustrates ordered supercells but publishes no machine-readable or
tabulated target fractional coordinates or occupancies.  None are inferred.

The NAL material carries Miura et al.'s complete `P63/m` CaMg2Al6O12 Rietveld
model solely as a topology proxy.  Every proxy site label says so, and its
elements and occupancies are preserved verbatim.  The primary proxy PDF
(American Mineralogist 85, 1799–1803; doi:10.2138/am-2000-11-1223) has SHA-256
`023aa07442ecd7b9ec42f37a59cef23196568eb4b7d36ec8d267c06cbf26f330`.

The CF material carries the complete deposited ideal NaAlSiO4 `Pbnm` model
from Qin et al. (2023), official MSA deposit AM-23-128432, solely as a topology
proxy.  It does not assert Mg substitution or an Al/Si distribution for the
target.  The deposit and CIF hashes are recorded in the material file.

## Source inconsistencies and limitations

- Table 1 prints a superscript `b` after the CF `B0=213.2` value, but provides
  no `b` footnote.  The number is clear; the superscript is left unresolved.
- The Figure 4 caption writes the two layer formulas with O12, while the prose
  gives charge-balanced, 28-atom O16 formulas.  The computation description
  and 84-atom supercell require O16, so the caption strings are documented as
  typographical errors.
- The source says pressure–volume differences between configurations are
  small but does not define the configuration aggregation used for Table 1.
- LDA underestimates volumes and the calculation has no thermal contribution;
  the authors caution that these account for much of the difference from
  experimental comparison data.

## Disposition of all nine same-DOI candidates

| LitCurate identifier | Candidate | Disposition | Reason |
|---|---|---|---|
| `litcurate_c42ee5edc21ddfe5` | This-study NaMg2Al5SiO12 NAL: `35.8`, `217.7`, `4.08` | **ACCEPT — production record** | Source-reported static-LDA BM3; phase, composition, cell basis, and Figure 2c curve verified. |
| `litcurate_3b528f9e781f555b` | Sanehira et al. experimental NAL: `36.932`, `189`, `4 fixed` | **REJECT as a Kawai-source record** | Table 1 is only a citation comparison.  Audit under Sanehira et al.'s primary DOI before any separate record. |
| `litcurate_d289c4c72b25a892` | Shinmei et al. experimental NAL: `36.69`, `176(2)`, `4.9(3)` | **REJECT as a Kawai-source record** | Citation-reported comparison, not a result of this paper. |
| `litcurate_b8fa7310ff9994db` | Vanpeteghem et al. experimental NAL: `37.05`, `214(2)`, `3.0(1)` | **REJECT as a Kawai-source record** | Citation-reported comparison, not a result of this paper. |
| `litcurate_f6c0b448f52c70d9` | Guignot and Andrault experimental NAL: `37.7`, `184(16)`, `4 fixed` | **REJECT as a Kawai-source record** | Citation-reported comparison, not a result of this paper. |
| `litcurate_527df255bad7681e` | This-study NaMg2Al5SiO12 CF: `35.2`, `213.2`, `4.12` | **ACCEPT — production record** | Source-reported static-LDA BM3; corrected calcium-ferrite structure, cell basis, and Figure 2d curve verified. |
| `litcurate_dba4bb2c50bd7a50` | Kawai and Tsuchiya (2010) theoretical Na3Al3Si3O12 CF: `35.3`, `196.5`, `4.5` | **REJECT as a Kawai-2012 source record** | Citation-reported comparison with a different composition and underlying paper; audit under the 2010 primary source. |
| `litcurate_74a5d6778933b7d6` | Guignot and Andrault experimental CF: `36.3`, `200(3)`, `4 fixed` | **REJECT as a Kawai-source record** | Citation-reported comparison, not a result of this paper. |
| `litcurate_6217113a1e9daf50` | Yutani et al. experimental CF: `36.1`, `241(3)`, `4 fixed` | **REJECT as a Kawai-source record** | Citation-reported comparison, not a result of this paper. |

Final production count for DOI `10.2138/am.2012.3915`: **2 records**.  No
same-DOI candidate is held: both source-reported records are accepted, and all
seven comparison rows are rejected from this DOI while remaining eligible for
separate primary-source audits.

## Reproduction assets

- `peritheos/data/materials/namg2al5sio12_nal.eosmat`
- `peritheos/data/materials/namg2al5sio12_cf.eosmat`
- `peritheos/data/datasets/namg2al5sio12-nal-kawai-2012-figure2c-vector-digitized.csv`
- `peritheos/data/datasets/namg2al5sio12-cf-kawai-2012-figure2d-vector-digitized.csv`
- `scripts/reproduce_kawai_2012_nal_cf.py`
- `tests/test_kawai_2012_nal_cf.py`

Run `python scripts/reproduce_kawai_2012_nal_cf.py` for the independent curve
and refit diagnostics.
