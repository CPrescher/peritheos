# Errandonea et al. (2006): laser-annealed fcc Ar and hcp coexistence

Primary source: D. Errandonea, R. Boehler, S. Japel, M. Mezouar, and
L. R. Benedetti, *Physical Review B* **73**, 092106 (28 March 2006),
[doi:10.1103/PhysRevB.73.092106](https://doi.org/10.1103/PhysRevB.73.092106).
The four-page publisher PDF was obtained from the
[APS full-text endpoint](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.73.092106/fulltext).
Identity, local PDF location, checksum, Zotero lookup, and integration details
are recorded in [the handoff manifest](../handoffs/argon-errandonea-2006.json).
The copyrighted paper and extracted figures are not committed to the repository.

## Published fcc equation and reference state

Page 092106-3 explicitly fits the **fcc** observations to third-order
Birch-Murnaghan, with V0=143(11) Å³, B0=6.5(5) GPa, B0'=5.1(3).
The record is `argon_fcc_errandonea_2006_bm3`, using the existing native/Python
BM3 implementation. No substitute equation or new native model is needed.
The independent check uses Eulerian strain f=((V0/V)^(2/3)-1)/2 and
P=3 B0 f (1+2f)^(5/2) [1+1.5(B0'-4)f].

The prose calls V0 an atomic volume, but the numerical value and Figure 5's
40–160 Å³ ordinate are **four-atom conventional fcc-cell volumes**.
143 Å³ per atom would imply 572 Å³ per fcc cell and be incompatible with the
displayed fcc observations and their comparison to Ross (1986).
We preserve the printed 143 and explicitly interpret it on the cell basis;
the atomic equivalent is 35.75(2.75) Å³. This is a normalization correction,
not a refit. All source parameter errors are retained, but their confidence
convention and covariance are not specified.

BM3 references P=0 at V0, an extrapolation rather than an ambient stable-solid
measurement. The measurements are at room temperature; Figure 4 explicitly
labels the room-temperature pattern **300 K**. Fe-foil laser annealing above
2000 K precedes diffraction collection. No thermal EOS is published here.
34.9–114 GPa is the explicitly identified diffraction span; 34.9 GPa is the
lowest new-study pattern identified in the text, not proof of a complete
numerical fit selection. The complete fit inputs are not tabulated.

## Pressure scale, phases, and supporting observations

Page 1 specifies ruby fluorescence with reference 17, Mao, Xu and Bell
(1986); pressures agreed with a W-gasket EOS from reference 18, Ruoff, Xia
and Xia (1992). The main scale is linked to `ruby_mao_1986`. No raw ruby
wavelengths or W lattice data are supplied, so pressure recalculation is
unavailable and the W cross-check remains only partially resolved.

The 34.9 GPa pattern indexes as fcc. Hcp reflections first appear at 49.6 GPa
after stress-relieving laser annealing; fcc and hcp coexist through 114 GPa.
The EOS is fitted to **fcc reflections even within coexistence**, not an
undifferentiated mixture. Hcp volumes are reported equal to fcc volumes
within experimental accuracy, but no hcp volume table, independent hcp
coefficients, or separate hcp equation is provided. Hcp c/a stays close to
1.633; neither individual lattice measurements nor errors are given.

`argon-errandonea-2006-phase-observations.csv` preserves the three explicit
room-temperature pressure/phase observations and the heated 49.6 GPa
comparison. The heater reached 2250 K; Ar had a thermal gradient with no
single numerical temperature. Hcp peaks disappeared during heating.
This asset is supporting evidence, not P–V points or input to Wittlinger's
historical hcp fit. It remains standalone because the material dataset
schema requires a linked EOS; inventing an hcp EOS link would be misleading.
Wittlinger remains hcp and `not_reproduced`, unchanged by this addition.

The approximately 0.3 value at 114 GPa is an intensity-derived relative
amount. Figure 5's axis is Ihcp/(Ihcp+Ifcc), whereas prose uses “hcp/fcc
ratio.” Its caption says fcc(220), while the text specifies the nonoverlapping
fcc(200) reflection. These inconsistencies preclude treating it as a precise
thermodynamic phase fraction. The estimated completion near 300 GPa is an
extrapolation, not a measurement or an EOS-validity limit.

## Observations and numerical diagnostic

The bundled `argon_fcc_errandonea_2006_figure5` dataset contains eight isolated
filled-circle centers from Figure 5. Open Ross (1986) circles, the Anderson
and Swenson diamond, the fitted curve, and the inset are excluded. Overlapping
new-study marks near 50 GPa are deliberately omitted: this is a **partial**
recovery, not all published data. The 895×668 embedded figure bitmap defines
x=117 at 0 GPa, x=874 at 120 GPa, y=578 at 40 Å³, and y=11 at 160 Å³.
The CSV retains source pixel centers. The stated two-pixel bounds correspond
to 0.3170 GPa and 0.4233 Å³; these are analyst digitization bounds, not
experimental uncertainties or fit weights. Digitized endpoints need not
equal the rounded pressures named in the prose.

Run `python -m scripts.reproduce_argon_errandonea_2006` to regenerate
[the diagnostic report](../data/argon-errandonea-2006-reproduction.json).
Independent Eulerian-strain pressures agree with the native implementation
within 1.3e-13 GPa over 48–143 Å³. Atomic/cell basis invariance and native
inversion are tested. At 34.9, 55, and 114 GPa, the published equation gives
67.0484, 59.5782, and 48.8343 Å³ per conventional cell.

The published coefficients yield **1.917 GPa RMS** against the recovered
eight observations. An equal-pressure-weight diagnostic refit reduces this
to about 0.886 GPa, but V0 hits the imposed 300 Å³ bound and the other
parameters become poorly constrained. Numerical solver convergence does not
validate that fit. No alternative coefficients are installed. Rounded
published parameters, sparse high-pressure points, digitization, and missing
weights/covariance prevent an exact reconstruction of the original fit.
Consequently, primary-source transcription and equation evaluation are
verified, but the fit remains explicitly **`not_reproduced`**.

## Integration and validation

The Studio adapter can use generic `pressure_gpa` / `volume_a3` columns on
the conventional-cell basis. Digitization bounds must not be presented as
measured error bars. Phase observations require a supporting-study view,
not an hcp fit or mixed-phase curve. The integration task owns Studio,
common generated ledgers, and import of the verified PDF to Zotero C18.

Validation: `python -m pytest tests/test_argon_errandonea_2006.py
tests/test_argon.py -q` passes 7 tests, covering schema/resource hashes,
published parameters, native equation agreement, normalization, inversion,
round trip, source phase distinctions, and existing argon regressions.

## Combined-data test (26 September 2026)

The user requested fitting all the data after questioning whether the authors
also used the older measurements. The complete Ross (1986) Table I is now
available from the verified primary journal PDF. Its 42 measurements are
used directly instead of redigitizing overlapping open circles. Together
with the eight recoverable Errandonea points, this gives a 50-point union.
This is not a claim to recover Errandonea's full raw data or exact Figure 5
selection. Ross row 17 is printed as 247 kbar, although its volume lies near
the adjacent 347 kbar row; retain 24.7 GPa and test omission separately.

Ross's Table I explicitly reports **298 K**. Anderson and Swenson (1975)
report cryogenic argon isotherms at **4.2–77 K**, not room temperature.
Their diamond in Errandonea's Figure 5 is approximately (0 GPa, 149.63 Å³),
consistent with the original 4.2 K zero-pressure molar volume of 22.56
cm³/mol, but the 2006 caption does not specify its temperature. Pooling that
point into a 300 K isotherm is physically inconsistent without a thermal
correction. It is included only as the requested all-marker sensitivity.
The 298/300 K union also neglects a 2 K thermal correction, explicitly.

Run `python -m scripts.fit_argon_errandonea_combined` to regenerate the
[report](../data/argon-errandonea-2006-combined-fit.json) and
[comparison plot](../data/argon-errandonea-2006-combined-fit.png).
The [input snapshot](../data/argon-errandonea-2006-combined-inputs.json)
retains all Ross table fields, primary PDF hashes, units, and the diamond's
pixel coordinates. Source pressure uncertainties are preserved but not used
as cross-study weights, because the Errandonea uncertainties and common
covariance are unavailable. Both objectives give each observation equal
weight, not each study equal weight. Bounds and solver status are reported.

| Data / objective | N | V0 (Å³/cell) | K0 (GPa) | K0′ | Pressure RMS (GPa) |
|---|---:|---:|---:|---:|---:|
| Published parameters | — | 143(11) | 6.5(5) | 5.1(3) | — |
| Ross + Errandonea / pressure residuals | 50 | 141.693 | 6.1975 | 5.2591 | 1.934 |
| Above + cryogenic diamond / pressure residuals | 51 | 142.633 | 5.9787 | 5.2903 | 1.916 |
| Ross + Errandonea, omit Ross row 17 / pressure residuals | 49 | 138.833 | 7.3277 | 5.0186 | 1.203 |
| Ross + Errandonea / volume residuals | 50 | 184.444 | 0.9025 | 9.1532 | 2.069 |
| Above + cryogenic diamond / volume residuals | 51 | 150.282 | 4.5186 | 5.5162 | 1.935 |

The 50-point pressure-residual fit lies inside each published parameter-error
interval. The 51-point fit is close, but K0 is just below the printed lower
interval endpoint of 6.0 GPa. This agreement supports the plausibility of a
combined-data interpretation; it does not establish the authors' selection
or weighting. The substantially different volume-residual solution and the
sensitivity to the suspicious Ross row prevent a unique inference.
No source pressure was corrected, no catalog parameters were replaced, and
the original fit remains unresolved. The safe result is a clearly labeled
**diagnostic combined fit**, not a newly validated 300 K EOS.

Following the user's source-quality decision, the **preferred diagnostic now
excludes Ross Table I row 17** (247(13) kbar, 4.047 Å, 9.98 cm³/mol).
Its pressure is suspected to be a printing error because the neighboring
near-identical volume is reported at 347 kbar. The intended pressure is not
established: no correction to 34.7 or 42.7 GPa is applied. The archived source
row remains unchanged, and the report marks it excluded with this reason.
The unfiltered fits above remain explicit sensitivity comparisons.

The preferred 49-point, equal-pressure-weight fit therefore gives
**V0=138.833 Å³/cell, K0=7.3277 GPa, K0'=5.0186**, with pressure RMS
**1.203 GPa**. The cryogenic diamond remains excluded as well. These are
diagnostic coefficients; the published catalog equation is unchanged.
