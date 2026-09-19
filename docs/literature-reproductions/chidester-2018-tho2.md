# Chidester et al. (2018): experimental ThO2 thermal equations of state

Two published experimental records are accepted, as separate materials:
`thorianite_chidester_2018_bm3_linear_thermal` and
`tho2_cotunnite_chidester_2018_bm3_linear_thermal`. No computed EOS or
literature-comparison parameter set is imported. Published coefficients are
preserved; diagnostic refits do not create additional executable records.

## Primary sources and retrieval

- Chidester, B.A., Pardo, O.S., Fischer, R.A., Thompson, E.C., Heinz, D.L.,
  Prescher, C., Prakapenka, V.B., and Campbell, A.J. (2018), *High-pressure
  phase behavior and equations of state of ThO2 polymorphs*, American
  Mineralogist **103**, 749–756,
  [DOI 10.2138/am-2018-6212](https://doi.org/10.2138/am-2018-6212).
  The publication-of-record PDF was recovered from Zotero item `NU5RF7W2`,
  attachment `33PKK488`, outside Methods → EOS Library (`JT8V6LUL`). Its
  journal pagination, DOI, Table 1, equations, figure and final supplement
  reference were checked. No Zotero library changes were made.
- [Official deposit AM-18-56212](http://www.minsocam.org/MSA/AmMin/TOC/2018/May2018_data/AM-18-56212.zip),
  linked by the [May 2018 supplements page](http://www.minsocam.org/MSA/AmMin/TOC/2018/May2018_data/May2018_data.html).
  HTTPS failed; HTTP supplied the publisher's ZIP and two complete CSV tables.
  Original ZIP and CSV bytes, retrieval URL and SHA-256 hashes are bundled in
  `peritheos/data/datasets/chidester_2018_sources/`. The publisher states that
  pure datasets are not copyrightable and other supplemental content is
  copyrighted by the authors. The article PDF is checksummed, not redistributed.
- The [author's publication page](https://raf.scholars.harvard.edu/publications/high-pressure-phase-behavior-and-equations-state-tho2-polymorphs)
  and [Caltech record](https://authors.library.caltech.edu/records/x6z0p-c8e51)
  independently identify the final article; Caltech is metadata-only.
  DOI/title searches for corrections and errata, these primary landing pages,
  the official deposit index and Crossref metadata were checked on 2026-09-19.
  No correction was located; Crossref has no update/relation entries. The
  publisher's GSW landing page could not be fetched, so this is not a claim
  that every publisher correction channel was accessible.
- Idiri, M., Le Bihan, T., Heathman, S., and Rebizant, J. (2004), *Behavior of
  actinide dioxides under pressure: UO2 and ThO2*, Physical Review B **70**,
  014113, [publisher full text](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.70.014113/fulltext).
  This primary experimental source supplies crystallography only.

Article PDF SHA-256:
`015178940c8c453a52a7b6842376c153240b57ba35e364de5ccb5e494b0ceb70`.
Idiri PDF SHA-256:
`e404759786d45f2b3951ed7ae196ea76d9ff403375c099cb5176e7256b91b775`.

## Candidate inventory and exact model

The two Table 1 rows labelled “This study”, method LH-DAC, are the entire
source-owned EOS inventory. The remaining rows quote earlier measurements
(Idiri, Dancausse, Olsen, Clausen, Macedo) or calculations (Olsen, Perry, Wang,
Li, Boettger, Kanchana, Boudjemline, Kelly, Shein, Sevik, Harding). They are
comparison literature, not Chidester fits. Raman-mode slopes and mode
Grüneisen parameters in Table 2 do not define another thermodynamic EOS.

Equations (1)–(3), p.754, specify

\[
P=3K_0f(1+2f)^{5/2}\{1+\tfrac32 f(K'_0-4)\}
  +C(T-300),\qquad
f=\tfrac12[(V_0/V)^{2/3}-1],\quad C=\alpha K_T.
\]

The coefficient **product** C is constant, not expansivity alone. The source
calls this “Mie–Grüneisen type”, but no Debye temperature, atom-count thermal
energy normalization, Grüneisen parameter or volume-dependent thermal term is
present. Existing `BM3` with fixed `K0_prime=4` plus `LinearThermalPressure`
is exact (equivalently BM2 for the reference isotherm). No model changes are
needed. Reference pressure is zero at the 300 K reference volume.

| Quantity | Thorianite | Cotunnite-type ThO2 | Source/provenance |
|---|---:|---:|---|
| V0, cm³/mol ThO2 | 26.379(7), fixed | 24.75(6), fitted | Table 1; p.754 text |
| K0, GPa | 204(2), fitted | 190(3), fitted | Table 1 |
| K0 prime | 4, fixed | 4, fixed | Bold Table 1; p.754 |
| C, GPa/K | 0.0035(3), fitted | 0.0037(4), fitted | Table 1; Equation 3 |
| Tr, K | 300, fixed | 300, fixed | Equations 1–3 |

The thorianite V0 error is an independent ambient measurement error and is
retained despite V0 being fixed in the regression. Fixed K0 prime has no
reported uncertainty; null is used, not zero. Parameter and observation error
confidence conventions, covariance, objective, weights and exact numerical
fitting sequence are not specified. No covariance is invented.

## Volume and diffraction identity

Both phases have Z=4 ThO2 formula units (Th4O8) in the conventional cell.
Table 1 volumes are **cm³ per mole of ThO2**, not Å³, nor molar volume per mole
of atoms. The supplement provides both unit-cell Å³ and molar volumes, which
independently confirms this normalization. Conversion is
`V_cell = V_molar * 4e24 / 6.02214076e23`, and the same factor applies to
V0 errors. Rounded duplicate volume columns agree within 0.10 Å³. All
original rounded quantities remain present; no column is silently replaced
by a conversion of another one.

Thorianite uses Fm-3m, No.225, Th 4a (0,0,0), O 8c (1/4,1/4,1/4), full
occupancies, following Idiri's introduction. The structural lattice is
Chidester's measured ambient a=5.5958(5) Å; its cube and the converted Table 1
V0 differ only through published rounding.

Cotunnite uses Pnma, No.62, Idiri p.014113-2's experimental 36 GPa cell
(a,b,c)=(5.898(8),3.600(3),6.862(7)) Å and three full 4c sites:
Th (0.261,1/4,0.111), O (0.163,1/4,0.449), O (0.12,3/4,0.27).
The second oxygen is preserved as printed, as a symmetry-equivalent member
of its 4c orbit, not moved to y=1/4 without transforming x and z. The
compressed structural cell volume (~145.7 Å³) intentionally differs from the
EOS's extrapolated zero-pressure V0. Each structure has exactly four Th and
eight O atoms; no fallback or borrowed isostructural compound is needed.

## Observations, exclusions and pressure calibration

All **41** Supplemental Table 1 rows and **59** Supplemental Table 2 rows
are bundled, including sample IDs, standard identities, P/T/volume errors,
lattice lengths and errors, both volume bases, and the cotunnite axis ratios.
Normalized headers and delimiters make the tables executable. Tests compare
every original CSV field against the original publisher bytes. Added
`included_in_fit` flags encode only the documented selection.

The fluorite fit uses all 21 high-temperature rows and the eight ambient/low-P
300 K rows, excluding 12 room-temperature rows at P≥15 GPa. All excluded
rows remain available up to 54.6 GPa. This is essential: a volume anomaly
starts near 15–16 GPa, well before the high-pressure structure is fully
resolved. The source does not exclude the high-temperature phase-resolved
volumes where both structures coexist. They are separate phase observations,
not mixture-average volumes.

The resulting thorianite marginal fit envelope is 0–28 GPa and 300–2135 K;
its room-temperature branch is constrained only below 15 GPa. Cotunnite has
59 rows at 16.6–62 GPa and 1018–2553 K. There are **no usable room-temperature
cotunnite volume observations in this study**; its 300 K reference curve,
V0 and K0 are inferred from hot data. The marginal envelopes do not guarantee
phase stability at every pressure/temperature combination. Heating kinetics,
coexistence and structural distortion restrict extrapolation; no >62 GPa or
room-temperature cotunnite experimental validation is claimed.

B12 uses Ar with Ross et al. (1986) pressure calibration and a separately
stated 3% accuracy. B5 and B25 use KBr; B70 uses KCl, both on Dewaele et al.
(2012), Physical Review B 85, 214105. Only the exact KCl calibration already
has a bundled executable link. Pt was a laser absorber, not the reported
pressure standard. Diamond Raman (Akahama–Kawamura 2007) belongs to the
separate Raman phase-transition experiments and does not calibrate these
EOS observations.

Sample temperatures were corrected by -3% for axial gradients; the alkali
halide temperature was estimated at the midpoint between sample surface and
diamond. The supplement reports sample temperatures and calculated pressures,
not raw Ar/KBr/KCl volumes or marker temperatures. Original observation-level
pressure recalculation is therefore unavailable, even though standard
identities and source citations are known. No missing readings are synthesized.

## Independent reproduction and refit

Run `uv run python scripts/reproduce_chidester_2018_tho2.py --check`.
The committed numerical ledger was generated in the locked Python 3.9.6,
NumPy 2.0.2, SciPy 1.13.1 environment.
The [machine-readable result](../data/chidester-2018-tho2-reproduction.json)
contains every selection, independent source-equation/native comparison,
measured-state benchmark, and weighting diagnostic. Source equations are
implemented directly in molar units independently of Peritheos. Native
pressures agree within 2e-13 GPa across all selected observations.

Four independent, off-reference measured states reproduce the source curves
within twice the printed pressure errors (a numerical tolerance, not a
claimed 95% interval):

| Phase / sample | T, K | Observed P, GPa | Published EOS P, GPa |
|---|---:|---:|---:|
| Thorianite B5_209 | 1734 | 18.7(8) | 17.416224 |
| Thorianite B5_253 | 1986 | 28(1) | 27.868357 |
| Cotunnite B25_446 | 1781 | 17.4(9) | 17.402518 |
| Cotunnite B70_562 | 2480 | 62(1) | 61.213837 |

Published pressure RMS residuals are 0.910071 and 0.967702 GPa respectively,
consistent with the scatter displayed in Figure 7. No fitted parameter or
reference-state identity is counted as an independent observation benchmark.

Each diagnostic jointly fits the source-free parameters to the selected
molar-volume rows, holding K0 prime and Tr fixed and also holding thorianite
V0 fixed. Three explicitly chosen objectives test the unreported weighting:
unweighted pressure residuals; pressure-error weights; and frozen propagated
P/V/T error weights evaluated at the published coefficients, adding Ar's
3% accuracy in quadrature. The latter treats scale error diagonally as a
sensitivity assumption, not a claim of independent calibrant errors.

| Phase / diagnostic | V0, cm³/mol | K0, GPa | C, GPa/K | P RMS, GPa |
|---|---:|---:|---:|---:|
| Thorianite unweighted | 26.379 fixed | 200.095161 | 0.003733790 | 0.900684 |
| Thorianite P errors | 26.379 fixed | 207.730742 | 0.003216832 | 0.933504 |
| Thorianite propagated errors | 26.379 fixed | 203.127708 | 0.003661393 | 0.916805 |
| Cotunnite unweighted | 24.795462 | 189.566704 | 0.003401732 | 0.956805 |
| Cotunnite P errors | 24.769131 | 189.600471 | 0.003596966 | 0.961826 |
| Cotunnite propagated errors | 24.763289 | 191.039762 | 0.003524437 | 0.959574 |

Every varied coefficient of the propagated-error diagnostics is within its
printed parameter error. This supports numerical similarity; neither the
original weighting nor statistical parity is recovered. The ledger records
`similar`, not `parity`. The unweighted thorianite K0 is nearly two printed
errors lower and demonstrates the significance of weighting. The published
records remain unchanged; no new refit is promoted to a production record.

One minor prose discrepancy does not affect execution: p.754 calls the
cotunnite phase 6.2% denser at ambient conditions. The printed V0 values imply
a 6.18% volume decrease, or a 6.58% density increase. This is an extrapolated
contrast, not independent ambient cotunnite evidence.

## Acceptance boundary

Both source-owned experimental fits satisfy the equation, normalization,
structure, observation and numerical-reproduction gates. No candidate is
withheld. The remaining limits are unreported fit weights/confidence/covariance,
missing raw calibrant readings, and lack of measured room-temperature
cotunnite volumes. These limits are recorded in each card and the refit ledger.
The tests cover source transcriptions/checksums, cell contents, normalization,
fixed/fitted parameters, thermal evaluation/inversion, validity envelopes,
JSON Schema, metadata-preserving Python and Rust interchange, and refits.
