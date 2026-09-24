# Maltby, Hammer and Wilhelmsen (2024): fcc argon Helmholtz EOS

Status: **not_reproduced; experimental implementation, no bundled executable EOS record**.
This is an independent EOS publication, not merely a supporting measurement study.
The literal published equations are implemented in Python and Rust, but numerical
reproduction of the reported 70 K sample state is unresolved. Do not present this
candidate as a validated pressure standard or substitute a Vinet/BM fit.

## Primary sources and PDF identity

Tage W. Maltby, Morten Hammer, Øivind Wilhelmsen, “Equation of State for Solid
Argon Valid for Temperatures up to 300 K and Pressures up to 16 GPa,” *Journal
of Physical and Chemical Reference Data* **53**(4), 043102, published online
19 December 2024, DOI [10.1063/5.0237497](https://doi.org/10.1063/5.0237497).

The publisher page and publisher full text were readable through the user's
existing ResearchGate institutional access. A stable, legitimately open PDF
was then retrieved from NTNU's NVA repository:

- [NTNU/NVA record](https://nva.sikt.no/registration/0198cc531207-be3c0dc3-6af4-4953-97be-777a88a705a8),
  formerly [handle 11250/3175369](https://hdl.handle.net/11250/3175369).
- [Public file-link API](https://api.nva.unit.no/publication/0198cc531207-be3c0dc3-6af4-4953-97be-777a88a705a8/filelink/b4a4e90d-77d9-4854-a2cb-af09e365233e).
  This returns an authorized expiring download link. No authentication bypass.
- Accepted author manuscript, 23 pages, dated 3 January 2025; repository marks
  `AcceptedVersion`. Its title, three authors, model, and Tables 1, 3 and 8 match
  the paper; the NVA publication metadata explicitly links the same DOI.
  Do not call this file the publisher version of record.
- Saved as `/Users/clemens/Documents/Argon-EOS-Papers/maltby-2024/author-manuscript.pdf`,
  1,292,435 bytes, SHA256
  `da12e25c1cfd36917876fc138dcb71b197df3c2ee5c512b56b11c4684214bc4d`.
- Official [supplement](https://aip.figshare.com/articles/journal_contribution/Supplementary_material/27846567),
  [download](https://ndownloader.figshare.com/files/50603193), CC BY 4.0,
  DOI 10.60893/figshare.jpr.27846567.v1; saved as
  `/Users/clemens/Documents/Argon-EOS-Papers/maltby-2024/supplement.pdf`, 215,281 bytes,
  SHA256 `08986ed607d464ed16565b6007c0f83ce080f326daa84ddd5260b5c92e086229`.

No copyrighted PDF is committed. Bibliographic metadata, exact retrieval routes,
file identities and source coefficients are in `curation/argon/maltby-2024.json`.
The initial read-only Zotero search found no matching item. The coordinator has
since imported and hash-verified the accepted manuscript in Methods / EOS Library
(C18): article `CGV58DF9`, attachment `IDGZPNCM`. Supplement attachment is pending
at this handoff. The coordinator owns all Zotero writes; this task made none.

## Equation and unit audit

The phase is monatomic **fcc Ar**, four atoms per conventional cell. The paper
notes partial fcc-to-hcp deformation at pressures as low as 4 GPa; this model
contains no hcp branch or phase-transition prediction. Wittlinger remains hcp
and is untouched.

The molar Helmholtz energy consists of six terms:

1. Static fcc coordination-sphere sum, Eqs. (12)--(16), using the Buckingham
   pair potential and its explicit density-dependent three-body multiplier.
2. Separate zero-point energy, Eq. (4), Table 1; it is not the Debye zero-point
   term and must not be counted again in the vibrational contributions.
3. Debye external vibrations, Eqs. (7)--(10). **The paper's D3 integral is one
   third of the conventional normalized D3 used elsewhere in Peritheos.**
   The Debye temperature depends on both T and v.
4. Three Einstein terms, Eq. (11), with integrated volume laws and q_i=1.
5. Temperature-dependent anharmonic energy, Eq. (18), including theta_D,T(T).
6. Volume-dependent anharmonic correction, Eq. (19).

All coefficients are literal Tables 1 and 3 values. No optimization was used to
repair the mismatch. The Buckingham sigma is the outer physical zero of Eq. (15),
computed as 3.384249109109384 Å; the unphysical inner root is excluded.

The implementation uses molar volume in J/bar/mol (one unit = 10 cm3/mol),
pressure GPa, and molar energy J/mol. For the repository's conventional-cell
volumes, `v_cm3_mol = V_cell_A3 * 0.602214076 / 4`. Pressure is `-dA/dv`;
the derivative includes the effective potential's explicit v dependence and
the moving tail-integration cutoff. The tail integral is evaluated analytically.
The independent script instead uses SI units and direct quadrature of the
potential tail and phonon density of states.

The parameter v0=22.56 cm3/mol is a volume-law reference, not a measured
ambient-pressure cell or an assumed 300 K zero-pressure state. The absolute
pressure includes zero-point and temperature-dependent terms. The optional
pressure-increment API subtracts 300 K only to define an increment; total
pressure never subtracts it. T=0 is supported as the limiting cold state.

Eqs. (21)--(23) align the solid Gibbs energy/entropy to Tegeler's fluid EOS at
83.8058 K, 68.891 kPa, with melting entropy 14.3 J/mol/K. These additive energy
and entropy offsets are not evaluated here. Therefore the unshifted A cannot
be compared directly to source A, g or s. They do not explain a pressure or
volume mismatch. This candidate does not calculate solid-fluid equilibria.

## Cutoff and reproduction findings

The article specifies r_cut=r_L but does not give a numerical L or r_cut. The
supplement's FCC shell table extends to (r/r_NN)^2=64; that is evidence for the
available shell list, not evidence that the authors used that cutoff. Both
constructors require an explicit `shell_cutoff_squared` argument. The audit
uses 64 and tests 16, 36, 100, 144 and 256 separately. Shells are generated from
integer fcc coordinates. A complete follow-up comparison finds that Table SI.1
incorrectly includes 48 neighbours at squared distance 14, where fcc geometry
has none; every other entry through 64 matches. Using that literal table
worsens the discrepancy. See the [rounding and refit audit](argon-maltby-2024-refit.md).

At the Table 8 state T=70 K, v=23.97 cm3/mol, the candidate gives:

| Quantity | Table 8 | Literal equations, cutoff squared 64 |
|---|---:|---:|
| Pressure (MPa) | 1 | -4.22598 |
| Volume at 1 MPa (cm3/mol) | 23.97 | 23.89042753 |
| Cv (J/mol/K) | 22.93 | 22.91770 |
| Cp (J/mol/K) | 30.35 | 30.39589 |
| Expansivity (mK^-1) | 1.684 | 1.697439 |
| Thermal Grüneisen parameter | 2.745 | 2.746207 |
| Isothermal compressibility (GPa^-1) | 0.6412 | 0.646484 |
| Isentropic compressibility (GPa^-1) | 0.4844 | 0.487432 |

Table 8 describes either (T,v) or (T,p) input. Its 0.01 cm3/mol volume rounding
cannot account for the 0.07957 cm3/mol difference (0.332%). Across the larger
cutoffs the recovered volume remains about 23.8905 cm3/mol. This is an
unresolved reproduction discrepancy; it is not proof that the publication is
wrong. Unprinted coefficient precision, implementation conventions or an error
in this transcription remain possible. No coefficient is altered to conceal it.

Native and Python pressure agree with a separately coded SI quadrature and
numerical free-energy derivative within 8e-9 GPa for the three diagnostic states
(23.97 cm3/mol,70 K), (18,300), and (12,0). This checks transcription and
thermodynamic differentiation; **it is not independent validation of the
published numerical parameterization**. Full results and cutoff sensitivity:
`docs/data/argon-maltby-2024-reproduction.json`.

## Experiments, fit selection, pressure scales and uncertainty

The study is **hybrid**: regression to experimental thermodynamic properties,
with zero-point coefficients fitted to coupled-cluster calculations of
Schwerdtfeger et al. It reports no new experimental P-V-T table. Its primary
experimental sources and counts are transcribed from Table 2 in the curation
manifest, and Table 8 values are explicitly **theory check values**, not
measurements. Tables 4--6 summarize deviations and data counts, not raw rows.

Reported model validity is solid fcc Ar up to 300 K and 16 GPa. This does not
mean the whole rectangle is a stable solid. Primary data include sublimation
properties at 4--83 K, sublimation pressure 36.803--83.804 K, room-temperature
compression at 293 K, Dewaele data 10--300 K and melting pressures up to
260.646 K. Secondary comparisons extend beyond the model's advertised limits.

The fit uses the mean squared property residual normalized by each measured
property's uncertainty (Eq. 20), with a Nelder--Mead search. Coexistence fits use
Gibbs-energy differences and propagated experimental T/P uncertainty (Eq. 24).
No parameter standard errors/covariance are published. The full weighted,
selected primary raw observations are not supplied in the main paper or SI;
therefore an independent refit of the complete objective cannot be claimed.

There is no single newly established pressure scale: literature datasets retain
source-specific scales. The candidate does not recalibrate them. Dewaele's
already bundled observations retain their Au/ruby provenance (see the existing
argon-fcc audit). Maltby does not supply raw calibrant observations or sufficient
instructions to harmonize all the older measurements onto a new common scale.

The reproducible script compares every existing Dewaele supplementary row with
finite P in (0,16] GPa and T in [0,300] K: 130 rows, repeats and all runs retained.
This is a **diagnostic comparison subset**, not Maltby's 38-row primary selection.
RMS pressure residual is 0.48719 GPa; mean absolute relative pressure deviation
is 4.82986%. Original observations remain in their original checked dataset;
no theoretical sampled curve is stored as an observation.

Reported estimated model uncertainties: sublimation molar volume 0.03%; Cp 1%
from 20--50 K, 2% above, 6% below; isentropic compressibility 4%, isothermal 3%
(as abstract/Section 4.6; conclusion reverses the ordering); expansivity 2%
above 30 K, 5% below; melting pressure 0.4%, sublimation pressure 5%. The
all-data pressure deviation is 5%. These property estimates are not parameter
errors and are not adopted as validation tolerances for the Table 8 calculation.

## Execution and integration

Python: `peritheos.eos.experimental.Maltby2024Published(64)`.
Rust: `peritheos::experimental::maltby_2024::Maltby2024Published::new(64)`.
Both provide scalar pressure, bulk modulus, stable volume inversion, energy and
300 K pressure increments. Native code is portable and does not require Python.
Volume inversion searches 8--26 cm3/mol for a positive-bulk-modulus root; this is
an explicit numerical interval, not a claim of physical validity everywhere.

Run `PYTHONPATH=. python scripts/reproduce_maltby_2024_argon.py`,
`pytest tests/test_maltby_2024_argon.py tests/test_argon.py`, and
`cargo test -p peritheos experimental::maltby_2024 --lib`.
No EOSMAT dispatch/schema entry or Studio executable record is added while the
source reproduction is unresolved. Studio should expose the study's pending
status and source PDFs if it has a supporting-study view, with no selectable
validated pressure curve and no measured-point label on Table 8.

## Follow-up rounding and refit audit

The [constrained-refit study](argon-maltby-2024-refit.md) checks all printed
coefficient rounding intervals, all occupied cutoffs through squared distance
256, the complete SI shell table, five held-out experimental runs and external
Ono observations. It improves pressure errors on a source-informed subset but
worsens low-pressure and external checks; neither the published candidate nor
the diagnostic refit is promoted to validated EOSMAT.
