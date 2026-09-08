# Fratanduono et al. (2020) 298 K Cu isotherm

The dedicated published record is `copper_fratanduono_2020_vinet3_298k`,
on the `copper` material card. It implements `Vinet3` (`vinet_3`),
using the final [primary article](https://doi.org/10.1103/PhysRevLett.124.015701),
Table I (page 015701-3), and the [official supplemental material](https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.124.015701),
Section S4, Eq. (2), page S2. Both source PDFs were inspected, including visual
checks of the equation and coefficient table; their SHA-256 hashes are stored
in the record's primary-source provenance.

With $X=\rho_0/\rho=V/V_0$ and $y=1-X^{1/3}$, the exact published form is

$$P=3K_0\frac{y}{X^{2/3}}\exp(\eta y+\beta y^2+\psi y^3).$$

| Central 298 K coefficient | Published value | Printed error |
|---|---:|---:|
| $\rho_0$ (g/cm³) | 8.939 | not stated |
| $K_0$ (GPa) | 133.6 | ±0.8 |
| $\eta$ | 6.29 | ±0.8 |
| $\beta$ | 2.06 | ±0.4 |
| $\psi$ | 1.65 | ±0.6 |

The API uses volume like the other equilibrium EOS classes. Conversion with
63.546 g/mol, four Cu atoms per conventional fcc cell, and exact Avogadro
constant gives $V_0=47.218085048721$ Å³. This is an explicit mass convention,
not a newly fitted volume. The record fixes V0 and leaves its uncertainty
unknown. It preserves the printed coefficient errors with unknown confidence
level and null covariance. In particular, eta is not K0-prime;
$K'_0=1+2\eta/3$ follows from differentiating this form. Ordinary Vinet is
recovered only when beta and psi vanish.

Table I and S4 explicitly identify 298 K, although the adjoining main-text
prose loosely says 300 K. The record follows the table and equation section.
The ramp experiments reach 2.30 TPa and diffraction establishes fcc structure
to 1.15 TPa; these are not exact fit endpoints for the reduced 298 K isotherm.
The record therefore supplies no numerical isothermal pressure interval.

Supplemental Section S5 and Table S1 contain the **reduced isentrope**.
They do not supply 298 K isothermal observations and are not registered as
fit input for this record. The numerical 298 K fit dataset, weights, covariance,
and confidence convention for Table I errors are unavailable in the inspected
sources. Separate upper/lower isotherm coefficient rows are envelope fits;
they are neither independent data nor errors on the central coefficients.
The refit ledger correctly classifies this executable parameterization as
`not_refittable`. No pseudo-data, pressure errors, weights or covariance are
created.

The [Shen–Smith reconstruction](shen-smith-2026-cu-refits.md) now obtains its
reference volume and pressures from this executable catalog record. Its
same-experiment/run pairing, phase/run masks, Pt first/last rows, and unweighted
fixed-V0 regressions are preserved. Equation tests independently transcribe
the density form, differentiate pressure to check bulk modulus, check the
ordinary-Vinet limit, and round-trip volumes. Record tests cover schema and
export round trips and verify same-run pressure reuse. Analytical checkpoints
are implementation checks, not a refit dataset.

## Attempt to reconstruct the thermal reduction

**A refit of the published 298 K EOS is not possible with the current information.**
Fratanduono's explicit Grüneisen discussion concerns reduction of measured
longitudinal stress to the principal isentrope. It does not explicitly specify
the subsequent Debye isentrope-to-298 K calculation used below. Applying
Kraus Eq. (16), its Debye caloric model, and its reference Debye temperature
is our inference from the cited predecessor, not a confirmed reproduction
of Fratanduono's thermal reduction. The candidate points are not recovered
original 298 K observations; their pressure agreement does not establish
coefficient parity or a superior EOS. The ledger remains `not_refittable`
and the catalog remains `parameterization_only`.

The missing numerical 298 K table does **not** mean that no reduction can be
attempted. We followed the citation to [Kraus et al. (2016)](https://doi.org/10.1103/PhysRevB.93.134105)
and checked the final primary PDF. Its Section III.B.4 Eq. (16) explicitly gives

$$P_{298}(\rho)=P_S(\rho)-\gamma(\rho)\rho
[e_{\rm th}(\rho,T_S)-e_{\rm th}(\rho,298)].$$

The thermal energy follows the Debye integral. Section III.B.3 gives
$\theta_0=343.5$ K and describes a quasiharmonic Debye solid, without
anharmonic or electronic heat-capacity contributions. Section III.B.1 Eq. (12)
gives $\gamma=1.41+(2.0-1.41)(\rho_0/\rho)^{13.6}$.
The exponent 13.6 is the *Grüneisen-law* eta, unrelated to Vinet3 eta.
Integrating $d\ln T_S/d\ln\rho=\gamma$ gives

$$F=(\rho/\rho_0)^{1.41}
\exp\left[\frac{0.59}{13.6}\left(1-(\rho_0/\rho)^{13.6}\right)\right],
\qquad T_S=298F,\quad\theta=343.5F.$$

Fratanduono explicitly adopts Kraus's Grüneisen law only for
$\rho_0/\rho\geq0.64$. Seven compressed Table S1 rows lie in this interval,
from 23.6 to 164.9 GPa along the source isentrope. The calculated thermal
correction increases from 0.383 to 1.760 GPa. Carrying the predecessor's
Debye temperature and caloric assumptions into 2020 remains an explicit
inheritance assumption; those quantities are not separately printed in 2020.

At higher compression Fratanduono determines gamma from the Hugoniot and
isentrope. We extracted the gray dashed **Our Hugoniot Fit** from the vector
paths of Supplemental Figure S2, using all five horizontal and six vertical
major ticks for affine calibration. The implementation then applies Kraus
Eq. (11):

$$\gamma=\frac{P_H-P_S}{\rho(e_H-e_S)},\qquad
e_H-e_0=\tfrac12P_H(1/\rho_0-1/\rho),\qquad
e_S-e_0=\int_{\rho_0}^{\rho}\frac{P_S(r)}{r^2}\,dr.$$

This supplies a **candidate** reduction for 46 additional table rows.
Whether the plotted Hugoniot is exactly the corrected internal input used
in the authors' iteration is unconfirmed. Its plotted extent ends at
22.548 g/cm³, so the final supported table density is 22.521 g/cm³,
corresponding to 1,248.3 GPa on the isentrope. The remaining 46 compressed
rows, through 2,331.8 GPa, are left unreconstructed. No high-pressure
Hugoniot or gamma extrapolation is used.

The source's ambient Table S1 entry is 8.938 g/cm³, whereas the specified
reference is 8.939. We preserve the printed row but exclude it from fitting;
thermodynamic integration starts at zero pressure and 8.939 g/cm³.
Table S1 pressure errors are not assigned to the reconstructed 298 K points:
they describe the isentrope, and the propagated thermal-model errors and
cross-row covariance are unavailable.

### Numerical outcome

The candidate 53-row isotherm differs from the published 298 K Vinet3 curve
by **1.562 GPa RMS**, but coefficient parity is **not** recovered:

| Parameter | Published central 298 K | Seven-row inherited-Debye diagnostic | 53-row candidate diagnostic |
|---|---:|---:|---:|
| K0 (GPa) | 133.6 ± 0.8 | 131.877 | 154.815 |
| eta | 6.29 ± 0.8 | 8.600 | 3.907 |
| beta | 2.06 ± 0.4 | -31.137 | 14.285 |
| psi | 1.65 ± 0.6 | 126.940 | -18.392 |

Both diagnostics minimize unweighted pressure residuals with fixed rho0.
These coefficients describe exploratory partial-range fits and are **not new
catalog EOS records**. Their strong variation illustrates that close pressure
curves do not imply recovery of the four correlated shape coefficients.
No coefficient covariance is inferred.

A useful control is fitting the **original isentrope** before any thermal
correction: unweighted fitting gives K0=106.748 GPa, and weighting by the
printed isentrope pressure errors gives K0=147.069 GPa, versus Table I's
138.9 GPa. Neither protocol reproduces its full coefficient set. The exact
fitting grid, objective, and weights cannot be established from the supplied
table alone.

### Figure labeling check

Figure S2's solid blue curve is labeled *NIF reduced Isotherm* in its legend,
while the caption calls it the isentrope. Its vector-extracted coordinates
match the numerical Table S1 isentrope to 0.863 GPa RMS over the compressed
rows. We therefore do not relabel that curve as independent isothermal data.
Subtracting the published isotherm fit from Table S1 even gives negative
apparent thermal corrections for 16 rows, reinforcing that subtraction of
separately approximated curves cannot recover the original correction.

### Reproducible artifacts

Run `python -m scripts.reproduce_fratanduono_2020_cu` or add `--check` to
verify the saved JSON. The optional `--extract-supplement path/to/PDF` uses
pdfplumber to re-extract the exact hash-checked source PDF.

- [Numerical reconstruction and fit report](../data/fratanduono-2020-cu-thermal-reduction.json)
- [53 qualified reconstructed 298 K points](../data/fratanduono-2020-cu-reconstructed-298k.csv)
- [All 100 original isentrope rows](../data/fratanduono-2020-cu-table-s1-isentrope.csv)
- [Digitized Hugoniot](../data/fratanduono-2020-cu-figure-s2-hugoniot.csv)
- [Digitized blue curve](../data/fratanduono-2020-cu-figure-s2-blue.csv)
- [Extraction provenance and axis calibration](../data/fratanduono-2020-cu-extraction.json)

These resources belong to the reconstruction audit, not the production
record's fit datasets. The published coefficients and the Shen–Smith pressure
reconstruction remain unchanged. An exact full-range refit still needs the
unclipped corrected Hugoniot or the authors' high-density gamma branch and
original 298 K fitting inputs/protocol.
