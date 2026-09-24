# Ross et al. (1986): supporting study, EOS not reproduced

Requested publication: M. Ross, H. K. Mao, P. M. Bell and J. A. Xu,
*The equation of state of dense argon: A comparison of shock and static studies*,
J. Chem. Phys. **85**(2), 1028-1033 (15 July 1986),
[doi:10.1063/1.451346](https://doi.org/10.1063/1.451346).

**Outcome: supporting study pending primary full text, `not_reproduced`.**
No executable Ross EOS is registered. This does not establish that the journal
paper lacks a reproducible EOS: its full text was not accessible in this audit.
The paper's abstract confirms static room-temperature measurements to 800 kbar
(80 GPa), comparison with a theoretical reduction of shock data, and a proposed
pressure-standard extrapolation to 3-4 Mbar. The latter is not measured coverage.

## Source identity and remaining PDF blocker

The full text retrieved is the **distinct July 1985 conference precursor**,
[UCRL-93030 / CONF-850736-45](https://www.osti.gov/biblio/5471606), prepared
for the APS Spokane meeting of July 22-25, 1985. Its cover, title, authors,
equations and figures were visually checked after Poppler rendering. The
seven-page PDF is not the six-page 1986 journal article. The journal has 22
references on its publisher page; the precursor has 12. Do not assign the
journal DOI to this PDF or claim that the requested journal PDF was obtained.

The AIP browser page explicitly reported no current access and offered purchase.
The Crossref version-of-record PDF URL returned HTTP 403. OSTI journal record
5817597 has no full-text link and its purl returned 404. OpenAlex provided only
the closed publisher location; Semantic Scholar reported CLOSED without a PDF
URL. Local Zotero searches for the title and argon found no Ross item. Focused
author/repository searches located the precursor only; the historical
`mao.gl.ciw.edu` author website timed out. No purchase, author contact or Zotero
write was performed. The exact journal PDF request remains open.

Related PDF, outside git:
`/Users/clemens/Documents/Peritheos-sources/argon-ross-1986/precursor.pdf`.
Origin: <https://www.osti.gov/servlets/purl/5471606>.
SHA256: `fbc357ed1cd980fd495f26b22932adb3540d3d36d9c10349b32899fc42800302`.
Size: 194349 bytes; seven pages. The [handoff manifest](../data/argon-ross-1986-handoff.json)
records journal metadata, source identity, access blockers and import guidance.

## What the precursor actually supplies

The following findings refer to **UCRL-93030**, not unverified journal content:

| Quantity | Verified convention |
| --- | --- |
| Static phase and temperature | fcc solid, 293 K (introduction; Figs. 1-3); no hcp fit |
| Static measurements | DAC diffraction and ruby pressure; headline maximum 800 kbar |
| Figure 1 volume | cm3/mol of Ar atoms; not atomic A3 or conventional-cell A3 |
| Package normalization | Four atoms per fcc conventional cell: `Vcell = Vmolar * 4e24 / NA` |
| Pressure conversion | 1 kbar = 0.1 GPa; 1 Mbar = 100 GPa |
| Ruby scale, Eq. (1) | `P(Mbar)=19.04/7.665*((1+delta_lambda/lambda0)^7.665-1)`; lambda0 at one bar |
| Pressure recalculation | Original ruby wavelengths absent; reported pressures retained |
| Potential, Eq. (2) | exponential-six, epsilon/kB=122 K, r*=3.85 A, alpha=13.2 and 13.0 |
| Theoretical method | Monte Carlo pressure and energy; pair correlation functions matter |
| Shock selection | Liquid shock data below 400 kbar used to constrain the potential; higher-temperature electronic excitation discussed separately |
| Reference state | No verified journal P0/V0/K0 set; precursor isotherm is at 293 K, not a zero-K cold curve |
| Uncertainty | No individual errors recovered; alpha alternatives and an approximately 10% discussion are not parameter standard deviations |

The pair potential is

```text
phi(r)/kB = 122/(alpha-6) *
            [6*exp(alpha*(1-r/3.85)) - alpha*(3.85/r)^6].
```

`precursor_exp6_energy_kelvin` implements only this verified pair energy in the
audit script, on separations at least 2 A. Its minimum is -122 K at 3.85 A.
This is **not a pressure EOS**. Replacing the Monte Carlo calculation by an
unrelaxed static lattice sum, a Birch-Murnaghan fit or a Vinet fit would not
reproduce the published finite-temperature calculation. No such replacement
is made, and no fit residual or parameter covariance is claimed.

Fig. 1 separately labels a liquid Hugoniot with calculated temperatures
(4633, 12000 and 16600 K). These are not 293 K solid observations.
Fig. 2 gives calculated pair correlations; Fig. 3 gives calculated 293 K
isotherms to 5 Mbar. Neither figure supplies new experimental P-V points.

## Recoverable static observations

`argon_ross_1985_precursor_figure1_subset` retains **19 clearly resolved DAC
triangle symbols** from precursor Fig. 1, roughly 1.81-57.66 GPa. It is an
explicitly incomplete digitization, not the complete 80 GPa dataset and not a
journal-table transcription. Overlapping high-pressure triangles were omitted;
the subset is unsuitable for claiming a reproduction of the original fit.
Calculated lines, shock symbols and figure-2 example states are excluded.

The CSV preserves selected pixel coordinates, source molar volume, converted
cell volume, pressure and 293 K. The affine calibration is reproducible from
the metadata: render PDF page 4 at 240 dpi, crop `(580,170,1280,1140)`, use
`x=60,620` for `V=8,19 cm3/mol` and `y=920,125` for `P=0,80 GPa`.
Manual coordinate bounds of +/-10 pixels correspond to approximately
0.20 cm3/mol and 1.01 GPa. These conservative reading bounds are **not**
published measurement uncertainties or one-sigma errors. Stored decimals
preserve conversion arithmetic and do not indicate experimental precision.

The dataset has `used_by_eos_records=[]`; its observations were not used to fit
either Dewaele or Ono. The `argon_fcc.supporting_studies` extension identifies
the requested 1986 paper and explicitly relates the separate precursor dataset.
Wittlinger's hcp identity and `not_reproduced` status are unchanged.

## Audit and Studio integration

Run `python -m scripts.reproduce_argon_ross` to regenerate
[the audit report](../data/argon-ross-1986-reproduction.json). It checks the
coordinate conversions and audits the precursor pair energy; it deliberately
reports no EOS fit and a null pressure residual. Tests cover source/version
separation, checksum, 293 K, conventional-cell normalization, absent fabricated
uncertainties, metadata round-trip and the pair-potential minimum/derivative.

Studio should list Ross as a supporting study with a pending journal-source
badge, offer the distinctly labeled precursor subset as historical comparison
observations, and show no Ross EOS curve. It must not place the precursor
dataset in another paper's fitted-data selection. All local paths and the
related-PDF SHA256 are in the handoff manifest; imports are serialized by the
integration task. Obtaining the actual journal PDF remains a required follow-up.
