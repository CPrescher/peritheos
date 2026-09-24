# Barker (1987): verified theoretical study and graphical data

J. A. Barker, *High pressure equation of state for solid argon from interatomic
potentials*, **The Journal of Chemical Physics 86**(3), 1509–1511 (1 February
1987), [doi:10.1063/1.452187](https://doi.org/10.1063/1.452187).

The user supplied the matching full paper on 24 September 2026. All three article
pages were visually inspected. The PDF has four pages including its publisher
cover. The earlier full-text access blocker is **resolved**. The paper remains
**theoretical**, **not_reproduced**, and a **supporting-study/data outcome**, because
it does not print a self-contained numerical EOS or the potential coefficients.
The `peritheos/data/studies/argon-barker-1987.json` now includes
checksummed numerical assets rather than an abstract-only record.

## Actual calculation and reference state

Section II, p. 1509, uses Monte Carlo with the Barker–Fisher–Watts pair potential
and Axilrod–Teller three-body interaction (BFW+AT), referring to Barker, Fisher
and Watts (1971), [doi:10.1080/00268977100101821](https://doi.org/10.1080/00268977100101821),
for the method. For the Aziz–Chen and Koide–Meath–Allnatt comparisons it assumes
the same thermal pressure as BFW at equal density. These are not independent
Monte Carlo simulations of every displayed alternative. Ross's effective exp-6
potential is a separate comparison, not Barker's model.

Figures 4–5 display the 298 K isotherm. Figure 3 displays 0 K results up to
20 kbar (2 GPa), compared with Anderson and Swenson (1975). Figure 3's caption
incorrectly labels that experimental source as reference 15; Section II and the
bibliography correctly identify it as reference **16**. Figures 4–5 use Ross,
Mao, Bell and Xu (1986), [doi:10.1063/1.451346](https://doi.org/10.1063/1.451346),
for experimental comparisons. There are no new solid compression measurements.

The three-page article contains no displayed mathematical equation, potential
coefficient table, explicit V0/K0 fit, simulation size, convergence criterion,
or cutoff/tail-correction recipe. Exact reproduction requires its cited method
and potential sources and a numerical Monte Carlo implementation. No static
lattice sum, Vinet fit or interpolation is substituted as a callable EOS here.
This limitation does not establish that the original calculation is irreproducible.

The solid structure is not explicitly called fcc or hcp in this paper. Solid
argon is verified, but a specific phase is not assigned from expectation. The
reference pressure and zero-pressure volume remain unspecified. The 298 K and
0 K curves are distinct isotherms, not a fitted rectangular P–T domain.

## Units, coverage and pressure scale

Axes in Figures 3–5 use kbar and **cm³/mol**, not atomic or conventional-cell
volume. Pressure in GPa is kbar/10. Atomic volume in Å³ is molar cm³/mol times
10²⁴/N_A, with exact SI N_A=6.02214076e23. A four-atom conventional-cell conversion
would require another factor of four, but is not applied to this phase-unassigned
supporting study. The current fcc material's cell convention is unchanged.

Figure 4 spans 5–20 cm³/mol and 0–900 kbar; Figure 5 enlarges 8–22 cm³/mol and
0–450 kbar. These are **plot axes**, not declared model validity or fit bounds.
Section I describes experimental coverage to 800 kbar (80 GPa).
Section II identifies ruby fluorescence and cites Mao et al. (1978, 1979),
DOIs 10.1063/1.325277 and 10.1063/1.1135966. It allows an overall pressure-scale
uncertainty possibly as large as 10%. This is not a one-sigma uncertainty for
every point. Raw ruby shifts and a pointwise scale/recalibration recipe are absent.

## Numerical assets and honest diagnostics

`peritheos/data/datasets/argon-barker-1987-table1-liquid.csv`
transcribes all seven rows at **100 K and 27.04 cm³/mol, liquid argon**. The two
source quantities are dimensionless pV/NkT and U/NkT. One row is labeled
Experiment, and six rows are theoretical BFW/AC/KMA with/without AT. The caption
says the tabulated values include quantum corrections; no uncertainty columns
or explicit source for the experimental row are supplied. The table must not be
shown as a solid isotherm or counted as seven experimental measurements.

`peritheos/data/datasets/argon-barker-1987-figure5-selected.csv`
contains **21 samples of the solid BFW+AT theoretical curve and 11 selected,
visually separable Ross experimental circles**. This is deliberately incomplete:
merged or ambiguous symbols were omitted, and the figure covers only the enlarged
low-pressure portion of Figure 4. No numerical values were fabricated for the
undigitized high-pressure region or the 0 K comparison.

Digitization uses a 300 dpi rendering of PDF page 4, cropped with Poppler
`-x 120 -y 120 -W 1080 -H 1120`. In that crop, x=150/958 maps to 8/22 cm³/mol,
and y=918/114 maps to 0/450 kbar. Original pixel coordinates are retained for every
point. Selected curve centers were inspected against the black line, and a colored
overlay was visually reviewed. A conservative ±3-pixel picking bound corresponds
to ±0.05198 cm³/mol and ±0.16791 GPa. These are digitization bounds, **not** measured
error bars or Monte Carlo uncertainty. Six printed decimal places retain conversion
reproducibility and do not imply that precision in the original figure.

Run `python -m scripts.reproduce_argon_barker_1987` to regenerate the
[diagnostic report](../data/argon-barker-1987-reproduction.json). It converts Table I
to dimensional pressure/energy and compares the selected Ross circles with linear
interpolation of the digitized BFW+AT curve. Two circles fall outside the sampled
curve domain and are excluded rather than extrapolated. The remaining nine
comparisons are graphical diagnostics only. They do not independently evaluate the
potential, perform a fit, validate an EOS, or reproduce the Monte Carlo calculation.

## PDF provenance and prior access history

Stable PDF: `/Users/clemens/Documents/Peritheos/argon-sources/barker-1987/barker-1987.pdf`.
Original user file: `/Users/clemens/Downloads/1.452187.pdf`.
SHA256: `d8a43b6186bdfbeaa00b8113e5275efe396bf44bab2fc4664fa6a68c2d6a6b98`.
The title, author, DOI and printed pages establish identity. Provenance is
**user supplied**; the user's acquisition URL was not supplied. The DOI and
[publisher PDF locator](https://pubs.aip.org/aip/jcp/article-pdf/86/3/1509/18963166/1509_1_online.pdf)
are bibliographic links, not a claimed download origin. The copyrighted PDF and
rendered images remain outside git.

Earlier publisher/IBM/repository/Zotero checks found no accessible full text.
The exact ResearchGate publication 254001099 offered only Request full-text,
including in the signed-in session; no author request was sent. These findings
remain historical, and no longer block reading the supplied paper.

## Studio and Zotero integration

Render a supporting theoretical study with separate **published theoretical curve**,
**Ross experimental comparison (digitized subset)** and **liquid Table I** labels.
Do not register a selectable EOS or display independent reproduction validation.
The existing material loader does not automatically index `data/studies/`; Studio
integration must expose the supporting data explicitly. Preserve the Ross origin
of every experimental mark and do not duplicate these as new Barker measurements.

The integration task owns attaching the supplied PDF to Zotero. The handoff records
its stable path/hash and bibliographic metadata. No Zotero write is performed here.
Wittlinger's hcp identity and not_reproduced status are untouched.
