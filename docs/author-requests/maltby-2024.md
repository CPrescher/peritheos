# Draft request: Maltby 2024 solid argon reproducibility

Status: **draft, not sent**. Prepared 2026-09-24.
Intended corresponding author: Øivind Wilhelmsen,
`oivind.wilhelmsen@ntnu.no`, as listed on the
[publisher article page](https://doi.org/10.1063/5.0237497).
No contact or permission to redistribute additional materials has been obtained.

Subject: Reproducing the solid argon EOS in JPCRD 53, 043102 (2024)

Dear Professor Wilhelmsen and coauthors,

We are implementing the Helmholtz-energy EOS from your paper for a scientific
EOS library. We have the accepted manuscript and official supplement and would
appreciate clarification of a few reproducibility details.

Using the printed Tables 1 and 3 coefficients, an exact fcc coordination-shell
sum through (r/rNN)^2=64, and the full volume derivative of Eqs. (12)--(16), we
obtain -4.22598 MPa at 70 K and 23.97 cm3/mol. Solving for 1 MPa gives
23.89042753 cm3/mol, whereas Table 8 specifies 23.97 cm3/mol. Separate SI-unit
quadrature and numerical differentiation agree with our analytic derivative,
but an interpretation or transcription error on our side remains possible.

Could you share or clarify:

1. The full-precision final coefficients and numerical coordination-shell cutoff
   used to produce Table 8.
2. A reference implementation or a short executable calculation of that table,
   including the treatment of the explicit v dependence in the effective pair
   potential and differentiation of the long-range cutoff.
3. Table SI.1 lists 48 neighbours at (r/rNN)^2=14. Exact fcc enumeration gives
   no neighbours at this distance. Is this a table typo, and which neighbour
   list was used for the published calculations? Including the printed entry
   increases our volume discrepancy.
4. The selected primary observations, especially the 38 Dewaele rows, their
   weights/uncertainties, fixed/free parameters and any staged fitting settings.
   The caloric, expansivity and compressibility observations would allow us to
   check that a pressure-only fit does not degrade the other properties.

We have not changed your published coefficients or attributed our exploratory
refits to your paper. A public code/data deposit would be ideal; otherwise,
please let us know any attribution or redistribution conditions for material
you are able to share.

Thank you for your help.
