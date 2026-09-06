# Jackson et al. (2006): single-crystal elasticity of (Mg0.94Fe0.06)O

Primary source: J. M. Jackson, S. V. Sinogeikin, S. D. Jacobsen, H. J. Reichmann, S. J. Mackwell, and J. D. Bass, *Journal of Geophysical Research: Solid Earth* **111**, B09203 (2006), DOI [10.1029/2005JB004052](https://doi.org/10.1029/2005JB004052).

## Audit result

The primary paper reports Brillouin velocities and pressure-dependent single-crystal elastic constants to 20 GPa. Its third-order finite-strain analysis parameterizes elastic moduli as functions of pressure; it is not a hydrostatic P-V EOS fit.

- **REJECT — LitCurate source row 533.** K0S = 163(3) GPa and K0S' = 3.9(2) are adiabatic aggregate-elasticity coefficients fitted to all Brillouin observations. Although the paper calls the method a third-order finite-strain EOS, the fitted observables are acoustic velocities and elastic constants, not volumes. The sample's ambient lattice parameter is inherited from Jacobsen et al. (2002), and combining it with the acoustic fit as if it were an independently fitted P-V reference state would change the source model.
- **REJECT — LitCurate source row 534.** K0S = 163(3) GPa and K0S' = 3.8(2) are the hydrostatic-subset sensitivity fit from the same sample and experiment. It is neither an independent material record nor a P-V EOS.
- **REJECT — LitCurate citation row 535.** K0T' approximately 5.5(2) summarizes static-compression results for three more Fe-rich compositions from Jacobsen et al. (2002). Those values are not source results of Jackson et al. (2006) and must be assessed in their own primary paper.

## Reproduction status

No executable record was added. Peritheos volume EOS require a source-defined volume reference state and a pressure-volume relation; this paper instead constrains the adiabatic elastic tensor and its pressure derivatives.
