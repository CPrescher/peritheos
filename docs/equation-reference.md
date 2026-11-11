# Equation reference

This reference defines the equations evaluated by Peritheos. Constructor parameter
names are shown in code font; mathematical symbols use the corresponding
subscripts. The definitions here follow the implementation and the numerical
reference cases described under [Validation](validation.md). Original model
publications are collected under [References](references.md).

Each model family has a separate page. Use the category inventories to compare
parameters, then open a family for its full equations and qualifications.

## Notation and reference state

For an isothermal equation of state (EOS), $V_0$ is the volume at zero model
pressure and $K_0$ is the isothermal bulk modulus there. Peritheos uses

\[
K_T=-V\left(\frac{\partial P}{\partial V}\right)_T,
\qquad
K_0'=\left(\frac{\partial K_T}{\partial P}\right)_{P=0},
\qquad
K_0''=\left(\frac{\partial^2K_T}{\partial P^2}\right)_{P=0}.
\]

Thus $K_0'$ is dimensionless and $K_0''$ has inverse-pressure units. Except
for Holzapfel and thermal models, $V$ and $V_0$ may use any internally
consistent volume unit. Public pressure and modulus values are documented in
GPa. Every model requires $V>0$, with additional domains noted below.

## Isothermal equations

[Compare isothermal models and parameters](models/isothermal.md).

- <span id="birch-murnaghan-family" data-equation-target="../equations/isothermal/birch-murnaghan/"></span><span id="finite-pressure-bm2-reference" data-equation-target="../equations/isothermal/birch-murnaghan/#finite-pressure-bm2-reference"></span><span id="standard-zero-pressure-forms" data-equation-target="../equations/isothermal/birch-murnaghan/#standard-zero-pressure-forms"></span>[Birch-Murnaghan family](equations/isothermal/birch-murnaghan.md)

- <span id="natural-strain-family" data-equation-target="../equations/isothermal/natural-strain/"></span>[Natural-strain family](equations/isothermal/natural-strain.md)

- <span id="murnaghan" data-equation-target="../equations/isothermal/murnaghan/"></span>[Murnaghan](equations/isothermal/murnaghan.md)

- <span id="vinet" data-equation-target="../equations/isothermal/vinet/"></span>[Vinet](equations/isothermal/vinet.md)

- <span id="baonza-pseudospinodal" data-equation-target="../equations/isothermal/baonza/"></span>[Baonza pseudospinodal](equations/isothermal/baonza.md)

- <span id="generalized-rydberg-stacey" data-equation-target="../equations/isothermal/rydberg-stacey/"></span>[Generalized Rydberg-Stacey](equations/isothermal/rydberg-stacey.md)

- <span id="morse-and-sun-morse" data-equation-target="../equations/isothermal/morse/"></span>[Morse and Sun-Morse](equations/isothermal/morse.md)

- <span id="modified-tait" data-equation-target="../equations/isothermal/modified-tait/"></span>[Modified Tait](equations/isothermal/modified-tait.md)

- <span id="holzapfel" data-equation-target="../equations/isothermal/holzapfel/"></span>[Holzapfel](equations/isothermal/holzapfel.md)

## Thermal equations

[Compare thermal models and parameters](models/thermal.md).

Thermal models use $V$ in `J bar^-1 mol^-1`, $T$ in K, energy in
`J mol^-1`, and pressure in GPa. Their total pressure is

\[
P(V,T)=P_{\mathrm{ref}}(V)+\Delta P_{\mathrm{th}}(V,T),
\qquad
\Delta P_{\mathrm{th}}(V,T_r)=0.
\]

This reference-isotherm form also applies to the double-Debye Helmholtz models
when their optional `Tr` is set. With `Tr` omitted, those models instead retain
their original absolute simulated free-energy formulation.

Because energy divided by the public molar-volume unit produces bar, the factor
$10^{-4}$ converts thermal pressure to GPa.

- <span id="debye-plus-quadratic-thermal-pressure" data-equation-target="../equations/thermal/debye-quadratic/"></span>[Debye plus quadratic thermal pressure](equations/thermal/debye-quadratic.md)

- <span id="double-debye-helmholtz" data-equation-target="../equations/thermal/double-debye/"></span><span id="benedict-et-al-diamond-example" data-equation-target="../equations/thermal/double-debye/#benedict-et-al-diamond-example"></span><span id="logarithmic-moment-double-debye-variant" data-equation-target="../equations/thermal/double-debye/#logarithmic-moment-double-debye-variant"></span>[Double-Debye Helmholtz](equations/thermal/double-debye.md)

- <span id="temperature-dependent-reference-state" data-equation-target="../equations/thermal/reference-state/"></span>[Temperature-dependent reference state](equations/thermal/reference-state.md)

- <span id="mie-gruneisen-debye-and-einstein" data-equation-target="../equations/thermal/mie-gruneisen/"></span><span id="symbols-parameters-and-units" data-equation-target="../equations/thermal/mie-gruneisen/#symbols-parameters-and-units"></span><span id="gruneisen-parameter-and-debye-temperature" data-equation-target="../equations/thermal/mie-gruneisen/#gruneisen-parameter-and-debye-temperature"></span><span id="debye-energy-and-heat-capacity" data-equation-target="../equations/thermal/mie-gruneisen/#debye-energy-and-heat-capacity"></span><span id="thermal-pressure-and-its-reference-state" data-equation-target="../equations/thermal/mie-gruneisen/#thermal-pressure-and-its-reference-state"></span><span id="thermoelastic-derivatives" data-equation-target="../equations/thermal/mie-gruneisen/#thermoelastic-derivatives"></span><span id="entropy-free-energy-and-normalization" data-equation-target="../equations/thermal/mie-gruneisen/#entropy-free-energy-and-normalization"></span><span id="worked-mgd-calculation" data-equation-target="../equations/thermal/mie-gruneisen/#worked-mgd-calculation"></span><span id="einstein-alternative" data-equation-target="../equations/thermal/mie-gruneisen/#einstein-alternative"></span>[Mie-Gruneisen Debye and Einstein](equations/thermal/mie-gruneisen.md)

- <span id="variable-exponent-debye-temperature-law" data-equation-target="../equations/thermal/variable-exponent-debye/"></span>[Variable-exponent Debye-temperature law](equations/thermal/variable-exponent-debye.md)

- <span id="tange-2009-mgo-thermal-model" data-equation-target="../equations/thermal/tange-2009/"></span>[Tange 2009 MgO thermal model](equations/thermal/tange-2009.md)

- <span id="asymptotic-power-law-debye-with-t2-excess" data-equation-target="../equations/thermal/asymptotic-debye-excess/"></span><span id="asymptotic-power-law-debye-with-t2-excess_1" data-equation-target="../equations/thermal/asymptotic-debye-excess/"></span>[Asymptotic-power-law Debye with T2 excess](equations/thermal/asymptotic-debye-excess.md)

- <span id="sound-velocity-quasi-debye-helmholtz-eos" data-equation-target="../equations/thermal/sound-velocity-debye/"></span><span id="sound-velocity-quasi-debye-helmholtz-eos_1" data-equation-target="../equations/thermal/sound-velocity-debye/"></span>[Sound-velocity quasi-Debye Helmholtz EOS](equations/thermal/sound-velocity-debye.md)

- <span id="linear-thermal-pressure" data-equation-target="../equations/thermal/linear-thermal-pressure/"></span>[Linear thermal pressure](equations/thermal/linear-thermal-pressure.md)

- <span id="second-order-temperature-compression-thermal-pressure" data-equation-target="../equations/thermal/second-order-taylor/"></span>[Second-order temperature-compression thermal pressure](equations/thermal/second-order-taylor.md)

- <span id="logarithmic-volume-linear-thermal-pressure" data-equation-target="../equations/thermal/log-volume-thermal-pressure/"></span>[Logarithmic-volume linear thermal pressure](equations/thermal/log-volume-thermal-pressure.md)

- <span id="holland-powell-thermal-pressure" data-equation-target="../equations/thermal/holland-powell/"></span><span id="thermal-modified-tait" data-equation-target="../equations/thermal/holland-powell/"></span>[Holland--Powell thermal pressure](equations/thermal/holland-powell.md)

- <span id="dewaele-2006-hcp-fe-thermal-pressure-scale" data-equation-target="../equations/thermal/dewaele-2006/"></span>[Dewaele 2006 hcp-Fe thermal pressure scale](equations/thermal/dewaele-2006.md)

- <span id="dorogokupets-oganov-2007-four-oscillator-model" data-equation-target="../equations/thermal/dorogokupets-oganov-2007/"></span>[Dorogokupets--Oganov 2007 four-oscillator model](equations/thermal/dorogokupets-oganov-2007.md)

- <span id="multi-oscillator-gruneisen-thermal-pressure" data-equation-target="../equations/thermal/multi-oscillator/"></span><span id="paper-versus-spreadsheet" data-equation-target="../equations/thermal/multi-oscillator/#paper-versus-spreadsheet"></span>[Multi-oscillator Gruneisen thermal pressure](equations/thermal/multi-oscillator.md)

- <span id="monatomic-debye-and-anharmonic-helmholtz-energy" data-equation-target="../equations/thermal/debye-anharmonic/"></span>[Monatomic Debye and anharmonic Helmholtz energy](equations/thermal/debye-anharmonic.md)

## Numerical inversion

`volume(P)` and `volume(P, T)` solve the applicable pressure equation with a
bracketed scalar root for every broadcast point. `temperature(P, V)` solves

\[
\Delta P_{\mathrm{th}}(V,T)=P-P_{\mathrm{ref}}(V)
\]

on the positive-temperature branch nearest $T_r$. The distinct two-volume
DAC inversion and forward `volume_with_dac_confinement(P_cold,T,f_dac)` solve
are documented under
[Diamond-anvil-cell thermal-pressure contribution](dac-thermal-pressure.md).
