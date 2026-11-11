# Sound-velocity quasi-Debye Helmholtz EOS

<nav class="eos-context" aria-label="Equation location">
<a href="../../../models/thermal/">Thermal equations</a>
<span aria-hidden="true">›</span>
<span>Sound-velocity quasi-Debye Helmholtz EOS</span>
</nav>

[Notation and reference state](../../equation-reference.md#notation-and-reference-state) · [Units](../../units.md) · [Numerical inversion](../../equation-reference.md#numerical-inversion)

`SoundVelocityDebyeHelmholtz` implements the quasi-Debye construction used by
Luo et al. (2023), Equations (1)--(6) and Appendix A. Independent linear
room-temperature regressions

\[
C_l^{\rm fit}=a_l+b_l\rho,\qquad C_s^{\rm fit}=a_s+b_s\rho
\]

define $\tau=C_l^{\rm fit}/C_s^{\rm fit}$ and
$\mu=(\tau^2-2)/[2(\tau^2-1)]$. The wrapped cold curve supplies
$B_s=-V\,\partial P_c/\partial V$; the model velocities are then

\[
G=\frac{3(1-2\mu)}{2(1+\mu)}B_s,\quad
C_l^2=\frac{B_s+4G/3}{\rho},\quad C_s^2=\frac{G}{\rho},
\]

\[
\frac{3}{C_{\rm eff}^3}=\frac{1}{C_l^3}+\frac{2}{C_s^3},\qquad
\Theta=\frac{\hbar}{k_B}\left(6\pi^2\frac{N}{V}\right)^{1/3}C_{\rm eff}.
\]

With $E_{\rm ph}=3nRTD_3(\Theta/T)$ and
$\gamma=-\partial\ln\Theta/\partial\ln V$, the absolute vibrational pressure
is $P_{\rm ph}=\gamma E_{\rm ph}/V$. The implementation exposes vibrational
energy, entropy, Helmholtz energy, heat capacity, longitudinal/shear/effective
velocities, cold compression energy, and the Rankine--Hugoniot temperature
solver. It omits zero-point energy, following Luo's integral of heat capacity
from zero temperature. Constructor volumes are molar J bar$^{-1}$ mol$^{-1}$,
velocities are km/s, and density is g cm$^{-3}$.

The four velocity-regression coefficients are required inputs. Luo et al. name
Kono et al. (2010) as their source but do not print the fitted coefficients, so
the catalog's published Luo pressure scale continues to execute the Appendix-B
Taylor surrogate. Supplying coefficients from another experiment creates an
explicit sensitivity model, not a source-author reconstruction.

---

**End of Sound-velocity quasi-Debye Helmholtz EOS.**

[All thermal models](../../models/thermal.md) · [Equation index](../../equation-reference.md)
