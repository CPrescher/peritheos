# Myhill (2022): anisotropic periclase EOS

Primary source: R. Myhill, “An anisotropic equation of state for high-pressure, high-temperature applications,” *Geophysical Journal International* **231**, 230–242 (2022), <https://doi.org/10.1093/gji/ggac180>.

Myhill simultaneously optimized a new anisotropic periclase model against the ab-initio isentropic stiffness data of Karki et al. (2000). Table 2 gives the scalar volumetric parameters `V0=11.33682 cm3/mol`, `K0=160.20484 GPa`, `K0'=4.17945`, `Theta0=767.09770 K`, `gamma0=1.53632`, and `q0=1.01212`, plus 18 anisotropic coefficients. The scalar part uses the finite-strain thermal formulation of Stixrude and Lithgow-Bertelloni (2005); it is not identified as a standalone BM3 isotherm. Peritheos currently lacks the complete SLB2005 anisotropic tensor EOS and cannot faithfully execute the published coefficient set. Mapping only the first three values to BM3 would discard the fitted thermal and tensor terms and change the reported model.

| Source row | Candidate | Disposition | Rationale |
|---:|---|---|---|
| 1152 | `litcurate_14aa941c91c0120d` | HOLD | Genuine source-owned simultaneous refit, but the complete SLB2005 anisotropic model is unsupported; the LitCurate `unknown` model must not be coerced to BM3. |

Result: **0 production records and 1 held source row**. This paper is a concrete implementation candidate for a future anisotropic/SLB finite-strain EOS family; its full Table 2 coefficient set should be preserved together when that model exists.
