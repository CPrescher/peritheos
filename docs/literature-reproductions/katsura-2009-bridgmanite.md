# Katsura et al. (2009): MgSiO3 bridgmanite P-V-T audit

## Primary source and corrigendum

This audit covers the MgSiO3 bridgmanite candidates traceable to Katsura et
al., *P-V-T relations of MgSiO3 perovskite determined by in situ X-ray
diffraction using a large-volume high-pressure apparatus*, *Geophysical
Research Letters* **36**, L01305 (2009),
[doi:10.1029/2008GL035658](https://doi.org/10.1029/2008GL035658). LitCurate was
used only for candidate discovery; the primary article and formal correction
are the scientific authorities.

The correction, *Geophysical Research Letters* **36**, L16309,
[doi:10.1029/2009GL039318](https://doi.org/10.1029/2009GL039318), is essential.
It states that the original Table 1 contained wrong numbers because of an
editorial mistake and replaces the final column with bridgmanite `V/V0`. The
corrected table also changes the printed pressure uncertainties. The correction
explicitly says that the original EOS parameters do not change. Peritheos
therefore bundles all 89 rows from the corrected table and does not preserve the
superseded numerical table as fit data.

## Material and experiment

The phase is pure MgSiO3 orthorhombic perovskite (bridgmanite), not the generic
`(Mg,Fe)SiO3` assigned by a later Japanese review. Below 30 GPa, the starting
material was a powdered forsterite plus MgO mixture; above 30 GPa it was mostly
a sintered perovskite plus MgO mixture converted in a separate run. The source
used tungsten-carbide (WC) and sintered-diamond (SD) Kawai-type multi-anvil
assemblies at SPring-8 BL04B1, with energy-dispersive X-ray diffraction.

The measured ambient lattice is `a0=4.7769(2)`, `b0=4.9298(2)`, and
`c0=6.8956(3) A`. Their product is `V0=162.3855988669 A^3` for the conventional
Pbnm cell with `Z=4`; independent-error propagation gives `0.0118124861 A^3`.
This exact pure-MgSiO3/Pbnm identity matches the existing `bridgmanite`
material.

## Pressure scale and corrected data

MgO was the internal pressure marker. The authors calculated each pressure
with both the Matsui et al. (2000) and Speziale et al. (2001) MgO scales, then
reported their average. The corrected Table 1 footnote says the pressure errors
include the difference between those scales. The complete bundled dataset
contains 36 WC and 53 SD observations, spanning 18.9-52.6 GPa and 300-2500 K,
with corrected MgO `V/V0`, mean pressure, and bridgmanite `V/V0` values and
uncertainties.

Both exact historical thermal MgO implementations are not bundled, so the
reduced pressure cannot yet be recalculated from first principles. The
published averaged coordinate and its scale-difference uncertainty are
preserved rather than recast as a one-standard-deviation experimental error.

## Source parameterizations

Section 4.1 fits the 300 K reference isotherm with a third-order
Birch--Murnaghan equation. It explicitly reports two fits to the same data:

| Role | V0 (A^3/cell) | K0 (GPa) | K0' | Peritheos record |
| --- | ---: | ---: | ---: | --- |
| Preferred fit | 162.3855988669, ambient lattice fixed | 256(2), fitted | 3.8(2), fitted | `bridgmanite_katsura_2009_bm3` |
| Brillouin-comparison sensitivity | same ambient lattice fixed | 253, assumed | 4.1(2), fitted | `bridgmanite_katsura_2009_bm3_k0_fixed` |

The second row is source-generated, not a LitCurate reconstruction. It remains
nonpreferred and is explicitly marked as a same-data sensitivity fit.

The paper also defines a complete Mie--Gruneisen--Debye fit with fixed
`theta0=1030 K`, `gamma0=2.6(1)`, and `q=1.7(1)`, and separately reports
`(dK_T/dT)_P=-0.035(2) GPa/K`, Anderson--Gruneisen `delta_T=6.5(5)`, and an
ambient thermal-expansion law. These quantities are documented on the
preferred static record but are not combined into a misleading extra thermal
record: the LitCurate candidates request only the two BM3 reference fits, and
the source's additional high-temperature formulation is not a single exact
match for Peritheos's current thermal representation.

## Numerical reproduction

Run:

```bash
.venv/bin/python scripts/reproduce_katsura_2009_bridgmanite.py
```

Nine corrected rows were collected at 300-308 K. Against those rounded values,
the preferred published curve has a `0.320208 GPa` pressure RMSE and a
`0.569935 GPa` maximum absolute pressure residual. With `K0'=3.8` fixed, an
independent unweighted conditional fit gives `K0=255.507510 GPa`, inside the
published 2 GPa interval.

Allowing both strongly correlated coefficients to vary against the four-decimal
`V/V0` values gives `K0=261.437273 GPa` and `K0'=3.404892`, with a
`0.246790 GPa` pressure RMSE. This rounded-table diagnostic is recorded rather
than hidden, but it does not replace the source fit: the article does not
publish fit weights, covariance, or unrounded corrected observations, and the
formal correction independently confirms that the EOS coefficients are
unchanged.

For the sensitivity fit, fixing `K0=253 GPa` and weighting by the corrected
pressure uncertainties gives `K0'=4.128233`, only 0.14 of the reported 0.2
interval from `4.1`. Its rounded published curve has a `0.469978 GPa` pressure
RMSE. Together, the direct curve checks and conditional refits establish parity
for both source-generated parameterizations without substituting a new fit.

## LitCurate disposition

| Candidate | Quantity represented | Decision |
| --- | --- | --- |
| `litcurate_210dfe29f03ef11e` | Katsura preferred `K0=256(2) GPa`, `K0'=3.8(2)` | **Accepted** as `bridgmanite_katsura_2009_bm3`. |
| `litcurate_5d0e3a79a7f88494` | Katsura same-data fixed-`K0=253 GPa`, `K0'=4.1(2)` | **Accepted / nonpreferred sensitivity** as `bridgmanite_katsura_2009_bm3_k0_fixed`. |
| `litcurate_9bfd5533244f8b1e` | Fiquet et al. (2000) `K0=253(9) GPa`, `K0'=3.9(2)` comparison | **Rejected as a Katsura-source record / citation trace only.** The underlying primary paper must control its own data and pressure scale. |
| `litcurate_9ecfcfc2c0c6cee3` | Sinogeikin et al. (2004) adiabatic `K_S0=253(3) GPa` comparison | **Rejected as incomplete and not a volumetric Katsura EOS.** It is an acoustic result from another primary paper. |
| `litcurate_f2803f2c9b31c424` | Jackson et al. (2005) adiabatic `K'_S0=3.7(3)` comparison | **Rejected as incomplete and not a volumetric Katsura EOS.** It is an acoustic result from another primary paper. |
| `litcurate_3b7821d637ffbb33` | Later Japanese review repeats Katsura preferred fit and labels composition `(Mg,Fe)SiO3` | **Rejected as a duplicate citation trace.** The primary source establishes pure MgSiO3. |
| `litcurate_5c0608da3c3eecdb` | Later Japanese review repeats Katsura fixed-K0 sensitivity fit | **Rejected as a duplicate citation trace.** It maps to the accepted primary-source sensitivity record. |

The exact production count from this audit is therefore **two** EOS records.

## Zotero status

The Zotero local API and connector were running during the audit. An exact DOI
search for `10.1029/2008GL035658` returned no existing item. This scoped audit
did not import the paper.
