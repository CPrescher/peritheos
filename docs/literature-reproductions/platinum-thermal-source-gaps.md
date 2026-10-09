# Platinum thermal scales: Matsui (2009) and Yokoo (2009) source gaps

## Investigation begun 2026-10-08, updated 2026-10-09

The [Fei (2007) thermal-scale audit](fei-2007-platinum.md) adds the complete
published Vinet-MGD pressure equation to the canonical Pt catalog. Both
Matsui and Yokoo full texts were supplied on 2026-10-09; the resolved access
gaps and remaining model dependencies are described below.

### Matsui (2009): recovered pressure model, original-fit gaps remain

The primary Matsui paper and both supplied Tsuchiya-Kawamura documents were
read and visually checked. The [Matsui audit](matsui-2009-platinum.md)
registers `platinum_matsui_2009_vinet_mgd_electronic` alongside the unchanged
`platinum_matsui_2009_vinet_300k`. The recovered 51 electronic pressure nodes
close the missing numeric pressure term. The 2003 erratum corrects only Au.
Raw Pt pressure at 300 K is subtracted once, using the four-atom fcc cell
conversion and integrated-q Debye law.

Linear interpolation between the upstream 100 K nodes is an explicit
numerical implementation choice; the authors supply no interpolation rule.
No extrapolation outside the table is allowed. All 35 Table III outputs agree
within 0.0071 GPa and 12 Table I calculated pressures within 0.0086 GPa.
These verify the rounded pressure parameterization, not hidden precision or
the original optimization. Static validation, derived output, assessed
expansion and shock observations remain separate. Marsh inputs, original
Hugoniot/expansion sampling, weights and the joint-fit objective are still
missing; original author-fit reproduction remains unavailable.

### Yokoo et al. (2009): supplied primary, accepted Vinet fit and thermal dependency

The primary PDF supplied on 2026-10-09 has now been inspected visually in
full. The [Yokoo audit](yokoo-2009-platinum.md) preserves Table IV, all 168 Pt
Table V states, the absolute cold/phonon/electronic decomposition, and the
recovered Holmes shock and paired Dewaele comparison inputs.

The published fitted 300 K Vinet branch is registered as
`platinum_yokoo_2009_vinet_300k`, with K0=276.4 GPa, K0'=5.48 and the
explicit density-derived V0=60.55 Å³ per four-atom fcc cell. It differs from
the full thermal model's 300 K Table V curve by up to 6.0093 GPa over the
printed grid and is labeled as a separate fitted reference isotherm.
Sakai (2018)'s named Pt Vinet equation reference now resolves to this branch;
its original paired Pt/Re volumes remain missing, preventing historical
pressure recalculation.

Both metals now have bounded executable
[derived PVT reconstructions](yokoo-2009-pvt-library.md), using the recovered
electronic-pressure table and separate inferred 0 K cold volumes. Au agrees
with unmarked Table III outputs within 0.0211 GPa; Pt agrees with unmarked
Table V outputs within 0.00594 GPa using an explicitly derived residual
pressure correction. The Pt correction is stored separately from electronic
pressure and has no established physical interpretation. These records
reconstruct calculated outputs, not the original author joint fit.
Both have scientific status `not_reproduced` and require explicit acceptance
for canonical loading. The [pressure-convention audit](yokoo-2009-pressure-conventions.md)
shows that the printed analytical PVT remains unreproduced even after the
stated last-digit rounding allowances, independently of fitting procedures
or caloric inputs.

## Next source inputs

The shared Tsuchiya-Kawamura pressure table and Au-only erratum are now
recovered. They support Matsui's pressure model and the derived Yokoo PVT
records. Electronic **energy** is required for the original shock/caloric
calculation, not bounded PVT evaluation. Exact source cold-volume/phonon
normalization and electronic interpolation remain unresolved.
Original shock/ambient-property selections, numerical input rows, weights
and objective details are needed for independent author-fit replay.
Validated parameterizations and derived source output remain distinct from
joint-fit reproduction.
