# Yokoo 2009 Au evidence needed to reopen the PVT investigation

Reviewed 2026-10-12. **The numerical investigation remains closed as of
2026-10-09, pending author inputs.** The printed analytical Au PVT remains
unreproduced. This follow-up concerns pressure evaluation for
[Yokoo et al. (2009)](https://doi.org/10.1103/PhysRevB.80.104114), not recovery
of the original optimization.

The current `gold_yokoo_2009_pvt_reconstruction` record is nondefault,
requires explicit selection and acceptance, and retains
`scientific_validation.status="not_reproduced"`. Its altered phonon
coefficients reconstruct derived Table III outputs; they are not recovered
author parameters. In particular, gamma0 changes from 2.96 to 2.923264661.
The RMS/maximum residuals of 0.008923/0.021133 GPa on 156 unmarked states
do not establish published-equation reproduction. See the
[Au reconstruction report](yokoo-2009-gold.md) and
[archived reconstruction](../data/yokoo-2009-gold-thermal-reconstruction.json).

## Established mismatch tests

The [closed convention audit](yokoo-2009-pressure-conventions.md) and its
[numerical evidence](../data/yokoo-2009-pressure-conventions.json) establish:

- At fixed volume, `P(V,T) - P(V,0)` cancels the cold curve, any volume-fixed
  zero-point pressure, and any constant electronic reference offset.
  Printed Au thermal increments differ from Table III by RMS **0.119080 GPa**
  and maximum **0.209408 GPa** over 135 unmarked warm states.
- Conditional half-last-digit phonon adjustments still leave a maximum
  difference of **0.124643 GPa**. At `V/V0=1`, a and b cancel; the 1500 K
  interval gap remains **0.09125 GPa**, including the audit's density,
  electronic-pressure and warm-minus-cold table rounding bounds. These
  bounds are hypotheses about printed digits, not source uncertainties.
- Every compared temperature is an exact electronic-table node.
  Between-node interpolation cannot explain these residuals. Independent
  Debye integration agrees within **10^-10 GPa**.
- Substituting Vc for ambient V0 in the phonon law leaves a maximum
  increment difference of **0.490218 GPa**. The previously tested independent
  Pytheos electronic polynomial leaves **0.209111 GPa** and is not author
  confirmation. Neither alternative is adopted.

## Pressure evaluation evidence checklist

- [ ] **Adopted phonon parameters and definitions.** Obtain the actual
  gamma0, a, b and theta0 used for Table III, with unrounded values, units,
  gamma(V), theta(V), Debye-energy definition, atom/mole normalization and
  reference volume/temperature. Clarify any zero-point or subtraction term.
  The specification must address the ambient 1500 K mismatch and compressed
  isochores. Extra digits confined to the tested rounding bounds are
  insufficient; a/b changes alone cannot resolve the ambient test.
- [ ] **Adopted electronic-pressure routine.** Obtain the author code,
  spreadsheet formulas, or a complete numerical prescription for Pel(V,T),
  including the actual nodes or coefficients, units, volume dependence,
  interpolation and domain/extrapolation rules, pressure reference and
  any subtraction. Identify whether the corrected Tsuchiya–Kawamura table
  was used. A proposed difference must affect the tested node values or
  their volume dependence to address the increment mismatch; a constant
  offset or between-node interpolation choice alone cannot do so.
- [ ] **Vc separately from V0 for absolute pressures.** Obtain the adopted
  0 K cold reference Vc and ambient 300 K V0, their cell/molar/density
  conversion, and the cold BM3 coefficients actually used. Confirm that
  Table III ratios use ambient V0 and specify how each pressure term uses
  volume. The reconstructed Vc/V0 is a fit to source outputs, not an
  independently recovered input. Vc is needed for absolute-pressure replay,
  but cannot repair the fixed-volume increment discrepancy.
- [ ] **Any source correction and author output checkpoints.** Identify
  exact corrections to printed equations, Table II coefficients or Table
  III pressures, with author/publisher provenance and superseded values.
  Supply pressure-term checkpoints at 0 and 1500 K at `V/V0=1`, and at 0
  and 3000 K at `V/V0=0.6`, preferably the complete generated Table III.
  Preserve rounding and first-liquid annotations. Numerical checkpoints
  without the adopted evaluation prescription do not recover the analytical
  model; a fit to those outputs does not satisfy this checklist.

The **51 corrected Tsuchiya–Kawamura electronic-pressure nodes are already
recovered**. The known 2003 upstream erratum changes Au nodes at 3100–3900 K,
outside Table III's range. Missing electronic energy, entropy, heat capacity,
original fitting weights or covariance are not PVT evaluation obstacles;
caloric and original-fit recovery are separate tasks.

## Evidence review before any new calculation

Reopening requires provenance-confirmed author code/spreadsheet, an equivalent
complete pressure specification, or a source correction that addresses the
established mismatch. First record its origin/version, relevant formulas and
parameters, affected checkpoints, completeness and relationship to the four
items above. Assess those conditions before running any new pressure
calculation. Partial evidence remains a documented gap; independent software
or a successful refit alone does not justify reopening.

On 2026-10-12, the accessible
[APS article record](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.80.104114)
showed no correction or supplemental-material link. Targeted DOI/title
searches for corrections, supplements, code and spreadsheets, including APS
and Ehime University sources, identified no new author artifact. The article
text on the landing page requires subscription access; this limited check
does not prove that unpublished or unindexed author material is absent.
**No newly identified evidence satisfies the reopening conditions.**

This follow-up adds one report only. No EOS calculations, parameter changes,
Pt work, outreach or global regeneration were performed. The initially clean
checkout acquired concurrent changes to `gold.eosmat` and the shared
`pressure-calibrations.json`; those changes belong to other work and were
preserved. This report generates no shared files requiring reconciliation.
Any later reconciliation of the concurrently changed calibration file belongs
with its owning follow-up work.
