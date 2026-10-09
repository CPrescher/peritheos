# Yokoo PVT: audit of the remaining thermal-pressure discrepancy

The published analytical Au/Pt PVT parameterizations remain unreproduced.
The separate table-output reconstructions are retained as **diagnostics**,
with `scientific_validation.status="not_reproduced"`, nondefault selection
and explicit acceptance required for canonical material loading. Their
small fitted residuals do not validate them as replacements for the source
pressure standard. No published coefficients or source observations have
been changed.

## Direct source checks

Yokoo (2009), Eqs. 1–12, Tables II–V and their footnotes were rechecked
against the primary PDF. The cold BM3 uses Vc at 0 K; the phonon law uses
ambient V0. The 51 Tsuchiya–Kawamura (2002) Table I electronic-pressure
nodes were visually rechecked. Electronic pressure is absolute relative to
0 K. The supplied erratum affects Au alone.

The [machine-readable audit](../data/yokoo-2009-pressure-conventions.json)
compares each source pressure with the 0 K pressure at the **same volume**:

`P(V,T) − P(V,0) = Pph(V,T) − Pph(V,0) + Pel(T) − Pel(0)`.

Consequently, the entire cold curve and its volume normalization cancel.
A volume-fixed zero-point term cancels too. A constant electronic reference
offset, including subtraction of Pel(300 K), cancels. Changing those terms
cannot resolve a discrepancy in these increments.

The implementation evaluates the printed Debye equations using SI molar
volume, four atoms per fcc cell, and Gauss–Legendre quadrature. Independent
adaptive integration agrees within 10^-10 GPa. Every compared temperature
is an actual electronic-table node, so between-node interpolation cannot
explain the observed differences.

| Metal | Unmarked warm states | Printed increment RMS (GPa) | Maximum (GPa) |
|---|---:|---:|---:|
| Au | 135 | 0.119080 | 0.209408 |
| Pt | 145 | 0.286919 | 0.380666 |

These compare increments and therefore differ slightly from full-pressure
residuals using an inferred cold-volume fit. All 162 Au and 168 Pt source
states remain archived, with first-liquid states excluded from the reported
unmarked metrics and no blank Au cells fabricated.

## Conditional last-digit rounding check

The audit permits gamma0/a ±0.005, b ±0.05 and theta0 ±0.5 K, based on their
printed digits. These are an explicit rounding hypothesis, **not** source
uncertainties or a recovery of unrounded author coefficients. Source a/b
uncertainties are kept distinct from this numerical test.

Even the best joint adjustment within those bounds leaves maximum increment
differences of 0.124643 GPa for Au and 0.296241 GPa for Pt.

At V/V0=1, a/b cancel identically from both gamma and theta. Direct bounds
on gamma0/theta0 then provide a check independent of the optimizer. Including
the half-last-digit electronic-pressure interval (±0.005 GPa), density
interval (±0.005 Mg/m³) and warm-minus-cold table interval (±0.01 GPa), Pt
at 1500 K still has a **0.27557 GPa gap between the permitted intervals**.
Au at 1500 K retains a 0.09125 GPa gap. These are incompatibility checks
under the stated rounding assumptions, not experimental accuracy claims.

## Other conventions examined

Substituting cold-volume normalization into the phonon law, contrary to
Eq. 10's stated ambient reference, leaves maximum increment differences
of 0.4902 GPa for Au and 0.3504 GPa for Pt. It is not adopted.

An independent [Pytheos implementation](https://github.com/SHDShim/pytheos/tree/b524b3439dfacfee8fa46a3e0718d50520ac37c5)
was inspected at the pinned commit. Its Tsuchiya electronic-pressure
polynomial also fails to close the increment discrepancy when combined with
the printed phonon equations. This tests the electronic convention alone,
not Pytheos' separate 300 K BM3 reference curve. Its Au class explicitly
offers a changed room-temperature K0′ to improve table agreement.
Independent software is useful comparison evidence, not author confirmation;
none of its coefficients are substituted into the published Peritheos records.

## Remaining evidence needed

**Investigation closed for now on 2026-10-09, pending author inputs.** This
closes the current Au/Pt PVT audit with an unresolved reproduction outcome;
the published analytical EOS remains unreproduced. The derived records stay
nondefault, explicit-selection diagnostics with status `not_reproduced`.

Reopen the investigation when the authors' code or spreadsheet used to
generate Tables III and V becomes available, or an equivalent complete
specification or source correction resolves the pressure discrepancy. The
required pressure-evaluation specification is:

- Exact adopted phonon parameters (gamma0, a, b, theta0), their definitions
  and volume/reference conventions.
- The adopted electronic-pressure function or numerical table, including
  volume dependence, interpolation and pressure-reference conventions.
- Cold reference volume Vc at 0 K separately from ambient V0 at 300 K,
  to reproduce absolute pressure.
- Clarification of any correction to the printed equations, parameters or
  Tables III/V.

The existing electronic-pressure table is recovered. Extra digits alone
within the tested rounding bounds do not resolve the mismatch, and cold
reference volumes cannot fix the thermal-increment discrepancy. The cause
remains unidentified. No author outreach has been sent in this investigation.

The available printed phonon parameters and recovered electronic table do
not establish the exact adopted analytical thermal-pressure implementation.
Resolving the discrepancy requires the author calculation/code or enough
unrounded phonon coefficients and the actual adopted electronic-pressure
routine to explain the source outputs. The current evidence does not
identify a unique cause or establish a source error.

Electronic energy, entropy, heat capacity and original fitting weights are
**not PVT evaluation dependencies**. They are not reasons to reject a
complete pressure parameterization. Here the obstacle is a directly tested
pressure discrepancy. The empirical Pt correction and refined Au phonon
coefficients remain diagnostic reconstructions, not recovered author inputs.

Reproduce the audit with:

```bash
.venv/bin/python -m scripts.audit_yokoo_2009_pressure_conventions --check
```
