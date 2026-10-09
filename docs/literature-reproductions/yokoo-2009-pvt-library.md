# Yokoo (2009): diagnostic Au and Pt PVT reconstructions

Two **derived**, nondefault thermal records expose the bounded pressure
reconstructions through the normal Peritheos catalog:

- `gold_yokoo_2009_pvt_reconstruction`
- `platinum_yokoo_2009_pvt_reconstruction`

Both have scientific validation status **`not_reproduced`** and
`catalog_access="explicit_selection"`. Numerical checks of these diagnostics
do not validate the published analytical Yokoo EOS. Canonical material loading
requires `require_primary_validation=False` to accept them explicitly. The
[pressure-convention audit](yokoo-2009-pressure-conventions.md) records the
remaining discrepancy without fitted empirical pressure corrections.

The Au/Pt published-PVT investigation was **closed for now on 2026-10-09,
pending author inputs**. The
[closure and reopening requirements](yokoo-2009-pressure-conventions.md#remaining-evidence-needed) identify the
source calculation or complete specification needed to resolve the mismatch.

They preserve the [published Au branches](yokoo-2009-gold.md) and
[published Pt branch](yokoo-2009-platinum.md). Their provenance explicitly
credits Yokoo et al. for the equations and calculated tables, Tsuchiya and
Kawamura for the electronic calculations, and the Peritheos reconstruction
for derived coefficients and numerical interpolation choices. They are not
the original author optimization or independently validated experimental
parameter estimates.

```python
from peritheos import get_eos_record

pt = get_eos_record("platinum_yokoo_2009_pvt_reconstruction")
p = pt.pressure(48.44, 1500)       # Å³/four-atom cell, K -> GPa
v = pt.volume(p, 1500)            # -> 48.44 Å³/cell
t = pt.eos.temperature(p, 48.44 * pt.volume_scale)  # -> 1500 K
```

This worked state has V/V0=0.8 and gives **118.330845 GPa**, compared with
**118.33 GPa** in Table V. Internally, `volume_scale=N_A*10^-25/4` converts
cell Å³ to J/bar/mol of atoms. Thus 48.44 Å³ is about 0.7293 J/bar/mol,
or 7.293 cm³/mol. Input pressure is always GPa.

## Pressure equation and reference volumes

`AsymptoticDebyeTabulatedPressure` evaluates

`P(V,T)=Pc_BM3(V;Vc)+gamma(V) E_D(V,T)/Vm + Pel(T) + deltaP_reconstruction(T)`.

The library converts `gamma E_D/Vm` to GPa using its molar J/bar/mol
convention. `Pc` uses the **cold 0 K Vc**. The phonon ratio is instead
`r=V/V0`, with ambient 300 K `V0=Vc/cold_volume_ratio`:

`gamma(r)=gamma0 [1+a(r^b−1)]`,

`theta(r)=theta0 r^[-(1−a)gamma0] exp[-(gamma(r)−gamma0)/b]`,

`E_D=3 R T D3(theta/T)` per mole of metal atoms.

The phonon energy is zero at 0 K without an added zero-point term. Both
phonon and electronic pressure are **absolute**; neither is subtracted at
300 K in the total equation. For DAC confinement, the thermal pressure
*increment* subtracts the full 300 K contribution. Its inverse uses the
full modeled 300 K isotherm, not the 0 K cold curve.

All 51 electronic pressure nodes are retained as immutable constructor
configuration. Au uses the corrected column; the erratum changes no Pt
values. Linear interpolation is explicit and is not an author-verified
prescription. Evaluations and inversions are restricted to **0.6≤V/V0≤1**
and **0≤T≤3000 K**, without extrapolation. These numerical bounds do not
certify solid-phase stability throughout the rectangle. The source's six
Au and two Pt first-liquid markers are retained in provenance and excluded
from the primary reconstructions. Au's six blank cells remain absent from
the source CSV; no observations are manufactured there.

## Au reconstruction

The [Au reconstruction report](../data/yokoo-2009-gold-thermal-reconstruction.json)
fixes K0=180 GPa and theta0=170 K and refines Vc/V0, K0′, gamma0, a and b
against 156 unmarked Table III outputs, with equal pressure weights.
`deltaP_reconstruction` is identically zero. The diagnostic coefficients are

| Parameter | Au |
|---|---:|
| Ambient V0 (Å³/cell) | 67.72 |
| Cold Vc (Å³/cell) | 67.054862283 |
| Cold K0 (GPa) | 180 |
| Cold K0′ | 5.610130908 |
| gamma0 | 2.923264661 |
| a | 0.447532583 |
| b | 4.215495747 |
| theta0 (K) | 170 |

The RMS/max differences are **0.008923/0.021133 GPa** on 156 unmarked
states and **0.009076/0.018672 GPa** on 74 withheld states from alternate
complete isochores. The six first-liquid states are assessed separately.
The gamma0 shift exceeds decimal rounding, and the exact author's
normalization remains unresolved. These are derived reconstruction
coefficients, not substituted published values.

## Pt reconstruction and explicit residual correction

The printed cold branch, phonon coefficients and recovered electronic nodes
do not reproduce Table V exactly. Even refining the same five quantities as
for Au leaves RMS/max discrepancies of **0.08398/0.17244 GPa**. That
diagnostic remains archived; it is not silently discarded.

The diagnostic pressure reconstruction fits Vc to the selected 0 K rows,
fixing K0=288.4 GPa and K0′=5.05. It then fixes gamma0=2.63 and theta0=230 K,
refines a/b, and profiles a separate pressure correction at each source
temperature: `deltaP(T_j)=mean[P_table−Pc−Pph−Pel]` over the selected
volume states. Its 0 K value is fixed to zero. All stages use equal pressure
weights; withheld and first-liquid states enter neither the cold fit nor
the correction averages.

| Parameter | Pt |
|---|---:|
| Ambient V0 (Å³/cell) | 60.55 |
| Cold Vc (Å³/cell) | 60.175801510 |
| Cold K0 (GPa) | 288.4 |
| Cold K0′ | 5.05 |
| gamma0 | 2.63 |
| a | 0.385056641 |
| b | 5.189776556 |
| theta0 (K) | 230 |

The separate residual table is approximately:

| T (K) | Derived deltaP (GPa) |
|---:|---:|
| 0 | 0 |
| 300 | 0.118625 |
| 500 | 0.186512 |
| 1000 | 0.287514 |
| 1500 | 0.313714 |
| 2000 | 0.281877 |
| 2500 | 0.238428 |
| 3000 | 0.200061 |

This correction is interpolated linearly, stored separately from the raw
electronic nodes, and **not identified as electronic pressure**. Its
physical origin is unestablished. The additional flexibility reconstructs
the source pressure outputs; it does not establish new phonon/electronic
physics or recover electronic energies.

The [Pt reconstruction report](../data/yokoo-2009-platinum-thermal-reconstruction.json)
archives the complete objective, coefficients, correction nodes, residuals,
input hashes, multiple starts, density-conversion sensitivity and the less
accurate diagnostics. On 166 unmarked states, RMS/max differences are
**0.002919/0.005940 GPa**. The 79 withheld states give
**0.003224/0.005992 GPa**. The two first-liquid markers remain withheld.
These are interpolation checks on derived source outputs, not independent
experimental validation. Agreement is close to table precision, with some
residuals slightly larger than its 0.005 GPa half-step.

## Library verification and limits

Run:

```bash
.venv/bin/python -m scripts.fit_yokoo_2009_gold_thermal --check
.venv/bin/python -m scripts.fit_yokoo_2009_platinum_thermal --check
.venv/bin/python -m scripts.register_yokoo_2009_pvt --check
```

The [library validation report](../data/yokoo-2009-pvt-library-validation.json)
compares all 330 source states with independent SI pressure calculations.
Maximum pressure differences are below 10^-9 GPa; volume inversion errors
are below 10^-8 Å³/cell and temperature inversion errors below 10^-6 K.
Explicitly accepted canonical material round trips preserve the configurations.
The legacy snapshot importer rejects the unvalidated composition. Tests
also check intermediate-state pressure monotonicity, integrated gamma/theta
consistency, bounded inversion failures and DAC confinement round trips.

Parameter errors and covariance remain null: residuals of rounded model
outputs are not experimental uncertainties. Entropy, heat capacity and the
original shock-temperature calculation are outside this pressure-only
reconstruction. Both records are marked **diagnostic source-output reconstructions**
in the refit ledger, with `published_analytical_pvt_reproduced=False`; the
published 300 K records retain their own statuses.
