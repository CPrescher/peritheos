# Kubo et al. (2000) aluminous bridgmanite compression audit

## Scope and primary source

This audit covers LitCurate candidates `litcurate_1b9444bb980ed86d` and
`litcurate_7fa232488f133567`, both discovered under DOI
`10.2183/pjab.76.103`. Scientific authority comes from Kubo, Yagi, Ono, and
Akaogi, “Compressibility of Mg0.9Al0.2Si0.9O3 perovskite,” *Proceedings of the
Japan Academy, Series B* 76, 103-107 (2000). The final article and PDF are
[openly available from J-STAGE](https://www.jstage.jst.go.jp/article/pjab1977/76/8/76_8_103/_article).
The downloaded publisher PDF had SHA-256
`70fd8aac56e920c21442ed6884b24024e839f34b2b561f7f0b55aab178a3c214`.

## Material and experiment

The sintered specimen came from synthesis run en90G3. A glass containing 90
mol% MgSiO3 and 10 mol% Al2O3 was held at 27 GPa and 2000 degrees C for 83
minutes in a 6-8 multianvil apparatus. Microfocus diffraction showed a
single-phase orthorhombic perovskite with space group `Pbnm`; EPMA confirmed
`Mg0.9Al0.2Si0.9O3` within analytical error. The recovered grains were smaller
than 1 micrometre.

The conventional Pbnm cell has `Z=4`, inferred from the named structure and
the approximately 163.6 A3 cell. The source does not print `Z`, fractional
coordinates, or Al occupancies. Peritheos consequently uses the same pure
MgSiO3 positional-topology proxy as its existing bridgmanite material and
labels it explicitly. The paper's later equal-site Al discussion is
hypothetical and is not encoded as a measured site distribution.

Two room-temperature DAC decompression runs used a 4:1 methanol-ethanol
pressure medium. Mao, Xu, and Bell's 1986 ruby scale is the adopted pressure
coordinate. Five ku-series pressures were also calculated from the Anderson,
Isaak, and Yamamoto 1989 gold EOS; these are cross-checks, not fit pressures.
The exact observed range is 0-8.71 GPa, although the article summarizes the
study as extending to 9 GPa and says the samples were initially compressed to
9-10 GPa.

## Primary observations

Table I contains ten complete rows: eight finite-pressure diffraction patterns
and one atmospheric reference from each run. The CSV transcribes every ruby
pressure, diagnostic Au pressure, lattice constant, conventional-cell volume,
axial ratio, and printed last-digit error. No digitization was necessary.

The two atmospheric cells are mutually consistent:

| Run | a (A) | b (A) | c (A) | V (A3) |
|---|---:|---:|---:|---:|
| ap107 | 4.774(4) | 4.941(4) | 6.935(6) | 163.6(4) |
| ku117 | 4.776(5) | 4.939(5) | 6.933(8) | 163.5(5) |

The article says `a/a0`, `b/b0`, `c/c0`, and `V/V0` were fixed to unity at
atmospheric pressure. Thus the published reduction used each run's separate
ambient reference rather than fitting one common `V0`. An executable EOS needs
one volume reference; the catalog adopts the more precise ap value of
163.6(4) A3 and retains both values in the dataset.

## Equation and parameterizations

The paper names a third-order Birch-Murnaghan EOS but does not print its
algebraic form. Peritheos uses the standard Eulerian convention:

\[
P(V)=\frac{3K_0}{2}
\left[\left(\frac{V_0}{V}\right)^{7/3}
-\left(\frac{V_0}{V}\right)^{5/3}\right]
\left\{1+\frac{3}{4}(K'_0-4)
\left[\left(\frac{V_0}{V}\right)^{2/3}-1\right]\right\}.
\]

Both source-reported reductions are retained:

| Catalog status | V0 (A3) | K0 (GPa) | K'0 | Source interpretation |
|---|---:|---:|---:|---|
| Default/adopted | 163.6(4) | 225.5(1.2) | 4 fixed | Abstract result and Figure 2 curve |
| Secondary/free | 163.6(4) | 215.4(4.4) | 7.2(1.4) | Valid fit, derivative explicitly described as inaccurate |

The experiment reached less than 4% compression. The authors therefore state
that the free derivative is not accurately determined and adopt `K'0=4`,
consistent with MgSiO3 perovskite. The abstract reports only this fixed fit,
Figure 2 plots it, and Table II rounds it to `226(1) GPa`. The unconstrained
fit is not a failed extraction, but it must not be presented as the preferred
physical result. Because the source calls both reductions third order, the
fixed record remains BM3 with an explicit fixed parameter; it is algebraically
the BM2 special case.

The source does not state a confidence convention for its EOS `+/-` values and
does not publish a covariance matrix. No confidence level or covariance is
invented.

## Numerical reproduction

Direct evaluation against the ten printed volumes, with the separate ap and ku
atmospheric normalizations described by the source, gives:

| Record | Pressure RMSE (GPa) | Maximum absolute residual (GPa) |
|---|---:|---:|
| Fixed `K'0=4` | 0.0834714544 | 0.1501695655 |
| Free `K'0` | 0.0591814632 | 0.1270385963 |

An unweighted pressure-space refit of the eight finite-pressure rows gives
`K0=226.86340115 GPa` for fixed `K'0=4`. This is 1.14 printed standard errors
above the published value. The volume column is printed to only 0.1 A3; using
products of the more finely printed lattice constants gives
`K0=224.83957212 GPa`, 0.55 printed standard errors below the published value.
The bracket supports coefficient parity within Table I's rounding.

The unconstrained refit gives `K0=215.61948074 GPa` and
`K'0=7.53074835`, respectively 0.05 and 0.24 published standard errors from
the printed coefficients. This numerical agreement does not remove the
source's scientific warning about the strong derivative-modulus tradeoff.

Run the deterministic check with:

```bash
.venv/bin/python scripts/reproduce_kubo_2000_aluminous_bridgmanite.py
```

## Calibration and remaining limitations

The ruby calibration is executable in Peritheos as `ruby_mao_1986`, and the
diagnostic gold reference exists as `gold_anderson_1989_bm3_1`. Exact raw-data
pressure recalculation is impossible because Kubo et al. print neither ruby R1
shifts nor Au lattice parameters. They report ruby and Au pressures agreeing
within error and a less than 0.08 GPa reflection-to-reflection Au pressure
spread, supporting hydrostatic conditions.

The final inclusion decision is:

- accept `mg09al02si09o3_bridgmanite_kubo_2000_bm3_1` as the default adopted
  fixed-derivative record;
- accept `mg09al02si09o3_bridgmanite_kubo_2000_bm3_2` as a non-default,
  explicitly poorly constrained free-derivative record; and
- do not treat the paper's comparison-table EOS values for other compositions
  as new Kubo et al. measurements.

No exact-composition duplicate was present. The existing
`mg095al010si095o3_bridgmanite` material is a distinct 5 mol% Al2O3 specimen
from Daniel et al. (2004).
