# Schoelmerich et al. (2020) stishovite audit

## Decision

No production EOS record is retained from Schoelmerich et al. (2020). The
paper reports a 300 K third-order Birch-Murnaghan fit with
`V0=46.5 A3` fixed, `K0=307(4) GPa`, and `K0'=4.66(15)`, but the
five shock-corrected 300 K fit states are shown only graphically.

Digitizing markers from the authors' fitted figure is not independent evidence
for reproducing the coefficients plotted in that same figure. The former
`stishovite_schoelmerich_2020_shock_300k_bm3` record and its derived
Figure 3 dataset were therefore removed from the production catalog.

## Retained primary data

The six rows of primary Table 1 remain packaged as the audit fixture
`peritheos/data/datasets/stishovite-schoelmerich-2020-table1-shock.csv`
under the article's CC BY 4.0 license. The fixture is not registered as an EOS
fit dataset. It preserves one ambient state and five shock experiments,
including particle and shock velocity, Rankine-Hugoniot pressure and energy,
cell volume, density, and the available lattice parameters.

These rows are Hugoniot states, not 300 K isotherm fit observations. Pressure
and internal energy are Rankine-Hugoniot reductions rather than independent
measurements. LCLS-233 has no diffraction-refined lattice parameters; its
volume and density are obtained from velocimetry.

The rounded Table 1 velocities reproduce the reported pressures within
0.84 GPa, comfortably inside the printed 2-15 GPa uncertainties. Volumes also
reproduce the reported densities within 0.026 g/cm3.

## Source limitations

The article and supplement do not provide:

- numerical shock-corrected 300 K pressure states;
- corrected-state uncertainties or covariance;
- the complete isentropic elastic and heat-capacity inputs;
- the exact EosFit objective and weights; or
- an explicit statement about including the velocimetry-only LCLS-233 state.

The supplement gives `gamma0=1.35` and `q=2.65`, while the first author's
2021 dissertation prints `q=2.6`. Literature-completed versions of the
thermal correction reproduce the reported Hugoniot-temperature scale but do
not recover the published BM3 coefficients, so they do not repair the missing
source protocol.

Table 1 also labels `E-E0` as kJ/mol, although the numerical values satisfy
the printed Rankine-Hugoniot equation in MJ/kg (equivalently kJ/g). The bundled
CSV preserves the literal values and records this source-unit defect.

## Reproduce the retained checks

```bash
.venv/bin/python scripts/audit_schoelmerich_2020_stishovite.py --check
```

The machine-readable result is
[`docs/data/schoelmerich-2020-stishovite-audit.json`](../data/schoelmerich-2020-stishovite-audit.json).
