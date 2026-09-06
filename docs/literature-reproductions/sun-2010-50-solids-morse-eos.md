# Sun et al. (2010): Morse equations for 50 solids

## Decision

This paper is accepted as a source of explicit, executable legacy benchmark
parameterizations. It defines the ordinary three-dimensional Morse EOS (MRS3)
and the Sun Jiu-Xun--Morse EOS for `n=3` (SMS3) and `n=4` (SMS4), then prints
all three fitted parameter pairs for 50 solids in Table 2. The paper concludes
that SMS4 gives the lowest mean pressure error of the tested universal EOS
families, so SMS4 is promoted first; SMS3 and MRS3 are retained as explicit
same-data model sensitivities where included.

These records are deliberately labelled *legacy compression aggregates*. The
source pooled older compression compilations over broad pressure ranges and did
not preserve a phase-resolved observation table, row-wise pressure calibration,
uncertainties, or regression weights. The curves are consequently useful for
model comparison and historical benchmarking, not as modern phase-pure pressure
standards. No structure is inferred from a chemical label.

The Brass row is held. The paper gives neither alloy composition nor a
reproducible formula/molar mass, so it cannot define a Peritheos material or
volume basis without invention.

## Primary source

- Jiu-Xun Sun, Qiang Wu, Yang Guo, and Ling-Cang Cai, “Two Universal Equations
  of State for Solids,” *Zeitschrift für Naturforschung A* **65**, 34–44
  (2010), DOI `10.1515/zna-2010-1-202`.
- Checked locations: equations (3)–(15), Tables 1–4, Figures 1–6, and the
  discussion and conclusion.
- Audited PDF SHA-256:
  `57b6c4ff8b59dcc3734aa9507de0f62610839363a8e5ae790ed8f71e6fa66db4`.
- Exact Table 1/2 transcription:
  `peritheos/data/datasets/sun-2010-table2-morse-eos-parameters.csv`, SHA-256
  `e71613962b7df52559b8633546a7513b902ce1b5d8d4a7666e0476464123ca68`.

## Equations

For `X=(V/V0)^(1/n)`, the ordinary Morse form in equation (6) is

`P = n K0 / (beta X^(n-1)) [exp(2 beta (1-X)) - exp(beta (1-X))]`,

where `n=3` and `beta=n(K0'-2)/3+1=K0'-1`. The Sun-Morse form in equation
(8) is

`P = n K0 / alpha [exp(2 alpha (1-X)) - exp(alpha (1-X))]`,

with `alpha=(n K0'+1)/3`; the paper evaluates `n=3` and `n=4`. Equation (15)
provides the analytic Sun-Morse bulk modulus. The native implementations recover
`P(V0)=0`, `K(V0)=K0`, and `(dK/dP)(V0)=K0'` numerically for every one of the
150 printed curves.

## Data lineage and limitations

The source says that normalized compression observations come from Kennedy and
Keeler (1972), except for solid hydrogen, tungsten, and NaCl. Those row-level
observations are not republished. Table 1 does provide the exact fitted pressure
range and mean pressure error for every solid/model combination, while Table 2
provides `V0`, `K0`, and `K0'`. The stored CSV is therefore a parameter and fit-
diagnostic table, not primary P–V observations. Published molar volumes are
converted to one-formula-unit volumes with the exact Avogadro constant:
`1 cm3/mol = 1.6605390671738466 A3/formula unit`.

No coefficient uncertainties, covariance, confidence convention, or exact
pressure-scale equations are stated. Those fields remain null/unresolved. The
fit ranges are preserved exactly, but they must not be interpreted as verified
single-phase stability ranges.

## Reproduction

Run:

```bash
uv run python scripts/reproduce_sun_2010_morse_eos.py
uv run pytest -q tests/test_sun_morse.py tests/test_sun_2010_morse_records.py
```

The reproduction checks the source table checksum, all model reference-state
identities, the analytic bulk modulus, and the link from every promoted record
back to its exact Table 2 row.
