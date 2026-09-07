# AIP Handbook Cu Vinet refit

## Scope and source reading

This note records a limited transcription of the Cu column in the *American
Institute of Physics Handbook*, section 4d, Table 4d-12, p. 4-100. The table
reports relative volume, `V/V0`, as a function of pressure. The 25 populated Cu
cells from 15 to 340 kbar were manually checked against the attached printed
copy on 2026-09-07. The 5 and 10 kbar Cu cells are blank and are not treated as
observations.

The 15 kbar cell is printed as `09900`. This is a printing error, not a scan or
OCR artifact, and is transcribed as `0.9900`. The remaining printed Cu values
were confirmed as read. Only this small factual excerpt is stored in
`aip-handbook-table4d12-cu.csv`; the handbook PDF is not part of the repository.

| Pressure (kbar) | `V/V0` | Pressure (kbar) | `V/V0` | Pressure (kbar) | `V/V0` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 15 | 0.9900 | 20 | 0.986 | 25 | 0.983 |
| 30 | 0.980 | 35 | 0.977 | 40 | 0.974 |
| 45 | 0.971 | 50 | 0.968 | 60 | 0.962 |
| 70 | 0.956 | 80 | 0.951 | 90 | 0.945 |
| 100 | 0.940 | 120 | 0.930 | 140 | 0.921 |
| 160 | 0.912 | 180 | 0.904 | 200 | 0.896 |
| 220 | 0.889 | 240 | 0.881 | 260 | 0.874 |
| 280 | 0.868 | 300 | 0.861 | 320 | 0.855 |
| 340 | 0.849 |  |  |  |  |

## Diagnostic fit

The reproduction script fits the two free Vinet parameters while fixing
`V0 = 11.81473546294192 A^3/formula unit`, the value associated with the Cu
curve in the Sun et al. catalog transcription. Since the observations are
ratios `V/V0`, the chosen numerical `V0` cancels from the fitted `K0` and
`K0_prime`.

Following equation (20) of Sun et al. (2010), the minimized objective is

`sum((P_model - P_table)^8)`.

The deterministic multi-start Nelder-Mead calculation uses the 25 populated
finite-pressure rows. The identity point `P=0, V/V0=1` is not counted as an
observation because every fixed-`V0` Vinet curve satisfies it exactly.

| Curve | `K0` (GPa) | `K0_prime` | RMSE (GPa) | MAE (GPa) | Max abs. residual (GPa) | L8 objective (GPa^8) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Diagnostic refit | 139.23431 | 5.001127 | 0.070154 | 0.062704 | 0.112737 | 8.92820e-08 |
| Sun et al. published Vinet | 140.95 | 4.798 | 0.083664 | 0.070923 | 0.154252 | 9.73909e-07 |

The diagnostic result does not reproduce the published coefficients exactly:
`K0` is 1.22% lower and `K0_prime` is 4.23% higher. Rounded table values and
unreported numerical details are plausible causes. The published curve remains
a source-author record; the corrected cell does not authorize silently changing
its coefficients.

## Catalog implication

This dataset can support a separate Peritheos refit EOS in principle. Such a
record should use an identifier such as `cu_sun_2010_low_vn_refit`, set
`record_kind` to `refit`, link to the published Vinet record with
`derived_from_record`, and carry this checked dataset plus the exact `m=8` fit
provenance. It must not replace or masquerade as `cu_sun_2010_low_vn`.

Promotion is deliberately left for a coherent catalog update, where the dataset
can be embedded or checksummed in the material document and the generated global
refit ledger can be updated atomically. The present commit documents and tests
the evidence and calculation without mixing it into unrelated catalog changes.

## Reproduction

```console
uv run --frozen python scripts/reproduce_aip_handbook_cu_vinet_refit.py
uv run --frozen pytest -q tests/test_aip_handbook_cu_vinet_refit.py
```
