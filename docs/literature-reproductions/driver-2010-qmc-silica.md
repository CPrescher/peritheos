# Driver et al. (2010): QMC silica equations of state

## Audit result

Three direct QMC Vinet parameterizations are accepted: quartz, stishovite, and alpha-PbO2-type SiO2 (seifertite). LitCurate found the latter two but omitted the equally direct quartz QMC column in the same primary table. Two experimental comparison columns are citation traces and are not imported under the 2010 DOI.

Primary source: [author manuscript (arXiv 1001.2066)](https://arxiv.org/pdf/1001.2066), checked against DOI [10.1073/pnas.0912130107](https://doi.org/10.1073/pnas.0912130107). Downloaded manuscript SHA-256: `bfe2710909b34e8bbfead19adac49750529f0fe0dcfb177b7bb19db7e6e8519d`.

The Equations of State section explicitly says Vinet equations were fitted to the Helmholtz-free-energy isotherms. Table 1 reports volumes in bohr3 per SiO2 and parenthesized one-standard-deviation uncertainties. Conversion uses `1 bohr = 0.529177210903 A`, followed by crystallographic `Z=3`, `2`, and `4` for quartz, stishovite, and alpha-PbO2-type conventional cells.

## Same-DOI row disposition

| LitCurate ID | Primary candidate | Disposition | Peritheos record / reason |
|---|---|---|---|
| not emitted by LitCurate | quartz QMC, Table 1: 254(2), 32(6), 7(1) | accepted after exhaustive same-table audit | `alpha_quartz_driver_2010_qmc_300k_vinet` |
| `litcurate_fc269c4a99bde326` | stishovite QMC: 159.0(4), 305(20), 3.7(6) | accepted | `sio2_stv_andr_driver_2010_qmc_300k_vinet` |
| `litcurate_c4d862798e2cba40` | stishovite experimental range | rejected as citation row | Values belong to Refs. 14-15; the cited primary studies must be audited under their own DOIs. |
| `litcurate_3bebcc8155727eba` | alpha-PbO2 QMC: 154.8(1), 329(4), 4.1(1) | accepted | `seifertite_driver_2010_qmc_300k_vinet` |
| `litcurate_3181bf31d9a5246c` | alpha-PbO2 experimental comparison | rejected as citation row | Values belong to Ref. 1 and the table does not even identify its EOS family. |

The individual QMC E(V) observations and fit covariance are not tabulated. Figure 1 displays QMC statistical envelopes and the 300 K curves. Numerical validation therefore preserves the table values and uncertainties and checks executable inversion/round trips inside each plotted pressure envelope; it does not fabricate digitized fit data or label a coefficient replay as a refit.

## Reproduction

Run:

```bash
uv run python scripts/reproduce_zhang_driver_thermoelastic_eos.py
```

## Zotero-ready metadata

```bibtex
@article{DriverEtAl2010QMC,
  author = {Driver, K. P. and Cohen, R. E. and Wu, Zhigang and Militzer, B. and L{\'o}pez R{\'i}os, P. and Towler, M. D. and Needs, R. J. and Wilkins, J. W.},
  title = {Quantum Monte Carlo Computations of Phase Stability, Equations of State, and Elasticity of High-Pressure Silica},
  journal = {Proceedings of the National Academy of Sciences},
  year = {2010},
  volume = {107},
  pages = {9519--9524},
  doi = {10.1073/pnas.0912130107}
}
```
