# Yang et al. (2015): (Mg0.92Fe0.08)O ferropericlase

## Primary-source audit

- Citation: J. Yang, X. Tong, J.-F. Lin, T. Okuchi, and N. Tomioka, “Elasticity of Ferropericlase across the Spin Crossover in the Earth’s Lower Mantle,” *Scientific Reports* **5**, 17188 (2015), DOI [`10.1038/srep17188`](https://doi.org/10.1038/srep17188).
- Publisher article: <https://www.nature.com/articles/srep17188>
- Official supplement: <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fsrep17188/MediaObjects/41598_2015_BFsrep17188_MOESM1_ESM.pdf>
- Local primary-PDF audit SHA-256: `d141f132653bf8eb0b10bf7f0e7d63c2d9a8e8f90eeb46f924683be8a81f5c4f`
- Local supplementary-PDF audit SHA-256: `ddb7c64834ae8ba46547080c35874c535efce0b46681f26f3a0323f8ff225efc`
- Audited 2026-09-05.

The Results give high-spin `KT0=152.5(2.4) GPa`, `KT0'=4.1(0.2)`, low-spin `KT0=161.6(7.1) GPa`, and fixed low-spin `KT0'=4`. The Methods give the directly measured ambient cubic lattice parameter `a=4.1996(4) Å`; therefore the conventional-cell high-spin reference volume is `a^3=74.06683401593601 Å^3` (four formula units). Supplementary Equations 1–5 establish that these are pure-spin reference curves coupled through the low-spin fraction, not two measured run splits.

## Exhaustive LitCurate disposition

| Candidate | Source claim | Decision | Reason |
|---|---|---|---|
| `litcurate_681529622576473b` | high-spin BM3 `K0=152.5`, `K0'=4.1` | **accepted** as `mg092fe008o_yang_2015_hs_bm3_reference` | The Methods ambient cell measurement completes the branch; composition, volume basis, uncertainties, and model are primary-source explicit. |
| `litcurate_0cb0d7a43e4cda51` | low-spin BM3 `K0=161.6`, fixed `K0'=4` | **held, incomplete** | No independent low-spin `V0` is printed. Reusing the ambient high-spin volume would silently erase the source model’s spin-volume collapse and is not defensible. |

This yields one production record from two same-DOI candidates. The low-spin row remains documented but is neither padded nor treated as a pseudo-run.

## Reproduction

Run `uv run python scripts/reproduce_yang_2015_ferropericlase.py`. The script checks `V0=a^3`, exact source-coefficient transcription, three forward BM3 pressure checkpoints, and inverse volume recovery. The article publishes no machine-readable P–V table; numerical validation is consequently an analytical execution check of the complete published branch rather than a refit.

The XRD experiment used Au powder as pressure calibrant and Ne medium; BLS/ISS used ruby spheres. Because row-wise Au or ruby observables are absent, pressure recalculation is marked unresolved.
