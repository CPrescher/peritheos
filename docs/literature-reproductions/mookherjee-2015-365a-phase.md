# Mookherjee et al. (2015): MgSi(OH)6 3.65 Å phase

## Disposition

Primary source: Mookherjee, Speziale, Marquardt, Jahn, Wunder, Koch-Müller, and Liermann (2015), *Equation of state and elasticity of the 3.65 Å phase: Implications for the X-discontinuity*, American Mineralogist 100, 2199–2208, DOI [10.2138/am-2015-5312](https://doi.org/10.2138/am-2015-5312).

Net production is **three records** in `mgsioh6_365a_phase`:

1. the paper's preferred experimental BM3;
2. its explicitly tabulated experimental BM4 sensitivity fit to the same observations; and
3. its static PAW-GGA model-crystal BM4.

These are scientifically distinct parameterizations, not three experimental runs. The source deposits observations from three DAC experiments, but reports only pooled experimental BM3/BM4 fits. No run-specific EOS is claimed or created. The BM4 experiment is deliberately non-default, and the theoretical record is explicitly a 0 K model crystal.

## Candidate-by-candidate audit

| LitCurate row | LitCurate interpretation | Primary-source disposition |
|---|---|---|
| `litcurate_b98adada3a08d56f` | MgSi(OH)6 BM3, V0=194.52 Å³, K0=83 GPa, K′=4.9 | **Accepted.** Preferred experimental EOS. Table 2 gives errors 0.02 Å³, 1 GPa, and 0.2, with V0 fixed. |
| `litcurate_292be0594883e2b5` | MgSi(OH)6 BM4, V0=194.52 Å³, K0=77 GPa; derivatives incomplete | **Accepted after recovery.** Equations 1–3 and Table 2 identify BM4 with K′=7.9(8) and K″=-0.7(2) GPa⁻¹. This is a sensitivity fit on the experimental dataset, not the preferred curve. |
| `litcurate_e091013c8442ddc8` | theoretical MgSi(OH)6 BM4, V0=202.02 Å³, K0=80 GPa, K′=3.4 | **Accepted after recovery.** Table 2 gives K″=-0.05 GPa⁻¹. It is a static PAW-GGA P2₁ model-crystal result, not experimental data. |
| `litcurate_efccbe0ac02295bc` | MgO, K0=160.2 GPa | **Citation trace only.** This is the periclase comparison in Table 4, attributed there to Karki et al. (2000); it is not a result or refit of this paper. |
| `litcurate_d87b5d28e5c194af` | SiO2, K0=294 GPa | **Citation trace only.** This is a silica comparison in Table 4, attributed to Karki et al. (1997); it is not a source EOS here. |
| `litcurate_f7574e211806046d` | MgSiO3, K0=256.7 GPa | **Citation trace only.** This is a bridgmanite/perovskite comparison in Table 4, attributed to Karki et al. (2001); it is not a source EOS here. |
| `litcurate_7a13e5afc23a7fe5` | MgSi2H2O6 phase D, K0=166 GPa | **Citation trace only.** This is the phase-D comparison in Table 4, attributed to Rosa et al. (2013); it is not a result of Mookherjee et al. (2015). |

Thus all seven same-DOI discovery rows are disposed exactly once: three source records and four citation traces.

## Material identity and crystal basis

The primary source does **not** describe MgSiO3. It identifies the 3.65 Å phase as stoichiometric `MgSi(OH)6`, a hydroxide perovskite with a vacant conventional perovskite A site. Mg and Si occupy octahedral units linked by hydroxyl hydrogen bonds. The experimental material has mean Mg/Si=1.02±0.06 from 43 crystals and 34±3 wt.% water, supporting the nominal formula. The phase is monoclinic P2₁ (space-group 4), primitive, Z=2, with 28 atoms per cell. The name is derived from its intense (002) reflection at about 3.65 Å.

The 2015 powder work used Le Bail refinement with the model of Wunder et al. (2012); it did not refine atom coordinates. Accordingly, the material card uses the high-precision experimental ambient lattice fixed in the EOS fit—`a=5.1131`, `b=5.1898`, `c=7.3303 Å`, `β=90.03°`, volume 194.5166 Å³—and clearly labels the fractional coordinates as a positional-topology proxy from COD 9016723/Wunder et al. (2012), DOI 10.2138/am.2012.4022. The fully occupied general sites yield exactly Mg2Si2H12O12 (28 atoms). No sample-specific site populations are invented.

The source's in-cell ambient measurement, 5.116(5), 5.192(5), 7.339(10) Å, β=89.9(1)°, V=194.94(38) Å³, is retained in the observation table but is not substituted for fixed V0.

## Primary observations

The official MSA deposit `AM-15-105312` contains two PDFs. Supplemental Table 1 gives P and unit-cell V; Supplemental Table 2 gives full lattice parameters. The EOS dataset bundles every one of the **95** P–V rows from Supplemental Table 1:

- one initial ambient observation;
- experiment I: 38 compression plus 12 decompression observations;
- experiment II: 19 compression plus 2 decompression observations; and
- experiment III: 23 compression observations.

This totals 80 explicitly labelled compression rows, 14 decompression rows, and one initial ambient row. The pressure range is 0–41 GPa. The two nominal 0 GPa rows have no printed pressure uncertainty, represented by blank CSV fields rather than invented numbers. Duplicate pressure values are genuine repeated observations and are preserved.

The deposit ZIP is available at [the official MSA data URL](https://msaweb.org/MSA/AmMin/TOC/2015/Oct2015_data/AM-15-105312.zip). The transcribed CSV SHA-256 is `85484dad4825baeea039b488f42f17ab42b3a5de341f4406c97b46ea24dd748e`; the source Supplemental Table 1 PDF SHA-256 is `aa12c93568e4157db69306ba292389d3ea5d2461c5394b6cdd5e5bdbc8abbe96`.

The computational energy-volume/P-V grid is plotted as green diamonds in Figure 3a but is not numerically deposited. Because Table 2 fully specifies the theoretical BM4 and the paper defines its equation convention, the theory record is retained as a source parameterization with plot confirmation. A coarse raster digitization would add false precision and is therefore not bundled.

## Equation convention and recovered coefficients

Equations 1–3 define Eulerian strain

`f = 0.5 [(V/V0)^(-2/3) - 1]`

and normalized pressure

`F = P / [3 f (1 + 2f)^(5/2)]`.

The fourth-order relation is

`F = K0 + (3K0/2) [(K′-4)f + {K0 K″ + (K′-4)(K′-3) + 35/9} f²]`.

For the paper's third-order relation, the entire coefficient multiplying `f²` is set to zero. This is the BM3/BM4 convention already implemented by Peritheos.

Table 2 gives:

| Status | V0 (Å³/cell) | K0 (GPa) | K′ | K″ (GPa⁻¹) |
|---|---:|---:|---:|---:|
| experimental BM3, preferred | 194.52(2), fixed | 83(1) | 4.9(2) | implicit BM3 truncation |
| experimental BM4, sensitivity | 194.52(2), fixed | 77(2) | 7.9(8) | -0.7(2) |
| static PAW-GGA BM4 | 202.02 | 80 | 3.4 | -0.05 |

The abstract reports K′=4.9±0.1 for BM3, while detailed Table 2 reports 4.9(2). The table value is stored and the inconsistency is explicit in record provenance. The source does not state a confidence convention for these parenthetical errors and publishes no covariance, so the error-confidence and covariance fields remain null. A fixed V0 still carries the external measurement uncertainty 0.02 Å³; it is not misrepresented as a freely fitted standard error.

## Independent numerical check

`scripts/reproduce_mookherjee_2015_365a_phase.py` evaluates Equations 1–3 directly, without calling Peritheos. Against all 95 official observations:

| Published curve | pressure RMSE (GPa) | maximum absolute residual (GPa) |
|---|---:|---:|
| preferred experimental BM3 | 0.419743 | 1.026109 |
| experimental BM4 sensitivity | 0.397045 | 1.877289 |

The source does not disclose the exact fit row mask, objective, weights, or covariance. To avoid reverse-engineering an undocumented choice, the reproduction adds a declared diagnostic orthogonal-distance regression of the 80 labelled compression rows, fixing V0=194.52 and using both printed P and V uncertainties. It obtains:

| Diagnostic | K0 (GPa) | K′ | K″ (GPa⁻¹) |
|---|---:|---:|---:|
| BM3 ODR | 83.3580 | 4.92820 | — |
| BM4 ODR | 78.0880 | 7.31930 | -0.54439 |

Every diagnostic coefficient falls within the corresponding printed parameter uncertainty of the source result. The BM4 source curve also has slightly lower all-row pressure RMSE than BM3, but this does not make it the preferred record: the authors explicitly state that the data are well represented by BM3 and foreground that parameterization in the abstract. The higher-order terms are strongly covariant and must not be extrapolated beyond the measured range.

The theoretical record is independently evaluated from the printed coefficients and source equation. Its calculated pressures at V=140, 160, 180, 200, and 220 Å³ are 52.9900, 27.4066, 11.2043, 0.8178, and -5.8930 GPa, respectively. These checkpoints establish executable equation identity; they are not substitutes for unpublished raw DFT energies.

## Experimental and theoretical scope

The sample was synthesized at 10 GPa and 425 °C for 77 h. Three symmetric DAC runs used 300 μm culets, Re gaskets, and Ne medium. Ruby fluorescence on the Mao et al. (1986) scale determined pressure in all runs; selected experiment-I images also included Au as a cross-check. Measurements were ambient-temperature synchrotron angle-dispersive XRD at PETRA III P02.2. Pressure and volume uncertainties include the authors' treatment of nonhydrostatic stress inferred from peak broadening/microstrain. Row-wise ruby wavelengths and Au volumes are absent, so pressure-scale recalculation is not possible.

The GGA calculation used the P2₁ primitive model crystal, PAW in VASP, a 1000 eV cutoff, and a 2×2×2 Monkhorst-Pack mesh; energies were converged within 5 meV/atom. Its `temperature_ref=0 K` records a static-lattice calculation, while the experimental records use nominal ambient 300 K. The GGA V0 is 3.8% larger than experimental V0, exactly as the paper notes; the two bases are not conflated.

The experimental validity envelope 0–41 GPa describes observation coverage, not equilibrium phase stability. The authors discuss formation from the hydrous 10 Å phase near 9–10 GPa and metastable retention below that range. The theoretical range follows the plotted Figure 3 calculation envelope and is likewise not a phase-stability assertion.

## Files and checks

- Material: `peritheos/data/materials/mgsioh6_365a_phase.eosmat`
- Official P–V transcription: `peritheos/data/datasets/mgsioh6-365a-mookherjee-2015-supplement-table1-pv.csv`
- Reproduction: `scripts/reproduce_mookherjee_2015_365a_phase.py`
- Focused tests: `tests/test_mookherjee_2015_365a_phase.py`

No aggregate inventory, global ledger, manifest, Zotero library, branch, commit, or remote was changed as part of this audit.
