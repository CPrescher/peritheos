# Shock Hugoniot candidate catalog

Research cut-off: 2026-09-04

## Decision summary

This review found **two high-confidence linear Hugoniot parameterizations** that are scientifically ready for implementation review, **four conditional candidates** that need an explicit phase/domain decision or a reviewed local regression, and several records that should remain on hold. The initial research pass made no production change. Following review, the phase-restricted B1 MgO relation and the reproducible B1 NiO regression were added on 2026-09-04; transformed and phase-pooled candidates remain excluded.

The strongest records are:

1. **B1 MgO, 14--133 GPa**: a directly published, single-phase linear relation with one-standard-deviation coefficient errors.
2. **Liquid Cu, above completion of melting**: a directly published liquid-branch fit with a covariance matrix; it requires a liquid-Cu material/phase representation rather than attaching the curve to the existing fcc document.

The complete repository inventory is in [`docs/data/hugoniot-material-inventory.csv`](../data/hugoniot-material-inventory.csv). It contains 141 material documents, 224 equilibrium records (171 isothermal and 53 thermal), and two production Hugoniot records, for 226 records total. Phase identity is taken from the explicit phase field where present and otherwise from the material name and crystal structure. Regenerate it with `python scripts/generate_hugoniot_material_inventory.py`.

## Interpretation and acceptance rules

A linear shock relation

\[
U_s=c_0+s u_p
\]

plus an initial state defines one Rankine--Hugoniot path. It is not a general thermodynamic equation of state and it does not supply arbitrary-temperature or off-Hugoniot pressures. A phase change does not automatically inherit the precursor's relation. A relation may be reused on more than one branch only when the primary source itself demonstrates that one fit spans those branches; that reuse must still be recorded as such.

For the domains below, pressure and density were either transcribed from the primary table or derived from the quoted linear relation and the jump conditions. Derived limits are labeled. “Experimental coverage” means the extrema of included measurements, not a claimed phase boundary. “Phase stability” means the authors supplied a phase or melt boundary; approximate intervals must not be stored as exact transitions.

No SESAME table was accessed, digitized, fitted, or redistributed. A paper's comparison to a named SESAME table is not treated as permission to copy it.

## Prioritized catalog

| Priority | Material and represented branch | Coefficients, km/s | Recommended status | Main reason |
|---|---|---|---|---|
| A1 | MgO, B1/periclase, untransformed | `c0=6.87±0.10`, `s=1.24±0.04` | Safe after metadata review | Directly quoted, single-phase interval and stated 1σ errors |
| A2 | Cu, liquid, transformed from ambient fcc | `c0=4.272±0.077`, `s=1.413±0.015` | Safe after adding/approving a liquid phase representation | Directly quoted liquid fit with covariance; strong phase diagnosis |
| B1 | NiO, rhombohedral B1, untransformed | `c0=5.355±0.053`, `s=1.214±0.030` | Conditional derived record | Primary tabulated observations permit a transparent OLS fit, but velocity errors are not reported row by row |
| B2 | MgSiO3, post-perovskite, transformed from enstatite | `c0=4.47±0.81`, `s=1.46±0.18` | Conditional derived record; no extrapolation | Five same-precursor points with errors, but the phase labels are model assignments and the fit is narrow/highly correlated |
| B3 | LiF, B1 solid, untransformed | `c0=5.144±0.010`, `s=1.355±0.004` | Conditional quoted record | Excellent modern fit, but coefficients are global across solid and liquid data rather than a solid-only regression |
| B4 | Ag, fcc solid, untransformed | `c0=3.21±0.02`, `s=1.62±0.02` | Conditional quoted record | Existing fcc phase can host a restricted branch, but the coefficients were fit globally across phase changes |
| C1 | LiF, liquid, transformed | same as B3 | Hold pending liquid-phase representation and boundary policy | Source supports a smooth global relation; mixed region must be excluded |
| C2 | Ag, bcc and liquid transformed branches | same as B4 | Hold pending new phase documents and boundary policy | Phase intervals come from companion diffraction/sound work, not separate linear regressions |
| D1 | Pt, phase unresolved over fit range | `c0=3.64±0.05`, `s=1.541±0.027` | Do not implement as phase-specific | Good kinematic fit, but no defensible single-phase domain has been established |

## Candidate records

### A1 — MgO, B1/periclase principal branch

- **Repository target:** `mgo`; cubic B1, `Fm-3m`, four formula units per conventional cell.
- **Path and branch:** principal loading; untransformed B1 final states from an ambient B1 precursor.
- **Initial state:** nominal ambient pressure and room temperature; measured low-porosity polycrystal density `rho0 = 3.562±0.006 g/cm3`. The paper does not quote a more precise starting temperature in the fit statement. For a four-formula-unit basis, the density corresponds to `V0 = 75.156 A3`; this value is derived, not a static-EOS reference volume.
- **Published relation:** `Us = (6.87±0.10) + (1.24±0.04) up`, with both velocities in km/s. Parentheses are explicitly one-standard-deviation errors.
- **Domain:** `P = 14--133 GPa` as stated by the authors. Derived from the fit, this is approximately `up = 0.523--3.377 km/s` and `rho = 3.828--5.128 g/cm3`. These are experimental coverage limits, not phase boundaries. Do not extend the line to the elastic precursor or beyond 133 GPa without another source.
- **Source:** Duffy and Ahrens, *J. Geophys. Res.* **100** B1, 529--542 (1995), [doi:10.1029/94JB02065](https://doi.org/10.1029/94JB02065). The relation, range, and uncertainty convention are printed in the abstract on p. 529 and as the MgO entry in Table 4 on p. 533; sample density is reported on p. 530 and individual observations in Table 3 on p. 533. An [author-hosted primary PDF](https://duffy.princeton.edu/document/27) is publicly readable.
- **Provenance:** quoted coefficients; no new regression.
- **Access and rights:** the author-hosted article is accessible, but the article states AGU copyright and no reusable data license was found. Transcribing the two published coefficients, uncertainty convention, and derived limits with citation is low risk; do not republish the article's observation tables or figures as a bundled dataset without a rights check.
- **Residual ambiguity:** exact initial laboratory temperature is not numerically reported. If the schema requires `298.15 K`, mark it as a project normalization of “room temperature,” not a quoted measurement.
- **Recommendation:** implemented as `mgo_b1_duffy_ahrens_1995_hugoniot_5` after provenance review.

### A2 — Cu liquid principal transformed branch

- **Repository target:** no liquid-Cu phase document exists. The current `copper` document represents ambient fcc Cu and should be used only as the precursor link, not as the represented liquid phase.
- **Path and branch:** principal loading; transformed liquid branch from ambient fcc Cu.
- **Initial state:** `P0` approximately zero, room temperature (not numerically specified), `rho0 = 8.930±0.003 g/cm3`; a four-formula-unit precursor basis gives derived `V0 = 47.266 A3`.
- **Published relation:** `Us = 4.272 + 1.413 up`, km/s. Table II reports variance/covariance terms `var(c0)=5.964e-3 (km/s)^2`, `var(s)=2.315e-4`, and `cov(c0,s)=-1.116e-3 km/s`, giving `sigma(c0)=0.0772 km/s`, `sigma(s)=0.0152`, correlation about `-0.95`.
- **Domain:** the fit selects shock-melted data above about `265 GPa` and includes literature measurements to about `2000 GPa`. The new high-precision experiments cover `623--1132 GPa`. For a conservative first record, use `P = 265--1200 GPa`, derived `up = 3.314--8.357 km/s` and `rho = 14.177--18.592 g/cm3`; label 265 GPa as an approximate melt-completion boundary and 1200 GPa as a conservative data-quality cap, not a physical upper boundary. The paper discusses caution at `up` near `10 km/s` (about 1.6 TPa under the fit), so the highest legacy points should not silently enlarge the operational domain.
- **Source:** McCoy, Knudson, and Root, *Phys. Rev. B* **96**, 174109 (2017), [doi:10.1103/PhysRevB.96.174109](https://doi.org/10.1103/PhysRevB.96.174109). Fit selection and Monte Carlo procedure are in Sec. III A of the [accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevB.96.174109), approximately pp. 10--12; coefficients and covariance are in Table II, manuscript p. 21.
- **Provenance:** quoted coefficients and covariance. The authors use weighted least squares and `10^6` Monte Carlo realizations of the experimental uncertainties.
- **Access and rights:** the accepted manuscript is publicly readable but carries APS copyright and no open table-data license. Coefficients/covariance may be transcribed with citation. Do not bundle article tables or any referenced SESAME 3325 values.
- **Residual ambiguity:** the exact numeric initial temperature is absent, and the appropriate Peritheos representation of a liquid phase must be agreed. The paper also develops sound-speed/Grüneisen information; that off-Hugoniot thermodynamic analysis must not be collapsed into the linear path record.
- **Recommendation:** scientifically safe once a liquid-phase material representation and conservative upper bound are approved.

### B1 — NiO rhombohedral B1 principal branch

- **Repository target:** `nickel_oxide`; antiferromagnetic rhombohedral B1, `R-3m`, three formula units per represented hexagonal cell.
- **Path and branch:** principal loading; untransformed final/plastic B1 branch. The elastic precursor measurements are a separate branch and were not included in the regression.
- **Initial state:** ambient pressure and room temperature; measured sample bulk density `rho0 = 6.781±0.003 g/cm3` (the separately reported X-ray density is `6.808±0.004 g/cm3`). The conservation calculation must use the bulk density. The three-formula-unit basis gives derived `V0 = 54.873 A3`.
- **New relation:** `Us = (5.355±0.053) + (1.214±0.030) up`, km/s; covariance matrix `[[2.821e-3, -1.419e-3], [-1.419e-3, 8.739e-4]]`, correlation `-0.904`. Errors are formal residual-scaled 1σ OLS errors, not source-reported experimental coefficient uncertainties.
- **Domain:** eight final-state observations, `up = 0.431--2.547 km/s`, `P = 17.7--147.6 GPa`, and reported `V/V0 = 0.930--0.702`, corresponding to `rho = 7.292--9.660 g/cm3`. These are included-data extrema. The source reports no discontinuous volume collapse, but this is not an in-situ structural phase boundary.
- **Source:** Noguchi et al., *J. Phys. Chem. Solids* **60**, 509--514 (1999), [doi:10.1016/S0022-3697(98)00296-0](https://doi.org/10.1016/S0022-3697(98)00296-0), Table 1, journal p. 510. The local primary transcription is `peritheos/data/datasets/nickel-oxide-noguchi-1999-table1-shock.csv`. A shorter, publicly readable [primary conference paper](https://www.jstage.jst.go.jp/article/jshpreview1992/7/0/7_0_832/_pdf) covers the experiment to 132 GPa on pp. 832--834.
- **Provenance and fit:** newly fit in this review by ordinary least squares of final `Us` on final `up`, all eight rows, equal weights, `n=8`, `RMSE(Us)=0.0556 km/s`. The paper says pressure and volume errors are below about 2% but does not print row-wise `Us`/`up` uncertainties, so an errors-in-variables covariance cannot be justified.
- **Access and rights:** the Elsevier version is restricted; the J-STAGE conference version is readable but no explicit CC/data license was found. The two derived coefficients can be distributed with provenance. Redistribution status of the already bundled full table should be reviewed separately.
- **Residual ambiguity:** structural persistence is inferred from kinematics/no volume discontinuity rather than direct diffraction, and exact initial temperature is not stated.
- **Recommendation:** implemented as the explicitly derived record `nickel_oxide_noguchi_1999_linear_hugoniot_2` after native Peritheos fit reproduction, retaining the limitations above.

### B2 — MgSiO3 post-perovskite transformed branch

- **Repository target:** `mgsio3_post_perovskite`, `Cmcm`, four formula units; precursor link `orthoenstatite`.
- **Path and branch:** principal loading; transformed post-perovskite final states from ambient enstatite.
- **Initial state:** nominal ambient pressure and room temperature; `rho0 = 3.199±0.009 g/cm3` for the enstatite shots. The four-formula-unit operational basis gives derived `V0 = 208.436 A3`.
- **New relation:** `Us = (4.47±0.81) + (1.46±0.18) up`, km/s. Absolute input-uncertainty covariance is `[[0.6543, -0.1481], [-0.1481, 0.03357]]`, correlation `-0.9993`. The large uncertainty and correlation are material facts, not rounding noise.
- **Domain:** five same-precursor rows assigned to post-perovskite: `up = 4.26--4.90 km/s`, `P = 142.6--181.6 GPa`, `rho = 5.34--5.54 g/cm3`. These are observation extrema only. A sixth post-perovskite-assigned row starts from an `80% perovskite + 20% majorite` mixture at `rho0=3.989 g/cm3`; it was excluded because a linear Hugoniot cannot mix initial states.
- **Source:** Mosenfelder et al., *J. Geophys. Res.* **114**, B01203 (2009), [doi:10.1029/2008JB005900](https://doi.org/10.1029/2008JB005900), Table 2 and Sec. 4.1. The [primary Wiley article](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JB005900) renders the table. The local primary transcription is `peritheos/data/datasets/mgsio3-post-perovskite-mosenfelder-2009-table2-shock.csv`.
- **Provenance and fit:** newly fit in this review using orthogonal distance regression with the five printed `sigma(Us)` and `sigma(up)` pairs, treated as absolute one-sigma inputs. The uncertainty convention should be reconfirmed from the article before implementation. The source's thermodynamic phase assignment was selected only when competing assignments differed by at least `chi^2=1` (described as at least 90% confidence); it is not in-situ diffraction.
- **Access and rights:** the primary article/table is readable, but no explicit reusable data license was found. Derived coefficients and row identifiers are safe to cite; confirm the rights basis before redistributing the complete table outside the repository.
- **Residual ambiguity:** phase labels are model-dependent; the narrow `up` span makes the intercept poorly constrained; the stated PPv interval is not a measured stability boundary.
- **Recommendation:** conditional and non-extrapolatable. Require independent fit reproduction plus review of the phase-assignment language.

### B3/C1 — LiF [100], B1 solid and liquid branches

- **Repository target:** `lif_b1` for the solid branch. No liquid-LiF document exists.
- **Path and branch:** principal loading. Solid B1 is untransformed; liquid LiF is transformed from the B1 precursor. The partial-melt interval must be omitted.
- **Initial state:** ambient pressure, room temperature, `[100]` UV-grade single crystal; `rho0 = 2.640±0.002 g/cm3`. A four-formula-unit basis gives derived `V0 = 65.262 A3`.
- **Published relation:** `Us = (5.144±0.010) + (1.355±0.004) up`, km/s.
- **Fit method:** global weighted orthogonal-distance regression of the authors' measurements plus selected prior measurements; Appendix A documents exclusions for unknown errors and one outlier. This is not an independently fit solid line and liquid line. The source explicitly reports one smooth relation with no velocity discontinuity, so branch-limited reuse is supportable if that provenance is retained.
- **Solid domain:** the study identifies solid behavior through approximately `182 GPa`; derived upper state `up≈5.483 km/s`, `rho≈4.682 g/cm3`. The global historical dataset reaches low pressure, but the paper does not print one clean numerical lower calibration bound for every included source. Store the lowest audited included datum, not `0 GPa`, when the observation ledger is assembled. Mark 182 GPa as an approximate onset boundary.
- **Mixed/liquid domains:** melting occurs approximately `182--195 GPa`; no single-phase record should cover this interval. A liquid branch can begin at about `195 GPa` (`up≈5.725 km/s`, `rho≈4.746 g/cm3`) and is directly constrained to `231.1 GPa`. Table II's highest measured state is `up=6.333 km/s`, `rho=4.88 g/cm3`; small differences from values derived from the global fit are expected residuals, not transcription errors.
- **Source:** Hawreliak et al., *Phys. Rev. B* **107**, 014104 (2023), [doi:10.1103/PhysRevB.107.014104](https://doi.org/10.1103/PhysRevB.107.014104). Relation: Eq. (3), accepted-manuscript p. 6; fit protocol: Appendix A, pp. 11--12; new states: Table II, p. 32. The [APS accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevB.107.014104) is publicly readable.
- **Access and rights:** public accepted manuscript, APS copyright, no explicit open table-data license. Coefficient transcription is acceptable; do not copy the compiled historical observation table into Peritheos without source-by-source provenance and rights review.
- **Residual ambiguity:** phase boundary values are intervals inferred from sound speed, not exact thermodynamic boundaries. The fit pools phases. The crystal orientation is part of the provenance even though the hydrodynamic branch is commonly used as a window standard.
- **Recommendation:** solid branch is conditional; liquid branch remains on hold until a liquid phase representation is approved. Never implement the 182--195 GPa mixed region as either phase.

### B4/C2 — Ag fcc, bcc, and liquid principal branches

- **Repository target:** `silver` represents fcc Ag. No bcc-Ag or liquid-Ag documents exist.
- **Path and branches:** principal loading from ambient fcc Ag. The fcc interval is untransformed; bcc and liquid intervals are transformed.
- **Initial state:** ambient pressure, room temperature, polycrystalline `99.95%` Ag; `rho0 = 10.50±0.01 g/cm3`. A four-formula-unit basis gives derived `V0 = 68.236 A3`.
- **Published relation:** `Us = (3.21±0.02) + (1.62±0.02) up`, km/s.
- **Fit and total data domain:** a linear fit to data over about `30--300 GPa`; Table II spans `up=0.682--3.34 km/s`, `P=30.9--301 GPa`, and `V/V0=0.842--0.611`. Some table states were inferred with the fit, so they must not be treated as independent observations in a new regression.
- **Phase-restricted domains:** companion diffraction places fcc-to-bcc near `150 GPa`; sound-speed work places the last clearly solid points near `171--172 GPa`, a mixed/melting region thereafter, and complete liquid near `218 GPa`. Conservative derived bounds are fcc `30.9--150 GPa` (`up≈0.682--2.140`, `rho≈12.471--15.452 g/cm3`), bcc about `150--171.5 GPa` (`up≈2.140--2.336`, `rho≈15.452--15.765`), omit the mixed interval, and liquid `218--301 GPa` (`up≈2.724--3.331`, `rho≈16.338--17.130`). All phase boundaries are approximate.
- **Primary sources:** kinematic fit and Table II: Wallace, Winey, and Gupta, *Phys. Rev. B* **104**, 014101 (2021), [doi:10.1103/PhysRevB.104.014101](https://doi.org/10.1103/PhysRevB.104.014101), Eq. (1), accepted-manuscript p. 8 and Table II p. 18, [accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevB.104.014101). Phase-boundary support: Briggs et al., *Phys. Rev. Lett.* **124**, 235701 (2020), [doi:10.1103/PhysRevLett.124.235701](https://doi.org/10.1103/PhysRevLett.124.235701), and Wallace et al., *Phys. Rev. B* **104**, 214106 (2021), [doi:10.1103/PhysRevB.104.214106](https://doi.org/10.1103/PhysRevB.104.214106), [accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevB.104.214106).
- **Provenance:** quoted global coefficients, then restricted to phase intervals using independent primary diagnostics. They are not separate phase regressions.
- **Access and rights:** APS accepted manuscripts are publicly readable but not openly licensed datasets. Coefficients and derived limits may be cited; tables/figures should not be bundled without a rights review.
- **Residual ambiguity:** the fcc--bcc and melt boundaries are not exact, and one fit spanning all phases does not demonstrate identical phase thermodynamics. The sound-speed paper also develops a liquid Mie--Gruneisen treatment; it is an off-Hugoniot construction and must remain separate.
- **Recommendation:** fcc can be considered after boundary policy review. Hold bcc/liquid until phase documents exist. Omit all mixed-region states.

### D1 — Pt, phase unresolved

- **Repository target:** `platinum`, ambient fcc.
- **Path and branch:** principal loading from ambient fcc Pt; final phase is not adequately established across the full relation, so neither `untransformed` nor a transformed phase label is currently defensible.
- **Initial state:** ambient pressure and room temperature; `rho0 = 21.43±0.03 g/cm3`; derived four-formula-unit `V0 = 60.466 A3`.
- **Published relation:** Hawreliak et al. reproduce the Holmes relation as `Us = (3.64±0.05) + (1.541±0.027) up`, km/s, citing the original source.
- **Domain:** Holmes et al. discuss Pt shock data from about `32--660 GPa`. The seven new Table III points bundled locally cover `P=218.9--659.3 GPa`; a local errors-in-variables check of their printed two-sigma velocity errors gives `c0=3.651±0.037`, `s=1.538±0.015`, corroborating but not replacing the quoted fit. The range is kinematic coverage only, not a single-phase range.
- **Primary source:** Holmes, Moriarty, Gathers, and Nellis, *J. Appl. Phys.* **66**, 2962--2967 (1989), [doi:10.1063/1.344177](https://doi.org/10.1063/1.344177), Table III, p. 2964. The local transcription is `peritheos/data/datasets/platinum-holmes-1989-table3-shock.csv`. The coefficient restatement is Hawreliak et al. (2023), Table IV, accepted-manuscript p. 34.
- **Provenance:** use the Holmes coefficients as quoted via the modern primary study only after checking the original article; the local seven-point fit is an audit calculation, not the preferred record.
- **Access and rights:** the AIP article is copyrighted and no open data license was found. Coefficient facts may be transcribed with citations; do not redistribute the complete table merely because a copy is readable online.
- **Blocking ambiguity:** no primary phase-resolved interval was verified. A relation spanning possible high-pressure solid/liquid changes cannot be attached to the fcc document as if it were phase-specific.
- **Recommendation:** hold.

## Screened but not recommended

- **Alpha quartz (`alpha_quartz`):** Knudson and Desjarlais report substantial curvature and use `Us = a + b up - c up exp(-d up)`, not a `LinearUsUpHugoniot`. Primary source: *Phys. Rev. Lett.* **103**, 225501 (2009), [doi:10.1103/PhysRevLett.103.225501](https://doi.org/10.1103/PhysRevLett.103.225501), [full text](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.103.225501/fulltext). Reject for this model class.
- **Legacy Al, Cu, and Ta standards:** Mitchell and Nellis, *J. Appl. Phys.* **52**, 3363--3374 (1981), [doi:10.1063/1.329160](https://doi.org/10.1063/1.329160), provide useful linear calibrant relations over broad ranges. Later sources quote, for example, Al `c0≈5.385`, `s≈1.339`, Cu `c0≈3.97`, `s≈1.48`, and Ta `c0≈3.32`, `s≈1.30`. They are not promoted here because the ranges span or approach phase changes, the exact material grade (for example alloy versus pure Al) matters, and the modern restatements combine multiple sources. Re-audit the original tables and phase diagnostics before proposing branch records.
- **Diamond (`diamond`):** high-pressure data are approximately linear only over selected mixed/melt intervals. Without a clean phase-specific interval and matching phase document, a straight-line record would encode the wrong boundary semantics.
- **Iron polymorphs:** common impedance-matching relations pool bcc/fcc/hcp and/or melt regimes. No global “iron” line should be assigned to one repository polymorph.
- **Precompressed paths:** no repository material in this pass yielded a primary, phase-specific *linear* `Us-up` relation with a fully reported positive `P0`, `T0`, and `rho0` suitable for Peritheos. Precompression papers often report individual Hugoniot points or use a thermodynamic EOS to calculate paths rather than calibrating a new straight line. Do not relabel those data as a principal relation or reuse an ambient `rho0`.

## Data-access and licensing policy

This section is a conservative provenance/redistribution screen, not legal advice. The recommended implementation artifact should contain only the small set of published coefficients, uncertainty/covariance metadata, bibliographic locators, fit-selection notes, and derived domain values needed to execute the model. Public readability is not the same as an open redistribution license.

For newly fit records, the observation file must have a documented rights basis in addition to scientific provenance. If redistribution is unclear, retain a row manifest/checksum and fit result without adding the article's complete numeric table. Never copy a SESAME table or reconstruct it from another publication.

## Implementation gate

Before any production record is added:

1. A second reviewer should verify every quoted number against the cited page/table and confirm whether coefficient errors are 1σ, 2σ, or another convention.
2. Reproduce NiO and MgSiO3 fits with `fit_linear_us_up`, preserving row selection, absolute-error convention, covariance, residuals, and a checksum of the input file.
3. Decide how liquids are represented as material/phase documents. Never attach liquid Cu, liquid LiF, or liquid Ag to an ambient crystal document without an explicit represented-phase field.
4. Store observed-domain and phase-boundary semantics separately. Approximate melt intervals must not become exact endpoints.
5. Record room-temperature normalization explicitly wherever the paper does not print a numeric `T0`.
6. Run the volume-basis consistency check using precursor `rho0`, molar mass, and formula-unit count; do not reuse a static-EOS `V0` automatically.
7. Review redistribution rights independently of scientific validation.

The current production set contains only the **MgO B1** and clearly labeled derived **NiO B1** principal Hugoniots. Liquid Cu and MgSiO3 post-perovskite remain transformed-branch candidates. LiF and Ag need a branch-boundary/phase-document decision. Pt and the screened legacy/global relations remain excluded.
