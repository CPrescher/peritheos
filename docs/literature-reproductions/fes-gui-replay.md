# Morard FeS: actual GUI fitting comparison

Audit date: 2026-10-19. **The distributed 20211014 GUI reproduces the tested console fits closely and still exceeds the ±3 GPa envelope.** The published coefficients are accepted unchanged on parameter agreement; the residual discrepancy is informational.

[Original audit and author-input confirmation](fes-eosfit-direct.md) · [Archived model files, covariance, hashes and residuals](../data/fes-gui-2026-10-19/manifest.json)

The unmodified official Mac GUI was run on the confirmed author-used inputs. The accepted 167-row selection excludes all 11 rejected quenched states. Full MGD is selected, q-compromise off, θ₀=417 K, n=2 and q=1 fixed. The exact author 20210609 executable remains unavailable; this later 2021 build was already pinned in the version audit. These results concern the tested controls, and do not identify every original GUI setting.

## Results

| GUI fit | γ₀ ± conditional esd | Hot-row RMS (GPa) | Hot max | Hot rows outside ±3 |
|---|---:|---:|---:|---:|
| gamma-unweighted | 2.44 ± 0.03 | 1.63 | 5.18 | 12/146 |
| gamma-weighted | 2.46 ± 0.04 | 1.63 | 5.28 | 12/146 |
| staged-weighted | 2.46 ± 0.04 | 1.63 | 5.29 | 12/146 |

The weighted 167-row GUI gamma differs from the console by 0.00001; the staged gamma differs by 0.00002. These differences are much smaller than the conditional parameter errors. Full output precision remains in the saved models and manifest.

The GUI cold fit gives V₀ = 15.4 ± 0.5 cm³/mol, K₀ = 120 ± 30 GPa and K′ = 5.0 ± 0.7. The next GUI data load clears the active model, so the three saved cold tokens were inserted into a copy of the author EOS at Tr=300 K for the hot stage. This is a saved-precision handoff; the archived initial files show exactly what was loaded. No estimation or q-compromise switch was used. Gamma-only controls use the saved author cold model at Tr=298 K and retain all 167 rows; the staged hot fit uses 146 rows.

RMS and rowwise residuals below are native evaluations of GUI-saved models. The GUI screenshots independently preserve convergence, selected weights, refined parameters, esds and maximum residuals. GUI and native maxima agree within 0.02 GPa at displayed precision. Parameter covariance is the actual saved GUI matrix. Gamma errors condition on fixed cold coefficients and omit cold/calibrant uncertainty. Curves show labelled isotherms, while residuals include each row's measured temperature.

### gamma-unweighted

![gamma-unweighted curves and residuals](../data/fes-gui-2026-10-19/gamma-unweighted-curves.svg#only-light)
![gamma-unweighted curves and residuals](../data/fes-gui-2026-10-19/gamma-unweighted-curves-dark.svg#only-dark)

[GUI fit output screenshot](../data/fes-gui-2026-10-19/gamma-unweighted.png)

### gamma-weighted

![gamma-weighted curves and residuals](../data/fes-gui-2026-10-19/gamma-weighted-curves.svg#only-light)
![gamma-weighted curves and residuals](../data/fes-gui-2026-10-19/gamma-weighted-curves-dark.svg#only-dark)

[GUI fit output screenshot](../data/fes-gui-2026-10-19/gamma-weighted.png)

### staged-weighted

![staged-weighted curves and residuals](../data/fes-gui-2026-10-19/staged-weighted-curves.svg#only-light)
![staged-weighted curves and residuals](../data/fes-gui-2026-10-19/staged-weighted-curves-dark.svg#only-dark)

[GUI fit output screenshot](../data/fes-gui-2026-10-19/staged-weighted.png)

[Cold-fit GUI output](../data/fes-gui-2026-10-19/cold-weighted.png)

## Reproduction

Repository-only `uv run python -m scripts.audit_fes_gui_replay` recomputes the saved-model pressure comparisons and verifies the archived inputs. Add `--check` to check the manifest, page and plots. It does not rerun the GUI. To repeat the actual fits in EosFit: open the selected data and initial EOS, enable only the listed refinement flags and weights, click Fit, accept the converged solution and save it. For the staged case, first fit the cold dataset with thermal model None and V0/K0/K′ free, then load the hot dataset and archived thermal initial EOS, fix all but gamma0, enable all P/V/T weights and fit full MGD.

The GUI comparison is now completed for these controls. Remaining limits are the original exact build/workflow and the source-output discrepancy, not missing author-used input files. No further input-file request or replacement EOS is implied.
