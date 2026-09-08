# Ismailova et al. (2016) Fe-deficient bridgmanite audit

Audit date: 2026-09-08.

## Source and licensing

The primary source is Ismailova et al., “Stability of Fe,Al-bearing bridgmanite in the lower mantle and synthesis of pure Fe-bridgmanite,” *Science Advances* 2, e1600427 (2016), [doi:10.1126/sciadv.1600427](https://doi.org/10.1126/sciadv.1600427). The open [PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC4956391/) identifies `1600427_SM.pdf` as the official supplementary file. The [University of Chicago repository copy](https://knowledge.uchicago.edu/records/z7csb-zr884) is 2,696,968 bytes with MD5 `d772465eb207e08a511b987d73a729c9`, matching the supplementary-media metadata in PMC's full-text package; its SHA-256 is `5ab1c1c60e7d5486299cf5f0653d05b6b8bed6957a382e839ba292e2456ecd34`.

The article and supplement are distributed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). Peritheos does not redistribute the source PDF. It bundles a four-row factual CSV transcription of the EOS-relevant Table S2 columns with attribution and the source license recorded on the dataset. Commercial users of the broader MIT-licensed package should treat this dataset's explicit CC BY-NC metadata as the governing source-use notice.

## Published parameterization

Figure 3A and supplementary Table S3 report a 300 K Birch-Murnaghan fit with `V0=178.98(6) A^3`, `K0=190(4) GPa`, and `K0'=4` fixed. The executable record remains BM2; its implicit fixed derivative is now explicit in `implicit_parameters` and in the fit-provenance metadata. No parameter covariance or confidence convention is published.

The ambient volume is an extrapolated fit parameter, not an observed zero-pressure state. Figure 3A plots compression and decompression measurements from approximately 14 to 132 GPa, broader than the approximately 45–110 GPa synthesis interval described in the article.

## Lossless Table S2 checkpoint data

Supplementary Table S2 says that it contains crystallographic data “at selected pressures.” All four pressure, reported high-temperature-condition, lattice-parameter, and conventional-cell-volume states are transcribed in source-column order. Exact printed tokens are retained beside normalized values and expanded parenthetical uncertainties.

| P (GPa) | Reported high-T condition (K) | a (Å) | b (Å) | c (Å) | V (Å³) |
|---:|---:|---:|---:|---:|---:|
| 44.0(5) | — | 4.6332(17) | 4.815(5) | 6.654(3) | 148.44(18) |
| 67.0(5) | 2200(100) | 4.5516(3) | 4.7753(15) | 6.5497(4) | 142.36(5) |
| 107.0(5) | 2300(100) | 4.4364(6) | 4.7079(10) | 6.328(5) | 132.17(11) |
| 129.0(5) | 1835(100) | 4.3494(15) | 4.6354(14) | 6.296(6) | 126.93(13) |

These four columns are not labeled compression or decompression, and the article does not map them to individual Figure 3 markers. They are therefore checkpoint observations, not asserted source-fit inputs.

## Pressure calibration

Materials and Methods says that pressure was determined from neon lattice parameters using Fei et al. (2007), [doi:10.1073/pnas.0609013104](https://doi.org/10.1073/pnas.0609013104); supplementary Figure S5 repeats the neon-lattice assignment. The record now links the executable `neon_fcc_fei_2007_vinet_2` reference EOS. Recalculation remains blocked because neither the article nor Table S2 gives row-wise neon lattice parameters.

Ruby chips were present only below 15 GPa. The source does not say that ruby defines any Table S2 pressure, so no ruby calibration is assigned to these four rows.

## Fit reconstruction and result

Figure 3A says that the plotted compression and decompression data were fitted, but it gives no explicit exclusions. It also does not report fitting software, residual direction, effective-variance convention, numerical weights, scalar fit statistic, or row-level inclusion flags.

The bundled rows were nevertheless tested under two transparent diagnostics using `scripts/reproduce_ismailova_2016_fe_bridgmanite.py`:

| Diagnostic | V0 (Å³) | K0 (GPa) | Pressure RMSE (GPa) | Reduced chi-square |
|---|---:|---:|---:|---:|
| Published BM2 evaluated at printed Table S2 volumes | 178.98 | 190.0 | 4.0437115 | — |
| Four-row unweighted pressure refit | 172.0922008 | 233.7443976 | 2.2368751 | 10.0072203 |
| Four-row pressure/volume errors-in-variables refit | 172.1799536 | 234.4924103 | 1.2209359 | 21.1758931 |

For the published curve, the pressure residuals in Table S2 order are `+7.7131`, `+1.8149`, `-1.1859`, and `+1.1023 GPa`. The mismatch cannot be explained by the printed uncertainties, and the diagnostic fits do not reproduce either published free coefficient.

## Disposition

The earlier `supplement_available_not_bundled` label was incomplete. The record now has a checksummed `partial_published_table` dataset, an explicit fixed `K0'=4`, a resolved Fei-neon pressure-scale lineage, documented absence of exclusions and weights, and a deterministic diagnostic fit.

The source fit remains `not_refittable`. Exact parity requires the complete numerical Figure 3 P-V series, row-level compression/decompression and fit-inclusion labels, the actual objective/weights, and paired neon calibrant volumes. Digitizing the graph would add pseudo-precision and still would not recover those missing methodological choices.
