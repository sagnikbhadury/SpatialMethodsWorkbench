# Spatial Methods Workbench 0.4.0

- Re-targeted the three advanced network paths to the current public engines. The ISPat 2D and ISPat 3D paths now call the CRAN package ISPAT3D 0.3.1 (`ispat3d_fit_2d` for the section-wise planar fit and `ispat3d_fit` for the pooled volumetric fit); the GP-GHS path calls GPGHS 1.0.1 (`gpghs_graph`). The retired ISPAT 0.1.0, ISPAT3D 0.1.1, and GPGHS 0.1.0 APIs no longer exist upstream, so v0.3.1 could not run these analyses against the released packages.
- ISPat results now report the engine's own signed zone partial correlations from the full fitted covariances instead of re-inverting the matrices. The shared low-rank component is no longer displayed as a network because it has none. New controls expose factor rank, anchor budget, Vecchia neighbors, optimizer iterations, section handling, threads, and seed; the GP adjustment log and GP fallback counts are returned with every result.
- The mapped sample/section/volume ID is used as the ISPAT3D section identifier only when it looks like a section label (at most one distinct value per ten rows); otherwise rows are pooled into one section and the result says so. A control overrides the rule.
- GP-GHS selection now uses the shrinkage-factor rule of the article (0.5 by default, previously 0.15) with the AND rule, an optional minimum-detection filter, a reproducible seed, and the GPGHS edge table (score, sign, sign change, directed shrinkage factors) alongside the coefficient-field RMS weights.
- Readiness checks now require at least two zones with more complete rows than selected features before an ISPAT3D-backed path is offered.
- Result manifests record the engine package, version, and function that produced each result.
- The complete source distribution bundles the exact CRAN ISPAT3D 0.3.1 source and GPGHS 1.0.1; the standalone ISPAT package is retired from the bundle. The installer now needs only `gpboost` as a compiled engine dependency and no longer installs Stan.
- Added real adapter smoke tests that run when the engines are installed, including in GitHub Actions.
- Citations for the ISPat paths now include the ISPAT3D CRAN package and its archived software DOI; the GP-GHS citation names the GPGHS package version.

# Spatial Methods Workbench 0.3.1

- Added visible figures for every registered analysis pipeline and automatic navigation to the Figure tab after a successful run.
- Added ISPAT-style network diagrams for conditional-network, ISPAT, and ISPAT-3D analyses.
- Added a GP-GHS network diagram whose edge widths are derived from posterior mean coefficient-field RMS values, not kappa values.
- Added deterministic, ready-to-analyze synthetic datasets and expanded compatibility and figure tests.
- Corrected the GitHub Actions dependency list so the complete test and vignette workflow compiles with the network renderer.

## Funding acknowledgement

Research reported in this software was supported by the National Cancer Institute of the National Institutes of Health under Award Number **R37CA214955-01A1**. The content is solely the responsibility of the author and does not necessarily represent the official views of the National Institutes of Health.
