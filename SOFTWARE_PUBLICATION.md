# Software publication record

## Current public release

- Software: **SpatialMethodsWorkbench: Spatial Methods Workbench**
- Author: Sagnik Bhadury (University of Pennsylvania)
- Version: **v0.4.0**
- Release date: 2026-10-24
- DOI across all versions (concept DOI): <https://doi.org/10.5281/zenodo.21763606>
- Version DOI: <https://doi.org/10.5281/zenodo.23627817> (Zenodo record 23627817, created 2026-10-24 by the GitHub integration from release v0.4.0)
- GitHub repository: <https://github.com/sagnikbhadury/SpatialMethodsWorkbench>
- GitHub release: <https://github.com/sagnikbhadury/SpatialMethodsWorkbench/releases/tag/v0.4.0>
- Live application: <https://sagnikbhadury.shinyapps.io/spatial-methods-workbench/>
- Website showcase: <https://sagnikbhadury.github.io/spatial-methods-workbench/>
- Website guide: <https://sagnikbhadury.github.io/workbench-guide/>

Suggested citation:

> Bhadury, S. (2026). *SpatialMethodsWorkbench: Spatial Methods Workbench* (Version 0.4.0) [Computer software]. Zenodo. <https://doi.org/10.5281/zenodo.21763606>

Users must also cite the method-specific publications written to `CITATION.txt` in each exported result bundle.

## Engines executed by v0.4.0

| Workbench path | Package and version | Function called | Source of record |
|---|---|---|---|
| ISPat 2D spatially adjusted network | `ISPAT3D` 0.3.1 (CRAN) | `ispat3d_fit_2d()` | <https://cran.r-project.org/package=ISPAT3D>, GitHub commit `18626f0`, Zenodo 10.5281/zenodo.22798429 |
| Volumetric ISPat 3D network | `ISPAT3D` 0.3.1 (CRAN) | `ispat3d_fit()` | same |
| Spatially varying GP-horseshoe network | `GPGHS` 1.0.1 | `gpghs_graph()` | <https://github.com/sagnikbhadury/GP-GHS>, commit `c42319e` |

The standalone `ISPAT` package (CAVI/SVI with Stan) is retired from the Workbench as of v0.4.0. Versions v0.3.1 and earlier targeted engine interfaces that no longer exist upstream (ISPAT 0.1.0, ISPAT3D 0.1.1, GPGHS 0.1.0), so those releases cannot run the three advanced paths against the currently released packages.

## Public scope

The release contains 12 executable paths: spatial quality control, conditional networks, ISPat 2D, GP-GHS, ISPat 3D, spatial machine learning, spatial clustering, neural prediction, scalar-on-image screening, image-to-image screening, mediation, and tumor-shape PCA.

The software reports numerical/statistical outputs and reproducibility artifacts. It does not explain results, assign scientific or biological meaning, or make clinical conclusions. Study-specific interpretation is reserved for direct discussion or collaboration with Sagnik Bhadury.

Only already-public, release-scoped methods and repositories belong in the public application. Do not add non-public or confidential material.

## Release artifacts

- `README.md`: overview, capabilities, installation, citation, and limitations.
- `docs/USER_GUIDE.md`: restrictions, data preparation, installation, Docker, deployment, and responsible use.
- `docs/PIPELINE_REFERENCE.md`: complete function-by-function and workflow reference.
- `CITATION.cff` and `.zenodo.json`: machine-readable software citation and Zenodo metadata.
- `DESCRIPTION`: R dependency and version metadata.
- `vendor/`: exact engine sources and curated reference sources with `manifest.json`.
- `tests/testthat/`: regression, engine-adapter, and Shiny server tests.

The Zenodo v0.4.0 archive corresponds to GitHub release tag `v0.4.0`. The concept DOI resolves to the newest version; use the version-specific DOI displayed on the v0.4.0 record when citing that exact release.
