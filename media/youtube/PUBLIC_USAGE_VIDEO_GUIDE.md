# Public tutorial companion: Spatial Methods Workbench v0.2.1

## Purpose

Spatial Methods Workbench is a Shiny interface for twelve analysis paths spanning spatial autocorrelation, conditional networks, spatially adjusted networks, prediction, clustering, aligned-image screening, mediation, and contour shape analysis. The public tutorial teaches operation and reproducibility without automatically explaining the biological or clinical meaning of a result.

## Before using real data

1. Complete segmentation, registration, normalization, batch handling, and quality control upstream.
2. Document the observational unit and nesting within patients, specimens, slides, sites, times, or volumes.
3. Remove protected identifiers. Never upload protected health information to the public service.
4. Verify coordinate units, feature names, factor levels, missingness, transformations, and identifiers against a data dictionary.
5. Use a controlled local or institutional deployment for governed data or long Bayesian computations.

## Method-selection map

| Question | Workbench path | Minimum scientific output to preserve |
|---|---|---|
| Is one feature globally spatially autocorrelated? | Moran's I | Statistic, neighbor graph, permutation scheme, seed |
| Which features remain linearly associated conditionally? | Conditional-network screen | Covariance, precision, partial-correlation matrix, threshold |
| What network is shared versus specific to discrete tissue zones? | ISPat 2D | Spatial adjustment, continuous shared/zone matrices, edge summaries, fitted object |
| Which conditional relationships vary continuously over space? | GP-GHS | Basis and sampler settings, spatial edge summaries, selection rule, fitted object |
| What network structure occurs in a registered 3D volume? | ISPat 3D | Registration provenance, shared/zone volumetric matrices, fitted object |
| Can measured features predict a scalar or class outcome? | Linear/logistic prediction | Group-aware split, predictions, errors/probabilities, coefficients |
| Can records be partitioned using features and location? | Spatial clustering | Scaling, coordinate weight, k, assignments, centers, stability |
| Does a compact nonlinear predictor improve held-out performance? | Shallow neural prediction | Architecture, decay, seeds, grouped validation, baseline comparison |
| Which aligned image features screen as predictors of a scalar? | Scalar-image ridge | Mask/alignment, penalty, coefficients, held-out errors |
| Can one aligned image predict another through latent scores? | Image-image latent screen | Factor counts, loadings/scores, reconstructions, held-out error |
| Is there a model-based indirect association through a mediator? | Mediation | Component models, indirect/direct/total effects, bootstrap design |
| What geometric modes summarize processed contours? | Shape PCA | Correspondence rules, scores, loadings, explained variance |

## Run checklist

1. Upload a non-sensitive CSV or load the synthetic demonstration.
2. Map coordinates, zones, identifiers, outcome, exposure, mediator, and features.
3. Review the Readiness tab. A ready status confirms only the interface contract.
4. Select the workflow from the scientific estimand.
5. Record every control before computation; plan sensitivity analyses where relevant.
6. Acknowledge both the Workbench citation and the selected method publications.
7. Run, then audit sample counts, labels, parameters, warnings, convergence, and validation design.
8. Download the reproducibility ZIP and archive it with the frozen input checksum, data dictionary, preprocessing log, and protocol.

## Local installation

```bash
git clone https://github.com/sagnikbhadury/SpatialMethodsWorkbench.git
cd SpatialMethodsWorkbench
Rscript scripts/install_dependencies.R
Rscript tests/testthat.R
R -e "shiny::runApp('.')"
```

Install the public ISPAT, GP-GHS, and ISPAT-3D engines using the repository instructions before enabling their adapters. ISPAT requires a working Stan toolchain. Run an engine-specific synthetic smoke test and archive `sessionInfo()` before formal use.

## Citation and links

- Live application: https://sagnikbhadury.shinyapps.io/spatial-methods-workbench/
- Project website: https://sagnikbhadury.github.io/spatial-methods-workbench/
- Full guide: https://sagnikbhadury.github.io/workbench-guide/
- Source: https://github.com/sagnikbhadury/SpatialMethodsWorkbench
- Version 0.2.1 DOI: https://doi.org/10.5281/zenodo.21764196
- All versions: https://doi.org/10.5281/zenodo.21763606

Spatial Methods Workbench is distributed under the MIT License. Cite the versioned software DOI and every method-specific paper emitted by the selected workflow.

## Interpretation boundary

The application reports statistical results, uncertainty summaries, and fitted objects. It does not assign biological, clinical, or causal meaning. For study-specific method selection, interpretation, deployment, or collaboration, contact Sagnik Bhadury through the application.
