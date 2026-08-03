# Reproducible simulated-data vignettes

The Workbench includes five executable vignettes built from deterministic, deliberately small simulated datasets. They document how to prepare inputs, reproduce the app's column mapping and readiness checks, call every workflow family, inspect outputs, and audit a downloaded result bundle.

## Start here

1. [Reproduce the Shiny app workflow](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/blob/master/vignettes/app-walkthrough.Rmd)
2. [Run spatial diagnostics and network workflows](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/blob/master/vignettes/spatial-networks.Rmd)
3. [Run prediction, clustering, neural, and mediation workflows](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/blob/master/vignettes/prediction-clustering-mediation.Rmd)
4. [Run imaging, shape, and volumetric workflows](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/blob/master/vignettes/imaging-shape-3d.Rmd)
5. [Audit reproducibility bundles](https://github.com/sagnikbhadury/SpatialMethodsWorkbench/blob/master/vignettes/reproducibility-bundles.Rmd)

## Small demonstration files

- [Spatial cell/spot CSV — 72 rows](https://raw.githubusercontent.com/sagnikbhadury/SpatialMethodsWorkbench/master/inst/extdata/spatial_demo.csv)
- [Aligned image-feature CSV — 48 subjects](https://raw.githubusercontent.com/sagnikbhadury/SpatialMethodsWorkbench/master/inst/extdata/image_demo.csv)
- [Ordered contour CSV — 6 contours × 16 points](https://raw.githubusercontent.com/sagnikbhadury/SpatialMethodsWorkbench/master/inst/extdata/contour_demo.csv)
- [Volumetric spatial CSV — 72 rows across 2 volumes](https://raw.githubusercontent.com/sagnikbhadury/SpatialMethodsWorkbench/master/inst/extdata/spatial_3d_demo.csv)

Download a CSV, turn off **Use synthetic spatial demonstration**, upload the file, and follow the mapping table in the relevant vignette. All files can also be regenerated from their fixed seeds with `demo_spatial_data()`, `demo_image_data()`, `demo_contour_data()`, and `demo_spatial_3d_data()`.

The fast examples use 19 permutations, 25 bootstrap replicates, small neural networks, and compact feature matrices. Those settings verify software execution; they are not default recommendations for a scientific analysis. Public Bayesian engines are opt-in because their runtime depends on hardware and requested posterior sampling.

## How the public forks are represented

| Public repository or fork | Workbench representation | Reproduced in |
|---|---|---|
| `sagnikbhadury/ISPAT` | Exact public R engine adapter | Spatial/network vignette |
| `sagnikbhadury/GP-GHS` | Exact public R engine adapter | Spatial/network vignette |
| `sagnikbhadury/ISPAT-3D` | Exact public R engine adapter | Imaging/shape/3D vignette |
| `sagnikbhadury/AI-For-Beginners` | Source credited for the executable spatial AI workflow design; not imported as an R package | Prediction/clustering/mediation vignette |
| `sagnikbhadury/ML-For-Beginners` | Source credited for executable prediction and clustering adaptations; not imported as an R package | Prediction/clustering/mediation vignette |
| `sagnikbhadury/DeepUnderstandingOfDeepLearning` | Source credited for the executable shallow-neural adaptation; not imported as an R package | Prediction/clustering/mediation vignette |
| `sagnikbhadury/SV-NN` | Related full method/software cited alongside the implemented fast scalar-on-image screen | Imaging/shape/3D vignette |
| `sagnikbhadury/STCAR` | Related full method/software cited alongside the implemented fast scalar-on-image screen | Imaging/shape/3D vignette |
| `sagnikbhadury/SBLF` | Related full method/software cited alongside the implemented fast image-to-image screen | Imaging/shape/3D vignette |

Only the first three entries are directly called package engines. The remaining forks are represented through executable adaptations or explicitly labeled related-method workflows. This distinction prevents a screening implementation from being misreported as an exact reproduction of another repository's algorithm.

The vignettes describe execution and numerical output structure. They do not interpret a user's results or assign biological, clinical, or causal meaning.
