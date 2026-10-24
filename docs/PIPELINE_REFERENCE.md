# Function-by-function pipeline reference

This manual describes every executable analysis path in Spatial Methods Workbench v0.4.0. Read the [Usage, Restrictions, and Local Installation Guide](USER_GUIDE.md) first. The toolkit produces statistical outputs but does not explain results or assign scientific meaning; interpretation requires the study design and domain context.

## Common programmatic interface

### Load the toolkit functions

From the repository root, without starting Shiny:

```r
invisible(lapply(
  list.files("R", pattern = "[.]R$", full.names = TRUE),
  source
))
```

Create synthetic cell/spot-level data:

```r
d <- demo_spatial_data(200, seed = 2026)
```

Small reproducible vignette inputs for each data family are available directly:

```r
spatial <- demo_spatial_data(72, seed = 2026)
volume <- demo_spatial_3d_data(72, volumes = 2, seed = 2026)
image <- demo_image_data(48, features = 8, seed = 2026)
contours <- demo_contour_data(6, points = 16, seed = 2026)
```

The generated CSV versions are stored under `inst/extdata/`. See the [reproducible vignette index](VIGNETTES.md) for app mappings, every analysis family, fork provenance, and bundle auditing.

### `run_analysis()`

```r
run_analysis(id, data, mapping, features, params = list())
```

This is the recommended dispatcher used by the Shiny server.

- `id`: one registry ID listed below.
- `data`: a data frame with one row per unit required by the method.
- `mapping`: named list mapping roles such as `x`, `y`, `z`, `zone`, `id`, `outcome`, `exposure`, and `mediator` to column names. Unused roles may be `NULL`.
- `features`: character vector of selected numeric feature-column names.
- `params`: named list of method controls.

Before executing, the dispatcher calls `method_compatibility()`. It then appends `analysis_id`, UTC start time, elapsed seconds, controls, mappings, feature names, and method-specific citations to the result.

### Standard result object

Successful workflows return a list containing most of:

| Field | Contents |
|---|---|
| `method` | Human-readable engine/workflow name |
| `plot` | `ggplot2` figure produced by every registered workflow |
| `table` | Primary numerical output table |
| `summary` | Named compact numerical/administrative summary |
| `model` or `raw` | Fitted R object or exact advanced-engine result |
| `matrices` | Network/covariance/partial-correlation matrices when applicable |
| `notes` | Method boundaries and assumption reminders, not study-specific interpretation |
| `citations` | Workbench and method-specific references |
| `params`, `mapping`, `features` | Reproducibility inputs |

### Export functions

```r
result_manifest(result)
write_result_bundle(result, "analysis-bundle.zip")
```

`result_manifest()` produces JSON-ready provenance. `write_result_bundle()` writes the table, figure when available, manifest, R object, citation file, and collaboration contact. It does not include the uploaded dataset or automated scientific interpretation.

## 1. Spatial exploration and autocorrelation

**Registry ID:** `spatial_qc`  
**Function:** `run_spatial_qc(data, mapping, features, params)`

### Purpose

Display one numeric feature over x/y coordinates and calculate global Moran's I under a symmetric k-nearest-neighbor graph with a permutation p-value.

### Required input

- numeric x and y columns;
- one selected numeric feature;
- at least 10 complete observations with nonzero feature variance.

### Controls

- `feature`: selected feature name; defaults to the first selected feature;
- `neighbors`: k in the neighbor graph; UI range 3–20, default 8;
- `permutations`: permutation count; UI range 99–999, default 199;
- `seed`: reproducibility seed.

For more than 4,000 complete rows, the implementation samples 4,000 rows using the seed to control interactive memory/time.

### Output

- spatial point map;
- table containing complete observations, k, Moran's I, and permutation p-value;
- stored null permutation distribution in the internal test object only when called directly through lower-level helpers.

### Boundaries

The result depends on k, coordinate scaling, sampling, boundaries, and exchangeability. The p-value is descriptive unless the design justifies permutation exchangeability. Do not interpret spatial autocorrelation as a biological mechanism.

### Example

```r
r <- run_analysis(
  "spatial_qc", d,
  mapping = list(x = "x", y = "y"),
  features = "APC",
  params = list(feature = "APC", neighbors = 8, permutations = 199, seed = 1)
)
```

## 2. Region-aware conditional networks

**Registry ID:** `conditional_network`  
**Function:** `run_conditional_network(data, mapping, features, params)`

### Purpose

Estimate ridge-regularized precision matrices and convert them to partial correlations, either overall or separately within a mapped zone.

### Required input

- at least three numeric features;
- more complete observations than selected features plus two within every analyzed group;
- optional zone/region column.

### Controls

- `threshold`: absolute partial-correlation display/export threshold, default 0.15;
- `ridge`: covariance regularization multiplier, default 0.05;
- `by_zone`: if `TRUE`, estimate one matrix for each mapped zone.

### Output

- faceted partial-correlation network graph for the overall or zone-specific estimates;
- edge table with network, source, target, weight, and sign;
- all partial-correlation matrices in `matrices`.

### Boundaries

This is a fast exploratory ridge workflow, not an ISPat posterior analysis. Results can be sensitive to scaling, feature filtering, missingness, ridge strength, and threshold. Conditional association does not establish contact, signaling, or causation.

### Example

```r
features <- c("epithelial", "APC", "CTL", "Treg", "T_helper")
r <- run_analysis(
  "conditional_network", d,
  mapping = list(zone = "zone"), features = features,
  params = list(threshold = 0.15, ridge = 0.05, by_zone = TRUE)
)
```

## 3. ISPat 2D spatially adjusted networks

**Registry ID:** `ispat`  
**Function:** `run_ispat(data, mapping, features, params)`  
**Exact package and call:** `ISPAT3D` (CRAN, >= 0.3.0), `ISPAT3D::ispat3d_fit_2d()`

### Purpose

Estimate zone-specific conditional-association networks among numeric markers (for example log-transformed cell-type densities) after section-wise planar spatial adjustment. This is the matched planar fit of the current ISPAT3D package: a Matérn-3/2 Gaussian process is fitted per marker within each zone-section group with a Vecchia likelihood on spatially balanced anchors and predicted at every cell; the residual covariance of each zone is then modelled as `Σ_q = ΦΦᵀ + Λ_qΛ_qᵀ + Ψ_q` by Gaussian maximum likelihood; the full fitted zone covariance is inverted to signed partial correlations. The historical CAVI/SVI ISPAT package is not called.

### Required input

- numeric x/y coordinates;
- zone/region column with at least two zones, each containing more complete rows than selected markers;
- optional sample/section ID used as the section label (see Controls);
- at least three numeric marker features;
- installed `ISPAT3D` package (>= 0.3.0) with its compiled dependency `gpboost`.

The adapter forms `Y` as rows × markers, passes x and y with a zero third coordinate (the planar fit uses only x and y), and passes the zone labels and section labels as character vectors.

### Controls

- `rank`: shared and zone-specific factor rank (package default 5; capped at markers minus one by the engine);
- `anchor_fraction`, `anchor_min`, `anchor_max`: GP anchor budget per zone-section group (package defaults 0.10, 300, 5000; the count is capped at the group size, so small demonstrations use every cell);
- `gp_neighbors`: Vecchia neighbors (package default 15);
- `gp_maxit`: GP optimizer iterations (package default 30);
- `factor_maxit`: covariance likelihood iterations (package default 350);
- `use_sections`: `"auto"` uses the mapped ID as the section label when it has at most one distinct value per ten rows, `"id"` always uses it, `"none"` pools all rows into one section;
- `threads`: GPBoost threads per GP fit, 1 or 2;
- `seed`: anchor selection and GP seed (package default 2026);
- `threshold`: partial-correlation display/export threshold.

### Output

- one signed partial-correlation matrix per zone in `matrices`, taken directly from `fit$partial`; the shared low-rank component is not displayed because it has no separately interpreted network;
- faceted zone network graph and thresholded edge table;
- `summary` with zone count, active variables, factor rank, GP fit counts and fallbacks, and covariance convergence;
- `gp_log` with the per-variable GP status, anchor count, and timing;
- `engine` recording the ISPAT3D version and the function called;
- exact ISPAT3D fit in `raw` (`shared`, `full`, `partial`, `covariances`, `gp_log`, `active_variables`).

### Boundaries

Each zone-section group with fewer than ten cells, a constant marker, or a failed planar GP is mean-centered and recorded in `gp_log` with a non-`ok` status; the result notes report the count. Constant markers are excluded from the covariance likelihood and re-embedded as placeholders. No edge-wise p-values or intervals are produced. Zone labels are relative within a specimen. Large studies belong on controlled compute infrastructure; the public host is for demonstrations.

### Citation

Bhadury S, Peruzzi M, Acharyya S, et al. (2026). Informed spatially aware patterns for multiplexed immunofluorescence data. *Scientific Reports* 16(1). [doi:10.1038/s41598-026-35341-8](https://doi.org/10.1038/s41598-026-35341-8). Software: ISPAT3D R package 0.3.1, [CRAN](https://cran.r-project.org/package=ISPAT3D), archived at [doi:10.5281/zenodo.22798429](https://doi.org/10.5281/zenodo.22798429).

### Example

```r
r <- run_analysis(
  "ispat", d,
  mapping = list(x = "x", y = "y", zone = "zone", id = "sample_id"),
  features = c("epithelial", "APC", "CTL", "Treg"),
  params = list(rank = 5, anchor_min = 300, anchor_max = 5000,
                gp_neighbors = 15, gp_maxit = 30, factor_maxit = 350,
                use_sections = "auto", threads = 1, seed = 2026, threshold = 0.15)
)
names(r$matrices)          # zones
table(r$gp_log$status)     # GP adjustment status per zone-section-marker
```

## 4. GP–GHS spatially varying networks

**Registry ID:** `gpghs`  
**Function:** `run_gpghs(data, mapping, features, params)`  
**Exact package and call:** `GPGHS` (>= 1.0.0), `GPGHS::gpghs_graph()`

### Purpose

Call the public Gaussian-process group-horseshoe engine for spatially varying nodewise relationships and construct a selected adjacency network. Each variable is regressed on the others with spatially varying coefficients represented in a Hilbert-space Matérn basis; the amplitude of each coefficient function carries a group horseshoe prior; a predictor is selected when the posterior mean shrinkage factor of its group falls below `edge_threshold`, and the two directed decisions are combined by the symmetry rule.

### Required input

- numeric x/y coordinates;
- at least three numeric features;
- installed `GPGHS` package (>= 1.0.0).

### Controls

- `basis`: HSGP basis functions per dimension, UI range 2–8 (article: 5 in simulations, 2 for 100-gene zones);
- `nu`: Matérn smoothness 0.5, 1.5, or 2.5;
- `nmc`: MCMC iterations; `burn`: burn-in iterations; `thin`: thinning, programmatic default 4 (article: 3000/1000/5);
- `edge_threshold`: shrinkage-factor threshold for selecting a predictor, default 0.5 (the rule of the article);
- `symmetry`: `"AND"` (article) or `"OR"`;
- `min_detect`: minimum share of locations at which a variable must be nonzero to be estimable (0 keeps every non-constant column; 0.10 for genes in the article);
- `seed`: worker seed; `cores`: parallel workers.

### Output

- selected network graph whose topology follows the symmetric fitted adjacency and whose edge widths are the smaller directional RMS values from the posterior mean coefficient fields;
- selected edge table from `GPGHS::gpghs_edges()` with the RMS weight, GPGHS score (`1 - max` of the two directed shrinkage factors), sign label (positive, negative, or mixed over space), sign-change flag, directed shrinkage factors, and directed mean coefficients;
- `summary` with the number of estimable variables, selected edges, basis size, retained draws, and the edge threshold;
- exact GPGHS fit in `raw`, including `edge_maps`, `kappa`, `score`, and `settings`.

### Boundaries

Use sufficient MCMC iterations for formal work and inspect posterior behavior. Basis dimension, length scale, shrinkage, threshold, and symmetry rule can materially change selection. A coefficient field may change sign over space, so the global graph intentionally uses unsigned RMS width and does not claim a single positive or negative direction; use the sign column of the edge table and `GPGHS::plot_edge_map()` on the downloaded fit. The low UI defaults are for interactive feasibility, not universal production settings.

### Example

```r
r <- run_analysis(
  "gpghs", d,
  mapping = list(x = "x", y = "y"),
  features = c("epithelial", "APC", "CTL", "Treg"),
  params = list(basis = 4, nu = 1.5, nmc = 3000, burn = 1000, thin = 5,
                edge_threshold = 0.5, symmetry = "AND", min_detect = 0,
                seed = 1, cores = 1)
)
```

## 5. ISPat 3D volumetric networks

**Registry ID:** `ispat3d`  
**Function:** `run_ispat3d(data, mapping, features, params)`  
**Exact package and call:** `ISPAT3D` (CRAN, >= 0.3.0), `ISPAT3D::ispat3d_fit()`

### Purpose

Estimate zone-specific conditional-association networks from registered volumetric or serial-section multiplexed imaging. An anisotropic Matérn-3/2 Gaussian process is fitted per marker and zone over registered x, y, z (pooling all sections within a zone) with a Vecchia likelihood on spatially balanced anchors, predicted at every selected cell; the residual covariances are modelled by the same shared-plus-zone covariance likelihood as the 2D path, and the full zone covariances are inverted to signed partial correlations.

### Required input

- numeric x/y/z coordinates in one common spatial unit;
- zone/region column with at least two zones, each containing more complete rows than selected markers;
- optional section or volume ID used to balance anchor selection across sections (see `use_sections`);
- at least three numeric markers;
- installed `ISPAT3D` package (>= 0.3.0) with `gpboost`.

The adapter forms `Y` as rows × markers and `coords` as rows × (x, y, z).

### Controls

Identical to the ISPat 2D path: `rank`, `anchor_fraction`, `anchor_min`, `anchor_max`, `gp_neighbors`, `gp_maxit`, `factor_maxit`, `use_sections`, `threads`, `seed`, and the display `threshold`. The package defaults correspond to the manuscript analyses (rank 5, up to 5,000 anchors, 15 neighbors, 30 GP iterations).

### Output

- one signed partial-correlation matrix per zone in `matrices`, taken from `fit$partial`;
- faceted zone network graph and thresholded edge table;
- `summary`, `gp_log`, and `engine` as in the 2D path;
- complete ISPAT3D fit in `raw` and therefore in the R result object.

### Boundaries

Registration, axis orientation, z spacing, section loss, masks, and zone consistency must be resolved before fitting. The 3D fit pools sections within each zone; compare it with the matched 2D path (same rows, same covariance estimator) to see the effect of volumetric pooling, remembering that the comparison changes both coordinate dimension and section pooling. Large 3D studies belong on controlled compute infrastructure.

### Example

```r
r <- run_analysis(
  "ispat3d", volume_data,
  mapping = list(x = "x", y = "y", z = "z",
                 zone = "zone", id = "volume_id"),
  features = c("marker_1", "marker_2", "marker_3", "marker_4"),
  params = list(rank = 5, gp_neighbors = 15, gp_maxit = 30,
                use_sections = "auto", threads = 2, seed = 2026, threshold = 0.15)
)
```

## 6. Spatial machine-learning prediction

**Registry ID:** `spatial_ml`  
**Function:** `run_spatial_ml(data, mapping, features, params)`

### Purpose

Fit a transparent linear model for a numeric outcome or logistic model for a binary outcome using selected features and any mapped coordinates, evaluated on a random 25% held-out set.

### Required input

- outcome column;
- at least three selected numeric features;
- optional mapped x/y/z added as predictors;
- at least 30 complete rows;
- binary outcome for classification, or numeric outcome for regression.

### Controls

- `seed`: controls the 75/25 random split.

### Output

- regression: held-out observed/predicted plot, RMSE, MAE, R-squared, and residual table;
- classification: probability-density plot, accuracy, observed/predicted class, and probability;
- fitted `lm` or `glm` object in `model`.

### Boundaries

Random rows are not appropriate validation when observations share patients, slides, fields, sections, sites, or spatial neighborhoods. Use grouped/spatial/external validation outside this quick screen.

### Example

```r
r <- run_analysis(
  "spatial_ml", d,
  mapping = list(outcome = "disease", x = "x", y = "y"),
  features = c("epithelial", "APC", "CTL", "Treg"),
  params = list(seed = 2026)
)
```

## 7. Spatially weighted phenotype clustering

**Registry ID:** `spatial_clustering`  
**Function:** `run_spatial_clustering(data, mapping, features, params)`

### Purpose

Cluster standardized multivariate phenotypes jointly with standardized x/y coordinates, allowing the analyst to control coordinate influence.

### Required input

- x/y coordinates;
- at least three numeric features;
- at least 20 complete rows.

### Controls

- `clusters`: k, UI range 2–12;
- `spatial_weight`: multiplier applied to standardized coordinates, UI range 0–2; zero ignores location;
- `seed`: controls k-means starts.

The implementation uses 25 random starts and 100 maximum iterations.

### Output

- spatial cluster map;
- x/y/cluster table;
- fitted k-means object and between/total sum-of-squares ratio.

### Boundaries

K-means assumes Euclidean, roughly spherical partitions and depends on scaling, k, outliers, and spatial weight. Clusters are exploratory, not validated cell states, tumor regions, or clinical subtypes.

### Example

```r
r <- run_analysis(
  "spatial_clustering", d,
  mapping = list(x = "x", y = "y"),
  features = c("epithelial", "APC", "CTL", "Treg"),
  params = list(clusters = 4, spatial_weight = 0.5, seed = 2026)
)
```

## 8. Spatial shallow-neural prediction

**Registry ID:** `neural_prediction`  
**Function:** `run_neural_prediction(data, mapping, features, params)`  
**Engine:** `nnet`

### Purpose

Fit a standardized single-hidden-layer neural network for numeric regression or binary classification using selected features and optional coordinates, evaluated on a random 25% holdout.

### Required input

- numeric or binary outcome;
- at least three numeric features;
- optional x/y/z predictors;
- at least 40 complete rows;
- no more than 250 predictors.

### Controls

- `hidden_units`: UI range 2–30, default 6;
- `decay`: L2 weight decay, UI range 0–0.2, default 0.01;
- `seed`: split and initialization seed.

The fit uses a 500-iteration maximum and a 20,000-weight safety ceiling.

### Output

- held-out predictions/probabilities and figure;
- RMSE/MAE/R-squared for regression or accuracy for binary classification;
- fitted `nnet` model.

### Boundaries

This is an applied spatial adaptation of public AI/deep-learning curricula, not a Bayesian SV-NN fit. Tune only within a nested validation design. Random split and small data can yield unstable or optimistic estimates.

### Example

```r
r <- run_analysis(
  "neural_prediction", d,
  mapping = list(outcome = "disease", x = "x", y = "y"),
  features = c("epithelial", "APC", "CTL", "Treg"),
  params = list(hidden_units = 6, decay = 0.01, seed = 2026)
)
```

## 9. Wide-image scalar regression screen

**Registry ID:** `scalar_image`  
**Function:** `run_scalar_image(data, mapping, features, params)`

### Purpose

Predict a numeric subject-level scalar from aligned pixel, voxel, or image-derived columns using standardized ridge regression and display the largest coefficients.

### Required input

- one row per subject;
- numeric outcome;
- at least five selected aligned image-feature columns;
- at least 30 complete subjects;
- at most 5,000 selected image features.

### Controls

- `image_ridge`: positive ridge penalty, UI range 0.1–100, default 10;
- `seed`: random 75/25 split.

The solver switches to a dual formulation when features outnumber training subjects.

### Output

- all standardized feature coefficients ordered by absolute magnitude;
- top-30 coefficient bar plot;
- held-out RMSE and MAE;
- coefficient and scaling information in `model`.

### Boundaries

All feature columns must refer to identical registered locations/definitions across subjects. Coefficient ranking is a regularized screen, not posterior region selection. This is not the exact SV-NN or ST-CAR algorithm; their papers are included as related advanced methodology.

### Example

```r
pixel_names <- grep("^pixel_", names(image_data), value = TRUE)
r <- run_analysis(
  "scalar_image", image_data,
  mapping = list(outcome = "clinical_score"),
  features = pixel_names,
  params = list(image_ridge = 10, seed = 2026)
)
```

## 10. Latent image-to-image regression

**Registry ID:** `image_to_image`  
**Function:** `run_image_to_image(data, mapping, features, params)`

### Purpose

Compress aligned predictor and outcome images separately by PCA, regress outcome scores on predictor scores, reconstruct held-out outcome images, and report feature-wise prediction error.

### Required input

- one row per subject;
- at least three selected predictor-image columns prefixed `input__`;
- at least three selected outcome-image columns prefixed `output__`;
- at least 30 complete subjects;
- common registration, mask, resolution, and feature ordering.

### Controls

- `latent_factors`: UI range 2–20, constrained by sample size and both image dimensions;
- `seed`: random 75/25 split.

### Output

- held-out RMSE for every outcome-image feature;
- observed/predicted plot for the first outcome feature;
- overall RMSE;
- input PCA, outcome PCA, and latent regression coefficients in `model`.

### Boundaries

This is a fast low-rank screening workflow, not the full Bayesian SBLF posterior engine. Random split is exploratory. PCA can miss low-variance predictive structure and does not encode irregular spatial adjacency.

### Example

```r
paired <- grep("^(input__|output__)", names(paired_images), value = TRUE)
r <- run_analysis(
  "image_to_image", paired_images,
  mapping = list(), features = paired,
  params = list(latent_factors = 5, seed = 2026)
)
```

## 11. Mediation with bootstrap uncertainty

**Registry ID:** `mediation`  
**Function:** `run_mediation(data, mapping, features, params)`

### Purpose

Fit a linear mediator model and linear outcome model, estimate product-of-coefficients indirect, direct, and total effects, and bootstrap the indirect effect.

### Required input

- numeric exposure, mediator, and outcome;
- at least 30 complete rows.

### Controls

- `bootstrap`: UI range 100–1,000, default 300;
- `seed`: bootstrap seed.

### Output

- indirect, direct, total, and percentile 95% bootstrap bounds;
- bootstrap histogram;
- summary values.

### Boundaries

The function does not adjust for covariates, multilevel structure, exposure–mediator interaction, measurement error, or unmeasured confounding. Output is not automatically causal. A defensible causal use requires an explicit estimand, temporal ordering, consistency, positivity, correct specification, and exchangeability assumptions.

### Example

```r
r <- run_analysis(
  "mediation", d,
  mapping = list(exposure = "epithelial", mediator = "APC", outcome = "CTL"),
  features = character(),
  params = list(bootstrap = 500, seed = 2026)
)
```

## 12. Tumor contour / shape PCA

**Registry ID:** `shape_pca`  
**Function:** `run_shape_pca(data, mapping, features, params)`

### Purpose

Close each ordered contour, interpolate equally spaced boundary landmarks, center and size-normalize the coordinates, then perform PCA across contours.

### Required input

- contour ID;
- ordered numeric x/y boundary coordinates;
- at least five contours;
- at least six points for every contour.

### Controls

- `landmarks`: resampled landmarks, UI range 20–100, default 40.

### Output

- PC1/PC2 score table by contour ID;
- shape-space scatter plot with percent variance labels;
- full `prcomp` object.

### Boundaries

Input points must already follow the boundary. The implementation centers and size-normalizes but does not fully Procrustes-align rotation or establish homologous anatomical landmarks. Results depend on segmentation and contour ordering.

### Example

```r
r <- run_analysis(
  "shape_pca", contour_data,
  mapping = list(id = "contour_id", x = "x", y = "y"),
  features = character(),
  params = list(landmarks = 40)
)
```

## Helper and validation functions

### `analysis_registry()`

Returns the authoritative method catalog: label, short description, family, runtime class, engine, input needs, optional package, and runner name.

### `recommend_methods(data, mapping, features)`

Evaluates every registry entry with `method_compatibility()` and returns readiness, package availability, and a human-readable missing-input reason. Readiness is a technical gate, not a scientific recommendation.

### `validate_dataset(data, mapping)`

Checks basic table structure, dimensions, duplicate names, missingness, and mapped-column validity. It does not detect biological mislabeling, batch effects, spatial registration errors, leakage, confounding, or study-design problems.

### `sanitize_for_analysis(data, columns)`

Selects requested columns and retains complete cases. Therefore every current analysis is a complete-case workflow for its selected variables. Assess missingness mechanisms and potential selection bias before use.

### `partial_correlation(X, ridge)` and `matrix_edges(mat, ...)`

Low-level utilities for exploratory ridge precision matrices and thresholded edge extraction. Prefer the registered workflow unless you intentionally need a programmatic component.

### `moran_permutation(values, coords, ...)`

Low-level Moran statistic and permutation distribution. It uses the toolkit's symmetric k-nearest-neighbor matrix.

### `resample_contour(x, y, landmarks)`

Low-level boundary interpolation, centering, and size normalization for one ordered closed contour.

## Common failures and corrective action

| Message/condition | Corrective action |
|---|---|
| Select x/y/z/zone/ID/outcome | Map the required existing column in the sidebar |
| Select at least three numeric features | Add valid measured numeric columns; do not add identifiers merely to pass the gate |
| Engine is not installed | Install the exact public package locally or use a deployment that includes it |
| Not enough complete observations | Reduce justified features, resolve missingness upstream, or collect/use an adequate dataset; do not impute casually inside the app |
| Numeric outcome must vary | Verify coding, filtering, and training split |
| Binary classification only | Recode only if scientifically justified or use another validated multiclass workflow |
| Singular/failed factor or network fit | Check constant/collinear features, scale, sample size, zones, initialization, and engine diagnostics |
| Browser disconnect or timeout | Use local/HPC/queued deployment for the intensive workflow |

## Reporting checklist

When reporting a Workbench analysis, record:

- software version and DOI;
- method/publication citation from `CITATION.txt`;
- observational unit, cohort, inclusion/exclusion, and preprocessing;
- coordinate system, neighbor graph, zones, features, transformations, and missing-data handling;
- all controls, seeds, priors/tuning, iteration counts, thresholds, and symmetry rules;
- validation split unit and leakage controls;
- convergence/diagnostic and sensitivity analyses performed outside the app;
- limitations and the fact that automated interpretation was not provided.

For study-specific method selection or interpretation, contact [bhadury@umich.edu](mailto:bhadury@umich.edu?subject=Spatial%20Methods%20Workbench%20collaboration).
