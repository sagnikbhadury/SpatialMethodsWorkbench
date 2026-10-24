# GP-GHS: spatially varying conditional-dependence networks for spatial transcriptomics and multiplexed imaging

`GPGHS` is an R package that estimates a sparse conditional-dependence network among
`p` variables observed at `n` locations of a tissue section, with an edge strength that
is a function of location. The variables can be the expression of genes at the spots of
a spatial transcriptomics array (Visium) or the local abundances of cell types in a
multiplexed image (CODEX, imaging mass cytometry, cyclic immunofluorescence). The method
uses only the values and the coordinates, so the same function fits both kinds of data.

The method is the **Gaussian process group horseshoe (GP-GHS)**:

* each variable is regressed on all the others with spatially varying coefficients
  (neighborhood selection);
* every coefficient function is a Hilbert space basis approximation to a Gaussian
  process with a Matern-3/2 covariance;
* the amplitude of each coefficient function has a horseshoe prior with one local scale
  shared by all basis coefficients of the function (a *group*), so that the evidence for
  an edge is pooled over the tissue;
* a predictor is selected when the shrinkage factor of its group is below one half, and
  the two directed decisions of a pair are combined by the AND rule.

It accompanies the article *Spatially Varying Conditional-Dependence Networks for
Spatial Transcriptomics and Multiplexed Imaging* by Sagnik Bhadury, Jeremy T. Gaskins
and Arvind Rao (submitted). The `analysis/` folder holds the scripts that produce every
table and figure of the article and its supplement.

## Contents

1. [Installation](#installation)
2. [Quick start on simulated data](#quick-start-on-simulated-data)
3. [Workflow for spatial transcriptomics (Visium)](#workflow-for-spatial-transcriptomics-visium)
4. [Workflow for multiplexed immunofluorescence or CODEX images](#workflow-for-multiplexed-immunofluorescence-or-codex-images)
5. [Reading the output](#reading-the-output)
6. [Choosing the settings](#choosing-the-settings)
7. [Checks we recommend](#checks-we-recommend)
8. [Reproducing the article](#reproducing-the-article)
9. [Function reference](#function-reference)
10. [Citation and license](#citation-and-license)

## Installation

The package is pure R (no compiled code) and depends only on base R and the `parallel`
package. `ggplot2` and `igraph` are optional, for the plots.

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/GP-GHS")

# optional, for plot_edge_map(), plot_network() and plot_adj_heatmap()
install.packages(c("ggplot2", "igraph"))
```

Requirements: R 4.1 or later. The article's analyses were run with R 4.4.1 and 4.5.2
on Windows; nothing in the package is platform specific. One nodewise regression with
100 predictors and `m = 2` takes about 25 seconds; one with 14 predictors and `m = 4`
about five seconds. The regressions of a network are independent and run in parallel
(`n_cores`).

To check the installation:

```r
library(GPGHS)
sim <- simulate_svggm(n = 200, p = 6, pi = 0.4, seed = 1)
fit <- gpghs_graph(sim$Y, sim$coords, m = 2, nmc = 600, burn = 300, thin = 3, n_cores = 1)
fit
graph_metrics(fit$adj, sim$adj)
```

## Quick start on simulated data

`simulate_svggm()` draws data from the spatially varying Gaussian graphical model of
the article: a sparse graph, and at each location a Gaussian vector whose precision
matrix changes smoothly over the tissue, so that the true coefficient functions change
sign across the domain.

```r
library(GPGHS)
set.seed(1)
sim <- simulate_svggm(n = 600, p = 15, pi = 0.3, scenario = "varying", seed = 1)

fit <- gpghs_graph(sim$Y, sim$coords,
                   m = 5,                       # 25 basis coefficients per function
                   nmc = 2000, burn = 500, thin = 5,
                   n_cores = 4, seed = 1)
fit                                   # number of edges, residual SD, global scale
graph_metrics(fit$adj, sim$adj)       # TPR, FDP, F1, MCC against the truth
head(gpghs_edges(fit))                # selected pairs, scores and signs
plot_edge_map(fit, "V1", "V2")        # one coefficient function over the domain
plot_network(fit)                     # the graph (needs igraph)
```

With `n = 600` and `p = 15` the 15 regressions take about four minutes on one core
and about as long as a single regression with one core per regression.

## Workflow for spatial transcriptomics (Visium)

This is the workflow of Section 4 of the article: one network of 100 spatially variable
genes within each pathologist-annotated zone of a section, then a replication check
across sections. The inputs are a gene-by-spot count matrix, the spot coordinates and
(optionally) a zone label per spot. For 10x Genomics Visium data these come from
`Seurat`, `scanpy` (h5ad) or the `Space Ranger` output; the article's reader for the
`CELLxGENE` h5ad files of Valdeolivas et al. (2024) is `analysis/st_prepare.R`.

### Step 1: counts to log-normalized expression

```r
library(GPGHS)
# counts: genes x spots (integer matrix or sparse dgCMatrix)
# coords: spots x 2 (array coordinates or pixel coordinates)
# zone:   one label per spot (optional)

keep   <- colSums(counts) >= 500                       # quality filter used in the article
counts <- counts[, keep]; coords <- coords[keep, ]; zone <- zone[keep]
E <- st_normalize(counts, scale_to = 1e4)              # spots x genes, log1p(counts per 10,000)
```

### Step 2: choose the genes

A section measures far more genes than a network model can use. The article keeps
genes detected in at least 10% of spots, drops mitochondrial and ribosomal protein
genes, keeps the spatially variable genes (Moran's I on the six-nearest-neighbor graph,
false discovery rate 0.05) and ranks them by a Kruskal-Wallis test across the zones.
`st_select_genes()` does all three steps; without `zone` it ranks by Moran's I.

```r
tab   <- st_select_genes(E, coords, zone = zone, n_genes = 100,
                         min_detect = 0.10, fdr = 0.05, k = 6)
genes <- tab$gene[tab$selected]
head(tab)
```

Any other screen can replace this step; the network only needs a matrix of the chosen
genes. Note that the ranking uses differences *between* zones while the networks use
the variation *within* a zone (each gene is standardized within the zone), so the
selection does not by itself create dependence among the selected genes.

### Step 3: one network per zone

Fit a separate network in each zone. The zones are biologically distinct by
construction and are interleaved at the scale of single spots, which a smooth basis
over the whole section could not follow. A zone need not be a connected region. With
100 genes and a few hundred spots use `m = 2` (four basis coefficients per function:
a level, a gradient in each direction and their interaction), so that a regression has
396 coefficients; a gene that is detected in fewer than 10% of the zone's spots is
not estimable and is set to zero (`min_detect = 0.10`).

```r
fits <- list()
for (z in unique(zone)) {
  idx <- which(zone == z)
  if (length(idx) < 150) next                          # the article skips zones below 150 spots
  fits[[z]] <- gpghs_graph(E[idx, genes], coords[idx, ],
                           m = 2, nmc = 3000, burn = 1000, thin = 5,
                           min_detect = 0.10, n_cores = 18, seed = 1)
}
fits[["tumor"]]
edges <- gpghs_edges(fits[["tumor"]])                   # selected pairs with sign and score
table(edges$sign)                                      # positive / negative / mixed
mean(edges$sign_change)                                # share whose coefficient changes sign over the zone
plot_edge_map(fits[["tumor"]], "IGHA1", "IGKC")
plot_network(fits[["tumor"]], label_top = 10)
```

A zone of 700 spots and 100 genes takes about 15 minutes on 18 cores.

### Step 4: permutation control

A false edge can come from the marginal distribution of a gene alone (many zeros).
Permute every gene independently over the spots of the zone, which keeps every
marginal distribution and removes all dependence, and refit: every edge selected in
the permuted data is a false edge. The article selected eight pairs out of 29,502 this
way, against 995 in the original data.

```r
idx <- which(zone == "tumor")
Ep  <- apply(E[idx, genes], 2, sample)
null_fit <- gpghs_graph(Ep, coords[idx, ], m = 2, nmc = 3000, burn = 1000, thin = 5,
                        min_detect = 0.10, n_cores = 18, seed = 2)
sum(null_fit$adj) / 2
```

### Step 5: replication across sections

When the study has two sections of the same tissue (serial sections, or two slides of
the same patient), fit the *same* genes in the matching zone of the second section and
compare the two networks. `network_overlap()` reports the common edges, the number
expected if the networks were unrelated, their ratio, the overlap (common edges divided
by the smaller edge count), the Spearman correlation of the pair scores and the share
of common edges with the same sign. In the article replicate sections shared about four
times the chance number of edges, but the overlap was only about 0.28 because each
section detects a different part of the same network at this sample size: a selected
edge is evidence of dependence, the absence of an edge is not evidence of its absence.

```r
fit_rep <- gpghs_graph(E2[idx2, genes], coords2[idx2, ], m = 2, nmc = 3000, burn = 1000,
                       thin = 5, min_detect = 0.10, n_cores = 18, seed = 1)
network_overlap(fits[["tumor"]], fit_rep)
```

## Workflow for multiplexed immunofluorescence or CODEX images

This is the workflow of Section 5 of the article: one network of cell types per image.
The input is a table of phenotyped cells (x, y, cell type) for each image. The package
ships one CODEX image of colorectal cancer from Schurch et al. (2020), with the cell
phenotypes supplied with that data.

### Step 1: cells to local abundance fields

An image is a point pattern, and the model needs a value of each variable at common
locations. `quadrat_counts()` divides the bounding box of the cells into an
`n_grid x n_grid` grid, counts the cells of each type in every quadrat, takes the
square root (which stabilizes the variance of a count) and standardizes within the
image. A type with fewer than `min_cells` cells is not estimable and is set to zero.

```r
library(GPGHS)
cells <- read.csv(system.file("extdata", "codex_image_21A.csv", package = "GPGHS"))
head(cells)                                            # x, y, type
table(cells$type)

q <- quadrat_counts(cells$x, cells$y, cells$type, n_grid = 15, min_cells = 5)
dim(q$Y)                                               # 225 quadrats x cell types
q$estimable
```

Why counts and not kernel density estimates: a density surface evaluated on a grid
gives neighboring values computed from the same cells, so the grid carries far fewer
independent observations than values, and the regressions become interpolations (in
the article a density grid made GP-GHS select 99.9% of the pairs). Quadrat counts share
no cells between locations.

Why 15 x 15: there is no canonical grid for a point pattern. The quadrat must be large
enough to contain several cells of the common types (a median of about eight cells per
quadrat in the article) and small enough to resolve the structures of interest, such as
a lymphoid aggregate. Choose it before fitting; the article's networks changed little
on 20 x 20 and 10 x 10 grids (82% of the decisions agreed, prevalence correlation
0.91).

### Step 2: the network of one image

With 225 locations and 15 cell types use `m = 4` (16 basis coefficients per function,
at most 224 per regression).

```r
fit <- gpghs_graph(q$Y, q$coords, m = 4, nmc = 3000, burn = 1000, thin = 5,
                   n_cores = 7, seed = 1)
fit
gpghs_edges(fit)
plot_edge_map(fit, "CD163+ macros", "CD8+ T cells")
plot_network(fit)
```

One image takes about five seconds with one core per regression.

### Step 3: many images

Fit every image, then summarize each pair by its *prevalence*, the share of images in
which the edge is selected among the images in which both cell types are estimable.
Restrict every summary to the images in which the pair is estimable: counting an image
without a cell type as an image without the edge turns a difference in abundance into
an apparent difference in dependence (Section 5 of the article).

```r
fits <- lapply(split(cells_all, cells_all$image), function(d) {
  q <- quadrat_counts(d$x, d$y, d$type, n_grid = 15, min_cells = 5, types = cell_types)
  gpghs_graph(q$Y, q$coords, m = 4, nmc = 3000, burn = 1000, thin = 5, n_cores = 7,
              seed = 1, verbose = FALSE)
})
pair_table <- do.call(rbind, lapply(names(fits), function(img) {
  e <- gpghs_edges(fits[[img]], all = TRUE); e$image <- img; e
}))
prevalence <- aggregate(selected ~ from + to, pair_table, mean)       # estimable pairs only
prevalence[order(-prevalence$selected), ][1:10, ]
```

To compare two groups of patients, compare the per-patient summaries (for example the
share of a patient's images in which the edge is selected) with the patient as the
unit, and adjust over pairs; see `analysis/crc_network_summary.R` and
`analysis/crc_group_gee.R` for the mixed-model and GEE versions used in the article.

## Reading the output

`gpghs_graph()` returns an object of class `"gpghs"`:

| Element | Meaning |
|---|---|
| `adj` | `p x p` symmetric adjacency matrix after the AND rule (1 = edge) |
| `adj_raw` | directed decisions, row = response, column = predictor |
| `kappa` | posterior mean shrinkage factor of predictor `column` in the regression of `row`; near 1 the function is shrunk to zero, near 0 its leading component is determined by the data |
| `score` | `1 - max(kappa[s, t], kappa[t, s])`, a threshold-free ranking of the pairs |
| `edge_maps`, `edge_maps_sd` | for each response, the `n x (p - 1)` posterior mean and SD of every coefficient function at the locations |
| `edge_summaries` | one data frame per response with the shrinkage factors and the mean, range and RMS of each coefficient function |
| `sigma`, `tau` | posterior mean residual SD and global scale of each regression (a residual SD near 0 means the regression is saturated; the article's fits had medians 0.80 and 0.89) |
| `estimable` | which variables were used as responses |

`gpghs_edges(fit)` lists the selected pairs with their two directed shrinkage factors,
the score, the spatial means of both coefficient functions and the sign of the edge:
`positive` when both means are positive, `negative` when both are negative, `mixed`
otherwise; `sign_change` says whether a function is positive at more than 10% of the
locations and negative at more than 10%. `coefficient_map(fit, from, to)` returns one
function with its posterior SD at every location.

Two cautions from the article. The `p` regressions are fitted separately, so the
coefficient functions should not be assembled into a precision matrix; the target is the
zero pattern. The pointwise intervals of the coefficient functions are too narrow
(coverage 0.74 to 0.83 at a nominal 0.95 in the varying scenario), because the part of a
function outside the span of the basis contributes error but no posterior variance;
the functions describe where an association is positive and where it is negative, and
their signs and shapes are recovered well, but the intervals should not be used for
inference at a single location.

## Choosing the settings

| Setting | Default | Guidance |
|---|---|---|
| `m` | 5 | Basis functions per coordinate; a function has `m^2` coefficients and a regression `(p - 1) m^2`. Keep that number comparable to `n`: `m = 5` for 15 to 25 variables at 600 locations, `m = 4` for 15 cell types on a 15 x 15 grid, `m = 2` for 100 genes at 150 to 900 spots. The cost of a sweep grows as the cube of `(p - 1) m^2`. |
| `rho` | median distance | Range on the rescaled coordinates; the default is close to one for a filled domain. |
| `L_factor` | 1.5 | Boundary of the expansion domain, so that no basis function is forced to zero at the edge of the tissue. |
| `nmc`, `burn`, `thin` | 3000, 1000, 5 | 400 retained draws. The global scale mixes slowly (effective sample size about 30 to 50 in the article); use longer chains when `tau` itself is of interest. |
| `edge_threshold` | 0.5 | The conventional horseshoe threshold; `score` gives the threshold-free ranking. |
| `symmetry` | `"AND"` | An edge needs both directed decisions. |
| `min_detect` | 0 | Share of locations at which a variable must be nonzero to be a response (0.10 for genes in the article). |

## Checks we recommend

* **Saturation.** Look at `fit$sigma`: a median residual SD near zero means the
  regressions reproduce the response and the edges are interpolation, not dependence
  (this is what happens with density surfaces on a fine grid).
* **Permutation control.** Permute every variable independently over the locations
  and refit; the number of edges selected is the number of false edges the marginal
  distributions produce.
* **Replication.** With two sections or two images of the same tissue,
  `network_overlap()` measures how much of a network is found again. Expect a ratio to
  chance well above one and an overlap well below one; edges found in both can be
  trusted, edges found once are candidates.
* **Convergence.** Run `gp_group_horseshoe_sampler()` for a few representative
  regressions with several seeds or dispersed `init` values and compare the chains of
  the shrinkage factor, the RMS of the coefficient function, `sigma` and `tau` (the
  `posterior` package gives R-hat and effective sample sizes).
* **Sensitivity.** Refit after removing a large-scale trend from every variable
  (project each column on the basis and keep the residual), with a different
  `min_cells` or `min_detect`, and with a coarser and a finer grid.

## Reproducing the article

`analysis/` contains the scripts, in the order given in `analysis/README.md`: the
sampler check against numerical integration, the simulation study (Tables 1, S1, S2,
S3, Figures S1, S2), the Visium analysis of the index section and the other 13
sections (Tables 2, 3, S4, Figures 2, 3, 4, S3), the CODEX analysis with its sensitivity
checks (Table 4, Figures S4, S5, Sections S4.3 to S4.5) and the workflow figure. The
scripts source `analysis/gp_group_horseshoe_spatial_graph.R`, which is the same code as
the package functions. The Visium data are available from the CZ CELLxGENE collection of
Valdeolivas et al. (2024) and the CODEX cell table from Schurch et al. (2020); the
scripts document the files read.

## Function reference

| Function | Purpose |
|---|---|
| `gpghs_graph()` | Fit the network (main function); `gp_group_horseshoe_graph()` is an alias |
| `gpghs_edges()`, `coefficient_map()` | Selected edges with signs; one coefficient function |
| `plot_edge_map()`, `plot_network()`, `plot_adj_heatmap()` | Plots |
| `st_normalize()`, `morans_i()`, `st_select_genes()`, `spatial_knn_weights()` | Spatial transcriptomics preprocessing |
| `quadrat_counts()` | Multiplexed-image preprocessing |
| `network_overlap()` | Replication summary between two fits |
| `simulate_svggm()`, `scale_free_adjacency()`, `graph_metrics()` | Simulation and evaluation |
| `build_hsgp_basis()`, `matern_spectral_density()` | The basis and the spectral weights |
| `gp_group_horseshoe_sampler()`, `summarize_node()` | One nodewise regression (lower level) |

Every function has a help page (`?gpghs_graph`).

## Citation and license

Bhadury, S., Gaskins, J. T. and Rao, A. (2026). Spatially varying conditional-dependence
networks for spatial transcriptomics and multiplexed imaging. Submitted.

The method builds on the grouped horseshoe of Xu et al. (2016), the Hilbert space
Gaussian process approximation of Solin and Sarkka (2020) and Riutort-Mayol et al.
(2023), and the shrinkage factor of Piironen and Vehtari (2017).

MIT license (see `LICENSE`). The example image in `inst/extdata` is from Schurch, C. M.
et al. (2020), Coordinated cellular neighborhoods orchestrate antitumoral immunity at the
colorectal cancer invasive front, *Cell* 182, 1341-1359.
