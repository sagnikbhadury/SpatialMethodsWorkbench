# GP-GHS: Gaussian Process Group Horseshoe for spatial Cell-Cell Interactions

An R package for estimating spatially varying cell-cell interaction
networks in tissue imaging data. Each cell type's expression is regressed on
every other cell type via a Hilbert Space Gaussian Process (HSGP) basis
expansion with Matern smoothness, combined with a group horseshoe prior
across neighbors for edge sparsity and a global half-Cauchy shrinkage
parameter.

## Important: bugs fixed while packaging

The original script could not run as committed. While restructuring it into
a package, four issues were fixed (all documented inline in the source and
in `?gp_group_horseshoe_sampler` / `?build_hsgp_basis` / `?matern32_kernel`):

1. `matern_spectral_density()` and `matern32_kernel()` were called but never
   defined anywhere in the repository. Standard, well-established closed
   forms are now supplied in `R/kernels.R` (the HSGP spectral density of
   Riutort-Mayol et al. 2020, and the standard Matern-3/2 kernel).
2. `build_hsgp_basis()` referenced an undefined `lam_y` when combining
   per-axis basis eigenvalues; fixed by reusing `lam_x` (the 1D eigenvalue
   formula depends only on mode index and domain width, not the coordinate
   axis).
3. `gp_group_horseshoe_sampler()` sampled `MASS::mvrnorm(10, ...)` instead
   of `MASS::mvrnorm(1, ...)` for the `theta` update, producing a vector 10x
   too long.
4. The `tau_sq` update referenced an undefined `global_norm`; now
   accumulated per the full-conditional derivation already documented in
   the source comments.

If you had different implementations of any of these in mind, they are
isolated to `R/kernels.R` (items 1) and `R/hsgp-basis.R` / `R/sampler.R`
(items 2-4) and can be swapped out independently.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/GP-GHS")
```

`plot_all_edges()` and `plot_adj_heatmap()` additionally require the
`tidyr` and `reshape2` packages respectively; `ggplot2` is used by all
plotting functions when available, with a base R fallback for
`plot_edge_map()`.

## Usage

```r
library(GPGHS)

result <- gp_group_horseshoe_graph(
  Y = expression_matrix,   # spots/cells x cell types, normalized
  coords = coords,         # n x 2 spatial coordinates
  m = 5, nu = 1.5,
  nmc = 3000, burn = 1000, thin = 5,
  symmetry = "AND"
)

print(result$adj)
plot_adj_heatmap(result)
plot_edge_map(result, "CT1", "CT3")
```

See `?gp_group_horseshoe_graph` for full documentation.
