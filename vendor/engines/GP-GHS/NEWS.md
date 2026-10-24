# GPGHS 1.0.1

* `gpghs_graph()` with `n_cores > 1` now works when the package is installed in a
  library that is not on the workers' default path (the library paths of the
  session are propagated to the workers).

# GPGHS 1.0.0

The package now contains the code of the article "Spatially Varying
Conditional-Dependence Networks for Spatial Transcriptomics and Multiplexed
Imaging" (Bhadury, Gaskins and Rao), as audited and rerun on 3 October 2026.

* The coefficient prior is scaled by the residual variance, as in the grouped
  horseshoe of Xu et al. (2016), and the full conditionals of the local scales,
  the global scale and the residual variance follow from that.
* The selection rule uses the shrinkage factor of Piironen and Vehtari (2017),
  which involves the global scale and the design; a group is summarized by its
  least-shrunk coefficient and selected when the posterior mean of that factor
  is below one half. The earlier `1 / (1 + lambda_j^2)` rule is gone.
* New high-level interface: `gpghs_graph()` (returns a `"gpghs"` object),
  `gpghs_edges()`, `coefficient_map()`, `network_overlap()`.
* New preprocessing: `st_normalize()`, `morans_i()`, `st_select_genes()` for
  spatial transcriptomics; `quadrat_counts()` for multiplexed images.
* New simulator of the spatially varying Gaussian graphical model of the
  article (`simulate_svggm()`) and `graph_metrics()`.
* Plotting: `plot_edge_map()`, `plot_network()` (igraph), `plot_adj_heatmap()`.
* The exploratory `gp_node()` backfitting function and the `matern32_kernel()`
  helper of version 0.1.0 were removed; `gp_group_horseshoe_graph()` is kept as
  an alias of `gpghs_graph()`.
* The `analysis/` folder holds the scripts that produce every table and figure
  of the article and its supplement.
