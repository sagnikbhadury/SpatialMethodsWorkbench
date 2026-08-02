analysis_registry <- function() {
  list(
    spatial_qc = list(
      label = "Spatial exploration & autocorrelation",
      short = "Map a numeric feature and test whether nearby observations have similar values.",
      family = "Spatial diagnostics", runtime = "Fast", engine = "Built-in validated workflow",
      needs = c("x", "y", "numeric_feature"), runner = "run_spatial_qc"
    ),
    conditional_network = list(
      label = "Region-aware conditional network",
      short = "Estimate regularized partial-correlation networks overall or within tissue regions.",
      family = "Network estimation", runtime = "Fast", engine = "Built-in exploratory workflow",
      needs = c("numeric_features"), runner = "run_conditional_network"
    ),
    ispat = list(
      label = "ISPat 2D Bayesian network",
      short = "Estimate shared and region-specific networks after spatial adjustment.",
      family = "Network estimation", runtime = "Intensive", engine = "ISPat",
      needs = c("x", "y", "zone", "numeric_features"), package = "ISPAT", runner = "run_ispat"
    ),
    gpghs = list(
      label = "Spatially varying GP–horseshoe network",
      short = "Estimate spatially varying edge surfaces with group-horseshoe shrinkage.",
      family = "Network estimation", runtime = "Intensive", engine = "GPGHS",
      needs = c("x", "y", "numeric_features"), package = "GPGHS", runner = "run_gpghs"
    ),
    ispat3d = list(
      label = "Volumetric ISPat 3D network",
      short = "Estimate shared and zone-specific networks in 3D or serial-section tissue data.",
      family = "Volumetric analysis", runtime = "Intensive", engine = "ISPAT3D",
      needs = c("x", "y", "z", "zone", "numeric_features"), package = "ISPAT3D", runner = "run_ispat3d"
    ),
    spatial_ml = list(
      label = "Spatial machine-learning prediction",
      short = "Predict a selected outcome with spatial coordinates and measured features using held-out validation.",
      family = "Prediction", runtime = "Fast", engine = "Built-in ML workflow",
      needs = c("outcome", "numeric_features"), runner = "run_spatial_ml"
    ),
    mediation = list(
      label = "Mediation with sensitivity diagnostics",
      short = "Estimate a model-based indirect effect with bootstrap uncertainty.",
      family = "Causal analysis", runtime = "Moderate", engine = "Built-in mediation workflow",
      needs = c("exposure", "mediator", "outcome"), runner = "run_mediation"
    ),
    shape_pca = list(
      label = "Tumor contour / shape PCA",
      short = "Align landmark coordinates and summarize shape variation with principal components.",
      family = "Imaging", runtime = "Fast", engine = "R implementation of public shape workflow",
      needs = c("id", "x", "y"), runner = "run_shape_pca"
    )
  )
}

method_compatibility <- function(method, data, mapping, features) {
  needs <- method$needs
  reason <- character()
  for (field in intersect(needs, c("x", "y", "z", "zone", "id", "outcome", "exposure", "mediator"))) {
    if (is.null(mapping[[field]]) || !mapping[[field]] %in% names(data)) reason <- c(reason, paste("Select", field))
  }
  if ("numeric_feature" %in% needs && length(features) < 1L) reason <- c(reason, "Select a numeric feature")
  if ("numeric_features" %in% needs && length(features) < 3L) reason <- c(reason, "Select at least three numeric features")
  package_ready <- is.null(method$package) || requireNamespace(method$package, quietly = TRUE)
  if (!package_ready) reason <- c(reason, paste(method$package, "engine is not installed on this server"))
  list(compatible = !length(reason), package_ready = package_ready, reason = paste(reason, collapse = "; "))
}

recommend_methods <- function(data, mapping, features) {
  registry <- analysis_registry()
  lapply(registry, method_compatibility, data = data, mapping = mapping, features = features)
}

workbench_citation <- function() {
  paste0(
    "Bhadury, S. (2026). Spatial Methods Workbench (Version 0.1.1) ",
    "[Computer software]. https://doi.org/10.5281/zenodo.21763607"
  )
}

method_citations <- function(id) {
  citations <- c(workbench_citation())
  if (id %in% c("ispat", "ispat3d")) {
    citations <- c(citations,
      "Bhadury, S. et al. (2026). ISPat: a flexible framework for Bayesian analysis of spatially variable tumor microenvironments. Scientific Reports. https://doi.org/10.1038/s41598-026-35341-8")
  }
  if (id == "gpghs") {
    citations <- c(citations,
      "Bhadury, S. GP-GHS software. https://github.com/sagnikbhadury/GP-GHS")
  }
  unique(citations)
}
