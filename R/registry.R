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
    spatial_clustering = list(
      label = "Spatially weighted phenotype clustering",
      short = "Discover multivariate phenotypes while controlling how strongly tissue coordinates influence the clusters.",
      family = "Spatial AI & ML", runtime = "Fast", engine = "Built-in k-means adaptation",
      needs = c("x", "y", "numeric_features"), runner = "run_spatial_clustering"
    ),
    neural_prediction = list(
      label = "Spatial shallow-neural prediction",
      short = "Fit a single-hidden-layer network to measured features and optional spatial coordinates with held-out validation.",
      family = "Spatial AI & ML", runtime = "Moderate", engine = "nnet shallow neural network",
      needs = c("outcome", "numeric_features"), package = "nnet", runner = "run_neural_prediction"
    ),
    scalar_image = list(
      label = "Wide-image scalar regression screen",
      short = "Predict a numeric subject-level outcome from selected pixel, voxel, or image-derived columns using regularized regression.",
      family = "Imaging regression", runtime = "Fast", engine = "Built-in ridge screening workflow",
      needs = c("outcome", "image_features"), runner = "run_scalar_image"
    ),
    image_to_image = list(
      label = "Latent image-to-image regression",
      short = "Link aligned predictor and outcome images through low-dimensional factors using input__ and output__ feature groups.",
      family = "Imaging regression", runtime = "Moderate", engine = "Built-in latent-factor screening workflow",
      needs = c("paired_image_features"), runner = "run_image_to_image"
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
  if ("image_features" %in% needs && length(features) < 5L) reason <- c(reason, "Select at least five numeric image features")
  if ("paired_image_features" %in% needs) {
    if (sum(startsWith(features, "input__")) < 3L || sum(startsWith(features, "output__")) < 3L) {
      reason <- c(reason, "Select at least three input__ and three output__ image columns")
    }
  }
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
    "Bhadury, S. (2026). SpatialMethodsWorkbench: Spatial Methods Workbench (Version 0.2.1) ",
    "[Computer software]. https://github.com/sagnikbhadury/SpatialMethodsWorkbench/releases/tag/v0.2.1"
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
  if (id == "spatial_clustering") {
    citations <- c(citations,
      "Microsoft. ML for Beginners curriculum. https://github.com/microsoft/ML-For-Beginners")
  }
  if (id == "neural_prediction") {
    citations <- c(citations,
      "Microsoft. AI for Beginners curriculum. https://github.com/microsoft/AI-For-Beginners",
      "Cohen, M. X. A deep understanding of deep learning. https://github.com/mikexcohen/DeepUnderstandingOfDeepLearning",
      "Venables, W. N. and Ripley, B. D. (2002). Modern Applied Statistics with S. Springer.")
  }
  if (id == "scalar_image") {
    citations <- c(citations,
      "Wu, B., Wu, K., and Kang, J. (2025). Bayesian Scalar-on-Image Regression with a Spatially Varying Single-layer Neural Network Prior. JMLR 26(116):1-38. https://jmlr.org/papers/v26/22-0246.html",
      "Xu, Y. et al. (2025). Bayesian Image Regression with Soft-thresholded Conditional Autoregressive Prior. ICLR 2025. https://openreview.net/forum?id=rnL3OafDdw")
  }
  if (id == "image_to_image") {
    citations <- c(citations,
      "Guo, C., Kang, J., and Johnson, T. D. (2022). A Spatial Bayesian Latent Factor Model for Image-on-Image Regression. Biometrics 78(1):72-84. https://doi.org/10.1111/biom.13420")
  }
  unique(citations)
}
