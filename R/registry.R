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
      label = "Spatially varying GP-horseshoe network",
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
    "[Computer software]. https://doi.org/10.5281/zenodo.21764196"
  )
}

method_citations <- function(id) {
  citations <- c(workbench_citation())
  if (id == "spatial_qc") {
    citations <- c(citations,
      "Moran, P. A. P. (1950). Notes on continuous stochastic phenomena. Biometrika 37(1-2):17-23. https://doi.org/10.1093/biomet/37.1-2.17")
  }
  if (id == "conditional_network") {
    citations <- c(citations,
      "Hoerl, A. E. and Kennard, R. W. (1970). Ridge Regression: Biased Estimation for Nonorthogonal Problems. Technometrics 12(1):55-67. https://doi.org/10.1080/00401706.1970.10488634")
  }
  if (id == "ispat") {
    citations <- c(citations,
      "Bhadury, S. et al. (2026). Informed spatially aware patterns for multiplexed immunofluorescence data. Scientific Reports 16(1). https://doi.org/10.1038/s41598-026-35341-8",
      "Bhadury, S. ISPAT software. https://github.com/sagnikbhadury/ISPAT")
  }
  if (id == "gpghs") {
    citations <- c(citations,
      "Bhadury, S., Gaskins, J. T., and Rao, A. (2026). Spatially Varying Graphical Models for Cell-Cell Interaction Networks in Multiplexed Tissue Imaging. https://doi.org/10.64898/2026.04.01.715977",
      "Bhadury, S. GP-GHS software. https://github.com/sagnikbhadury/GP-GHS")
  }
  if (id == "ispat3d") {
    citations <- c(citations,
      "Bhadury, S. (2026). ISPAT-3D software. https://github.com/sagnikbhadury/ISPAT-3D")
  }
  if (id == "spatial_ml") {
    citations <- c(citations,
      "Roberts, D. R. et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography 40(8):913-929. https://doi.org/10.1111/ecog.02881",
      "Microsoft. ML for Beginners curriculum. Upstream: https://github.com/microsoft/ML-For-Beginners; development fork: https://github.com/sagnikbhadury/ML-For-Beginners")
  }
  if (id == "spatial_clustering") {
    citations <- c(citations,
      "Chavent, M., Kuentz-Simonet, V., Labenne, A., and Saracco, J. (2018). ClustGeo: an R package for hierarchical clustering with spatial constraints. Computational Statistics 33(4):1799-1822. https://doi.org/10.1007/s00180-018-0791-1",
      "Microsoft. ML for Beginners curriculum. Upstream: https://github.com/microsoft/ML-For-Beginners; development fork: https://github.com/sagnikbhadury/ML-For-Beginners")
  }
  if (id == "neural_prediction") {
    citations <- c(citations,
      "Microsoft. AI for Beginners curriculum. Upstream: https://github.com/microsoft/AI-For-Beginners; development fork: https://github.com/sagnikbhadury/AI-For-Beginners",
      "Cohen, M. X. A Deep Understanding of Deep Learning. Upstream: https://github.com/mikexcohen/DeepUnderstandingOfDeepLearning; development fork: https://github.com/sagnikbhadury/DeepUnderstandingOfDeepLearning",
      "Venables, W. N. and Ripley, B. D. (2002). Modern Applied Statistics with S. Springer. https://doi.org/10.1007/978-0-387-21706-2",
      "Roberts, D. R. et al. (2017). Cross-validation strategies for structured data. Ecography 40(8):913-929. https://doi.org/10.1111/ecog.02881")
  }
  if (id == "scalar_image") {
    citations <- c(citations,
      "Hoerl, A. E. and Kennard, R. W. (1970). Ridge Regression: Biased Estimation for Nonorthogonal Problems. Technometrics 12(1):55-67. https://doi.org/10.1080/00401706.1970.10488634",
      "Wu, B., Wu, K., and Kang, J. (2025). Bayesian Scalar-on-Image Regression with a Spatially Varying Single-layer Neural Network Prior. JMLR 26(116):1-38. https://www.jmlr.org/papers/v26/22-0246.html. Software upstream: https://github.com/benwu233/SV-NN; development fork: https://github.com/sagnikbhadury/SV-NN",
      "Xu, Y. and Kang, J. (2025). Bayesian Image Regression with Soft-thresholded Conditional Autoregressive Prior. ICLR 2025. https://proceedings.iclr.cc/paper_files/paper/2025/hash/f418594e90047a10f4c158f70d6701cc-Abstract-Conference.html. Software upstream: https://github.com/yuliangxu/STCAR; development fork: https://github.com/sagnikbhadury/STCAR")
  }
  if (id == "image_to_image") {
    citations <- c(citations,
      "Guo, C., Kang, J., and Johnson, T. D. (2022). A Spatial Bayesian Latent Factor Model for Image-on-Image Regression. Biometrics 78(1):72-84. https://doi.org/10.1111/biom.13420. Software upstream: https://github.com/umich-biostatistics/SBLF; development fork: https://github.com/sagnikbhadury/SBLF")
  }
  if (id == "mediation") {
    citations <- c(citations,
      "Imai, K., Keele, L., and Tingley, D. (2010). A general approach to causal mediation analysis. Psychological Methods 15(4):309-334. https://doi.org/10.1037/a0020761")
  }
  if (id == "shape_pca") {
    citations <- c(citations,
      "Dryden, I. L. and Mardia, K. V. (2016). Statistical Shape Analysis, with Applications in R (2nd ed.). Wiley. https://doi.org/10.1002/9781119072492")
  }
  unique(citations)
}
