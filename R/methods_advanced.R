network_result <- function(method, matrices, notes = character(), raw = NULL, threshold = .15) {
  pcor <- lapply(matrices, function(Sigma) {
    Sigma <- as.matrix(Sigma); precision <- solve(Sigma + diag(1e-6, nrow(Sigma)))
    out <- -precision / sqrt(outer(diag(precision), diag(precision))); diag(out) <- 1; out
  })
  edges <- do.call(rbind, Map(function(mat, name) matrix_edges(mat, name, threshold), pcor, names(pcor)))
  first <- pcor[[1]]; heat <- expand.grid(row = rownames(first), column = colnames(first), stringsAsFactors = FALSE); heat$weight <- as.vector(first)
  plot <- ggplot2::ggplot(heat, ggplot2::aes(column, row, fill = weight)) + ggplot2::geom_tile(color = "white", linewidth = .3) +
    ggplot2::scale_fill_gradient2(low = "#d46245", mid = "#fffdf8", high = "#176b68", limits = c(-1, 1)) +
    ggplot2::coord_equal() + ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste(names(pcor)[1], "partial correlations"), x = NULL, y = NULL)
  list(method = method, plot = plot, table = edges, matrices = pcor, raw = raw,
       summary = c(networks = length(pcor), visible_edges = nrow(edges), threshold = threshold), notes = notes)
}

run_ispat <- function(data, mapping, features, params) {
  if (!requireNamespace("ISPAT", quietly = TRUE)) stop("The ISPat engine is not installed on this deployment.")
  clean <- sanitize_for_analysis(data, c(mapping$x, mapping$y, mapping$zone, features))
  S <- cbind(clean[[mapping$x]], clean[[mapping$y]], as.integer(factor(clean[[mapping$zone]])))
  Y <- t(as.matrix(clean[, features, drop = FALSE])); rownames(Y) <- features
  fit <- ISPAT::ISPAT(Y, S, ncores = params$cores %||% 1L, RefPrior = diag(length(features)),
                      use_ref = FALSE, Kernel = params$kernel %||% "Matern",
                      sGLM_method = params$spatial_fit %||% "VB", VB_MSFA = TRUE,
                      MSFA_method = params$factor_fit %||% "CAVI")
  matrices <- c(list(Shared = fit$Shared_Net), stats::setNames(fit$Cluster_Net, levels(factor(clean[[mapping$zone]]))))
  matrices <- lapply(matrices, function(x) { dimnames(x) <- list(features, features); x })
  network_result("ISPat 2D Bayesian network", matrices,
    c("Networks are conditional-association summaries from the selected public ISPat engine.", "They do not by themselves establish molecular signaling or causality."),
    fit, params$threshold %||% .15)
}

run_gpghs <- function(data, mapping, features, params) {
  if (!requireNamespace("GPGHS", quietly = TRUE)) stop("The GP–GHS engine is not installed on this deployment.")
  clean <- sanitize_for_analysis(data, c(mapping$x, mapping$y, features))
  Y <- as.matrix(clean[, features, drop = FALSE]); coords <- as.matrix(clean[, c(mapping$x, mapping$y), drop = FALSE])
  fit <- GPGHS::gp_group_horseshoe_graph(Y, coords, m = params$basis %||% 4L, nu = params$nu %||% 1.5,
    nmc = params$nmc %||% 1200L, burn = params$burn %||% 400L, thin = params$thin %||% 4L,
    edge_threshold = params$threshold %||% .15, symmetry = params$symmetry %||% "AND", n_cores = params$cores %||% 1L)
  adj <- fit$adj; storage.mode(adj) <- "numeric"
  edges <- matrix_edges(adj, "Selected edges", .5)
  plot_data <- expand.grid(row = rownames(adj), column = colnames(adj), stringsAsFactors = FALSE); plot_data$value <- as.vector(adj)
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(column, row, fill = value)) + ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "#fffdf8", high = "#176b68") + ggplot2::coord_equal() + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) + ggplot2::labs(title = "GP–GHS selected network", x = NULL, y = NULL)
  list(method = "Spatially varying GP–horseshoe network", plot = plot, table = edges, raw = fit,
       summary = c(cell_types = ncol(Y), selected_edges = nrow(edges)),
       notes = c("The adjacency matrix summarizes edges selected by the configured posterior shrinkage rule.", "Inspect spatial edge maps before drawing region-specific biological conclusions."))
}

run_ispat3d <- function(data, mapping, features, params) {
  if (!requireNamespace("ISPAT3D", quietly = TRUE)) stop("The ISPat 3D engine is not installed on this deployment.")
  clean <- sanitize_for_analysis(data, c(mapping$x, mapping$y, mapping$z, mapping$zone, mapping$id, features))
  S <- cbind(clean[[mapping$x]], clean[[mapping$y]], clean[[mapping$z]], as.integer(factor(clean[[mapping$zone]])))
  Y <- t(as.matrix(clean[, features, drop = FALSE])); rownames(Y) <- features
  spot <- if (is.null(mapping$id)) rep("volume_1", nrow(clean)) else as.character(clean[[mapping$id]])
  fit <- ISPAT3D::ISPAT_3D(Y, S, spot, ncores = params$cores %||% 1L,
                           Kernel = params$kernel %||% "Matern", MSFA_method = params$factor_fit %||% "CAVI")
  first <- fit[[1]]; matrices <- c(list(Shared = first$Shared_Net), first$Zone_Nets)
  matrices <- lapply(matrices, function(x) { dimnames(x) <- list(features, features); x })
  network_result("Volumetric ISPat 3D network", matrices,
    c("The displayed view is for the first uploaded volume/spot; the downloadable result object contains all fitted volumes.",
      "Serial-section alignment, zone definitions, and coordinate units should be validated before fitting."),
    fit, params$threshold %||% .15)
}

run_analysis <- function(id, data, mapping, features, params = list()) {
  registry <- analysis_registry()
  if (!id %in% names(registry)) stop("Unknown analysis method.")
  compatibility <- method_compatibility(registry[[id]], data, mapping, features)
  if (!compatibility$compatible) stop(compatibility$reason)
  runner <- get(registry[[id]]$runner, mode = "function")
  started <- Sys.time()
  result <- runner(data, mapping, features, params)
  result$analysis_id <- id
  result$started_at <- format(started, tz = "UTC", usetz = TRUE)
  result$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  result$params <- params
  result$mapping <- mapping
  result$features <- features
  result$citations <- method_citations(id)
  result
}
