network_result <- function(method, matrices, notes = character(), raw = NULL, threshold = .15) {
  pcor <- lapply(matrices, function(Sigma) {
    Sigma <- as.matrix(Sigma); precision <- solve(Sigma + diag(1e-6, nrow(Sigma)))
    out <- -precision / sqrt(outer(diag(precision), diag(precision))); diag(out) <- 1; out
  })
  edges <- do.call(rbind, Map(function(mat, name) matrix_edges(mat, name, threshold), pcor, names(pcor)))
  plot <- network_figure(pcor, threshold = threshold, title = method,
                         top_n = 30L, signed = TRUE,
                         weight_label = "|Partial correlation|")
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
  if (!requireNamespace("GPGHS", quietly = TRUE)) stop("The GP-GHS engine is not installed on this deployment.")
  clean <- sanitize_for_analysis(data, c(mapping$x, mapping$y, features))
  Y <- as.matrix(clean[, features, drop = FALSE]); coords <- as.matrix(clean[, c(mapping$x, mapping$y), drop = FALSE])
  fit <- GPGHS::gp_group_horseshoe_graph(Y, coords, m = params$basis %||% 4L, nu = params$nu %||% 1.5,
    nmc = params$nmc %||% 1200L, burn = params$burn %||% 400L, thin = params$thin %||% 4L,
    edge_threshold = params$threshold %||% .15, symmetry = params$symmetry %||% "AND", n_cores = params$cores %||% 1L)
  adj <- fit$adj; storage.mode(adj) <- "numeric"
  coefficient_rms <- gpghs_coefficient_matrix(fit)
  edges <- matrix_edges(adj, "Selected edges", .5)
  edge_index <- cbind(match(edges$source, rownames(coefficient_rms)),
                      match(edges$target, colnames(coefficient_rms)))
  edges$weight <- coefficient_rms[edge_index]
  edges$sign <- rep("Spatially varying", nrow(edges))
  plot_weights <- coefficient_rms
  plot_weights[adj != 0 & plot_weights == 0] <- .Machine$double.eps
  plot <- network_figure(list(`Selected spatial network` = plot_weights), threshold = 0,
                         title = "GP-GHS selected spatial network", top_n = 30L,
                         signed = FALSE, weight_label = "Coefficient-field RMS")
  list(method = "Spatially varying GP-horseshoe network", plot = plot, table = edges, raw = fit,
       summary = c(cell_types = ncol(Y), selected_edges = nrow(edges)),
       notes = c("The graph topology is the fitted symmetric adjacency; edge width is the smaller directional RMS of the posterior mean coefficient fields.",
                 "A coefficient field may change sign over space, so the global graph intentionally does not assign a positive or negative edge direction.",
                 "Inspect spatial edge maps before drawing region-specific biological conclusions."))
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
