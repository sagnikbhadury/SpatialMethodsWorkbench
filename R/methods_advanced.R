network_result <- function(method, matrices, notes = character(), raw = NULL, threshold = .15,
                           precomputed = FALSE, summary_extra = NULL) {
  pcor <- if (isTRUE(precomputed)) {
    lapply(matrices, as.matrix)
  } else {
    lapply(matrices, function(Sigma) {
      Sigma <- as.matrix(Sigma); precision <- solve(Sigma + diag(1e-6, nrow(Sigma)))
      out <- -precision / sqrt(outer(diag(precision), diag(precision))); diag(out) <- 1; out
    })
  }
  edges <- do.call(rbind, Map(function(mat, name) matrix_edges(mat, name, threshold), pcor, names(pcor)))
  plot <- network_figure(pcor, threshold = threshold, title = method,
                         top_n = 30L, signed = TRUE,
                         weight_label = "|Partial correlation|")
  list(method = method, plot = plot, table = edges, matrices = pcor, raw = raw,
       summary = c(networks = length(pcor), visible_edges = nrow(edges), threshold = threshold, summary_extra),
       notes = notes)
}

engine_version <- function(package, minimum) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s engine is not installed on this deployment.", package))
  }
  installed <- utils::packageVersion(package)
  if (installed < minimum) {
    stop(sprintf("%s %s is installed but the Workbench requires %s >= %s.", package, installed, package, minimum))
  }
  as.character(installed)
}

engine_sections <- function(clean, mapping, params) {
  n <- nrow(clean)
  rule <- params$use_sections %||% "auto"
  id <- mapping$id
  has_id <- !is.null(id) && id %in% names(clean)
  values <- if (has_id) as.character(clean[[id]]) else character()
  k <- length(unique(values))
  plausible <- has_id && k >= 2L && k <= floor(n / 10)
  use_id <- has_id && (identical(rule, "id") || (identical(rule, "auto") && plausible))
  if (!use_id) {
    note <- if (!has_id) {
      "No sample/section/volume ID was mapped, so all rows were treated as one section."
    } else if (identical(rule, "none")) {
      "All rows were pooled into one section by the selected control."
    } else {
      sprintf("The mapped ID column '%s' has %d distinct values among %d rows and was not treated as a section identifier; all rows were pooled into one section. Choose 'Use mapped ID' to override.", id, k, n)
    }
    return(list(sections = rep(1L, n), note = note, count = 1L))
  }
  list(sections = values, count = k,
       note = sprintf("The mapped ID column '%s' identified %d sections used for anchor balancing%s.", id, k,
                      if (isTRUE(params$planar)) " and for the separate planar GP fits" else ""))
}

ispat3d_engine_args <- function(params) {
  anchor_min <- as.integer(params$anchor_min %||% 300L)
  anchor_max <- max(as.integer(params$anchor_max %||% 5000L), anchor_min)
  list(
    rank = as.integer(params$rank %||% 5L),
    anchor_fraction = as.numeric(params$anchor_fraction %||% 0.10),
    anchor_min = anchor_min,
    anchor_max = anchor_max,
    neighbors = as.integer(params$gp_neighbors %||% 15L),
    gp_maxit = as.integer(params$gp_maxit %||% 30L),
    factor_maxit = as.integer(params$factor_maxit %||% 350L),
    seed = as.integer(params$seed %||% 2026L),
    threads = as.integer(min(max(as.integer(params$threads %||% 1L), 1L), 2L))
  )
}

ispat3d_network_result <- function(fit, method, features, params, notes, engine_function) {
  matrices <- lapply(fit$partial, function(x) { x <- as.matrix(x); dimnames(x) <- list(features, features); x })
  statuses <- as.character(fit$gp_log$status)
  summary_extra <- c(
    zones = length(fit$partial),
    active_variables = length(fit$active_variables),
    factor_rank = if (is.null(fit$loadings_shared)) NA_real_ else ncol(as.matrix(fit$loadings_shared)),
    gp_fits = length(statuses),
    gp_fits_ok = sum(statuses == "ok"),
    gp_fits_fallback = sum(statuses != "ok"),
    covariance_converged = as.numeric(isTRUE(fit$convergence == 0))
  )
  if (length(fit$inactive_variables)) {
    notes <- c(notes, sprintf("Constant variables excluded from the covariance likelihood and re-embedded as placeholders: %s. Do not interpret pairs involving them.",
                              paste(fit$inactive_variables, collapse = ", ")))
  }
  if (any(statuses != "ok")) {
    notes <- c(notes, sprintf("%d of %d GP adjustments did not complete a Vecchia fit and were mean-centered instead (statuses: %s). Inspect raw$gp_log before interpreting those zones.",
                              sum(statuses != "ok"), length(statuses), paste(unique(statuses[statuses != "ok"]), collapse = "; ")))
  }
  result <- network_result(method, matrices, notes, fit, params$threshold %||% .15,
                           precomputed = TRUE, summary_extra = summary_extra)
  result$engine <- list(package = "ISPAT3D", version = as.character(utils::packageVersion("ISPAT3D")),
                        call = engine_function, method = fit$method)
  result$gp_log <- fit$gp_log
  result
}

run_ispat <- function(data, mapping, features, params) {
  engine_version("ISPAT3D", "0.3.0")
  clean <- sanitize_for_analysis(data, c(mapping$x, mapping$y, mapping$zone, mapping$id, features))
  Y <- as.matrix(clean[, features, drop = FALSE]); colnames(Y) <- features
  coords <- cbind(clean[[mapping$x]], clean[[mapping$y]], 0)
  zones <- as.character(clean[[mapping$zone]])
  sections <- engine_sections(clean, mapping, c(params, list(planar = TRUE)))
  args <- ispat3d_engine_args(params)
  fit <- do.call(ISPAT3D::ispat3d_fit_2d,
                 c(list(Y = Y, coords = coords, zones = zones, sections = sections$sections), args))
  ispat3d_network_result(fit, "ISPat 2D spatially adjusted network", features, params,
    c("Networks are signed partial correlations from the full fitted zone covariances of the ISPAT3D planar fit (ispat3d_fit_2d): a Matern-3/2 Vecchia GP per variable within each zone-section group, then a shared-plus-zone factor covariance by Gaussian maximum likelihood.",
      "The shared low-rank component has no separately interpreted network, so only zone networks are displayed.",
      "Edges are conditional associations; they do not by themselves establish molecular signaling or causality. No edge-wise p-values or intervals are produced.",
      sections$note),
    "ISPAT3D::ispat3d_fit_2d")
}

run_gpghs <- function(data, mapping, features, params) {
  engine_version("GPGHS", "1.0.0")
  clean <- sanitize_for_analysis(data, c(mapping$x, mapping$y, features))
  Y <- as.matrix(clean[, features, drop = FALSE]); coords <- as.matrix(clean[, c(mapping$x, mapping$y), drop = FALSE])
  settings <- list(
    m = as.integer(params$basis %||% 4L), nu = as.numeric(params$nu %||% 1.5),
    nmc = as.integer(params$nmc %||% 1200L), burn = as.integer(params$burn %||% 400L), thin = as.integer(params$thin %||% 4L),
    edge_threshold = as.numeric(params$edge_threshold %||% .5), symmetry = params$symmetry %||% "AND",
    min_detect = as.numeric(params$min_detect %||% 0), n_cores = as.integer(params$cores %||% 1L),
    seed = as.integer(params$seed %||% 2026L))
  fit <- do.call(GPGHS::gpghs_graph, c(list(Y = Y, coords = coords), settings, list(verbose = FALSE)))
  adj <- fit$adj; storage.mode(adj) <- "numeric"
  coefficient_rms <- gpghs_coefficient_matrix(fit)
  selected <- GPGHS::gpghs_edges(fit)
  edges <- if (nrow(selected)) {
    index <- cbind(match(selected$from, rownames(coefficient_rms)), match(selected$to, colnames(coefficient_rms)))
    data.frame(network = "Selected spatial network", source = selected$from, target = selected$to,
               weight = coefficient_rms[index], score = selected$score,
               sign = selected$sign, sign_change = selected$sign_change,
               kappa_from_to = selected$kappa_from_to, kappa_to_from = selected$kappa_to_from,
               mean_from_to = selected$mean_from_to, mean_to_from = selected$mean_to_from,
               stringsAsFactors = FALSE, row.names = NULL)
  } else {
    data.frame(network = character(), source = character(), target = character(), weight = numeric(),
               score = numeric(), sign = character(), sign_change = logical(),
               kappa_from_to = numeric(), kappa_to_from = numeric(), mean_from_to = numeric(), mean_to_from = numeric())
  }
  plot_weights <- coefficient_rms
  plot_weights[adj != 0 & plot_weights == 0] <- .Machine$double.eps
  plot <- network_figure(list(`Selected spatial network` = plot_weights), threshold = 0,
                         title = "GP-GHS selected spatial network", top_n = 30L,
                         signed = FALSE, weight_label = "Coefficient-field RMS")
  list(method = "Spatially varying GP-horseshoe network", plot = plot, table = edges, raw = fit,
       engine = list(package = "GPGHS", version = as.character(utils::packageVersion("GPGHS")), call = "GPGHS::gpghs_graph"),
       summary = c(cell_types = ncol(Y), estimable_variables = sum(fit$estimable), selected_edges = nrow(edges),
                   basis_per_dimension = settings$m, retained_draws = floor((settings$nmc - settings$burn) / settings$thin),
                   edge_threshold = settings$edge_threshold),
       notes = c("A predictor is selected when the posterior mean shrinkage factor of its group falls below the edge threshold (0.5 is the rule of the GP-GHS article); the two directed decisions of a pair are combined by the chosen symmetry rule.",
                 "The graph topology is the fitted symmetric adjacency; edge width is the smaller directional RMS of the posterior mean coefficient fields.",
                 "The edge table reports the GPGHS score, sign (positive, negative, or mixed over space), and whether the coefficient field changes sign; a field may change sign over space, so the global graph does not assign a direction.",
                 "Inspect spatial edge maps (GPGHS::plot_edge_map on the downloaded fit) before drawing region-specific biological conclusions."))
}

run_ispat3d <- function(data, mapping, features, params) {
  engine_version("ISPAT3D", "0.3.0")
  clean <- sanitize_for_analysis(data, c(mapping$x, mapping$y, mapping$z, mapping$zone, mapping$id, features))
  Y <- as.matrix(clean[, features, drop = FALSE]); colnames(Y) <- features
  coords <- cbind(clean[[mapping$x]], clean[[mapping$y]], clean[[mapping$z]])
  zones <- as.character(clean[[mapping$zone]])
  sections <- engine_sections(clean, mapping, params)
  args <- ispat3d_engine_args(params)
  fit <- do.call(ISPAT3D::ispat3d_fit,
                 c(list(Y = Y, coords = coords, zones = zones, sections = sections$sections), args))
  ispat3d_network_result(fit, "Volumetric ISPat 3D network", features, params,
    c("Networks are signed partial correlations from the full fitted zone covariances of the ISPAT3D volumetric fit (ispat3d_fit): an anisotropic Matern-3/2 Vecchia GP per variable pooled within each zone over registered x, y, z, then a shared-plus-zone factor covariance by Gaussian maximum likelihood.",
      "The shared low-rank component has no separately interpreted network, so only zone networks are displayed.",
      "Serial-section registration, zone definitions, and a common coordinate unit across x, y, and z should be validated before fitting.",
      sections$note),
    "ISPAT3D::ispat3d_fit")
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
