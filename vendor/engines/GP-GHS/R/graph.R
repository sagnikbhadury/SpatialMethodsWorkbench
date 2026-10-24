#' Estimate a spatially varying conditional-dependence network
#'
#' Fits the nodewise regression of every variable on all the others with the
#' Gaussian process group horseshoe prior (GP-GHS) and combines the directed
#' decisions into a graph. This is the main function of the package: its
#' input is a matrix of `p` variables at `n` locations and the coordinates
#' of the locations, and nothing else about the measurement technology.
#'
#' Each column of `Y` is centered and scaled within the section. A column
#' with (near) zero variance, or one that is nonzero at fewer than
#' `min_detect` of the locations, is not estimable: it is not used as a
#' response and its column is set to zero, so that it contributes nothing as
#' a predictor (as in the gene and cell-type applications of the article).
#' The `p` regressions share no parameters and run in parallel.
#'
#' @param Y An `n x p` numeric matrix: variables (genes, cell types) in
#'   columns, locations in rows. Column names are used as node names.
#' @param coords An `n x 2` matrix of locations.
#' @param m Basis functions per coordinate; see [build_hsgp_basis()].
#' @param nu,rho,L_factor Matern smoothness, range and boundary; see
#'   [build_hsgp_basis()].
#' @param nmc,burn,thin MCMC settings (3,000 / 1,000 / 5 in the
#'   applications of the article, 2,000 / 500 / 5 in the simulations).
#' @param edge_threshold Shrinkage-factor threshold for selecting a predictor
#'   (0.5).
#' @param symmetry `"AND"` (an edge needs both directed decisions; the rule
#'   of the article) or `"OR"`.
#' @param min_detect Minimum share of locations at which a variable must be
#'   nonzero to be estimable (0.10 for genes in the article; set to 0 to
#'   keep every column with positive variance).
#' @param n_cores Number of parallel workers (default: physical cores minus
#'   one; 1 runs sequentially).
#' @param seed Optional integer; worker `s` uses `seed + s`.
#' @param verbose Print progress.
#' @return An object of class `"gpghs"`: a list with `adj` (the `p x p`
#'   symmetrized adjacency matrix), `adj_raw` (directed decisions, row =
#'   response), `kappa` (posterior mean shrinkage factor of predictor
#'   `column` in the regression of `row`), `score` (`1 - max` of the two
#'   directed factors), `edge_maps` (a list over responses of `n x (p - 1)`
#'   posterior mean coefficient functions), `edge_maps_sd`, `edge_summaries`,
#'   `sigma` and `tau` (posterior means per regression), `estimable`,
#'   `coords`, `phi`, `S_diag`, `hsgp_params` and `settings`.
#' @seealso [gpghs_edges()] for the selected edges with their signs,
#'   [coefficient_map()] for one coefficient function, [plot_edge_map()],
#'   [plot_network()], [network_overlap()].
#' @examples
#' sim <- simulate_svggm(n = 200, p = 5, pi = 0.4, seed = 2)
#' fit <- gpghs_graph(sim$Y, sim$coords, m = 2, nmc = 400, burn = 200, thin = 2,
#'                    n_cores = 1)
#' fit$adj
#' sim$adj
#' gpghs_edges(fit)
#' @export
gpghs_graph <- function(Y, coords, m = 5, nu = 1.5, rho = NULL, L_factor = 1.5,
                        nmc = 3000, burn = 1000, thin = 5,
                        edge_threshold = 0.5, symmetry = c("AND", "OR"),
                        min_detect = 0, n_cores = NULL, seed = NULL,
                        verbose = TRUE) {
  symmetry <- match.arg(symmetry)
  Y <- as.matrix(Y); coords <- as.matrix(coords)
  n <- nrow(Y); p <- ncol(Y)
  if (nrow(coords) != n) stop("Y and coords must have the same number of rows")
  if (p < 2) stop("Y needs at least two columns")
  node_names <- colnames(Y)
  if (is.null(node_names)) node_names <- paste0("V", seq_len(p))
  colnames(Y) <- node_names

  ## estimability and standardization within the section
  sds <- apply(Y, 2, stats::sd)
  det <- colMeans(Y != 0)
  estimable <- is.finite(sds) & sds > 1e-10 & det >= min_detect
  if (sum(estimable) < 2) stop("fewer than two estimable variables")
  Y_scaled <- scale(Y)
  Y_scaled[, !estimable] <- 0
  Y_scaled[!is.finite(Y_scaled)] <- 0
  Y_scaled <- unclass(Y_scaled)

  max_cores <- parallel::detectCores(logical = FALSE)
  if (is.null(n_cores)) n_cores <- max(1L, max_cores - 1L)
  n_cores <- min(n_cores, sum(estimable))

  if (verbose) cat("GP-GHS:", p, "variables (", sum(estimable), "estimable ) at", n,
                   "locations; m =", m, "; ", n_cores, "core(s)\n")
  hsgp   <- build_hsgp_basis(coords, m = m, nu = nu, rho = rho, L_factor = L_factor)
  phi    <- hsgp$phi
  S_diag <- hsgp$S_diag
  m2     <- hsgp$m2

  ## The worker function lives in a plain environment that holds everything it
  ## needs, including copies of the two package functions it calls, so that the
  ## parallel workers do not have to load the package (which fails when the
  ## package is installed in a library the workers do not see, or when a
  ## different version is installed than the one running).
  wenv <- new.env(parent = globalenv())
  wenv$sampler <- gp_group_horseshoe_sampler; environment(wenv$sampler) <- globalenv()
  wenv$summarize <- summarize_node;            environment(wenv$summarize) <- globalenv()
  for (nm in c("Y_scaled", "node_names", "phi", "S_diag", "m2", "p", "nmc", "burn", "thin",
               "edge_threshold", "seed"))
    assign(nm, get(nm, envir = environment()), envir = wenv)
  run_node <- function(s) {
    if (!is.null(seed)) set.seed(seed + s)
    y_s <- as.vector(Y_scaled[, s])
    X_s <- Y_scaled[, -s, drop = FALSE]
    q   <- ncol(X_s)
    fit  <- sampler(y = y_s, X_s = X_s, phi = phi, S_diag = S_diag,
                    nmc = nmc, burn = burn, thin = thin)
    summ <- summarize(fit, phi, S_diag, q, m2, node_names[-s], edge_threshold)
    list(s = s, nbr_idx = (1:p)[-s],
         edge_active = as.integer(summ$edge_summary$edge_active),
         kappa_mean = summ$kappa_mean, edge_maps = summ$edge_maps,
         edge_maps_sd = summ$edge_maps_sd, edge_summary = summ$edge_summary,
         sigma = summ$sigma_mean, tau = summ$tau_mean)
  }
  environment(run_node) <- wenv

  nodes <- which(estimable)
  if (n_cores == 1L) {
    node_results <- lapply(nodes, function(s) {
      if (verbose) cat("  node", s, ":", node_names[s], "\n")
      run_node(s)
    })
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    node_results <- parallel::parLapply(cl, nodes, run_node)
  }

  adj_matrix   <- matrix(0, p, p, dimnames = list(node_names, node_names))
  kappa_matrix <- matrix(NA_real_, p, p, dimnames = list(node_names, node_names))
  edge_maps <- edge_maps_sd <- edge_summaries <- stats::setNames(vector("list", p), node_names)
  sigma <- tau <- stats::setNames(rep(NA_real_, p), node_names)
  for (res in node_results) {
    s <- res$s
    adj_matrix[s, res$nbr_idx]   <- res$edge_active
    kappa_matrix[s, res$nbr_idx] <- res$kappa_mean
    edge_maps[[s]] <- res$edge_maps; edge_maps_sd[[s]] <- res$edge_maps_sd
    colnames(edge_maps[[s]]) <- colnames(edge_maps_sd[[s]]) <- node_names[-s]
    edge_summaries[[s]] <- res$edge_summary
    sigma[s] <- res$sigma; tau[s] <- res$tau
  }
  ## a non-estimable variable has no regression and no edges
  adj_matrix[!estimable, ] <- 0; adj_matrix[, !estimable] <- 0
  kmax <- pmax(kappa_matrix, t(kappa_matrix))
  score <- 1 - kmax; diag(score) <- NA
  if (symmetry == "AND") adj_sym <- adj_matrix * t(adj_matrix)
  else adj_sym <- ((adj_matrix + t(adj_matrix)) > 0) * 1
  diag(adj_sym) <- 0

  structure(list(adj = adj_sym, adj_raw = adj_matrix, kappa = kappa_matrix, score = score,
                 edge_maps = edge_maps, edge_maps_sd = edge_maps_sd,
                 edge_summaries = edge_summaries, sigma = sigma, tau = tau,
                 estimable = estimable, coords = coords, phi = phi, S_diag = S_diag,
                 hsgp_params = list(m = m, m2 = m2, nu = nu, rho = hsgp$rho, L = L_factor),
                 settings = list(nmc = nmc, burn = burn, thin = thin,
                                 edge_threshold = edge_threshold, symmetry = symmetry,
                                 min_detect = min_detect)),
            class = "gpghs")
}

#' @rdname gpghs_graph
#' @param ... Passed to [gpghs_graph()].
#' @export
gp_group_horseshoe_graph <- function(Y, coords, ...) gpghs_graph(Y, coords, ...)

#' @export
print.gpghs <- function(x, ...) {
  p <- ncol(x$adj); ne <- sum(x$adj[upper.tri(x$adj)])
  cat("GP-GHS network:", p, "variables (", sum(x$estimable), "estimable ),",
      nrow(x$coords), "locations, m =", x$hsgp_params$m, "\n")
  cat("  ", ne, "edges selected by the", x$settings$symmetry, "rule out of",
      sum(x$estimable) * (sum(x$estimable) - 1) / 2, "estimable pairs\n")
  cat("   median residual SD", round(stats::median(x$sigma, na.rm = TRUE), 3),
      "| median global scale tau", round(stats::median(x$tau, na.rm = TRUE), 3), "\n")
  invisible(x)
}

#' Selected edges with their signs
#'
#' Lists the pairs selected by the symmetrization rule, with the two
#' directed shrinkage factors, the pair score and the sign of the edge. As in
#' the article, an edge is `positive` when the spatial means of both
#' coefficient functions are positive, `negative` when both are negative and
#' `mixed` otherwise; `sign_change` says whether either function is positive
#' at more than `change_share` of the locations and negative at more than
#' `change_share`.
#'
#' @param fit A `"gpghs"` object.
#' @param all If `TRUE`, return every estimable pair with its score, not only
#'   the selected ones.
#' @param change_share Share of locations used to call a sign change (0.10).
#' @return A data frame with columns `from`, `to`, `selected`, `score`,
#'   `kappa_from_to`, `kappa_to_from`, `mean_from_to`, `mean_to_from`,
#'   `sign` and `sign_change`, sorted by score.
#' @export
gpghs_edges <- function(fit, all = FALSE, change_share = 0.10) {
  nodes <- rownames(fit$adj); p <- length(nodes)
  rows <- list()
  for (s in 1:(p - 1)) for (t in (s + 1):p) {
    if (!fit$estimable[s] || !fit$estimable[t]) next
    sel <- fit$adj[s, t] == 1
    if (!all && !sel) next
    f_st <- fit$edge_maps[[s]][, nodes[t]]; f_ts <- fit$edge_maps[[t]][, nodes[s]]
    m_st <- mean(f_st); m_ts <- mean(f_ts)
    sgn <- if (m_st > 0 && m_ts > 0) "positive" else if (m_st < 0 && m_ts < 0) "negative" else "mixed"
    flips <- function(f) mean(f > 0) > change_share && mean(f < 0) > change_share
    rows[[length(rows) + 1]] <- data.frame(
      from = nodes[s], to = nodes[t], selected = sel, score = fit$score[s, t],
      kappa_from_to = fit$kappa[s, t], kappa_to_from = fit$kappa[t, s],
      mean_from_to = m_st, mean_to_from = m_ts, sign = sgn,
      sign_change = flips(f_st) || flips(f_ts), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(-out$score), ]
}

#' Posterior mean of one coefficient function
#'
#' Returns the posterior mean and standard deviation of the coefficient of
#' predictor `to` in the regression of response `from`, at the `n` locations.
#'
#' @param fit A `"gpghs"` object.
#' @param from,to Node names (response and predictor).
#' @return A data frame with the coordinates, `mean` and `sd`.
#' @export
coefficient_map <- function(fit, from, to) {
  if (!from %in% rownames(fit$adj) || !to %in% rownames(fit$adj)) stop("unknown node name")
  if (!fit$estimable[from]) stop(from, " is not estimable")
  data.frame(x = fit$coords[, 1], y = fit$coords[, 2],
             mean = fit$edge_maps[[from]][, to], sd = fit$edge_maps_sd[[from]][, to])
}
