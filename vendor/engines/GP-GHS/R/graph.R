#' Recover spatial edge maps from a group horseshoe sampler fit
#'
#' Summarizes a \code{\link{gp_group_horseshoe_sampler}} fit for one node
#' into per-neighbor spatial edge-strength maps and group-level shrinkage
#' diagnostics.
#'
#' @param fit Output of \code{\link{gp_group_horseshoe_sampler}}.
#' @param phi n x m2 basis matrix.
#' @param S_diag m2 spectral densities.
#' @param q Number of neighbors.
#' @param m2 Basis dimension.
#' @param nbr_names Character vector of neighbor names.
#' @param edge_threshold Kappa cutoff for edge detection (currently informational;
#'   edge activity is thresholded at \code{kappa_mean < 0.5}, see Details).
#'
#' @return A list with \code{edge_maps} (n x q posterior mean spatial edge
#'   strength), \code{edge_maps_sd} (n x q posterior SD),
#'   \code{kappa_mean}/\code{kappa_sd} (group-level shrinkage,
#'   \eqn{\kappa_j = 1/(1+\lambda_j)}; near 1 = fully shrunk, near 0 =
#'   active edge), and \code{edge_summary} (a data frame of per-neighbor
#'   diagnostics).
#' @export
summarize_node <- function(fit, phi, S_diag, q, m2, nbr_names,
                           edge_threshold = 0.1) {
  n  <- nrow(phi)
  SS <- ncol(fit$theta_store)

  ## Group kappa: kappa_j = 1/(1 + lambda_j)
  ## kappa near 1 -> fully shrunk (no edge), near 0 -> active edge
  kappa_samples <- 1 / (1 + fit$lambda_store)   # q x SS
  kappa_mean    <- rowMeans(kappa_samples)
  kappa_sd      <- apply(kappa_samples, 1, sd)

  ## Recover spatial fields: beta_j(l) = phi %*% theta_j (per posterior sample)
  edge_maps    <- matrix(0, n, q)
  edge_maps_sd <- matrix(0, n, q)

  for (j in 1:q) {
    idx           <- ((j - 1) * m2 + 1):(j * m2)
    th_j_samples  <- fit$theta_store[idx, , drop = FALSE]   # m2 x SS
    field_samples <- phi %*% th_j_samples                    # n x SS

    edge_maps[, j]    <- rowMeans(field_samples)
    edge_maps_sd[, j] <- apply(field_samples, 1, sd)
  }

  ## Edge summary
  edge_active <- kappa_mean < 0.5 #threshold Kappa's

  edge_summary <- data.frame(
    neighbor    = nbr_names,
    kappa_mean  = round(kappa_mean, 4),
    kappa_sd    = round(kappa_sd, 4),
    edge_active = edge_active,
    map_max     = round(apply(edge_maps, 2, max), 4),
    map_min     = round(apply(edge_maps, 2, min), 4),
    map_range   = round(apply(edge_maps, 2, function(x) diff(range(x))), 4),
    map_sd      = round(apply(edge_maps, 2, sd), 4)
  )

  list(
    edge_maps    = edge_maps,
    edge_maps_sd = edge_maps_sd,
    kappa_mean   = kappa_mean,
    kappa_sd     = kappa_sd,
    edge_summary = edge_summary
  )
}

#' Estimate spatially varying cell-cell interaction networks
#'
#' Main entry point: for each cell type (column of \code{Y}), regresses its
#' expression on every other cell type's expression via
#' \code{\link{gp_group_horseshoe_sampler}} using an HSGP spectral basis
#' (\code{\link{build_hsgp_basis}}) built once and shared across all
#' nodewise regressions, recovers spatial edge maps
#' (\code{\link{summarize_node}}), and assembles a (optionally symmetrized)
#' cell-cell interaction adjacency matrix.
#'
#' @param Y n x p normalized expression matrix (spots/cells x cell types).
#' @param coords n x 2 spatial coordinate matrix.
#' @param m HSGP basis functions per dimension.
#' @param nu Matern smoothness (0.5, 1.5, 2.5).
#' @param rho GP length-scale (\code{NULL} = auto from median pairwise distance).
#' @param nmc,burn,thin MCMC settings for \code{\link{gp_group_horseshoe_sampler}}.
#' @param edge_threshold Passed through to \code{\link{summarize_node}}.
#' @param symmetry \code{"AND"} (conservative: both directions must be
#'   active) or \code{"OR"} (liberal: either direction active).
#' @param n_cores Number of parallel workers (default:
#'   \code{parallel::detectCores() - 1}). Set to 1 to run sequentially
#'   (useful for debugging).
#'
#' @return A list with \code{adj} (symmetrized adjacency matrix),
#'   \code{adj_raw} (directed adjacency before symmetrization),
#'   \code{kappa} (group-level shrinkage matrix), \code{edge_maps} (named
#'   list of per-node spatial edge maps), \code{edge_summaries}, \code{phi},
#'   \code{S_diag}, \code{coords}, and \code{hsgp_params}.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 600; p <- 12
#' coords <- matrix(runif(n * 2, 0, 10), n, 2)
#' Y <- matrix(rnorm(n * p), n, p)
#' colnames(Y) <- paste0("CT", 1:p)
#'
#' result <- gp_group_horseshoe_graph(Y, coords, m = 4, nu = 1.5,
#'                                     nmc = 3000, burn = 1000, thin = 5)
#' print(result$adj)
#' plot_adj_heatmap(result)
#' }
#'
#' @export
gp_group_horseshoe_graph <- function(Y, coords,
                                     m = 5, nu = 1.5, rho = NULL,
                                     nmc = 3000, burn = 1000, thin = 5,
                                     edge_threshold = 0.1,
                                     symmetry = "AND",
                                     n_cores = NULL) {

  n  <- nrow(Y)
  p  <- ncol(Y)
  cell_types <- colnames(Y)
  if (is.null(cell_types)) cell_types <- paste0("CT", 1:p)

  ## Resolve number of cores
  max_cores <- parallel::detectCores(logical = FALSE)
  if (is.null(n_cores)) n_cores <- max(1L, max_cores - 1L)
  n_cores <- min(n_cores, p)
  cat("Using", n_cores, "core(s) for", p, "nodewise regressions.\n\n")

  ## Build HSGP basis once — shared across all workers
  cat("Building HSGP basis (m =", m, ", m^2 =", m^2, ", nu =", nu, ")...\n")
  hsgp   <- build_hsgp_basis(coords, m = m, nu = nu, rho = rho)
  phi    <- hsgp$phi
  S_diag <- hsgp$S_diag
  m2     <- hsgp$m2
  cat("  Spectral density range: [", round(min(S_diag), 6),
      ",", round(max(S_diag), 4), "]\n\n")

  ## Scale Y once
  Y_scaled <- scale(Y)

  ## Per-node worker function
  run_node <- function(s) {
    y_s       <- as.vector(Y_scaled[, s])
    X_s       <- Y_scaled[, -s, drop = FALSE]
    nbr_names <- cell_types[-s]
    q         <- ncol(X_s)

    fit  <- gp_group_horseshoe_sampler(
      y = y_s, X_s = X_s, phi = phi, S_diag = S_diag,
      nmc = nmc, burn = burn, thin = thin
    )
    summ <- summarize_node(fit, phi, S_diag, q, m2, nbr_names, edge_threshold)

    list(
      s            = s,
      nbr_idx      = (1:p)[-s],
      edge_active  = as.integer(summ$edge_summary$edge_active),
      kappa_mean   = summ$kappa_mean,
      edge_maps    = summ$edge_maps,
      edge_summary = summ$edge_summary
    )
  }

  ## Run nodewise regressions
  if (n_cores == 1L) {
    node_results <- lapply(1:p, function(s) {
      cat("=== Node", s, ":", cell_types[s], "===\n")
      res <- run_node(s)
      cat("  Active edges:", sum(res$edge_active), "/", p - 1, "\n\n")
      res
    })

  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(cl, varlist = c(
      "Y_scaled", "cell_types", "phi", "S_diag", "m2", "p",
      "nmc", "burn", "thin", "edge_threshold",
      "gp_group_horseshoe_sampler", "summarize_node",
      "matern_spectral_density", "build_hsgp_basis"
    ), envir = environment())

    parallel::clusterEvalQ(cl, library(MASS))

    cat("Dispatching", p, "nodes to", n_cores, "workers...\n")
    node_results <- parallel::parLapply(cl, 1:p, run_node)
    cat("All nodes complete.\n\n")
  }

  ## Assemble results
  adj_matrix        <- matrix(0,  p, p, dimnames = list(cell_types, cell_types))
  kappa_matrix      <- matrix(NA, p, p, dimnames = list(cell_types, cell_types))
  edge_maps_list    <- vector("list", p)
  edge_summary_list <- vector("list", p)
  names(edge_maps_list) <- names(edge_summary_list) <- cell_types

  for (res in node_results) {
    s <- res$s
    adj_matrix[s, res$nbr_idx]   <- res$edge_active
    kappa_matrix[s, res$nbr_idx] <- res$kappa_mean
    edge_maps_list[[s]]          <- res$edge_maps
    edge_summary_list[[s]]       <- res$edge_summary
  }

  ## Symmetrize
  if (symmetry == "AND") {
    adj_sym <- (adj_matrix * t(adj_matrix))
  } else {
    adj_sym <- ((adj_matrix + t(adj_matrix)) > 0) * 1
  }
  diag(adj_sym) <- 0

  list(
    adj            = adj_sym,
    adj_raw        = adj_matrix,
    kappa          = kappa_matrix,
    edge_maps      = edge_maps_list,
    edge_summaries = edge_summary_list,
    phi            = phi,
    S_diag         = S_diag,
    coords         = coords,
    hsgp_params    = list(m = m, m2 = m2, nu = nu, rho = hsgp$rho)
  )
}
