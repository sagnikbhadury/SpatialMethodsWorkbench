#' Gibbs sampler for one nodewise regression with the group horseshoe prior
#'
#' Fits the regression of one response on `q` predictors with spatially
#' varying coefficients,
#' `y = sum_j Z_j theta_j + e`, `Z_j = diag(x_j) phi`, `e ~ N(0, sigma^2 I)`,
#' under the prior
#' `theta_{jk} | lambda_j, tau, sigma^2 ~ N(0, sigma^2 tau^2 lambda_j^2 S_k)`,
#' `lambda_j ~ C+(0, 1)`, `tau ~ C+(0, 1)`, `p(sigma^2) ~ 1 / sigma^2`
#' (equations (6) to (8) of the article). All full conditional
#' distributions are closed form once the half-Cauchy priors are written as
#' inverse-gamma mixtures (Makalic and Schmidt 2016); the updates are those of
#' Xu et al. (2016) with the spectral weights inserted. One sweep draws the
#' `q m^2` basis coefficients jointly through one Cholesky factorization.
#'
#' The sampler also records, at every retained draw, the shrinkage factor of
#' each predictor group: the Piironen and Vehtari (2017) factor
#' `kappa_{jk} = 1 / (1 + tau^2 lambda_j^2 S_k ||z_{jk}||^2)` of every
#' coefficient, summarized by the least-shrunk coefficient of the group
#' (`kappa_j = min_k kappa_{jk}`, equation (9) of the article). A value near
#' one means that the whole function is shrunk to zero.
#'
#' @param y Response vector of length `n` (centered and scaled).
#' @param X_s An `n x q` matrix of predictors (centered and scaled).
#' @param phi The `n x m^2` basis matrix from [build_hsgp_basis()].
#' @param S_diag The `m^2` spectral densities from [build_hsgp_basis()].
#' @param nmc,burn,thin Number of iterations, burn-in and thinning.
#' @param grouped `TRUE` (GP-GHS) gives one local scale per predictor
#'   group; `FALSE` gives one local scale per basis coefficient (the
#'   "HS-basis" ablation of the article).
#' @param a_sigma,b_sigma Shape and rate of an inverse-gamma prior on
#'   `sigma^2`; both zero (the default) is the Jeffreys prior.
#' @param init Optional list with `theta`, `sigma_sq`, `tau_sq` and
#'   `lambda_sq` to start from (for example the `last` element of a previous
#'   fit, or dispersed values for convergence checks).
#' @param verbose Print progress every 500 iterations.
#' @return A list with the retained draws: `theta_store` (`q m^2 x S`),
#'   `lambda_store` (`lambda_j^2`), `kappa_store` (`q x S`, the group
#'   shrinkage factors), `kappa_avg_store` (the block average of the
#'   coefficient factors), `tau_store` (`tau^2`), `sigma_store` (`sigma^2`),
#'   and `last`, the state at the final iteration.
#' @references
#' Xu, Z., Schmidt, D. F., Makalic, E., Qian, G. and Hopper, J. L. (2016).
#' Bayesian grouped horseshoe regression with application to additive
#' models. *Australasian Joint Conference on Artificial Intelligence*,
#' 229-240.
#'
#' Makalic, E. and Schmidt, D. F. (2016). A simple sampler for the horseshoe
#' estimator. *IEEE Signal Processing Letters*, 23, 179-182.
#'
#' Piironen, J. and Vehtari, A. (2017). Sparsity information and
#' regularization in the horseshoe and other shrinkage priors. *Electronic
#' Journal of Statistics*, 11, 5018-5051.
#' @examples
#' set.seed(1)
#' sim <- simulate_svggm(n = 150, p = 4, pi = 0.5, seed = 1)
#' hs  <- build_hsgp_basis(sim$coords, m = 2)
#' fit <- gp_group_horseshoe_sampler(sim$Y[, 1], sim$Y[, -1], hs$phi, hs$S_diag,
#'                                   nmc = 400, burn = 200, thin = 2)
#' rowMeans(fit$kappa_store)   # shrinkage factor of each predictor
#' @export
gp_group_horseshoe_sampler <- function(y, X_s, phi, S_diag,
                                       nmc = 3000, burn = 1000, thin = 5,
                                       grouped = TRUE,
                                       a_sigma = 0, b_sigma = 0,
                                       init = NULL, verbose = FALSE) {
  y   <- as.numeric(y)
  X_s <- as.matrix(X_s)
  n  <- length(y)
  q  <- ncol(X_s)
  m2 <- ncol(phi)
  P  <- q * m2
  if (nrow(X_s) != n || nrow(phi) != n) stop("y, X_s and phi must have the same number of rows")
  if (length(S_diag) != m2) stop("S_diag must have one entry per basis function")
  effective_ss <- floor((nmc - burn) / thin)
  if (effective_ss < 1) stop("nmc, burn and thin leave no retained draws")

  Z <- matrix(0, n, P)
  for (j in 1:q) {
    idx <- ((j - 1) * m2 + 1):(j * m2)
    Z[, idx] <- X_s[, j] * phi
  }
  ZtZ <- crossprod(Z)
  Zty <- as.vector(crossprod(Z, y))
  d_Z <- diag(ZtZ)

  block  <- rep(1:q, each = m2)
  S_rep  <- rep(S_diag, q)
  loc_id <- if (grouped) block else 1:P
  n_loc  <- if (grouped) q else P
  loc_sz <- if (grouped) rep(m2, q) else rep(1, P)
  clamp  <- function(x) pmin(pmax(x, 1e-100), 1e100)

  theta     <- rep(0, P)
  sigma_sq  <- stats::var(y)
  tau_sq    <- 1
  lambda_sq <- rep(1, n_loc)
  if (!is.null(init)) {
    theta <- init$theta; sigma_sq <- init$sigma_sq
    tau_sq <- init$tau_sq; lambda_sq <- init$lambda_sq
  }

  theta_store  <- matrix(0, P, effective_ss)
  lambda_store <- matrix(0, n_loc, effective_ss)
  kappa_store  <- matrix(0, q, effective_ss)
  kappa_avg_store <- matrix(0, q, effective_ss)
  tau_store    <- numeric(effective_ss)
  sigma_store  <- numeric(effective_ss)
  iter_saved <- 0

  for (iter in seq_len(nmc)) {
    ## theta | rest ~ N(A^{-1} Z'y, sigma^2 A^{-1}), A = Z'Z + D^{-1}
    prior_prec <- 1 / (tau_sq * lambda_sq[loc_id] * S_rep)
    A <- ZtZ
    diag(A) <- diag(A) + prior_prec
    R <- tryCatch(chol(A), error = function(e) {
      diag(A) <- diag(A) + 1e-8 * mean(diag(A)); chol(A)
    })
    theta_mu <- backsolve(R, backsolve(R, Zty, transpose = TRUE))
    theta    <- theta_mu + sqrt(sigma_sq) * backsolve(R, stats::rnorm(P))

    ## lambda_j^2 | rest
    w_norm <- as.vector(rowsum(theta^2 / S_rep, loc_id, reorder = TRUE))
    nu_g   <- 1 / stats::rgamma(n_loc, 1, rate = 1 + 1 / lambda_sq)
    lambda_sq <- clamp(1 / stats::rgamma(n_loc, (loc_sz + 1) / 2,
                                         rate = 1 / nu_g + w_norm / (2 * sigma_sq * tau_sq)))

    ## tau^2 | rest
    xi          <- 1 / stats::rgamma(1, 1, rate = 1 + 1 / tau_sq)
    global_norm <- sum(w_norm / lambda_sq)
    tau_sq <- clamp(1 / stats::rgamma(1, (P + 1) / 2,
                                      rate = 1 / xi + global_norm / (2 * sigma_sq)))

    ## sigma^2 | rest (the prior on theta carries sigma^2)
    E1 <- sum((y - Z %*% theta)^2)
    E2 <- global_norm / tau_sq
    sigma_sq <- 1 / stats::rgamma(1, (n + P) / 2 + a_sigma, rate = (E1 + E2) / 2 + b_sigma)

    if (verbose && iter %% 500 == 0)
      cat("    iter:", iter, "/ sigma:", round(sqrt(sigma_sq), 3),
          "/ tau:", signif(sqrt(tau_sq), 3), "\n")

    if (iter > burn && (iter - burn) %% thin == 0) {
      iter_saved <- iter_saved + 1
      theta_store[, iter_saved]  <- theta
      lambda_store[, iter_saved] <- lambda_sq
      tau_store[iter_saved]      <- tau_sq
      sigma_store[iter_saved]    <- sigma_sq
      k_coef <- 1 / (1 + tau_sq * lambda_sq[loc_id] * S_rep * d_Z)
      k_mat  <- matrix(k_coef, m2, q)
      kappa_store[, iter_saved]     <- apply(k_mat, 2, min)
      kappa_avg_store[, iter_saved] <- colMeans(k_mat)
    }
  }
  keep <- seq_len(iter_saved)
  list(theta_store = theta_store[, keep, drop = FALSE],
       lambda_store = lambda_store[, keep, drop = FALSE],
       kappa_store = kappa_store[, keep, drop = FALSE],
       kappa_avg_store = kappa_avg_store[, keep, drop = FALSE],
       tau_store = tau_store[keep],
       sigma_store = sigma_store[keep],
       last = list(theta = theta, sigma_sq = sigma_sq, tau_sq = tau_sq, lambda_sq = lambda_sq))
}

#' Summarize one nodewise fit
#'
#' Turns the retained draws of [gp_group_horseshoe_sampler()] into the
#' posterior mean shrinkage factor of every predictor, the posterior mean
#' and standard deviation of every coefficient function at the `n`
#' locations, and an edge summary table.
#'
#' @param fit Output of [gp_group_horseshoe_sampler()].
#' @param phi,S_diag Basis matrix and spectral densities used in the fit.
#' @param q Number of predictors.
#' @param m2 Number of basis functions (`m^2`).
#' @param nbr_names Names of the predictors.
#' @param edge_threshold A predictor is selected when its posterior mean
#'   shrinkage factor is below this value (0.5 in the article).
#' @return A list with `edge_maps` and `edge_maps_sd` (`n x q` matrices of
#'   posterior means and standard deviations of the coefficient functions),
#'   `kappa_mean`, `kappa_sd`, `kappa_avg_mean` and `edge_summary` (a data
#'   frame with one row per predictor).
#' @export
summarize_node <- function(fit, phi, S_diag, q, m2, nbr_names, edge_threshold = 0.5) {
  n  <- nrow(phi)
  kappa_samples <- fit$kappa_store
  kappa_mean    <- rowMeans(kappa_samples)
  kappa_sd      <- apply(kappa_samples, 1, stats::sd)
  edge_maps    <- matrix(0, n, q)
  edge_maps_sd <- matrix(0, n, q)
  for (j in 1:q) {
    idx           <- ((j - 1) * m2 + 1):(j * m2)
    field_samples <- phi %*% fit$theta_store[idx, , drop = FALSE]
    edge_maps[, j]    <- rowMeans(field_samples)
    edge_maps_sd[, j] <- apply(field_samples, 1, stats::sd)
  }
  edge_active <- kappa_mean < edge_threshold
  edge_summary <- data.frame(
    neighbor    = nbr_names,
    kappa_mean  = round(kappa_mean, 4),
    kappa_sd    = round(kappa_sd, 4),
    edge_active = edge_active,
    map_mean    = round(colMeans(edge_maps), 4),
    map_max     = round(apply(edge_maps, 2, max), 4),
    map_min     = round(apply(edge_maps, 2, min), 4),
    map_rms     = round(sqrt(colMeans(edge_maps^2)), 4),
    stringsAsFactors = FALSE)
  list(edge_maps = edge_maps, edge_maps_sd = edge_maps_sd,
       kappa_mean = kappa_mean, kappa_sd = kappa_sd,
       kappa_avg_mean = rowMeans(fit$kappa_avg_store),
       sigma_mean = mean(sqrt(fit$sigma_store)),
       tau_mean = mean(sqrt(fit$tau_store)),
       edge_summary = edge_summary)
}
