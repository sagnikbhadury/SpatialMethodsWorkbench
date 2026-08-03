#' Per-neighbor spatial GP edge-strength estimation (backfitting approximation)
#'
#' For each neighbor \code{j} of a node, fits a spatially varying
#' edge-strength function via a Gaussian process posterior with
#' heteroscedastic noise, using a backfitting approximation (residualizing
#' out the other neighbors via OLS) rather than the full group horseshoe
#' sampler in \code{\link{gp_group_horseshoe_sampler}}. This is a faster,
#' simpler alternative for exploratory use.
#'
#' @param y n-vector response (scaled).
#' @param X_s n x q matrix of neighbor expressions (scaled).
#' @param coords n x 2 spatial coordinate matrix.
#' @param rho Matern length-scale.
#' @param alpha GP marginal standard deviation.
#' @param sigma_sq_init Optional initial noise variance; if \code{NULL}
#'   (default), estimated from OLS residuals.
#' @param ci_majority Fraction of locations at which the pointwise 95\%
#'   credible interval must exclude zero for an edge to be marked active.
#'
#' @return A list with \code{beta_mean} (n x q posterior mean edge-strength
#'   maps), \code{beta_sd} (n x q posterior SD), \code{edge_active}
#'   (length-q 0/1 vector), and \code{sigma_sq}.
#' @export
gp_node <- function(y, X_s, coords, rho, alpha = 1,
                           sigma_sq_init = NULL, ci_majority = 0.5) {
  n  <- nrow(coords)
  q  <- ncol(X_s)

  ## Precompute kernel matrix once — shared across neighbors
  D   <- as.matrix(dist(coords))
  K   <- matern32_kernel(D, rho = rho, alpha = alpha)
  K   <- K + diag(1e-6, n)   # numerical jitter

  ## Estimate sigma_sq from residual variance (simple init)
  ## Use OLS residuals as a fast proxy
  if (is.null(sigma_sq_init)) {
    ols_fit  <- tryCatch(
      lm.fit(cbind(1, X_s), y),
      error = function(e) list(residuals = y)
    )
    sigma_sq <- max(var(ols_fit$residuals), 1e-4)
  } else {
    sigma_sq <- sigma_sq_init
  }

  beta_mean   <- matrix(0, n, q)
  beta_sd     <- matrix(0, n, q)
  edge_active <- integer(q)

  for (j in 1:q) {
    xj <- X_s[, j]   # n-vector

    ## Effective noise for this neighbor:
    ## y_j_eff = y / x_j treated as noisy obs of beta_j
    ## This is an approximation — marginalizing over other betas
    ## is intractable; we use the "backfitting" residual instead.
    ## Residual after removing other neighbors (OLS approximation):
    r_j <- tryCatch({
      if (q > 1) {
        Xother <- X_s[, -j, drop = FALSE]
        b_other <- lm.fit(Xother, y)$coefficients
        b_other[is.na(b_other)] <- 0
        y - Xother %*% b_other
      } else {
        y
      }
    }, error = function(e) y)

    ## Locations where xj is non-negligible
    ## GP posterior: beta_j | r_j ~ GP
    ## Observation model: r_j(l) = x_j(l) * beta_j(l) + noise
    ## => diag noise variance: sigma_sq / x_j(l)^2  (heteroscedastic)
    ## To avoid /0, add small floor
    xj_safe <- pmax(abs(xj), 1e-3) * sign(xj + 1e-10)

    ## Noise variance per location (heteroscedastic)
    noise_var <- sigma_sq / xj_safe^2   # n-vector

    ## GP posterior with heteroscedastic noise:
    ## Sigma_noise = diag(noise_var)
    ## K_post_inv  = (K + Sigma_noise)^{-1}
    Sigma_n    <- diag(noise_var)
    K_noise    <- K + Sigma_n

    ## Cholesky solve — O(n^3) — the main cost
    L <- tryCatch(
      chol(K_noise),
      error = function(e) chol(K_noise + diag(1e-4, n))
    )
    ## Posterior mean: K (K + Sigma_n)^{-1} r_j / x_j
    ## Effective observations: z_j = r_j / x_j
    z_j      <- r_j / xj_safe
    alpha_vec <- backsolve(L, forwardsolve(t(L), z_j))
    post_mean <- as.vector(K %*% alpha_vec)

    ## Posterior variance (diagonal): diag(K) - diag(K (K+Sn)^{-1} K)
    ## Efficient: v_k = L^{-1} K[,k], then var_k = K[k,k] - ||v_k||^2
    V <- backsolve(L, t(K))   # n x n
    post_var <- pmax(diag(K) - colSums(V^2), 1e-8)
    post_sd  <- sqrt(post_var)

    beta_mean[, j] <- post_mean
    beta_sd[, j]   <- post_sd

    ## Edge selection: 95% CI excludes 0 at majority of locations
    ci_lo <- post_mean - 1.96 * post_sd
    ci_hi <- post_mean + 1.96 * post_sd
    excludes_zero  <- (ci_lo > 0) | (ci_hi < 0)
    edge_active[j] <- as.integer(mean(excludes_zero) > ci_majority)
  }

  list(
    beta_mean   = beta_mean,    # n x q
    beta_sd     = beta_sd,      # n x q
    edge_active = edge_active,  # q-vector
    sigma_sq    = sigma_sq
  )
}
