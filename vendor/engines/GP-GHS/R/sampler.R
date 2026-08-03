#' Gibbs sampler: group horseshoe + GP spectral prior for one node's regression
#'
#' Fits the regression of one node's (cell type's) expression on its
#' neighbors, y = Z theta + eps, where each neighbor's coefficient function
#' is expanded in the HSGP basis (theta_{j,k} = delta_j * w_{j,k}), with a
#' group horseshoe prior across neighbors (for edge sparsity) and a global
#' half-Cauchy shrinkage parameter, via a Gibbs sampler.
#'
#' Model:
#' \preformatted{
#'   y = Z \%*\% theta + eps,  eps ~ N(0, sigma_sq * I)
#'   theta_\{j,k\} ~ N(0, lambda_j^2 * tau^2 * S_diag[k])
#'   lambda_j ~ Half-Cauchy(0,1)   [group horseshoe]
#'   tau      ~ Half-Cauchy(0,1)   [global shrinkage]
#'   sigma_sq ~ sparse Inverse-Gamma(0.1, 0.1)
#' }
#'
#' @param y n-vector response (scaled).
#' @param X_s n x q matrix of neighbor expressions (scaled).
#' @param phi n x m2 GP basis matrix (see \code{\link{build_hsgp_basis}}).
#' @param S_diag Length-m2 spectral densities (GP prior variances).
#' @param nmc,burn,thin MCMC settings.
#'
#' @return A list with \code{theta_store} (P x saved-iterations), and the
#'   corresponding \code{lambda_store}, \code{tau_store}, \code{sigma_store}.
#'
#' @section Fixes applied when packaging:
#' Two bugs in the original script prevented it from running at all and have
#' been corrected here:
#' \itemize{
#'   \item The theta update sampled \code{MASS::mvrnorm(10, ...)} instead of
#'     \code{MASS::mvrnorm(1, ...)}, producing a vector 10 times too long for
#'     the very next matrix multiplication to fail on.
#'   \item The \code{tau_sq} update referenced an undefined \code{global_norm}.
#'     Per the derivation comment directly above that line
#'     (\code{tau^\{-2\} ~ Gamma((q*m2+1)/2, 1/xi + sum_j ||theta_j||^2_\{S^\{-1\}\} / (2*lambda_j^2))}),
#'     it is now accumulated as \code{sum_j gp_weighted_norm_j / lambda_sq[j]}
#'     across the neighbor loop that computes \code{lambda_sq}.
#' }
#' @export
gp_group_horseshoe_sampler <- function(y, X_s, phi, S_diag,
                                       nmc = 3000, burn = 1000, thin = 5) {
  n  <- length(y)
  q  <- ncol(X_s)       # number of neighbors
  m2 <- ncol(phi)       # basis dimension
  P  <- q * m2          # total coefficients

  effective_ss <- floor((nmc - burn) / thin)

  ## --- Build blocked design matrix Z (n x P) ---
  ## Z[, ((j-1)*m2+1):(j*m2)] = x_j * phi
  Z <- matrix(0, n, P)
  for (j in 1:q) {
    idx <- ((j - 1) * m2 + 1):(j * m2)
    Z[, idx] <- X_s[, j] * phi
  }

  ZtZ <- crossprod(Z)   # P x P  (precompute)
  Zty <- crossprod(Z, y)

  ## --- Initialize ---
  theta     <- rep(0, P)
  sigma_sq  <- var(y)
  tau_sq    <- 1
  lambda_sq <- rep(1, q)
  nu_g      <- rep(1, q)
  xi        <- 1

  ## --- Storage ---
  theta_store  <- matrix(0, P, effective_ss)
  lambda_store <- matrix(0, q, effective_ss)
  tau_store    <- numeric(effective_ss)
  sigma_store  <- numeric(effective_ss)

  iter_saved <- 0

  for (iter in 1:nmc) {

    ## ----------------------------------------------------------
    ## 1. Update theta | y, lambda, tau, sigma_sq
    ##    Prior precision for theta_{j,k}: 1 / (lambda_sq[j] * tau_sq * S_diag[k])
    ##    D = blkdiag(lambda_j^2 * tau^2 * diag(S_diag))  for j = 1,...,q
    ## ----------------------------------------------------------
    prior_prec <- rep(0, P)
    for (j in 1:q) {
      idx <- ((j - 1) * m2 + 1):(j * m2)
      prior_prec[idx] <- 1 / (lambda_sq[j] * tau_sq * S_diag)
    }

    A     <- ZtZ / sigma_sq + diag(prior_prec)
    A_inv <- tryCatch(
      solve(A, tol = 1e-60),
      error = function(e) solve(A + diag(1e-8, P), tol = 1e-60)
    )
    theta_mu <- A_inv %*% Zty / sigma_sq
    # Draw exactly one P-dimensional sample (the original script's
    # mvrnorm(10, ...) drew ten stacked samples, giving a vector 10x too
    # long for the theta[idx] block indexing and Z %*% theta below).
    theta    <- as.vector(MASS::mvrnorm(1, mu = theta_mu, Sigma = A_inv))

    ## ----------------------------------------------------------
    ## 2. Update lambda_sq[j] | theta, tau
     ##    FC: lambda_j^{-2} ~ Gamma((m2+1)/2,  1/nu_j + ||theta_j||^2_{S^{-1}} / (2*tau^2))
    ##    where ||theta_j||^2_{S^{-1}} = sum_k theta_{j,k}^2 / S_diag[k]
    ##
    ##    global_norm accumulates sum_j ||theta_j||^2_{S^{-1}} / lambda_sq[j]
    ##    (using each neighbor's just-updated lambda_sq[j]) for use in the
    ##    tau_sq update in step 3.
    ## ----------------------------------------------------------
    global_norm <- 0
    for (j in 1:q) {
      idx  <- ((j - 1) * m2 + 1):(j * m2)
      th_j <- theta[idx]

      # Auxiliary variable update for Half-Cauchy prior on lambda_j
      nu_g[j]   <- 1 / rgamma(1, 1, rate = 1 + 1 / lambda_sq[j])

      # GP-adjusted weighted norm (spectral structure)
      gp_weighted_norm <- sum(th_j^2 / S_diag)

      # CORRECTED rate: no sigma_sq
      s_lambda     <- 1 / nu_g[j] + gp_weighted_norm / (2 * tau_sq)
      lambda_sq[j] <- 1 / rgamma(1, (m2 + 1) / 2, rate = s_lambda)

      global_norm <- global_norm + gp_weighted_norm / lambda_sq[j]
    }

    ## ----------------------------------------------------------
    ## 3. Update tau_sq | theta, lambda
    ##    FC: tau^{-2} ~ Gamma((q*m2+1)/2,  1/xi + sum_j ||theta_j||^2_{S^{-1}} / (2*lambda_j^2))
    ## ----------------------------------------------------------
    xi <- 1 / rgamma(1, 1, rate = 1 + 1 / tau_sq)

    s_tau  <- 1 / xi + global_norm / 2
    tau_sq <- 1 / rgamma(1, (P + 1) / 2, rate = s_tau)

    ## ----------------------------------------------------------
    ## 4. Update sigma_sq | y, theta
    ##    FC: sigma^{-2} ~ Gamma((n+0.2)/2,  (||y - Z*theta||^2 + 0.2) / 2)
    ## ----------------------------------------------------------
    resid    <- y - Z %*% theta
    E1       <- sum(resid^2)
    sigma_sq <- 1 / rgamma(1, (n + 0.2) / 2, rate = (E1 + 0.2) / 2)

    if (iter %% 500 == 0) cat("    iter:", iter, "/ sigma:", round(sqrt(sigma_sq), 3),
                              "/ tau:", round(sqrt(tau_sq), 4), "\n")

    ## --- Store ---
    if (iter > burn && (iter - burn) %% thin == 0) {
      iter_saved <- iter_saved + 1
      theta_store[, iter_saved]  <- theta
      lambda_store[, iter_saved] <- lambda_sq
      tau_store[iter_saved]      <- tau_sq
      sigma_store[iter_saved]    <- sigma_sq
    }
  }

  list(
    theta_store  = theta_store[, 1:iter_saved],
    lambda_store = lambda_store[, 1:iter_saved],
    tau_store    = tau_store[1:iter_saved],
    sigma_store  = sigma_store[1:iter_saved]
  )
}
