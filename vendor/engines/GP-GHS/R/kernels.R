## NOTE ON THIS FILE
## ------------------------------------------------------------------
## matern_spectral_density() and matern32_kernel() were referenced by
## build_hsgp_basis() and gp_node() in the original script but were never
## defined anywhere in this repository, so the code could not run as
## committed. The implementations below are standard, well-established
## closed forms (not something novel supplied by guesswork):
##
##   - matern32_kernel(): the standard Matern-3/2 covariance function,
##     K(d) = alpha^2 * (1 + sqrt(3) d / rho) * exp(-sqrt(3) d / rho).
##     This is the same functional form already used inline for the
##     Matern kernel elsewhere in the author's own ISPAT package.
##
##   - matern_spectral_density(): the closed-form spectral density of the
##     Matern-nu covariance in D spatial dimensions, as used in the
##     Hilbert Space Gaussian Process (HSGP) approximation of
##     Riutort-Mayol, Buerkner, Andersen, Solin & Vehtari (2020),
##     "Practical Hilbert Space Approximate Bayesian Gaussian Processes
##     for Probabilistic Programming". D is fixed at 2 here because
##     build_hsgp_basis() is documented as operating on 2D spatial
##     coordinates and passes omega_sq = lambda_x + lambda_y (the squared
##     eigenfrequency of the 2D tensor-product Laplacian eigenbasis),
##     which is exactly the argument this formula expects.
##
## If you had a different specific formula in mind when writing the
## original call sites, replace these two functions; nothing else in the
## package depends on their internals, only on their (omega_sq/d, nu, rho,
## alpha) -> numeric interface.
## ------------------------------------------------------------------

#' Matern-3/2 covariance kernel
#'
#' Standard isotropic Matern covariance function with smoothness
#' \eqn{\nu = 3/2}, applied elementwise to a distance (matrix or vector).
#'
#' @param D A distance matrix or vector (e.g. from \code{\link[stats]{dist}}).
#' @param rho Length-scale.
#' @param alpha Marginal standard deviation (so the kernel evaluates to
#'   \code{alpha^2} at distance 0).
#'
#' @return \eqn{K(d) = \alpha^2 (1 + \sqrt3\, d/\rho) \exp(-\sqrt3\, d/\rho)},
#'   with the same shape as \code{D}.
#' @export
matern32_kernel <- function(D, rho, alpha = 1) {
  r <- sqrt(3) * D / rho
  (alpha^2) * (1 + r) * exp(-r)
}

#' Spectral density of the Matern covariance (2D Hilbert-space GP approximation)
#'
#' Closed-form spectral density of an isotropic Matern-\eqn{\nu} covariance
#' function in 2 spatial dimensions, used as the prior variance on each
#' basis-function weight in the Hilbert Space Gaussian Process (HSGP)
#' approximation of Riutort-Mayol et al. (2020).
#'
#' @param omega_sq Squared eigenfrequency
#'   (\eqn{\|\omega\|^2 = \omega_x^2 + \omega_y^2}) at which to evaluate the
#'   spectral density.
#' @param nu Matern smoothness parameter.
#' @param rho Length-scale.
#' @param alpha Marginal standard deviation of the GP.
#'
#' @return The spectral density \eqn{S(\omega)}, a non-negative number (or
#'   vector, if \code{omega_sq} is a vector).
#' @export
matern_spectral_density <- function(omega_sq, nu, rho, alpha = 1) {
  D <- 2   # build_hsgp_basis() operates on 2D spatial coordinates
  const <- (2^D * pi^(D / 2) * gamma(nu + D / 2) * (2 * nu)^nu) / (gamma(nu) * rho^(2 * nu))
  (alpha^2) * const * (2 * nu / rho^2 + omega_sq)^(-(nu + D / 2))
}
