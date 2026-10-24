#' Spectral density of the Matern covariance in two dimensions
#'
#' Evaluates the spectral density of the isotropic Matern covariance with
#' smoothness `nu`, range `rho` and marginal standard deviation `alpha` at
#' squared frequency `omega_sq`, for a two-dimensional domain. The Hilbert
#' space approximation of a Gaussian process gives each basis coefficient a
#' prior variance equal to this density evaluated at the eigenfrequency of
#' its basis function (Solin and Sarkka 2020; Riutort-Mayol et al. 2023).
#'
#' @param omega_sq Squared frequency (a numeric vector).
#' @param nu Matern smoothness; 1.5 (the default) is used throughout the
#'   article.
#' @param rho Range (length-scale) on the rescaled coordinates.
#' @param alpha Marginal standard deviation.
#' @return A numeric vector of spectral densities, the same length as
#'   `omega_sq`.
#' @references
#' Solin, A. and Sarkka, S. (2020). Hilbert space methods for reduced-rank
#' Gaussian process regression. *Statistics and Computing*, 30, 419-446.
#'
#' Riutort-Mayol, G., Burkner, P.-C., Andersen, M. R., Solin, A. and
#' Vehtari, A. (2023). Practical Hilbert space approximate Bayesian Gaussian
#' processes for probabilistic programming. *Statistics and Computing*, 33, 17.
#' @examples
#' matern_spectral_density(c(0.5, 2, 10), nu = 1.5, rho = 1)
#' @export
matern_spectral_density <- function(omega_sq, nu = 1.5, rho = 1, alpha = 1) {
  d <- 2
  C_num <- (2^d) * (pi^(d / 2)) * gamma(nu + d / 2) * (2 * nu)^nu
  C_den <- gamma(nu) * rho^(2 * nu)
  alpha^2 * (C_num / C_den) * (2 * nu / rho^2 + omega_sq)^(-(nu + d / 2))
}
