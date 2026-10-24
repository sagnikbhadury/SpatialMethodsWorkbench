#' Hilbert space Gaussian process basis on a two-dimensional domain
#'
#' Rescales the coordinates linearly to `[-1, 1]` on each axis and evaluates
#' the `m^2` tensor-product eigenfunctions of the Laplacian on the larger
#' square `[-L, L]^2`, together with the Matern spectral density at each
#' eigenfrequency. The basis is the one in equation (4) of the article and
#' the spectral weights are equation (5).
#'
#' @param coords An `n x 2` matrix of locations (any units).
#' @param m Number of basis functions per coordinate, so that a coefficient
#'   function has `m^2` basis coefficients. The article uses `m = 5` when
#'   `n` is large relative to `p`, and a smaller `m` (2 for 100 genes at a few
#'   hundred spots, 4 for 15 cell types on a 15 x 15 grid) so that the number
#'   of coefficients of a regression, `(p - 1) m^2`, stays comparable to `n`.
#' @param nu Matern smoothness (1.5 in the article).
#' @param rho Range on the rescaled coordinates; `NULL` (the default) uses
#'   the median pairwise distance between the rescaled locations.
#' @param alpha Marginal standard deviation of the Gaussian process prior
#'   (absorbed into the horseshoe scales by [gp_group_horseshoe_sampler()], so
#'   its value does not matter there).
#' @param L_factor Boundary `L` of the expansion domain on the rescaled
#'   coordinates (1.5 in the article).
#' @param verbose Print the range that was chosen.
#' @return A list with `phi` (the `n x m^2` basis matrix), `S_diag` (the
#'   `m^2` spectral densities), `omega_sq` (the squared eigenfrequencies),
#'   `m`, `m2`, `rho`, `nu` and `coords_norm` (the rescaled coordinates).
#' @examples
#' coords <- matrix(runif(200), 100, 2)
#' hs <- build_hsgp_basis(coords, m = 3)
#' dim(hs$phi)
#' @export
build_hsgp_basis <- function(coords, m = 5, nu = 1.5,
                             rho = NULL, alpha = 1, L_factor = 1.5,
                             verbose = FALSE) {
  coords <- as.matrix(coords)
  if (ncol(coords) != 2) stop("coords must have two columns")
  n  <- nrow(coords)
  m2 <- m^2
  coords_norm <- apply(coords, 2, function(x) {
    r <- range(x)
    if (diff(r) == 0) return(rep(0, length(x)))
    2 * (x - r[1]) / (r[2] - r[1]) - 1
  })
  L <- L_factor
  if (is.null(rho)) {
    D   <- as.matrix(stats::dist(coords_norm))
    rho <- stats::median(D[upper.tri(D)])
    if (verbose) cat("  Auto length-scale rho =", round(rho, 4), "\n")
  }
  phi_1d <- function(x, m) {
    sapply(1:m, function(j) sin(j * pi * (x + L) / (2 * L)) / sqrt(L))
  }
  lambda_1d <- function(m) sapply(1:m, function(j) (j * pi / (2 * L))^2)
  phi_x <- phi_1d(coords_norm[, 1], m)
  phi_y <- phi_1d(coords_norm[, 2], m)
  lam   <- lambda_1d(m)
  phi    <- matrix(0, n, m2)
  S_diag <- rep(0, m2)
  omega_sq_vec <- rep(0, m2)
  idx <- 1
  for (i in 1:m) {
    for (j in 1:m) {
      phi[, idx]        <- phi_x[, i] * phi_y[, j]
      omega_sq_vec[idx] <- lam[i] + lam[j]
      S_diag[idx]       <- matern_spectral_density(omega_sq_vec[idx], nu = nu,
                                                   rho = rho, alpha = alpha)
      idx <- idx + 1
    }
  }
  list(phi = phi, S_diag = S_diag, omega_sq = omega_sq_vec,
       m = m, m2 = m2, rho = rho, nu = nu, L = L, coords_norm = coords_norm)
}
