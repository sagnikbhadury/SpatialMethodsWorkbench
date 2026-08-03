#' Build a 2D Hilbert Space Gaussian Process (HSGP) basis
#'
#' Constructs the Laplacian eigenfunction basis and associated spectral
#' densities used to approximate a 2D Matern Gaussian process, following
#' Riutort-Mayol et al. (2020). Spatial coordinates are first normalized to
#' \eqn{[-1, 1]}; the basis is a tensor product of 1D sine eigenfunctions of
#' the Dirichlet Laplacian on \eqn{[-L, L]}.
#'
#' @param coords n x 2 matrix of spatial coordinates.
#' @param m Number of basis functions per dimension (\code{m^2} total in 2D).
#' @param nu Matern smoothness parameter.
#' @param rho Length-scale; if \code{NULL} (default), set to the median
#'   pairwise distance on the normalized coordinates.
#' @param alpha GP marginal standard deviation.
#' @param L_factor Domain boundary expansion factor (the normalized domain
#'   used for the eigenbasis is \eqn{[-L, L]} with \eqn{L = }\code{L_factor}).
#'
#' @return A list with \code{phi} (n x m^2 basis matrix), \code{S_diag}
#'   (length-m^2 spectral densities, the prior variance for each basis
#'   weight), \code{omega_sq} (length-m^2 squared eigenfrequencies),
#'   \code{m}, \code{m2}, \code{rho}, and \code{nu}.
#'
#' @export
build_hsgp_basis <- function(coords, m = 5, nu = 1.5,
                             rho = NULL, alpha = 1, L_factor = 1.5) {
  n  <- nrow(coords)
  m2 <- m^2

  # Normalize coordinates to [-1, 1]
  coords_norm <- apply(coords, 2, function(x) {
    r <- range(x)
    2 * (x - r[1]) / (r[2] - r[1]) - 1
  })

  L <- L_factor  # boundary for normalized domain [-L, L]

  # Default rho: median pairwise distance on normalized coords
  if (is.null(rho)) {
    D   <- as.matrix(dist(coords_norm))
    rho <- median(D[upper.tri(D)])
    cat("  Auto length-scale rho =", round(rho, 4), "\n")
  }

  # 1D Laplacian eigenfunctions and eigenvalues
  # phi_j(x) = (1/sqrt(L)) * sin(j*pi*(x + L) / (2L))
  # lambda_j = (j*pi / (2L))^2
  phi_1d <- function(x, m) {
    sapply(1:m, function(j) sin(j * pi * (x + L) / (2 * L)) / sqrt(L))
  }
  lambda_1d <- function(m) sapply(1:m, function(j) (j * pi / (2 * L))^2)

  phi_x <- phi_1d(coords_norm[, 1], m)   # n x m
  phi_y <- phi_1d(coords_norm[, 2], m)   # n x m
  lam_x <- lambda_1d(m)                  # m eigenvalues (x)
  # The eigenvalue formula lambda_j = (j*pi/(2L))^2 depends only on the mode
  # index j and the domain half-width L, not on which coordinate axis it is
  # applied to -- so the y-axis eigenvalues are identical to the x-axis ones.
  # (This line was missing from the original script, which referenced an
  # undefined `lam_y` below; see the repository README for details.)
  lam_y <- lam_x

  # 2D tensor product: phi_{ij}(l) = phi_i(x) * phi_j(y)
  # omega^2_{ij} = lambda_i + lambda_j  (2D Laplacian eigenvalue)
  phi    <- matrix(0, n, m2)
  S_diag <- rep(0, m2)
  omega_sq_vec <- rep(0, m2)

  idx <- 1
  for (i in 1:m) {
    for (j in 1:m) {
      phi[, idx]        <- phi_x[, i] * phi_y[, j]
      omega_sq_vec[idx] <- lam_x[i] + lam_y[j]
      S_diag[idx]       <- matern_spectral_density(omega_sq_vec[idx],
                                                   nu = nu, rho = rho,
                                                   alpha = alpha)
      idx <- idx + 1
    }
  }

  list(phi = phi, S_diag = S_diag, omega_sq = omega_sq_vec,
       m = m, m2 = m2, rho = rho, nu = nu)
}
