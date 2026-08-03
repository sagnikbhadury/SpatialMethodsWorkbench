#' Informed Spatially Aware Patterns for multiplexed immunofluorescence data
#'
#' Fits a Bayesian spatial mixed-effects model to each marker at each spatial
#' location of a multiplexed immunofluorescence (mIF) tissue slide, de-spatialises
#' the resulting gene/marker expression, and then recovers shared and
#' region-specific cell-type interaction networks via a multi-study factor
#' analysis (MSFA) of the de-spatialised residuals.
#'
#' The spatial step regresses each marker independently on a Gaussian process
#' with a user-chosen covariance kernel (squared-exponential/"RBF" or Matern
#' 3/2), fit either by full Bayesian sampling ("MLE", via \code{rstan::sampling})
#' or variational inference ("VB", via \code{rstan::vb}). The spatially adjusted
#' residuals are then passed to \code{\link{cavi_msfa}} or \code{\link{svi_msfa}}
#' to estimate a shared interaction network and cluster/region-specific
#' interaction networks.
#'
#' @param Y A G x N matrix of marker/gene expression (G markers, N spatial
#'   locations/cells).
#' @param S An N x 3 matrix of spatial coordinates: the first two columns are
#'   the x/y spatial locations, the third column is an integer cluster or
#'   region annotation for each location.
#' @param ncores Number of cores to use for the parallel per-marker spatial fit.
#' @param RefPrior Reference prior covariance for the MSFA step (used when
#'   \code{use_ref = TRUE}); supply the squared actual covariance matrix.
#' @param use_ref Logical; whether to rescale the estimated networks using a
#'   graphical-lasso reference precision matrix.
#' @param Kernel Spatial covariance kernel: \code{"RBF"} (squared exponential)
#'   or \code{"Matern"} (Matern 3/2).
#' @param sGLM_method Estimation method for the per-marker spatial model:
#'   \code{"MLE"} (full sampling, recommended) or \code{"VB"} (variational
#'   inference, faster but less reliable).
#' @param VB_MSFA Logical; must be \code{TRUE}. Retained for interface
#'   compatibility with the original research code.
#' @param MSFA_method Multi-study factor analysis algorithm: \code{"CAVI"}
#'   (coordinate-ascent variational inference, recommended) or \code{"SVI"}
#'   (stochastic variational inference).
#'
#' @return A list with elements:
#' \describe{
#'   \item{Shared_Net}{G x G shared interaction network common to all clusters.}
#'   \item{Cluster_Net}{List of G x G cluster/region-specific interaction networks.}
#'   \item{glasso_refs}{Graphical lasso reference precision matrices per cluster.}
#'   \item{KS_c}{Cluster-specific spatial kernel matrices.}
#'   \item{VBfit_MSFA}{Raw output of the MSFA fit (\code{\link{cavi_msfa}} or
#'     \code{\link{svi_msfa}}).}
#'   \item{Z_est}{De-spatialised marker expression per cluster.}
#' }
#'
#' @examples
#' \dontrun{
#' mydata <- readRDS("Data - Analysis ready Marginally obtained Spatial intensities
#'                     of Cells on Tissue Slide_<diseasetype>_<patID>.rds")
#' rawdata <- mydata |> tidyr::pivot_wider(names_from = Cells, values_from = Intensities)
#' S <- as.matrix(rawdata[, c("X", "Y", "area_category")])
#' Y <- t(as.matrix(rawdata[, c("APC", "CTL", "Epi", "THelper", "Treg")]))
#'
#' fit <- ISPAT(Y, S, ncores = 4, RefPrior = diag(1, 5, 5), use_ref = FALSE,
#'              Kernel = "Matern", sGLM_method = "MLE",
#'              VB_MSFA = TRUE, MSFA_method = "CAVI")
#' }
#'
#' @export
ISPAT <- function(Y, S, ncores, RefPrior,
                   use_ref = c(TRUE, FALSE),
                   Kernel = c("RBF", "Matern"),
                   sGLM_method = c("MLE", "VB"),
                   VB_MSFA = c(TRUE, FALSE),
                   MSFA_method = c("CAVI", "SVI")) {

  Kernel <- match.arg(Kernel)
  sGLM_method <- match.arg(sGLM_method)
  MSFA_method <- match.arg(MSFA_method)

  ## Compile the two Gaussian spatial-GLM Stan models bundled with the package.
  ## Compilation happens here (at call time), not at package load, so that
  ## `library(ISPAT)` stays fast and does not require rstan at load time.
  stan_model_matern <- rstan::stan_model(
    file = system.file("stan", "sGLM_Matern.stan", package = "ISPAT")
  )
  stan_model_expo <- rstan::stan_model(
    file = system.file("stan", "sGLM_RBF.stan", package = "ISPAT")
  )

  ######   Global Parameters   ########
  N <- dim(Y)[2]  # Number of spatial locations
  G <- dim(Y)[1]  # Number of genes/markers
  C <- length(unique(S[, 3])) # Number of clusters
  Clusters <- unique(S[, 3])

  N_c <- numeric(C)
  Z_est <- vector("list", C)   ## latent gene expression matrix
  KS_c <- vector("list", C)    ## Cluster specific spatial kernels

  ### Cluster sizes ####
  for (c in seq_len(C)) {
    pos <- which(S[, 3] == Clusters[c])
    N_c[c] <- length(pos)
    Z_est[[c]] <- matrix(0, nrow = N_c[c], ncol = G)
    KS_c[[c]] <- matrix(0, N_c[c], N_c[c])
  }

  for (c in seq_len(C)) {
    pos <- which(S[, 3] == Clusters[c])

    ### length-scale heuristic ###
    a <- stats::dist(S[pos, -3])
    a_max <- log10(2 * max(a))
    a_min <- log10(min(a) / 2)
    a_seq <- seq(a_min, a_max, length.out = 10)
    lS <- 10^(a_seq[5])

    if (Kernel == "RBF") {
      S_loc_c <- S[pos, -3]
      for (i in seq_len(N_c[c])) {
        for (j in seq_len(i)) {
          dist_loc <- (S_loc_c[i, 1] - S_loc_c[j, 1])^2 + (S_loc_c[i, 2] - S_loc_c[j, 2])^2
          KS_c[[c]][i, j] <- exp(-dist_loc / (2 * lS^2))
          KS_c[[c]][j, i] <- KS_c[[c]][i, j]
        }
      }
    }

    if (Kernel == "Matern") {
      S_loc_c <- S[pos, -3]
      for (i in seq_len(N_c[c])) {
        for (j in seq_len(i)) {
          distS <- sqrt((S_loc_c[i, 1] - S_loc_c[j, 1])^2 + (S_loc_c[i, 2] - S_loc_c[j, 2])^2) / lS
          KS_c[[c]][i, j] <- (1 + distS * sqrt(3)) * exp(-distS * sqrt(3))
          KS_c[[c]][j, i] <- KS_c[[c]][i, j]
        }
      }
      if (!matrixcalc::is.positive.definite(KS_c[[c]])) {
        KS_c[[c]] <- as.matrix(Matrix::nearPD(KS_c[[c]])$mat)
      }
    }

    message("Fitting parallel spatial mixed effects model for cluster - ", c)
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    Z_est[[c]] <- suppressWarnings({
      foreach::foreach(
        my_count = seq_len(G),
        .packages = c("base", "Matrix", "rstan", "matrixcalc", "CVglasso"),
        .combine = cbind,
        .export = c("stan_model_matern", "stan_model_expo")
      ) %dopar% {
        ext_loop(my_count, c, N_c, pos, Y, S, Clusters, KS_c, lS,
                 stan_model_expo, stan_model_matern, sGLM_method, Kernel)
      }
    })
    parallel::stopCluster(cl)
    message("Finished with parallel spatial mixed effects model for cluster - ", c)
  }
  message("Finished with parallel spatial mixed effects model for all clusters")

  refs <- lapply(seq_len(C), function(c) {
    Glasso <- CVglasso::CVglasso(X = Z_est[[c]], nlam = 20, lam.min.ratio = 0.0000001,
                                  diagonal = FALSE, crit.cv = "BIC", maxit = 120,
                                  adjmaxit = NULL, K = 10, cores = 10)
    Glasso$Omega
  })

  ####   Multi-study Factor model on latent gene expression matrix   ####
  message("Fitting multi-study factor model")

  if (isTRUE(VB_MSFA) || identical(VB_MSFA, c(TRUE, FALSE))) {
    if (MSFA_method == "CAVI") {
      VBfit_MSFA <- cavi_msfa(Z_est, floor(2 * log(G)), rep(floor(2 * log(G)), C), scale = FALSE)
      Shared_Net <- tcrossprod(VBfit_MSFA$mean_phi)
      Cluster_Net <- lapply(seq_len(C), function(c) {
        tcrossprod(VBfit_MSFA$mean_lambda_s[[c]]) + tcrossprod(VBfit_MSFA$mean_phi)
      })
    } else {
      VBfit_MSFA <- svi_msfa(Z_est, floor(2 * log(G)), rep(floor(2 * log(G)), C), scale = FALSE)
      Shared_Net <- tcrossprod(VBfit_MSFA$mean_phi)
      Cluster_Net <- lapply(seq_len(C), function(c) {
        tcrossprod(VBfit_MSFA$mean_lambda_s[[c]])
      })
    }
  } else {
    stop("Please set VB_MSFA = TRUE and rerun")
  }

  if (isTRUE(use_ref) || identical(use_ref, TRUE)) {
    Cluster_Net <- lapply(seq_len(C), function(c) {
      (1 / sqrt(refs[[c]])) %*% Cluster_Net[[c]] %*% (1 / sqrt(refs[[c]]))
    })
  }

  message("Function finished - happy exploring!")
  list(Shared_Net = Shared_Net, Cluster_Net = Cluster_Net, glasso_refs = refs,
       KS_c = KS_c, VBfit_MSFA = VBfit_MSFA, Z_est = Z_est)
}
