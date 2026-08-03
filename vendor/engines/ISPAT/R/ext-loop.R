#' Fit the spatial GLM for one marker within one cluster
#'
#' Internal worker called in parallel, once per marker, by \code{\link{ISPAT}}.
#' Fits a spatial Gaussian mixed-effects model (via Stan) for a single marker's
#' expression across the spatial locations of one cluster, then returns the
#' de-spatialised (spatial-signal-removed) residual for that marker.
#'
#' @param count Marker/gene row index into \code{Y}.
#' @param c Cluster index.
#' @param N_c Vector of cluster sizes.
#' @param pos Indices (into \code{S}) of the locations belonging to cluster \code{c}.
#' @param Y Marker/gene expression matrix.
#' @param S Spatial coordinate + cluster-annotation matrix.
#' @param Clusters Vector of unique cluster labels.
#' @param KS_c List of cluster-specific spatial kernel matrices.
#' @param lS Length-scale parameter for the spatial kernel.
#' @param stan_model_expo Compiled Stan model for the RBF/exponential kernel.
#' @param stan_model_matern Compiled Stan model for the Matern kernel.
#' @param sGLM_method \code{"MLE"} (sampling) or \code{"VB"} (variational inference).
#' @param Kernel \code{"RBF"} or \code{"Matern"}.
#'
#' @return A one-column data frame of de-spatialised residuals for this marker.
#' @keywords internal
#' @noRd
ext_loop <- function(count, c, N_c, pos, Y, S, Clusters, KS_c, lS,
                      stan_model_expo, stan_model_matern, sGLM_method, Kernel) {

  mydata <- Y
  pos <- which(S[, 3] == Clusters[c])
  S_loc_c <- S[pos, -3]
  stan_Data <- list(N = N_c[c], y = as.vector(as.matrix(mydata)[count, pos]),
                     Sx = S_loc_c[, 1], Sy = S_loc_c[, 2], lS = lS)

  beta_est <- NA_real_
  sigmaS_est <- NA_real_
  sigmaEPS_est <- NA_real_

  if (Kernel == "RBF") {
    if (sGLM_method == "MLE") {
      MLE_RBF_G <- rstan::sampling(stan_model_expo, data = stan_Data, iter = 500,
                                    chains = 2, warmup = 200, thin = 1, verbose = TRUE)
      MLEs <- summary(MLE_RBF_G)$summary[, "mean"]
      beta_est <- MLEs[1]
      sigmaS_est <- MLEs[2]
      sigmaEPS_est <- MLEs[3]
    }
    if (sGLM_method == "VB") {
      message("We recommend 'sampling' from the posterior for sGLM!")
      VB_RBF_G <- rstan::vb(stan_model_expo, data = stan_Data, seed = 50,
                             output_samples = 5000, refresh = 0)
      VBests <- summary(VB_RBF_G)$summary[, "mean"]
      beta_est <- VBests[1]
      sigmaS_est <- VBests[2]
      sigmaEPS_est <- VBests[3]
    }
  }

  if (Kernel == "Matern") {
    if (sGLM_method == "MLE") {
      MLE_Matern_G <- rstan::sampling(stan_model_matern, data = stan_Data, iter = 500,
                                       chains = 2, warmup = 200, thin = 1, verbose = TRUE)
      MLEs <- summary(MLE_Matern_G)$summary[, "mean"]
      beta_est <- MLEs[1]
      sigmaS_est <- MLEs[2]
      sigmaEPS_est <- MLEs[3]
    }
    if (sGLM_method == "VB") {
      message("We recommend 'sampling' from the posterior for sGLM!")
      VB_Matern_G <- rstan::vb(stan_model_matern, data = stan_Data, seed = 50,
                                output_samples = 5000, refresh = 0)
      VBests <- summary(VB_Matern_G)$summary[, "mean"]
      beta_est <- VBests[1]
      sigmaS_est <- VBests[2]
      sigmaEPS_est <- VBests[3]
    }
  }

  if (sigmaS_est < 0.01 || sigmaEPS_est < 0.01) {
    V <- (sigmaS_est + 0.001) * KS_c[[c]] + (sigmaEPS_est + 0.001) * diag(N_c[c])
  } else {
    V <- sigmaS_est * KS_c[[c]] + sigmaEPS_est * diag(N_c[c])
  }

  if (sigmaEPS_est < 0.001) {
    storeZ <- ((sigmaEPS_est + 0.001) * solve(V)) %*% (as.matrix(Y)[count, pos] - beta_est)
  } else {
    storeZ <- (sigmaEPS_est * solve(V)) %*% (as.matrix(Y)[count, pos] - beta_est)
  }

  data.frame(storeZ)
}
