#' Simulate from a spatially varying Gaussian graphical model
#'
#' Generates the data of the simulation study of the article: `n` locations
#' uniform on `[0, 10]^2`, a Barabasi-Albert graph on `p` nodes with density
#' `pi`, and at each location an independent draw from
#' `N_p(0, Omega(l)^{-1})` with `Omega(l) = d I + W(l)`, where `W_{ss'}(l)` is
#' zero unless `(s, s')` is an edge. In the `varying` scenario
#' `W_{ss'}(l) = tanh(1.5 g_{ss'}(l))` with independent unit-variance
#' Matern-3/2 Gaussian processes `g` of range `range`, so the coefficient
#' functions change sign over the domain; in the `constant` scenario
#' `W_{ss'}` is a constant with absolute value uniform on `(0.5, 1)` and a
#' random sign. The constant `d` is 1.25 times the largest negative
#' eigenvalue of `W(l)` over the locations, so `Omega(l)` is positive
#' definite everywhere and the nodewise model is correctly specified with
#' `beta_{ss'}(l) = -W_{ss'}(l) / d`.
#'
#' @param n Number of locations.
#' @param p Number of variables.
#' @param pi Graph density (share of the `p (p - 1) / 2` pairs that are
#'   edges); 0 gives an empty graph.
#' @param scenario `"varying"` or `"constant"`.
#' @param range Matern range of the true coefficient functions on
#'   `[0, 10]^2` (3 in the article).
#' @param seed Optional seed.
#' @return A list with `Y` (standardized `n x p` data), `coords`, `adj` (the
#'   true adjacency matrix), `beta_true` (a list over responses of
#'   `n x (p - 1)` true coefficient functions on the standardized scale) and
#'   `d`.
#' @examples
#' sim <- simulate_svggm(n = 300, p = 8, pi = 0.3, seed = 1)
#' sum(sim$adj) / 2
#' @export
simulate_svggm <- function(n = 600, p = 15, pi = 0.3, scenario = c("varying", "constant"),
                           range = 3, seed = NULL) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)
  coords <- matrix(stats::runif(n * 2, 0, 10), n, 2)
  adj <- scale_free_adjacency(p, pi)
  edges <- which(upper.tri(adj) & adj == 1, arr.ind = TRUE)
  E <- nrow(edges)
  W <- array(0, c(p, p, n))
  if (E > 0) {
    if (scenario == "varying") {
      D <- as.matrix(stats::dist(coords))
      K <- (1 + sqrt(3) * D / range) * exp(-sqrt(3) * D / range)
      Lk <- t(chol(K + diag(1e-6, n)))
      Fld <- tanh(1.5 * (Lk %*% matrix(stats::rnorm(n * E), n, E)))
    } else {
      cc <- sample(c(-1, 1), E, replace = TRUE) * stats::runif(E, 0.5, 1)
      Fld <- matrix(cc, n, E, byrow = TRUE)
    }
    for (e in seq_len(E)) {
      W[edges[e, 1], edges[e, 2], ] <- Fld[, e]
      W[edges[e, 2], edges[e, 1], ] <- Fld[, e]
    }
    lam_min <- vapply(seq_len(n), function(i)
      min(eigen(W[, , i], symmetric = TRUE, only.values = TRUE)$values), numeric(1))
    d <- 1.25 * max(-lam_min)
  } else d <- 1
  Y <- matrix(0, n, p)
  for (i in seq_len(n)) {
    R <- chol(diag(d, p) + W[, , i])
    Y[i, ] <- backsolve(R, stats::rnorm(p))
  }
  sds <- apply(Y, 2, stats::sd)
  beta_true <- lapply(seq_len(p), function(s) {
    nb <- (1:p)[-s]
    sapply(nb, function(t) -W[s, t, ] / d * sds[t] / sds[s])
  })
  Ys <- scale(Y); Ys <- unclass(Ys); attr(Ys, "scaled:center") <- NULL; attr(Ys, "scaled:scale") <- NULL
  colnames(Ys) <- paste0("V", seq_len(p)); dimnames(adj) <- list(colnames(Ys), colnames(Ys))
  list(Y = Ys, coords = coords, adj = adj, beta_true = beta_true, d = d)
}

#' Barabasi-Albert graph with a given density
#'
#' A preferential-attachment tree on `p` nodes to which edges are added (or
#' from which edges are removed) at random until it has
#' `round(pi * p (p - 1) / 2)` edges.
#' @param p Number of nodes.
#' @param pi Target density.
#' @return A `p x p` binary adjacency matrix.
#' @export
scale_free_adjacency <- function(p, pi) {
  A <- matrix(0, p, p)
  if (pi == 0) return(A)
  A[1, 2] <- A[2, 1] <- 1
  degree <- c(1, 1, rep(0, p - 2))
  for (v in 3:p) {
    prob <- degree[1:(v - 1)] / sum(degree[1:(v - 1)])
    tgt  <- sample(1:(v - 1), 1, prob = prob)
    A[v, tgt] <- A[tgt, v] <- 1
    degree[v] <- degree[v] + 1; degree[tgt] <- degree[tgt] + 1
  }
  tgt_e <- round(pi * p * (p - 1) / 2)
  cur_e <- sum(A[upper.tri(A)])
  if (tgt_e > cur_e) {
    cands <- which(upper.tri(A) & A == 0, arr.ind = TRUE)
    idx   <- sample(nrow(cands), min(tgt_e - cur_e, nrow(cands)))
    for (k in idx) A[cands[k, 1], cands[k, 2]] <- A[cands[k, 2], cands[k, 1]] <- 1
  } else if (tgt_e < cur_e) {
    ons <- which(upper.tri(A) & A == 1, arr.ind = TRUE)
    idx <- sample(nrow(ons), cur_e - tgt_e)
    for (k in idx) A[ons[k, 1], ons[k, 2]] <- A[ons[k, 2], ons[k, 1]] <- 0
  }
  diag(A) <- 0
  A
}

#' Graph recovery metrics
#'
#' True positive rate, false discovery proportion, F1, Matthews correlation
#' coefficient and the number of false edges of an estimated graph against
#' the truth, over the `p (p - 1) / 2` pairs.
#' @param est,truth Binary adjacency matrices.
#' @return A named numeric vector.
#' @export
graph_metrics <- function(est, truth) {
  tv <- truth[upper.tri(truth)]; ev <- est[upper.tri(est)]
  TP <- sum(ev == 1 & tv == 1); FP <- sum(ev == 1 & tv == 0)
  TN <- sum(ev == 0 & tv == 0); FN <- sum(ev == 0 & tv == 1)
  den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  c(TPR = TP / max(TP + FN, 1), FDP = FP / max(TP + FP, 1),
    F1 = ifelse(2 * TP + FP + FN > 0, 2 * TP / (2 * TP + FP + FN), 0),
    MCC = ifelse(den > 0, (TP * TN - FP * FN) / den, 0), TP = TP, FP = FP)
}
