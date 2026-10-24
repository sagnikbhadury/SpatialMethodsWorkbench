#' Agreement between two estimated networks
#'
#' Compares the edge sets of two `"gpghs"` fits on the same variables (for
#' example two replicate tissue sections), as in the replication analysis of
#' the article: the number of common edges, the number expected if the two
#' networks were unrelated (the product of the two edge counts divided by
#' the number of pairs estimable in both), their ratio, the overlap (common
#' edges divided by the smaller edge count), the Spearman correlation of the
#' pair scores and the share of common edges with the same sign.
#'
#' @param fit1,fit2 Two `"gpghs"` objects with the same node names.
#' @return A one-row data frame.
#' @export
network_overlap <- function(fit1, fit2) {
  nodes <- intersect(rownames(fit1$adj), rownames(fit2$adj))
  if (length(nodes) < 2) stop("the fits share fewer than two variables")
  est <- fit1$estimable[nodes] & fit2$estimable[nodes]
  nodes <- nodes[est]
  A1 <- fit1$adj[nodes, nodes]; A2 <- fit2$adj[nodes, nodes]
  up <- upper.tri(A1)
  e1 <- sum(A1[up]); e2 <- sum(A2[up]); both <- sum(A1[up] == 1 & A2[up] == 1)
  n_pairs <- sum(up)
  expected <- e1 * e2 / n_pairs
  s1 <- fit1$score[nodes, nodes][up]; s2 <- fit2$score[nodes, nodes][up]
  sign_of <- function(fit, s, t) {
    m1 <- mean(fit$edge_maps[[s]][, t]); m2 <- mean(fit$edge_maps[[t]][, s])
    if (m1 > 0 && m2 > 0) 1 else if (m1 < 0 && m2 < 0) -1 else 0
  }
  common <- which(up & A1 == 1 & A2 == 1, arr.ind = TRUE)
  sign_agree <- if (nrow(common)) mean(apply(common, 1, function(ij)
    sign_of(fit1, nodes[ij[1]], nodes[ij[2]]) == sign_of(fit2, nodes[ij[1]], nodes[ij[2]]))) else NA
  data.frame(n_pairs = n_pairs, edges1 = e1, edges2 = e2, common = both,
             expected = expected, ratio = if (expected > 0) both / expected else NA,
             overlap = if (min(e1, e2) > 0) both / min(e1, e2) else NA,
             score_cor = suppressWarnings(stats::cor(s1, s2, method = "spearman", use = "complete.obs")),
             sign_agree = sign_agree)
}
