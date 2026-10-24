## ================================================================
## Spatial transcriptomics application, step 3: summaries and figures
##
## Input: st_prepared.rds (st_prepare.R) and fit_m<M>/st_fit.rds (st_fit.R).
## Output (in fit_m<M>/summary): zone table, edge list, gene degrees,
## overlap between zones, convergence diagnostics, summary.log; and the
## two figures of the application in manuscript/Image_s:
##   st_zone_networks.pdf   one gene network per zone (six panels)
##   st_coefficient_maps.png  two estimated coefficient functions in the largest zone
##
## Usage (R 4.4.1): Rscript st_network_summary.R
## ================================================================

suppressPackageStartupMessages({
  library(igraph); library(ggplot2); library(gridExtra); library(parallel); library(posterior)
})
source("gp_group_horseshoe_spatial_graph.R")

OUT_DIR <- Sys.getenv("ST_OUT_DIR", "st_results")
M_BASIS <- as.integer(Sys.getenv("ST_M_BASIS", "2"))
FIT_DIR <- file.path(OUT_DIR, sprintf("fit_m%d", M_BASIS))
SUM_DIR <- file.path(FIT_DIR, "summary")
FIG_DIR <- file.path("manuscript", "Image_s")
dir.create(SUM_DIR, showWarnings = FALSE)
sink(file.path(SUM_DIR, "summary.log"), split = TRUE)

prep <- readRDS(file.path(OUT_DIR, "st_prepared.rds"))
fit  <- readRDS(file.path(FIT_DIR, "st_fit.rds"))
fits <- fit$fits; genes <- fit$genes; p <- length(genes); K <- length(fits); set <- fit$settings
up <- upper.tri(matrix(0, p, p))
zone_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7")[seq_len(K)]
top_zone <- prep$gene_tab$top_zone[match(genes, prep$gene_tab$gene)]

## ---- 1. edges of every zone ---------------------------------------------------
edge_rows <- list(); zone_rows <- list()
for (f in fits) {
  ok <- outer(f$est, f$est) & up
  R <- suppressWarnings(cor(f$Y)); R[!is.finite(R)] <- NA
  ij <- which(ok, arr.ind = TRUE)
  bi <- f$beta_mean[ij]; bj <- t(f$beta_mean)[ij]
  ed <- data.frame(zone = f$zone, gene_i = genes[ij[, 1]], gene_j = genes[ij[, 2]],
                   edge_on = f$adj[ij], kappa_max = f$kappa_max[ij],
                   mean_ij = bi, mean_ji = bj,
                   rms_min = pmin(f$beta_rms[ij], t(f$beta_rms)[ij]),
                   share_pos_ij = f$share_pos[ij], share_pos_ji = t(f$share_pos)[ij],
                   marg_cor = R[ij])
  ed$sign <- ifelse(bi > 0 & bj > 0, "positive", ifelse(bi < 0 & bj < 0, "negative", "mixed"))
  ## a coefficient function varies in sign when 10% to 90% of the spots are positive
  ed$sign_varies <- (ed$share_pos_ij > 0.1 & ed$share_pos_ij < 0.9) |
                    (ed$share_pos_ji > 0.1 & ed$share_pos_ji < 0.9)
  edge_rows[[f$zone]] <- ed
  on <- ed[ed$edge_on == 1, ]
  zone_rows[[f$zone]] <- data.frame(
    zone = f$zone, name = f$name, n_spots = length(f$idx), n_genes = sum(f$est),
    n_coef = (p - 1) * ncol(f$phi), n_pairs = nrow(ed), n_edges = nrow(on),
    density = nrow(on) / nrow(ed),
    positive = mean(on$sign == "positive"), negative = mean(on$sign == "negative"),
    sign_varies = mean(on$sign_varies),
    marg_015 = mean(abs(ed$marg_cor) > 0.15), sel_given_marg = mean(ed$edge_on[abs(ed$marg_cor) > 0.15]),
    sel_low_cor = mean(abs(on$marg_cor) <= 0.15),
    median_abs_cor_edges = median(abs(on$marg_cor)),
    sigma_median = median(f$sigma, na.rm = TRUE), r2_median = median(f$r2, na.rm = TRUE),
    tau_median = median(f$tau, na.rm = TRUE))
}
edges <- do.call(rbind, edge_rows); zt <- do.call(rbind, zone_rows)
write.csv(edges, file.path(SUM_DIR, "pair_by_zone.csv"), row.names = FALSE)
write.csv(zt, file.path(SUM_DIR, "zone_summary.csv"), row.names = FALSE)
cat(sprintf("slide %s | m = %d | %d genes | %d zones\n", prep$slide, set$m_basis, p, K))
print(format(zt, digits = 3), row.names = FALSE)
cat("pooled: selected", sum(zt$n_edges), "of", sum(zt$n_pairs), "pairs (", round(sum(zt$n_edges) / sum(zt$n_pairs), 4),
    ") | marginal |r| > 0.15:", round(mean(abs(edges$marg_cor) > 0.15), 3),
    "| of those, selected:", round(mean(edges$edge_on[abs(edges$marg_cor) > 0.15]), 3), "\n")
on <- edges[edges$edge_on == 1, ]
cat("selected edges: positive", round(mean(on$sign == "positive"), 3), "| negative", round(mean(on$sign == "negative"), 3),
    "| a coefficient function changes sign over the zone:", round(mean(on$sign_varies), 3), "\n")
cat("residual sd over all regressions (quartiles):",
    round(quantile(unlist(lapply(fits, function(f) f$sigma)), c(.25, .5, .75), na.rm = TRUE), 3), "\n")

## ---- 2. recurrence of edges across zones -----------------------------------------
edges$pair <- paste(edges$gene_i, "--", edges$gene_j)
rec <- aggregate(cbind(n_zones = edge_on, n_est = 1) ~ pair + gene_i + gene_j, edges, sum)
rec$zones <- sapply(rec$pair, function(pp) paste(edges$zone[edges$pair == pp & edges$edge_on == 1], collapse = ""))
rec$signs <- sapply(rec$pair, function(pp) paste(substr(edges$sign[edges$pair == pp & edges$edge_on == 1], 1, 1), collapse = ""))
rec <- rec[order(-rec$n_zones, rec$pair), ]
write.csv(rec[rec$n_zones > 0, ], file.path(SUM_DIR, "edge_recurrence.csv"), row.names = FALSE)
cat("\npairs selected in k zones (k = 0..", K, "):", paste(table(factor(rec$n_zones, 0:K)), collapse = " "), "\n")
cat("pairs selected in at least one zone:", sum(rec$n_zones > 0), "| in exactly one:", sum(rec$n_zones == 1), "\n")
cat("most recurrent pairs:\n"); print(head(rec[, c("pair", "n_zones", "zones", "signs")], 25), row.names = FALSE)
jac <- matrix(NA, K, K, dimnames = list(paste0("z", 1:K), paste0("z", 1:K)))
for (a in 1:K) for (b in 1:K) {
  ok <- outer(fits[[a]]$est & fits[[b]]$est, fits[[a]]$est & fits[[b]]$est) & up
  x <- fits[[a]]$adj[ok] == 1; y <- fits[[b]]$adj[ok] == 1
  jac[a, b] <- sum(x & y) / max(sum(x | y), 1)
}
write.csv(round(jac, 3), file.path(SUM_DIR, "zone_jaccard.csv"))
cat("\nJaccard index of the edge sets between zones:\n"); print(round(jac, 2))

## ---- 3. degrees and hubs -----------------------------------------------------------
deg <- sapply(fits, function(f) rowSums(f$adj)); colnames(deg) <- paste0("zone", 1:K)
deg_df <- data.frame(gene = genes, top_zone = top_zone, deg, total = rowSums(deg))
write.csv(deg_df[order(-deg_df$total), ], file.path(SUM_DIR, "gene_degrees.csv"), row.names = FALSE)
cat("\nhubs (five largest degrees) per zone:\n")
for (k in 1:K) { o <- order(-deg[, k])[1:5]
  cat(sprintf("  zone %d %-36s %s\n", k, fits[[k]]$name, paste0(genes[o], " (", deg[o, k], ")", collapse = ", "))) }
cat("genes with no edge in any zone:", sum(rowSums(deg) == 0), "\n")

## ---- 4. network figure: one panel per zone ----------------------------------------
union_w <- Reduce(`+`, lapply(fits, function(f) f$adj))
gu <- graph_from_adjacency_matrix(union_w, mode = "undirected", weighted = TRUE, diag = FALSE)
set.seed(7)
L <- layout_with_fr(gu, weights = igraph::E(gu)$weight, niter = 2000)   # posterior also exports E
## the force-directed layout crowds the hubs at the centre; keep each gene's angle and
## replace its radius by the rank of the radius, so that the genes fill a disc evenly
L <- sweep(L, 2, apply(L, 2, median))
rad <- sqrt(rank(sqrt(rowSums(L^2))) / p); ang <- atan2(L[, 2], L[, 1])
L <- cbind(rad * cos(ang), rad * sin(ang))
halo_text <- function(x, y, labels, cex) {
  for (a in seq(0, 2 * pi, length.out = 13)[-1])
    text(x + 0.012 * cos(a), y + 0.012 * sin(a), labels, cex = cex, font = 2, col = "white")
  text(x, y, labels, cex = cex, font = 2, col = "black")
}
sign_cols <- c(positive = "#e34948", negative = "#2a78d6", mixed = "grey55")
cairo_pdf(file.path(FIG_DIR, "st_zone_networks.pdf"), width = 11, height = 9)
layout(matrix(c(1:6, 7, 7, 7), 3, 3, byrow = TRUE), heights = c(1, 1, 0.3))
par(mar = c(0.3, 0.3, 2.2, 0.3))
for (f in fits) {
  g <- graph_from_adjacency_matrix(f$adj, mode = "undirected", diag = FALSE)
  el <- as_edgelist(g, names = FALSE)
  bi <- f$beta_mean[el]; bj <- f$beta_mean[el[, 2:1, drop = FALSE]]
  sg <- ifelse(bi > 0 & bj > 0, "positive", ifelse(bi < 0 & bj < 0, "negative", "mixed"))
  d <- degree(g); lab <- rep(NA, p)
  hub <- order(-d)[seq_len(10)]; hub <- hub[d[hub] > 0]
  plot(g, layout = L, rescale = FALSE, xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05),
       vertex.size = ifelse(d > 0, 2 + 1.1 * sqrt(d), 1.2) * 1.5,
       vertex.color = ifelse(!f$est, NA, ifelse(d > 0, zone_cols[top_zone], "grey88")),
       vertex.frame.color = ifelse(!f$est, NA, ifelse(d > 0, "grey25", "grey75")), vertex.frame.width = 0.5,
       vertex.label = NA,
       edge.color = adjustcolor(sign_cols[sg], 0.7), edge.width = 0.8,
       edge.lty = c(positive = 1, negative = 5, mixed = 3)[sg])
  halo_text(L[hub, 1], L[hub, 2] + 0.065, genes[hub], cex = 0.78)
  title(main = sprintf("%s", f$name), cex.main = 1.15, line = 1.1)
  mtext(sprintf("%d spots, %d edges", length(f$idx), ecount(g)), side = 3, line = 0.05, cex = 0.72)
  box(col = "grey70")
}
par(mar = c(0, 0, 0, 0)); plot.new()
legend(0.5, 1.05, xjust = 0.5, horiz = TRUE, bty = "n", cex = 1, xpd = NA,
       legend = c("Positive association", "Negative association", "Mixed"),
       col = sign_cols, lty = c(1, 5, 3), lwd = 1.6)
legend(0.5, 0.78, xjust = 0.5, ncol = 3, bty = "n", cex = 1, pt.cex = 1.7, pch = 21, col = "grey25", xpd = NA,
       pt.bg = zone_cols, legend = sapply(fits, function(f) f$name),
       title = "Node color: zone in which the gene has its largest mean expression")
dev.off()

## ---- 5. zones on the tissue and two coefficient functions ---------------------------
obs <- prep$obs
obs$zone_name <- factor(prep$zone_tab$name[obs$zone], prep$zone_tab$name)
base <- theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 8),
        legend.key.size = unit(0.35, "cm"))
p_zone <- ggplot(obs, aes(x, -y, color = zone_name)) + geom_point(size = 0.75) +
  scale_color_manual(values = setNames(zone_cols, prep$zone_tab$name), name = NULL) + coord_equal() +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
  labs(x = NULL, y = NULL, title = "Zones annotated by the pathologist",
       subtitle = sprintf("%s spots of one Visium slide", format(nrow(obs), big.mark = ","))) +
  base + theme(legend.position = "bottom")
## examples from the largest zone: the selected edge with the smallest kappa, and the
## selected edge whose coefficient function has the largest spatial range relative to its RMS
z1 <- fits[[1]]; e1 <- edges[edges$zone == 1 & edges$edge_on == 1, ]
field_of <- function(f, s, t) { nb <- (1:p)[-s]; m2 <- ncol(f$phi)
  as.vector(f$phi %*% matrix(f$theta[[s]], m2, p - 1)[, match(t, nb)]) }
rng <- mapply(function(a, b) { v <- field_of(z1, match(a, genes), match(b, genes)); diff(range(v)) / sqrt(mean(v^2)) },
              e1$gene_i, e1$gene_j)
ex <- rbind(e1[which.min(e1$kappa_max), ], e1[which.max(rng), ])
map_plot <- function(row, tag) {
  s <- match(row$gene_i, genes); t <- match(row$gene_j, genes)
  d <- data.frame(z1$coords, beta = field_of(z1, s, t)); names(d)[1:2] <- c("x", "y")
  lim <- max(abs(d$beta))
  ggplot() + geom_point(data = obs[obs$zone != 1, ], aes(x, -y), color = "grey90", size = 0.6) +
    geom_point(data = d, aes(x, -y, color = beta), size = 0.85) +
    scale_color_gradient2(low = "#2a78d6", mid = "#f4f4f2", high = "#e34948", midpoint = 0,
                          limits = c(-lim, lim), name = "Coefficient") + coord_equal() +
    labs(x = NULL, y = NULL, title = sprintf("%s given %s", row$gene_i, row$gene_j),
         subtitle = sprintf("%s (%s)", z1$name, tag)) +
    base + theme(legend.position = "bottom", legend.key.width = unit(0.8, "cm"))
}
p_a <- map_plot(ex[1, ], "strongest edge"); p_b <- map_plot(ex[2, ], "most variable function")
ggsave(file.path(FIG_DIR, "st_coefficient_maps.png"), arrangeGrob(p_a, p_b, nrow = 1),
       width = 9, height = 4.8, dpi = 300)
cat("\ncoefficient maps drawn for:", ex$pair[1], "(kappa", round(ex$kappa_max[1], 3), ") and", ex$pair[2],
    "(range / RMS", round(max(rng), 2), ")\n")
saveRDS(list(layout = L, examples = ex), file.path(SUM_DIR, "figure_objects.rds"))

## ---- 6. convergence of representative fits ----------------------------------------
## The two example edges, both directions, four chains from dispersed starting values.
cl <- makeCluster(4); invisible(clusterEvalQ(cl, source("gp_group_horseshoe_spatial_graph.R")))
dg <- list()
for (k in 1:2) for (dirn in 1:2) {
  s <- match(if (dirn == 1) ex$gene_i[k] else ex$gene_j[k], genes)
  t <- match(if (dirn == 1) ex$gene_j[k] else ex$gene_i[k], genes)
  jj <- match(t, (1:p)[-s])
  ch <- parLapply(cl, 1:4, function(c, zf, s, jj, set, p) {
    set.seed(51000 + 10 * s + c); q <- p - 1; m2 <- ncol(zf$phi)
    init <- list(theta = rep(0, q * m2), sigma_sq = 1, tau_sq = c(1e-3, 1e-2, 1e-1, 1)[c],
                 lambda_sq = rep(c(0.1, 1, 10, 1)[c], q))
    S <- build_hsgp_basis(zf$coords, m = set$m_basis)$S_diag
    fit <- gp_group_horseshoe_sampler(zf$Y[, s], zf$Y[, -s], zf$phi, S, nmc = set$nmc, burn = set$burn,
                                      thin = set$thin, init = init)
    idx <- ((jj - 1) * m2 + 1):(jj * m2)
    b <- zf$phi %*% fit$theta_store[idx, , drop = FALSE]
    cbind(rms = sqrt(colMeans(b^2)), kappa = fit$kappa_store[jj, ],
          tau = sqrt(fit$tau_store), sigma = sqrt(fit$sigma_store))
  }, zf = z1, s = s, jj = jj, set = set, p = p)
  for (v in colnames(ch[[1]])) { m <- sapply(ch, function(z) z[, v])
    dg[[length(dg) + 1]] <- data.frame(response = genes[s], predictor = genes[t], quantity = v, mean = mean(m),
                                       rhat = rhat(m), ess_bulk = ess_bulk(m), ess_tail = ess_tail(m)) }
}
stopCluster(cl)
dg <- do.call(rbind, dg)
write.csv(dg, file.path(SUM_DIR, "representative_mcmc_diagnostics.csv"), row.names = FALSE)
cat("\nConvergence, representative regressions in zone 1 (4 chains, dispersed starts):\n")
print(format(dg, digits = 3), row.names = FALSE)
for (v in unique(dg$quantity)) cat(sprintf(" %-5s Rhat max %.3f | bulk ESS min %.0f\n", v,
                                           max(dg$rhat[dg$quantity == v]), min(dg$ess_bulk[dg$quantity == v])))
sink()
