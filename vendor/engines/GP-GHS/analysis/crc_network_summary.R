## Application summaries for the count-based CRC analysis (2026-10-03).
##
## Input: the fit written by crc_analysis.R with CRC_FEATURE=count.
## Output (in <FIT_DIR>/summary): per-image descriptives, edge prevalence
## among the images in which a pair is estimable, complete-case group
## comparisons, convergence diagnostics for representative fits, and the
## two application figures.
##
## Usage (R 4.4.1): CRC_OUT_DIR=crc_results_counts Rscript crc_network_summary.R

suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(ggplot2); library(posterior)
  library(parallel)
})
source("gp_group_horseshoe_spatial_graph.R")

FIT_DIR <- Sys.getenv("CRC_OUT_DIR", "crc_results_counts")
OUT_DIR <- file.path(FIT_DIR, "summary")
FIG_DIR <- file.path("manuscript", "Image_s")
dir.create(OUT_DIR, showWarnings = FALSE)
set <- readRDS(file.path(FIT_DIR, "run_settings.rds"))
N_GRID <- set$n_grid; M_BASIS <- set$m_basis; MIN_CELLS <- set$min_cells
stopifnot(set$feature == "count")

dat <- read.csv("CRC_data_cleaned_corrected.csv", stringsAsFactors = FALSE)
g   <- readRDS(file.path(FIT_DIR, "graph_list_final.rds"))
CT  <- g[[1]]$cell_types; p <- length(CT)
dat <- dat[dat$type %in% CT, ]
meta <- unique(dat[, c("img", "pathology", "patient")])
meta$group <- factor(ifelse(meta$pathology == 1, "CLR", "DII"), c("CLR", "DII"))
meta <- meta[match(names(g), meta$img), ]

count_matrix <- function(cells) {            # as in crc_analysis.R
  xr <- range(cells$X); yr <- range(cells$Y)
  xb <- seq(xr[1], xr[2], length.out = N_GRID + 1)
  yb <- seq(yr[1], yr[2], length.out = N_GRID + 1)
  grid <- expand.grid(X = (xb[-1] + xb[-(N_GRID + 1)]) / 2,
                      Y = (yb[-1] + yb[-(N_GRID + 1)]) / 2)
  Y <- matrix(0, nrow(grid), p, dimnames = list(NULL, CT))
  for (ct in CT) {
    pts <- cells[cells$type == ct, ]
    if (nrow(pts) < MIN_CELLS) next
    ix <- findInterval(pts$X, xb, all.inside = TRUE)
    iy <- findInterval(pts$Y, yb, all.inside = TRUE)
    Y[, ct] <- sqrt(as.vector(table(factor(ix, 1:N_GRID), factor(iy, 1:N_GRID))))
  }
  Y <- scale(Y); Y[!is.finite(Y)] <- 0
  list(Y = Y, coords = as.matrix(grid))
}

## ---- 1. per-image descriptives -------------------------------------------
est_mat <- t(sapply(g, function(x) !is.na(sapply(x$node_results, function(r) r$tau_mean))))
colnames(est_mat) <- CT
img_df <- data.frame(
  meta[, c("img", "patient", "group")],
  n_cells = as.numeric(table(dat$img)[names(g)]),
  n_types = rowSums(est_mat),
  n_pairs = choose(rowSums(est_mat), 2),
  n_edges = sapply(g, function(x) sum(x$adj) / 2))
img_df$density <- img_df$n_edges / img_df$n_pairs
sig <- unlist(lapply(g, function(x) sapply(x$node_results, function(r) r$sigma_mean)))
write.csv(img_df, file.path(OUT_DIR, "image_descriptives.csv"), row.names = FALSE)
cat("fits:", sum(est_mat), "of", length(est_mat), "\n")
cat("types per image (quartiles):", quantile(img_df$n_types, c(.25, .5, .75)), "\n")
cat("edges per image (quartiles):", quantile(img_df$n_edges, c(.25, .5, .75)), "\n")
cat("graph density (quartiles):", round(quantile(img_df$density, c(.25, .5, .75)), 3),
    "| pooled:", round(sum(img_df$n_edges) / sum(img_df$n_pairs), 3), "\n")
cat("residual sd (quartiles):", round(quantile(sig, c(.25, .5, .75), na.rm = TRUE), 3), "\n")
cat("density by group (LMM):\n")
print(coef(summary(lmer(density ~ group + (1 | patient), img_df))))

## ---- 2. pair-by-image table -----------------------------------------------
pairs <- t(combn(p, 2))
rows <- vector("list", length(g))
for (ii in seq_along(g)) {
  x <- g[[ii]]; img <- names(g)[ii]
  Yc <- count_matrix(dat[dat$img == img, ])$Y
  R  <- suppressWarnings(cor(Yc))
  rows[[ii]] <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(e) {
    i <- pairs[e, 1]; j <- pairs[e, 2]
    ok <- est_mat[ii, i] && est_mat[ii, j]
    ni <- x$node_results[[i]]; nj <- x$node_results[[j]]
    bi <- ni$edge_maps[, match(j, ni$nbr_idx)]; bj <- nj$edge_maps[, match(i, nj$nbr_idx)]
    data.frame(img = img, edge = paste(CT[i], "--", CT[j]), ct_i = CT[i], ct_j = CT[j],
               estimable = ok, edge_on = if (ok) x$adj[i, j] else NA,
               kappa_max = if (ok) x$kappa[i, j] else NA,
               rms_min = if (ok) min(sqrt(mean(bi^2)), sqrt(mean(bj^2))) else NA,
               mean_i = if (ok) mean(bi) else NA, mean_j = if (ok) mean(bj) else NA,
               marg_cor = if (ok) R[i, j] else NA)
  }))
}
pi_df <- do.call(rbind, rows)
pi_df <- merge(pi_df, meta[, c("img", "patient", "group")], by = "img")
write.csv(pi_df, file.path(OUT_DIR, "pair_by_image.csv"), row.names = FALSE)

## ---- 3. edge prevalence among estimable images -----------------------------
pe <- pi_df[pi_df$estimable, ]
prev <- do.call(rbind, lapply(split(pe, pe$edge), function(d) {
  on <- d[d$edge_on == 1, ]
  data.frame(edge = d$edge[1], ct_i = d$ct_i[1], ct_j = d$ct_j[1],
             n_img = nrow(d), n_CLR = sum(d$group == "CLR"), n_DII = sum(d$group == "DII"),
             prev = mean(d$edge_on),
             prev_CLR = mean(d$edge_on[d$group == "CLR"]),
             prev_DII = mean(d$edge_on[d$group == "DII"]),
             ## sign of the association in the images where the edge is selected
             share_pos = if (nrow(on)) mean(on$mean_i > 0 & on$mean_j > 0) else NA,
             share_neg = if (nrow(on)) mean(on$mean_i < 0 & on$mean_j < 0) else NA,
             median_rms = if (nrow(on)) median(on$rms_min) else NA,
             marg_prev = mean(abs(d$marg_cor) > 0.15),
             median_marg_cor = median(d$marg_cor))
}))
prev <- prev[order(-prev$prev), ]
write.csv(prev, file.path(OUT_DIR, "edge_prevalence.csv"), row.names = FALSE)
cat("\nTop edges by prevalence (pairs estimable in >= 50 images):\n")
print(head(prev[prev$n_img >= 50, c("edge", "n_img", "prev", "prev_CLR", "prev_DII",
                                    "share_pos", "share_neg", "median_rms",
                                    "marg_prev", "median_marg_cor")], 15),
      digits = 2, row.names = FALSE)
cat("pairs with prevalence > 0.5:", sum(prev$prev > 0.5 & prev$n_img >= 50),
    "| < 0.1:", sum(prev$prev < 0.1 & prev$n_img >= 50),
    "| of", sum(prev$n_img >= 50), "\n")
cat("pooled: selected", round(mean(pe$edge_on), 3), "| marginal |r| > 0.15:",
    round(mean(abs(pe$marg_cor) > 0.15), 3), "\n")
cat("of marginally correlated pair-images, share selected:",
    round(mean(pe$edge_on[abs(pe$marg_cor) > 0.15]), 3),
    "| of selected, share with |r| <= 0.15:",
    round(mean(abs(pe$marg_cor[pe$edge_on == 1]) <= 0.15), 3), "\n")
on <- pe[pe$edge_on == 1, ]
cat("selected pair-images: both directions positive", round(mean(on$mean_i > 0 & on$mean_j > 0), 3),
    "| both negative", round(mean(on$mean_i < 0 & on$mean_j < 0), 3), "\n")

## ---- 4. complete-case group comparisons -------------------------------------
cmp <- do.call(rbind, lapply(split(pe, pe$edge), function(d) {
  out <- data.frame(edge = d$edge[1], n_CLR = sum(d$group == "CLR"),
                    n_DII = sum(d$group == "DII"), prev = mean(d$edge_on),
                    logit_or = NA, logit_p = NA  ,
                    rms_delta = NA, rms_p = NA)
  if (out$n_CLR < 10 || out$n_DII < 10) return(out)
  ## Edge indicator: one proportion per patient (share of that patient's
  ## estimable images with the edge), compared by a Wilcoxon rank-sum test.
  ## A mixed logistic model was unstable here (singular fits).
  if (out$prev >= 0.05 && out$prev <= 0.95) {
    pp <- aggregate(edge_on ~ patient + group, d, mean)
    if (all(table(pp$group) >= 5)) {
      out$logit_or <- diff(tapply(pp$edge_on, pp$group, mean))
      out$logit_p  <- suppressWarnings(wilcox.test(edge_on ~ group, pp)$p.value)
    }
  }
  f <- tryCatch(suppressMessages(suppressWarnings(
    lmer(log(rms_min + 1e-4) ~ group + (1 | patient), d, REML = FALSE))), error = function(e) NULL)
  if (!is.null(f)) { cs <- coef(summary(f)); out$rms_delta <- cs["groupDII", 1]; out$rms_p <- cs["groupDII", "Pr(>|t|)"] }
  out
}))
cmp$logit_fdr <- p.adjust(cmp$logit_p, "BH"); cmp$rms_fdr <- p.adjust(cmp$rms_p, "BH")
write.csv(cmp[order(cmp$logit_p), ], file.path(OUT_DIR, "group_comparisons.csv"), row.names = FALSE)
cat("\nGroup comparison, complete case:\n edge indicator: tested", sum(!is.na(cmp$logit_p)),
    "| BH < 0.05:", sum(cmp$logit_fdr < 0.05, na.rm = TRUE), "| raw p < 0.05:",
    sum(cmp$logit_p < 0.05, na.rm = TRUE), "| min BH:", signif(min(cmp$logit_fdr, na.rm = TRUE), 3),
    "\n RMS magnitude: tested", sum(!is.na(cmp$rms_p)), "| BH < 0.05:",
    sum(cmp$rms_fdr < 0.05, na.rm = TRUE), "| raw p < 0.05:", sum(cmp$rms_p < 0.05, na.rm = TRUE),
    "| min BH:", signif(min(cmp$rms_fdr, na.rm = TRUE), 3), "\n")
print(head(cmp[order(cmp$logit_p), ], 6), digits = 3, row.names = FALSE)

## Signed companion of the magnitude summary (suggested by J. Gaskins): the
## directional spatial mean that is smaller in absolute value, compared
## between groups on the untransformed scale.
pe$signed_min <- ifelse(abs(pe$mean_i) <= abs(pe$mean_j), pe$mean_i, pe$mean_j)
sg <- do.call(rbind, lapply(split(pe, pe$edge), function(d) {
  out <- data.frame(edge = d$edge[1], signed_delta = NA, signed_p = NA)
  if (sum(d$group == "CLR") < 10 || sum(d$group == "DII") < 10) return(out)
  f <- tryCatch(suppressMessages(suppressWarnings(
    lmer(signed_min ~ group + (1 | patient), d, REML = FALSE))), error = function(e) NULL)
  if (!is.null(f)) { cs <- coef(summary(f)); out$signed_delta <- cs["groupDII", 1]; out$signed_p <- cs["groupDII", "Pr(>|t|)"] }
  out
}))
sg$signed_fdr <- p.adjust(sg$signed_p, "BH")
write.csv(sg[order(sg$signed_p), ], file.path(OUT_DIR, "group_comparisons_signed.csv"), row.names = FALSE)
cat(" signed mean: tested", sum(!is.na(sg$signed_p)), "| BH < 0.05:",
    sum(sg$signed_fdr < 0.05, na.rm = TRUE), "| raw p < 0.05:", sum(sg$signed_p < 0.05, na.rm = TRUE),
    "| min BH:", signif(min(sg$signed_fdr, na.rm = TRUE), 3), "\n")
cat(" median |signed mean| / RMS:", round(median(abs(pe$signed_min) / pe$rms_min, na.rm = TRUE), 3),
    "| two directions agree in sign:", round(mean(sign(pe$mean_i) == sign(pe$mean_j)), 3), "\n")

## estimability itself differs by group (abundance)
estl <- data.frame(meta[, c("patient", "group")], est_mat, check.names = FALSE)
cat("\nImages in which each type is estimable (CLR of 68, DII of 72):\n")
print(rbind(CLR = colSums(est_mat[meta$group == "CLR", ]), DII = colSums(est_mat[meta$group == "DII", ])))

## ---- 5. prevalence heat map --------------------------------------------------
ord <- CT
hm <- rbind(data.frame(x = prev$ct_i, y = prev$ct_j, v = prev$prev_CLR, n = prev$n_CLR, group = "CLR"),
            data.frame(x = prev$ct_j, y = prev$ct_i, v = prev$prev_CLR, n = prev$n_CLR, group = "CLR"),
            data.frame(x = prev$ct_i, y = prev$ct_j, v = prev$prev_DII, n = prev$n_DII, group = "DII"),
            data.frame(x = prev$ct_j, y = prev$ct_i, v = prev$prev_DII, n = prev$n_DII, group = "DII"))
hm$v[hm$n < 10] <- NA
hm$x <- factor(hm$x, ord); hm$y <- factor(hm$y, rev(ord))
p_hm <- ggplot(hm, aes(x, y, fill = v)) +
  geom_tile(color = "white", linewidth = 0.4) + facet_wrap(~ group) +
  scale_fill_gradient(low = "#eef4fd", high = "#0d366b", limits = c(0, 1),
                      na.value = "grey88", name = "Edge\nprevalence") +
  labs(x = NULL, y = NULL) + coord_equal() + theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(FIG_DIR, "crc_edge_prevalence.png"), p_hm, width = 10, height = 5, dpi = 300)

## ---- 6. representative coefficient maps --------------------------------------
top <- head(prev[prev$n_img >= 50, ], 2)
map_rows <- list(); rep_sel <- list()
for (k in seq_len(nrow(top))) {
  d <- pe[pe$edge == top$edge[k] & pe$edge_on == 1, ]
  img <- d$img[which.min(abs(d$rms_min - median(d$rms_min)))]
  x <- g[[img]]; i <- match(top$ct_i[k], CT); j <- match(top$ct_j[k], CT)
  cm <- count_matrix(dat[dat$img == img, ])
  ni <- x$node_results[[i]]
  lab <- paste0(top$ct_i[k], " and ", top$ct_j[k], " (image ", img, ")")
  cells <- dat[dat$img == img & dat$type %in% c(top$ct_i[k], top$ct_j[k]), ]
  map_rows[[k]] <- list(
    field = data.frame(pair = lab, X = cm$coords[, 1], Y = cm$coords[, 2],
                       beta = ni$edge_maps[, match(j, ni$nbr_idx)]),
    cells = data.frame(pair = lab, X = cells$X, Y = cells$Y,
                       type = ifelse(cells$type == top$ct_i[k], "first type", "second type")))
  rep_sel[[k]] <- data.frame(edge = top$edge[k], img = img, i = i, j = j)
}
lim <- max(sapply(map_rows, function(m) max(abs(m$field$beta))))
pretty_ct <- function(x) sub(" macros", " macrophages", sub("memory CD4+ T", "memory CD4+ T cells", x, fixed = TRUE), fixed = TRUE)
plots <- lapply(seq_along(map_rows), function(k) {
  fld <- map_rows[[k]]$field; cel <- map_rows[[k]]$cells
  cel$type <- factor(cel$type, c("first type", "second type"),
                     pretty_ct(c(top$ct_i[k], top$ct_j[k])))
  ggplot(fld, aes(X, Y)) +
    geom_tile(aes(fill = beta), width = diff(sort(unique(fld$X)))[1],
              height = diff(sort(unique(fld$Y)))[1]) +
    geom_point(data = cel, aes(shape = type), size = 0.8, alpha = 0.8, color = "grey10") +
    scale_shape_manual(values = c(16, 2), name = NULL) +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f4f4f2", high = "#e34948", midpoint = 0,
                         limits = c(-lim, lim), name = "Coefficient") +
    coord_equal() +
    labs(x = "x (pixels)", y = "y (pixels)",
         title = paste0("Image ", rep_sel[[k]]$img, ": ", pretty_ct(top$ct_i[k]),
                        " given ", pretty_ct(top$ct_j[k]))) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          plot.title = element_text(size = 10, face = "bold"))
})
ggsave(file.path(FIG_DIR, "crc_coefficient_maps.png"),
       gridExtra::arrangeGrob(grobs = plots, nrow = 1), width = 10.5, height = 5.2, dpi = 300)

## ---- 7. convergence of the representative fits --------------------------------
rep_sel <- do.call(rbind, rep_sel)
diag_rows <- list(); trace_rows <- list()
cl <- makeCluster(4)
invisible(clusterEvalQ(cl, source("gp_group_horseshoe_spatial_graph.R")))
for (k in seq_len(nrow(rep_sel))) {
  cm <- count_matrix(dat[dat$img == rep_sel$img[k], ])
  hs <- build_hsgp_basis(cm$coords, m = M_BASIS)
  for (dirn in 1:2) {
    s <- if (dirn == 1) rep_sel$i[k] else rep_sel$j[k]
    t <- if (dirn == 1) rep_sel$j[k] else rep_sel$i[k]
    jj <- match(t, (1:p)[-s])
    ch <- parLapply(cl, 1:4, function(c, Y, hs, s, jj, set) {
      set.seed(83021 + 10 * s + c)
      fit <- gp_group_horseshoe_sampler(Y[, s], Y[, -s], hs$phi, hs$S_diag,
                                        nmc = set$nmc, burn = set$burn, thin = set$thin)
      idx <- ((jj - 1) * hs$m2 + 1):(jj * hs$m2)
      b <- hs$phi %*% fit$theta_store[idx, , drop = FALSE]
      cbind(rms = sqrt(colMeans(b^2)), kappa = fit$kappa_store[jj, ],
            tau = sqrt(fit$tau_store), sigma = sqrt(fit$sigma_store))
    }, Y = cm$Y, hs = hs, s = s, jj = jj, set = set)
    reg_lab <- paste0("Image ", rep_sel$img[k], "\n", pretty_ct(CT[s]), " on\n", pretty_ct(CT[t]))
    for (v in colnames(ch[[1]])) {
      m <- sapply(ch, function(z) z[, v])
      trace_rows[[length(trace_rows) + 1]] <- data.frame(
        regression = reg_lab, quantity = v, draw = rep(seq_len(nrow(m)), ncol(m)),
        chain = factor(rep(seq_len(ncol(m)), each = nrow(m))), value = as.vector(m))
      diag_rows[[length(diag_rows) + 1]] <- data.frame(
        edge = rep_sel$edge[k], img = rep_sel$img[k], response = CT[s], quantity = v,
        mean = mean(m), rhat = rhat(m), ess_bulk = ess_bulk(m), ess_tail = ess_tail(m))
    }
  }
}
stopCluster(cl)
dg <- do.call(rbind, diag_rows)
write.csv(dg, file.path(OUT_DIR, "representative_mcmc_diagnostics.csv"), row.names = FALSE)

## trace plots of the same chains
tr <- do.call(rbind, trace_rows)
saveRDS(tr, file.path(OUT_DIR, "representative_mcmc_draws.rds"))
tr$regression <- factor(tr$regression, unique(tr$regression))
tr$quantity <- factor(tr$quantity, c("rms", "kappa", "tau", "sigma"),
                      c("RMS of the coefficient function", "Shrinkage factor",
                        "Global scale", "Residual SD"))
p_tr <- ggplot(tr, aes(draw, value, color = chain, group = chain)) +
  geom_line(linewidth = 0.3, alpha = 0.85) +
  facet_grid(quantity ~ regression, scales = "free_y", switch = "y") +
  scale_color_manual(values = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"), name = "Chain") +
  labs(x = "Retained draw", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_blank(), strip.placement = "outside",
        strip.text.x = element_text(size = 8.5), strip.text.y = element_text(size = 8.5))
ggsave(file.path(FIG_DIR, "crc_traceplots.png"), p_tr, width = 10.5, height = 8, dpi = 300)
cat("\nConvergence, representative fits (4 chains):\n")
print(dg, digits = 3, row.names = FALSE)
cat("Rhat range:", round(range(dg$rhat), 3), "| bulk ESS range:", round(range(dg$ess_bulk)),
    "| tail ESS range:", round(range(dg$ess_tail)), "\n")
