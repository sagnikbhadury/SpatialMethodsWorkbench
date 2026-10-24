## Sensitivity and model checks for the count-based CRC analysis (2026-10-03).
##   1. global scale tau by node, image and pathology group
##   2. residual checks of the primary fit (size and spatial correlation)
##   3. spatial detrending of every abundance field before fitting
##      (crc_analysis.R with CRC_DETREND=1), for the missing spatial intercept
##   4. rare-cell threshold of 10 instead of 5 cells (CRC_MIN_CELLS=10)
##
## Usage (R 4.4.1): Rscript crc_sensitivity.R

suppressPackageStartupMessages({ library(lme4); library(lmerTest) })
BASE <- "crc_results_counts"
ALT  <- c(detrend = "crc_results_counts_detrend", min10 = "crc_results_counts_min10")
OUT  <- file.path(BASE, "summary")

dat <- read.csv("CRC_data_cleaned_corrected.csv", stringsAsFactors = FALSE)
g0  <- readRDS(file.path(BASE, "graph_list_final.rds"))
set <- readRDS(file.path(BASE, "run_settings.rds")); NG <- set$n_grid
CT  <- g0[[1]]$cell_types; p <- length(CT)
dat <- dat[dat$type %in% CT, ]
meta <- unique(dat[, c("img", "pathology", "patient")])
meta$group <- factor(ifelse(meta$pathology == 1, "CLR", "DII"), c("CLR", "DII"))
meta <- meta[match(names(g0), meta$img), ]

## ---- 1. tau ------------------------------------------------------------------
tau <- do.call(rbind, lapply(names(g0), function(img) data.frame(
  img = img, response = CT,
  tau = sapply(g0[[img]]$node_results, function(r) r$tau_mean),
  sigma = sapply(g0[[img]]$node_results, function(r) r$sigma_mean))))
tau <- merge(tau, meta[, c("img", "patient", "group")], by = "img")
tau <- tau[!is.na(tau$tau), ]
write.csv(tau, file.path(OUT, "tau_by_image_node.csv"), row.names = FALSE)
cat("tau, posterior mean per regression (quartiles) by group:\n")
print(round(do.call(rbind, tapply(tau$tau, tau$group, quantile, c(.25, .5, .75))), 3))
ft <- lmer(log(tau) ~ group + (1 | patient) + (1 | response), tau)
cat("log tau, DII - CLR:\n"); print(round(coef(summary(ft))["groupDII", ], 4))
pt <- sapply(split(tau, tau$response), function(d) {
  f <- tryCatch(suppressMessages(suppressWarnings(lmer(log(tau) ~ group + (1 | patient), d))),
                error = function(e) NULL)
  if (is.null(f)) NA else coef(summary(f))["groupDII", "Pr(>|t|)"]
})
cat("per response type: BH < 0.05 in", sum(p.adjust(pt, "BH") < 0.05, na.rm = TRUE), "of",
    sum(!is.na(pt)), "| smallest BH p:", signif(min(p.adjust(pt, "BH"), na.rm = TRUE), 3), "\n")

## ---- 2. residuals of the primary fit -------------------------------------------
count_matrix <- function(cells, min_cells = set$min_cells) {
  xr <- range(cells$X); yr <- range(cells$Y)
  xb <- seq(xr[1], xr[2], length.out = NG + 1); yb <- seq(yr[1], yr[2], length.out = NG + 1)
  Y <- matrix(0, NG^2, p, dimnames = list(NULL, CT))
  for (ct in CT) {
    pts <- cells[cells$type == ct, ]
    if (nrow(pts) < min_cells) next
    Y[, ct] <- sqrt(as.vector(table(factor(findInterval(pts$X, xb, all.inside = TRUE), 1:NG),
                                    factor(findInterval(pts$Y, yb, all.inside = TRUE), 1:NG))))
  }
  Y <- scale(Y); Y[!is.finite(Y)] <- 0; Y
}
nb_cor <- function(r) {                      # correlation between adjacent quadrats
  R <- matrix(r, NG, NG)
  a <- c(R[-NG, ], R[, -NG]); b <- c(R[-1, ], R[, -1]); cor(a, b)
}
res <- do.call(rbind, lapply(names(g0), function(img) {
  Y <- count_matrix(dat[dat$img == img, ])
  do.call(rbind, lapply(g0[[img]]$node_results, function(nr) {
    if (is.na(nr$tau_mean)) return(NULL)
    y <- Y[, nr$s]; fit <- rowSums(Y[, nr$nbr_idx, drop = FALSE] * nr$edge_maps)
    data.frame(img = img, response = CT[nr$s], r2 = 1 - var(y - fit) / var(y),
               resid_nb_cor = nb_cor(y - fit), resp_nb_cor = nb_cor(y))
  }))
}))
write.csv(res, file.path(OUT, "residual_checks.csv"), row.names = FALSE)
cat("\nR^2 of the posterior mean fit (quartiles):", round(quantile(res$r2, c(.25, .5, .75)), 3), "\n")
cat("adjacent-quadrat correlation, response (quartiles):", round(quantile(res$resp_nb_cor, c(.25, .5, .75)), 3),
    "\nadjacent-quadrat correlation, residual (quartiles):", round(quantile(res$resid_nb_cor, c(.25, .5, .75)), 3),
    "| share above 0.2:", round(mean(res$resid_nb_cor > 0.2), 3), "\n")

## ---- 3-4. alternative fits ---------------------------------------------------------
est_of <- function(g) t(sapply(g, function(x) !is.na(sapply(x$node_results, function(r) r$tau_mean))))
up <- upper.tri(matrix(0, p, p))
for (nm in names(ALT)) {
  f <- file.path(ALT[[nm]], "graph_list_final.rds")
  if (!file.exists(f)) { cat("\n", nm, ": fit not found\n"); next }
  g1 <- readRDS(f); stopifnot(identical(names(g1), names(g0)))
  e0 <- est_of(g0); e1 <- est_of(g1)
  both <- a0 <- a1 <- agree <- 0; jac <- numeric(0)
  for (k in seq_along(g0)) {
    ok <- outer(e0[k, ] & e1[k, ], e0[k, ] & e1[k, ]) & up
    x0 <- g0[[k]]$adj[ok] == 1; x1 <- g1[[k]]$adj[ok] == 1
    both <- both + sum(ok); a0 <- a0 + sum(x0); a1 <- a1 + sum(x1); agree <- agree + sum(x0 == x1)
    if (sum(x0 | x1) > 0) jac <- c(jac, sum(x0 & x1) / sum(x0 | x1))
  }
  prev <- function(g, e) {
    num <- den <- matrix(0, p, p)
    for (k in seq_along(g)) { ok <- outer(e[k, ], e[k, ]); num <- num + g[[k]]$adj * ok; den <- den + ok }
    (num / den)[up & den >= 50]
  }
  den0 <- matrix(0, p, p); for (k in seq_along(g0)) den0 <- den0 + outer(e0[k, ] & e1[k, ], e0[k, ] & e1[k, ])
  pr0 <- { num <- matrix(0, p, p); for (k in seq_along(g0)) num <- num + g0[[k]]$adj * outer(e0[k, ] & e1[k, ], e0[k, ] & e1[k, ]); (num / den0)[up & den0 >= 50] }
  pr1 <- { num <- matrix(0, p, p); for (k in seq_along(g1)) num <- num + g1[[k]]$adj * outer(e0[k, ] & e1[k, ], e0[k, ] & e1[k, ]); (num / den0)[up & den0 >= 50] }
  cat("\n==", nm, "==\n estimable regressions:", sum(e1), "(primary:", sum(e0), ")\n",
      "pairs estimable in both:", both, "| selected primary:", round(a0 / both, 3),
      "| selected alternative:", round(a1 / both, 3), "| same decision:", round(agree / both, 3),
      "\n median Jaccard index per image:", round(median(jac), 3),
      "\n correlation of edge prevalence across pairs:", round(cor(pr0, pr1), 3),
      "| Spearman:", round(cor(pr0, pr1, method = "spearman"), 3), "\n")
}
