## ================================================================
## Grid-size sensitivity of the count-based CODEX analysis (2026-10-22).
##
## Raised by J. Gaskins (20 Oct 2026): how was the 15 x 15 grid chosen
## and would a finer grid help? The primary fit (crc_results_counts,
## 15 x 15, m = 4) is compared with refits on 20 x 20 and 10 x 10 grids
## (crc_analysis.R with CRC_N_GRID=20 / 10, CRC_FIT_ONLY=1, m = 4) on
## the pairs that are estimable in both. Same comparison as
## crc_sensitivity.R sections 3-4, plus the number of edges per image
## and the mean count per quadrat.
##
## Usage (R 4.4.1): Rscript crc_grid_sensitivity.R
## Output: crc_results_counts/summary/grid_sensitivity.log, grid_sensitivity.csv
## ================================================================
BASE <- "crc_results_counts"
ALT  <- c(grid20 = "crc_results_counts_grid20", grid10 = "crc_results_counts_grid10")
OUT  <- file.path(BASE, "summary")
g0  <- readRDS(file.path(BASE, "graph_list_final.rds"))
CT  <- g0[[1]]$cell_types; p <- length(CT)
est_of <- function(g) t(sapply(g, function(x) !is.na(sapply(x$node_results, function(r) r$tau_mean))))
up <- upper.tri(matrix(0, p, p))
dat <- read.csv("CRC_data_cleaned_corrected.csv", stringsAsFactors = FALSE)
dat <- dat[dat$type %in% CT, ]
n_cells <- table(dat$img)

sink(file.path(OUT, "grid_sensitivity.log"))
rows <- list()
e0 <- est_of(g0)
edges0 <- sapply(g0, function(g) sum(g$adj[up]))
cat("Primary 15 x 15: median edges per image", median(edges0),
    "| selected share", round(mean(sapply(seq_along(g0), function(k) {
      ok <- outer(e0[k, ], e0[k, ]) & up; mean(g0[[k]]$adj[ok]) })), 3), "\n")
cat("Mean cells per quadrat (median over images): 15x15",
    round(median(n_cells[names(g0)] / 225), 1), "| 20x20", round(median(n_cells[names(g0)] / 400), 1),
    "| 10x10", round(median(n_cells[names(g0)] / 100), 1), "\n")
for (nm in names(ALT)) {
  f <- file.path(ALT[[nm]], "graph_list_final.rds")
  if (!file.exists(f)) { cat("\n", nm, ": fit not found\n"); next }
  g1 <- readRDS(f); stopifnot(identical(names(g1), names(g0)))
  e1 <- est_of(g1)
  both <- a0 <- a1 <- agree <- 0; jac <- numeric(0)
  for (k in seq_along(g0)) {
    ok <- outer(e0[k, ] & e1[k, ], e0[k, ] & e1[k, ]) & up
    x0 <- g0[[k]]$adj[ok] == 1; x1 <- g1[[k]]$adj[ok] == 1
    both <- both + sum(ok); a0 <- a0 + sum(x0); a1 <- a1 + sum(x1); agree <- agree + sum(x0 == x1)
    if (sum(x0 | x1) > 0) jac <- c(jac, sum(x0 & x1) / sum(x0 | x1))
  }
  den0 <- matrix(0, p, p); num0 <- num1 <- matrix(0, p, p)
  for (k in seq_along(g0)) {
    ok <- outer(e0[k, ] & e1[k, ], e0[k, ] & e1[k, ])
    den0 <- den0 + ok; num0 <- num0 + g0[[k]]$adj * ok; num1 <- num1 + g1[[k]]$adj * ok
  }
  keep <- up & den0 >= 50
  pr0 <- (num0 / den0)[keep]; pr1 <- (num1 / den0)[keep]
  edges1 <- sapply(g1, function(g) sum(g$adj[up]))
  sig1 <- median(unlist(lapply(g1, function(g) sapply(g$node_results, function(r) r$sigma_mean))), na.rm = TRUE)
  cat("\n==", nm, "==\n estimable regressions:", sum(e1), "(primary:", sum(e0), ")\n",
      "median edges per image:", median(edges1), "(primary:", median(edges0), ")\n",
      "median residual SD:", round(sig1, 3), "\n",
      "pairs estimable in both:", both, "| selected primary:", round(a0 / both, 3),
      "| selected alternative:", round(a1 / both, 3), "| same decision:", round(agree / both, 3),
      "\n median Jaccard index per image:", round(median(jac), 3),
      "\n correlation of edge prevalence across pairs:", round(cor(pr0, pr1), 3),
      "| Spearman:", round(cor(pr0, pr1, method = "spearman"), 3), "\n")
  ## the ten most prevalent pairs of the primary analysis under the alternative grid
  top <- order(-pr0)[1:10]
  lab <- outer(CT, CT, paste, sep = " -- ")[keep][top]
  cat(" ten most prevalent primary pairs, prevalence primary vs", nm, ":\n")
  print(data.frame(pair = lab, primary = round(pr0[top], 2), alt = round(pr1[top], 2)), row.names = FALSE)
  rows[[nm]] <- data.frame(grid = nm, estimable = sum(e1), median_edges = median(edges1), median_sigma = sig1,
                           pairs_both = both, sel_primary = a0 / both, sel_alt = a1 / both,
                           same_decision = agree / both, median_jaccard = median(jac),
                           prev_cor = cor(pr0, pr1), prev_spearman = cor(pr0, pr1, method = "spearman"))
}
sink()
if (length(rows)) write.csv(do.call(rbind, rows), file.path(OUT, "grid_sensitivity.csv"), row.names = FALSE)
cat(readLines(file.path(OUT, "grid_sensitivity.log")), sep = "\n")
