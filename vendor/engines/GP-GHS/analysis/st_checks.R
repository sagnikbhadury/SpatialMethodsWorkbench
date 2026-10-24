## Negative control for the spatial transcriptomics application (2026-10-03).
## Every gene is permuted independently over the spots of its zone
## (st_fit.R with ST_PERMUTE=1 ST_TAG=_permuted). The permutation keeps the
## marginal distribution of each gene and removes all dependence between
## genes, so every pair selected in this run is a false edge.
##
## Usage (R 4.4.1): Rscript st_checks.R

OUT_DIR <- Sys.getenv("ST_OUT_DIR", "st_results")
f0 <- readRDS(file.path(OUT_DIR, "fit_m2", "st_fit.rds"))
fp <- readRDS(file.path(OUT_DIR, "fit_m2_permuted", "st_fit.rds"))
p <- length(f0$genes); up <- upper.tri(matrix(0, p, p))
sink(file.path(OUT_DIR, "fit_m2", "summary", "checks.log"), split = TRUE)
cat("== negative control: genes permuted independently within each zone ==\n")
tot <- 0; prs <- 0
for (f in fp$fits) {
  ok <- outer(f$est, f$est) & up; e <- sum(f$adj[ok]); tot <- tot + e; prs <- prs + sum(ok)
  cat(sprintf("  zone %d %-36s selected %2d of %d pairs | smallest kappa_max %.3f | median sigma %.2f\n",
              f$zone, f$name, e, sum(ok), min(f$kappa_max[ok]), median(f$sigma, na.rm = TRUE)))
}
cat(sprintf("  total: %d of %d pairs (%.5f) | primary analysis: %d\n", tot, prs, tot / prs,
            sum(sapply(f0$fits, function(f) sum(f$adj) / 2))))
sink()
