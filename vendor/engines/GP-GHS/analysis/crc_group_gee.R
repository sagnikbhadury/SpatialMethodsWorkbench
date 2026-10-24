## ================================================================
## Group comparison of edge presence by GEE (2026-10-22).
##
## Suggested by J. Gaskins (meeting of 20 Oct 2026): the logistic mixed
## model for the edge indicator gave singular fits for many pairs, so the
## main analysis compared per-patient proportions by a Wilcoxon test. A
## logistic GEE with an exchangeable working correlation within patient
## is less sensitive to the random-effect structure. This script fits it
## for every pair tested in the main analysis and reports the same
## summaries (BH over pairs), for the magnitude and the signed summary as
## well (Gaussian GEE). Pairs in which the edge is never or always
## selected in one group are complete separations and are not tested
## with the logistic GEE (the coefficient diverges).
##
## Usage (R 4.5.2, needs geepack):  Rscript crc_group_gee.R
## Input : crc_results_counts/summary/pair_by_image.csv
## Output: crc_results_counts/summary/group_comparisons_gee.csv, gee.log
## ================================================================
suppressPackageStartupMessages(library(geepack))
OUT <- "crc_results_counts/summary"
pe <- read.csv(file.path(OUT, "pair_by_image.csv"), stringsAsFactors = FALSE)
pe <- pe[pe$estimable, ]
pe$group <- factor(pe$group, c("CLR", "DII"))
pe$patient <- factor(pe$patient)
pe$signed_min <- ifelse(abs(pe$mean_i) <= abs(pe$mean_j), pe$mean_i, pe$mean_j)

gee_p <- function(formula, d, family) {
  d <- d[order(d$patient), ]
  f <- tryCatch(suppressWarnings(geeglm(formula, id = patient, data = d, family = family,
                                        corstr = "exchangeable")), error = function(e) NULL)
  if (is.null(f)) return(c(NA, NA))
  cs <- coef(summary(f))
  c(cs["groupDII", "Estimate"], cs["groupDII", "Pr(>|W|)"])
}

res <- do.call(rbind, lapply(split(pe, pe$edge), function(d) {
  out <- data.frame(edge = d$edge[1], n_CLR = sum(d$group == "CLR"), n_DII = sum(d$group == "DII"),
                    pat_CLR = length(unique(d$patient[d$group == "CLR"])),
                    pat_DII = length(unique(d$patient[d$group == "DII"])),
                    prev = mean(d$edge_on), prev_CLR = mean(d$edge_on[d$group == "CLR"]),
                    prev_DII = mean(d$edge_on[d$group == "DII"]),
                    edge_est = NA, edge_p = NA, rms_est = NA, rms_p = NA, signed_est = NA, signed_p = NA)
  if (out$n_CLR < 10 || out$n_DII < 10) return(out)
  separated <- out$prev_CLR %in% c(0, 1) || out$prev_DII %in% c(0, 1)
  if (out$prev >= 0.05 && out$prev <= 0.95 && !separated) {
    v <- gee_p(edge_on ~ group, d, binomial()); out$edge_est <- v[1]; out$edge_p <- v[2]
  }
  d$lrms <- log(d$rms_min + 1e-4)
  v <- gee_p(lrms ~ group, d, gaussian()); out$rms_est <- v[1]; out$rms_p <- v[2]
  v <- gee_p(signed_min ~ group, d, gaussian()); out$signed_est <- v[1]; out$signed_p <- v[2]
  out
}))
res$edge_fdr <- p.adjust(res$edge_p, "BH")
res$rms_fdr <- p.adjust(res$rms_p, "BH")
res$signed_fdr <- p.adjust(res$signed_p, "BH")
res <- res[order(res$edge_p), ]
write.csv(res, file.path(OUT, "group_comparisons_gee.csv"), row.names = FALSE)

sink(file.path(OUT, "gee.log"))
cat("Logistic GEE (exchangeable within patient), edge presence, DII vs CLR:\n")
cat(" pairs tested:", sum(!is.na(res$edge_p)), "| separated (not tested):",
    sum(res$n_CLR >= 10 & res$n_DII >= 10 & res$prev >= 0.05 & res$prev <= 0.95 & is.na(res$edge_p)),
    "| BH < 0.05:", sum(res$edge_fdr < 0.05, na.rm = TRUE),
    "| raw p < 0.05:", sum(res$edge_p < 0.05, na.rm = TRUE), "| min BH:", signif(min(res$edge_fdr, na.rm = TRUE), 3), "\n")
cat("Gaussian GEE, log RMS magnitude:\n")
cat(" pairs tested:", sum(!is.na(res$rms_p)), "| BH < 0.05:", sum(res$rms_fdr < 0.05, na.rm = TRUE),
    "| raw p < 0.05:", sum(res$rms_p < 0.05, na.rm = TRUE), "| min BH:", signif(min(res$rms_fdr, na.rm = TRUE), 3), "\n")
cat("Gaussian GEE, signed summary:\n")
cat(" pairs tested:", sum(!is.na(res$signed_p)), "| BH < 0.05:", sum(res$signed_fdr < 0.05, na.rm = TRUE),
    "| raw p < 0.05:", sum(res$signed_p < 0.05, na.rm = TRUE), "| min BH:", signif(min(res$signed_fdr, na.rm = TRUE), 3), "\n")
cat("\nSmallest p-values for edge presence:\n")
print(head(res[, c("edge", "n_CLR", "n_DII", "pat_CLR", "pat_DII", "prev_CLR", "prev_DII", "edge_est", "edge_p", "edge_fdr")], 8), digits = 3, row.names = FALSE)
cat("\nSmallest p-values for the signed summary:\n")
r2 <- res[order(res$signed_p), ]
print(head(r2[, c("edge", "n_CLR", "n_DII", "pat_CLR", "pat_DII", "signed_est", "signed_p", "signed_fdr", "rms_p", "rms_fdr")], 6), digits = 3, row.names = FALSE)
cat("\nB cells -- CD8+ T cells:\n")
print(res[res$edge == "B cells -- CD8+ T cells", ], digits = 3, row.names = FALSE)
sink()
cat(readLines(file.path(OUT, "gee.log")), sep = "\n")
