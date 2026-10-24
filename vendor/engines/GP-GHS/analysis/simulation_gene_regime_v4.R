## ================================================================
## Simulation in the regime of the gene application (2026-10-03).
##
## The main study (simulation_study_v4.R) has p <= 25 variables and m = 5.
## The spatial transcriptomics application has p = 100 genes, m = 2 and
## between 150 and 900 spots per zone, so that a regression has 396
## coefficients and often more coefficients than locations. This script
## simulates that regime from the same spatially varying Gaussian graphical
## model and adds one element: two independent datasets ("replicate
## sections") are drawn from the same graph and the same coefficient
## functions, and GP-GHS is fitted to each. The overlap of the two
## estimated graphs can then be read against their true sensitivity.
##
## Design: p = 100, density 0.05, varying scenario, n in {200, 400, 800},
## six datasets per n. Nothing is adjusted after the metrics are computed.
##
## Usage: Rscript simulation_gene_regime_v4.R
## ================================================================

Sys.setenv(SIM_SKIP_RUN = "1")
source("simulation_study_v4.R")
OUT_DIR <- "simulation_v4"
SEED_G  <- 20261301
P_GENE  <- 100; PI_GENE <- 0.05

## Fields and two replicate datasets from one spatially varying precision matrix
generate_pair <- function(adj, coords) {
  n <- nrow(coords); p <- ncol(adj)
  edges <- which(upper.tri(adj) & adj == 1, arr.ind = TRUE); E <- nrow(edges)
  D <- as.matrix(dist(coords))
  K <- (1 + sqrt(3) * D / RHO_TRUE) * exp(-sqrt(3) * D / RHO_TRUE)
  Fld <- tanh(1.5 * (t(chol(K + diag(1e-6, n))) %*% matrix(rnorm(n * E), n, E)))
  W <- array(0, c(p, p, n))
  for (e in seq_len(E)) { W[edges[e, 1], edges[e, 2], ] <- Fld[, e]; W[edges[e, 2], edges[e, 1], ] <- Fld[, e] }
  lam_min <- vapply(seq_len(n), function(i)
    min(eigen(W[, , i], symmetric = TRUE, only.values = TRUE)$values), numeric(1))
  d <- 1.25 * max(-lam_min)
  draw <- function() { Y <- matrix(0, n, p)
    for (i in seq_len(n)) Y[i, ] <- backsolve(chol(diag(d, p) + W[, , i]), rnorm(p)); scale(Y) }
  list(Y1 = draw(), Y2 = draw(), d = d, beta_rms = sqrt(mean(Fld^2)) / d)
}

run_gene_task <- function(task) {
  out_file <- file.path(OUT_DIR, "tasks_gene", sprintf("n%d_r%02d.rds", task$n, task$rep))
  if (file.exists(out_file)) return(readRDS(out_file))
  set.seed(task$seed)
  M_BASIS <<- 2; NMC <<- 3000; BURN <<- 1000; THIN <<- 5      # as in the gene application
  coords <- matrix(runif(task$n * 2, 0, 10), task$n, 2)
  adj <- generate_scale_free_adj(P_GENE, PI_GENE)
  dat <- generate_pair(adj, coords)
  t0 <- proc.time()["elapsed"]
  g1 <- run_horseshoe(dat$Y1, coords, "hsgp", TRUE); rt <- unname(proc.time()["elapsed"] - t0)
  g2 <- run_horseshoe(dat$Y2, coords, "hsgp", TRUE)
  others <- list(HS_const = run_horseshoe(dat$Y1, coords, "const", FALSE), Glasso = run_glasso(dat$Y1),
                 Nodewise_Lasso = run_nodewise_lasso(dat$Y1), Corr_Threshold = run_corr_threshold(dat$Y1))
  row <- function(nm, res) data.frame(n = task$n, rep = task$rep, method = nm,
                                      t(compute_metrics(res$adj, adj)), t(compute_auc(res$score, adj)))
  met <- rbind(row("GP_GHS", g1), row("GP_GHS_replicate", g2),
               do.call(rbind, lapply(names(others), function(nm) row(nm, others[[nm]]))))
  up <- upper.tri(adj); x <- g1$adj[up] == 1; y <- g2$adj[up] == 1; tr <- adj[up] == 1
  rep_df <- data.frame(n = task$n, rep = task$rep, n_true = sum(tr), d = dat$d, beta_rms = dat$beta_rms,
                       edges1 = sum(x), edges2 = sum(y), both = sum(x & y), both_true = sum(x & y & tr),
                       overlap = sum(x & y) / max(min(sum(x), sum(y)), 1),
                       expected = sum(x) * sum(y) / sum(up), runtime = rt)
  out <- list(met = met, rep = rep_df)
  saveRDS(out, out_file); out
}

dir.create(file.path(OUT_DIR, "tasks_gene"), showWarnings = FALSE)
tasks <- expand.grid(n = c(800, 400, 200), rep = 1:6)
tasks$seed <- SEED_G + seq_len(nrow(tasks))
cl <- makeCluster(min(as.integer(Sys.getenv("SIM_N_CORES", "18")), nrow(tasks)))
invisible(clusterEvalQ(cl, { Sys.setenv(SIM_SKIP_RUN = "1"); source("simulation_study_v4.R") }))
clusterExport(cl, c("generate_pair", "run_gene_task", "P_GENE", "PI_GENE"))
t0 <- proc.time()["elapsed"]
res <- clusterApplyLB(cl, split(tasks, seq_len(nrow(tasks))), function(tk)
  tryCatch(run_gene_task(tk), error = function(e) list(error = conditionMessage(e))))
stopCluster(cl)
bad <- vapply(res, function(r) !is.null(r$error), logical(1))
if (any(bad)) print(unique(unlist(lapply(res[bad], function(r) r$error))))
met <- do.call(rbind, lapply(res[!bad], function(r) r$met)); rp <- do.call(rbind, lapply(res[!bad], function(r) r$rep))
saveRDS(list(met = met, rep = rp), file.path(OUT_DIR, "simulation_gene_regime_raw_v4.rds"))
num <- c("F1", "MCC", "TPR", "FDP", "FP", "AUPRC")
tab <- aggregate(as.matrix(met[, num]) ~ n + method, met, mean, na.action = na.pass)
tab$n_rep <- aggregate(F1 ~ n + method, met, length)$F1
write.csv(tab[order(tab$n, tab$method), ], file.path(OUT_DIR, "simulation_gene_regime_table_v4.csv"), row.names = FALSE)
cat(sprintf("Done in %.1f min | p = %d, density %.2f, m = 2\n", (proc.time()["elapsed"] - t0) / 60, P_GENE, PI_GENE))
print(format(tab[order(tab$n, tab$method), ], digits = 3), row.names = FALSE)
cat("\nReplicate datasets from the same model (means by n):\n")
rt <- aggregate(cbind(n_true, d, beta_rms, edges1, edges2, both, both_true, overlap, expected, runtime) ~ n, rp, mean)
write.csv(rt, file.path(OUT_DIR, "simulation_gene_regime_replicates_v4.csv"), row.names = FALSE)
print(format(rt, digits = 3), row.names = FALSE)
