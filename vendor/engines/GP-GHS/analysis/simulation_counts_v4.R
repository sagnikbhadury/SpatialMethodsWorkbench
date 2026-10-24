## ================================================================
## Count-based simulation (2026-10-03): the whole application pipeline.
##
## Latent abundance fields follow the spatially varying Gaussian graphical
## model of simulation_study_v4.R (varying scenario) at the centres of a
## 15 x 15 grid of quadrats. The observed data are Poisson counts,
##   N_s(l) ~ Poisson( exp{a + b * Y_s(l)} ),
## with a chosen so that the mean count per quadrat is MU. The methods
## receive sqrt(N), standardised, exactly as in the CRC analysis
## (n = 225, m = 4, 3,000 iterations). The truth is the graph of the
## latent model. Nothing is adjusted after the metrics are computed.
##
## Usage: Rscript simulation_counts_v4.R
## ================================================================

Sys.setenv(SIM_SKIP_RUN = "1")
source("simulation_study_v4.R")

OUT_DIR  <- "simulation_v4"
N_GRID   <- 15
B_LOAD   <- 0.7                 # log-intensity per sd of the latent field
SEED_CNT <- 20261201

run_count_task <- function(task) {
  set.seed(task$seed)
  M_BASIS <<- 4; NMC <<- 3000; BURN <<- 1000; THIN <<- 5
  g <- (seq_len(N_GRID) - 0.5) / N_GRID * 10
  coords <- as.matrix(expand.grid(X = g, Y = g))
  adj <- generate_scale_free_adj(task$p, task$pi)
  lat <- generate_dataset(adj, coords, "varying")$Y        # standardised latent fields
  a   <- log(task$mu) - B_LOAD^2 / 2
  cnt <- matrix(rpois(length(lat), exp(a + B_LOAD * lat)), nrow(lat))
  Y   <- scale(sqrt(cnt)); Y[!is.finite(Y)] <- 0
  methods <- list(
    GP_GHS   = function() run_horseshoe(Y, coords, "hsgp", TRUE),
    HS_const = function() run_horseshoe(Y, coords, "const", FALSE),
    Glasso   = function() run_glasso(Y),
    Nodewise_Lasso = function() run_nodewise_lasso(Y),
    Corr_Threshold = function() run_corr_threshold(Y))
  do.call(rbind, lapply(names(methods), function(nm) {
    res <- methods[[nm]]()
    data.frame(p = task$p, mu = task$mu, pi = task$pi, rep = task$rep, method = nm,
               t(compute_metrics(res$adj, adj)), t(compute_auc(res$score, adj)),
               zero_share = mean(cnt == 0))
  }))
}

tasks <- expand.grid(p = 15, mu = c(1, 0.3), pi = c(0, 0.1, 0.3, 0.5), rep = 1:20)
tasks$seed <- SEED_CNT + seq_len(nrow(tasks))
cl <- makeCluster(as.integer(Sys.getenv("SIM_N_CORES", "20")))
invisible(clusterEvalQ(cl, {
  Sys.setenv(SIM_SKIP_RUN = "1"); source("simulation_study_v4.R") }))
clusterExport(cl, c("N_GRID", "B_LOAD", "run_count_task"))
res <- clusterApplyLB(cl, split(tasks, seq_len(nrow(tasks))), function(tk)
  tryCatch(run_count_task(tk), error = function(e) data.frame(error = conditionMessage(e))))
stopCluster(cl)
bad <- vapply(res, function(r) !is.null(r$error), logical(1))
if (any(bad)) print(unique(do.call(rbind, res[bad])))
raw <- do.call(rbind, res[!bad])
saveRDS(raw, file.path(OUT_DIR, "simulation_counts_raw_v4.rds"))
num <- c("F1", "MCC", "TPR", "FDP", "FP", "AUPRC", "zero_share")
tab <- aggregate(as.matrix(raw[, num]) ~ mu + pi + method, raw, mean, na.action = na.pass)
tab$n_rep <- aggregate(F1 ~ mu + pi + method, raw, length)$F1
tab <- tab[order(-tab$mu, tab$pi, tab$method), ]
write.csv(tab, file.path(OUT_DIR, "simulation_counts_table_v4.csv"), row.names = FALSE)
print(format(tab, digits = 3), row.names = FALSE)
cat("\nPooled over pi > 0:\n")
print(aggregate(as.matrix(raw[raw$pi > 0, num]) ~ mu + method, raw[raw$pi > 0, ], mean,
                na.action = na.pass), digits = 3)
