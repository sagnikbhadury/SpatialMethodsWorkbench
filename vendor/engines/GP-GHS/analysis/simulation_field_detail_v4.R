## ================================================================
## Coefficient-function recovery in detail (2026-10-22).
##
## Added for the JASA ACS revision after J. Gaskins asked for more
## empirical results on the estimated coefficient functions: the
## intervals undercover, but are the estimates directionally correct?
##
## Re-generates a subset of the v4 simulation datasets from their task
## seeds (same code path as simulation_study_v4.R), fits GP-GHS only,
## and records for every true directed edge:
##   cor          Pearson correlation of the posterior mean with the truth
##   coverage     share of locations with the truth inside mean +/- 1.96 sd
##   sign_all     share of locations where the posterior mean has the sign
##                of the truth
##   sign_strong  the same over the half of the locations where |truth| is
##                largest
##   mean_sign    1 if the spatial means of estimate and truth agree in sign
##   truth_flips  1 if the truth is positive at > 10% of locations and
##                negative at > 10% (changes sign over the domain)
##   est_flips    the same for the posterior mean
##   rmse, mean_sd, ratio  root-mean-square error, mean posterior sd and
##                their ratio (how much too narrow the intervals are)
##   selected     1 if the AND rule selects the pair
##
## Usage (R 4.5.2):  Rscript simulation_field_detail_v4.R
## Output: simulation_v4/field_detail_v4.csv (per edge) and
##         simulation_v4/field_detail_summary_v4.csv (by scenario and pi)
## ================================================================

Sys.setenv(SIM_SKIP_RUN = "1")
source("simulation_study_v4.R")
N_CORES_FD <- as.integer(Sys.getenv("FD_N_CORES", "10"))
OUT_DIR    <- Sys.getenv("SIM_OUT_DIR", "simulation_v4")

tasks <- make_tasks()
sel <- tasks[tasks$p == 15 & tasks$pi %in% c(0.1, 0.3, 0.5, 0.7) & tasks$rep <= 5, ]
cat(sprintf("Field detail | %d datasets | %d cores\n", nrow(sel), N_CORES_FD))

one_task <- function(task) {
  set.seed(task$seed)
  coords <- matrix(runif(N * 2, 0, 10), N, 2)
  adj <- generate_scale_free_adj(task$p, task$pi)
  dat <- generate_dataset(adj, coords, task$scenario)
  res <- run_horseshoe(dat$Y, coords, "hsgp", TRUE)
  p <- task$p; rows <- list()
  for (s in 1:p) {
    nb <- (1:p)[-s]
    for (pos in seq_along(nb)) {
      t <- nb[pos]
      if (adj[s, t] != 1) next
      bt <- dat$beta_true[[s]][, pos]; bh <- res$mean_list[[s]][, pos]
      sdj <- res$sd_list[[s]][, pos]
      strong <- abs(bt) >= median(abs(bt))
      rows[[length(rows) + 1]] <- data.frame(
        id = task$id, scenario = task$scenario, pi = task$pi, rep = task$rep,
        s = s, t = t,
        cor = if (sd(bt) > 1e-10 && sd(bh) > 1e-10) cor(bh, bt) else NA,
        coverage = mean(abs(bt - bh) <= 1.96 * sdj),
        sign_all = mean(sign(bh) == sign(bt)),
        sign_strong = mean(sign(bh[strong]) == sign(bt[strong])),
        mean_sign = as.integer(sign(mean(bh)) == sign(mean(bt))),
        truth_flips = as.integer(mean(bt > 0) > 0.1 && mean(bt < 0) > 0.1),
        est_flips = as.integer(mean(bh > 0) > 0.1 && mean(bh < 0) > 0.1),
        rmse = sqrt(mean((bt - bh)^2)), mean_sd = mean(sdj),
        rms_truth = sqrt(mean(bt^2)),
        selected = res$adj[s, t])
    }
  }
  do.call(rbind, rows)
}

task_list <- split(sel, seq_len(nrow(sel)))
cl <- makeCluster(min(N_CORES_FD, length(task_list)))
clusterEvalQ(cl, { suppressPackageStartupMessages(library(MASS)); source("gp_group_horseshoe_spatial_graph.R") })
clusterExport(cl, setdiff(ls(globalenv()), c("cl", "task_list")))
t0 <- proc.time()["elapsed"]
res <- clusterApplyLB(cl, task_list, function(tk)
  tryCatch(one_task(tk), error = function(e) data.frame(id = tk$id, error = conditionMessage(e))))
stopCluster(cl)
cat(sprintf("Done in %.1f min\n", (proc.time()["elapsed"] - t0) / 60))
ok <- vapply(res, function(r) is.null(r$error), logical(1))
if (any(!ok)) print(do.call(rbind, res[!ok]))
det <- do.call(rbind, res[ok])
det$ratio <- det$rmse / det$mean_sd
write.csv(det, file.path(OUT_DIR, "field_detail_v4.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(split(det, interaction(det$scenario, det$pi, drop = TRUE)), function(d) data.frame(
  scenario = d$scenario[1], pi = d$pi[1], n_edges = nrow(d), n_datasets = length(unique(d$id)),
  cor = mean(d$cor, na.rm = TRUE), coverage = mean(d$coverage),
  sign_all = mean(d$sign_all), sign_strong = mean(d$sign_strong),
  mean_sign = mean(d$mean_sign), truth_flips = mean(d$truth_flips),
  est_flips_given_truth_flips = mean(d$est_flips[d$truth_flips == 1]),
  est_flips_given_truth_const = if (any(d$truth_flips == 0)) mean(d$est_flips[d$truth_flips == 0]) else NA,
  rmse = mean(d$rmse), mean_sd = mean(d$mean_sd), ratio = median(d$ratio),
  rms_truth = mean(d$rms_truth), selected = mean(d$selected),
  cor_selected = mean(d$cor[d$selected == 1], na.rm = TRUE),
  sign_all_selected = mean(d$sign_all[d$selected == 1]),
  coverage_selected = mean(d$coverage[d$selected == 1]))))
rownames(summ) <- NULL
write.csv(summ, file.path(OUT_DIR, "field_detail_summary_v4.csv"), row.names = FALSE)
print(summ, digits = 3)
