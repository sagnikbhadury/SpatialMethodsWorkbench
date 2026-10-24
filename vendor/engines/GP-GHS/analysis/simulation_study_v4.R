## ================================================================
## Simulation Study v4: GP-GHS vs competitors (2026-10-03)
##
## Replaces simulation_study_v3.R (archive_pre_20261003/). Changes:
##   - Data come from a spatially varying Gaussian graphical model
##       Y(l) ~ N_p(0, Omega(l)^{-1}),  Omega(l) = d I + W(l),
##     so the true graph IS the conditional-independence graph and the
##     nodewise coefficients are beta_ss'(l) = -W_ss'(l) / d.
##   - Two scenarios: "varying" (sign-changing smooth W) and "constant"
##     (W constant over space, where GLasso and the lasso are well specified).
##   - m = 5 basis functions per coordinate, as in the imaging analysis.
##   - Three horseshoe fits that differ in one thing at a time:
##       HS_const   constant coefficients, scalar horseshoe
##       HS_basis   varying coefficients, one local scale per basis coefficient
##       GP_GHS     varying coefficients, one local scale per block
##     All three use the same rule: posterior mean kappa < 0.5, AND rule,
##     with kappa the shrinkage factor of the least-shrunk coefficient.
##   - Exact GP: posterior variance computed with the correct triangular
##     solve, coordinates on the same normalised scale as GP-GHS.
##   - Every method returns a continuous edge score, so AUPRC and AUROC
##     are reported besides the fixed operating point.
##   - Summaries are written exactly as computed. Nothing is adjusted.
##
## Usage:  Rscript simulation_study_v4.R            (full study)
##         SIM_PILOT=1 Rscript simulation_study_v4.R (one task per cell, timed)
## ================================================================

suppressPackageStartupMessages({
  library(MASS); library(parallel); library(glmnet); library(glasso)
})
source("gp_group_horseshoe_spatial_graph.R")

N         <- 600
M_BASIS   <- 5
NMC <- 2000; BURN <- 500; THIN <- 5
RHO_TRUE  <- 3          # Matern-3/2 range of the true fields on [0,10]^2
SEED_BASE <- 20261101
OUT_DIR   <- Sys.getenv("SIM_OUT_DIR", "simulation_v4")
N_CORES   <- as.integer(Sys.getenv("SIM_N_CORES", "20"))
PILOT     <- identical(Sys.getenv("SIM_PILOT", "0"), "1")
dir.create(file.path(OUT_DIR, "tasks"), recursive = TRUE, showWarnings = FALSE)


## ----------------------------------------------------------------
## 1. Scale-free adjacency (Barabasi-Albert tree, then edges added or
##    removed at random to reach round(edge_prob * p(p-1)/2) edges)
## ----------------------------------------------------------------
generate_scale_free_adj <- function(p, edge_prob) {
  A <- matrix(0, p, p)
  if (edge_prob == 0) return(A)
  A[1,2] <- A[2,1] <- 1
  degree <- c(1, 1, rep(0, p - 2))
  for (v in 3:p) {
    prob <- degree[1:(v-1)] / sum(degree[1:(v-1)])
    tgt  <- sample(1:(v-1), 1, prob = prob)
    A[v, tgt] <- A[tgt, v] <- 1
    degree[v] <- degree[v] + 1; degree[tgt] <- degree[tgt] + 1
  }
  max_e <- p*(p-1)/2; tgt_e <- round(edge_prob * max_e)
  cur_e <- sum(A[upper.tri(A)])
  if (tgt_e > cur_e) {
    cands <- which(upper.tri(A) & A == 0, arr.ind = TRUE)
    idx   <- sample(nrow(cands), min(tgt_e - cur_e, nrow(cands)))
    for (k in idx) { A[cands[k,1], cands[k,2]] <- A[cands[k,2], cands[k,1]] <- 1 }
  } else if (tgt_e < cur_e) {
    ons <- which(upper.tri(A) & A == 1, arr.ind = TRUE)
    idx <- sample(nrow(ons), cur_e - tgt_e)
    for (k in idx) { A[ons[k,1], ons[k,2]] <- A[ons[k,2], ons[k,1]] <- 0 }
  }
  diag(A) <- 0; A
}


## ----------------------------------------------------------------
## 2. Data: spatially varying Gaussian graphical model
##    Returns Y (columns standardised) and the true coefficient
##    fields on the standardised scale, beta_true[[s]] = n x (p-1).
## ----------------------------------------------------------------
generate_dataset <- function(adj, coords, scenario) {
  n <- nrow(coords); p <- ncol(adj)
  edges <- which(upper.tri(adj) & adj == 1, arr.ind = TRUE)
  E <- nrow(edges)
  W <- array(0, c(p, p, n))
  if (E > 0) {
    if (scenario == "varying") {
      D <- as.matrix(dist(coords))
      K <- (1 + sqrt(3) * D / RHO_TRUE) * exp(-sqrt(3) * D / RHO_TRUE)
      Lk <- t(chol(K + diag(1e-6, n)))
      Fld <- tanh(1.5 * (Lk %*% matrix(rnorm(n * E), n, E)))
    } else {
      cc <- sample(c(-1, 1), E, replace = TRUE) * runif(E, 0.5, 1)
      Fld <- matrix(cc, n, E, byrow = TRUE)
    }
    for (e in seq_len(E)) {
      W[edges[e, 1], edges[e, 2], ] <- Fld[, e]
      W[edges[e, 2], edges[e, 1], ] <- Fld[, e]
    }
    lam_min <- vapply(seq_len(n), function(i)
      min(eigen(W[, , i], symmetric = TRUE, only.values = TRUE)$values),
      numeric(1))
    d <- 1.25 * max(-lam_min)
  } else d <- 1
  Y <- matrix(0, n, p)
  for (i in seq_len(n)) {
    R <- chol(diag(d, p) + W[, , i])
    Y[i, ] <- backsolve(R, rnorm(p))
  }
  sds <- apply(Y, 2, sd)
  beta_true <- lapply(seq_len(p), function(s) {
    nb <- (1:p)[-s]
    sapply(nb, function(t) -W[s, t, ] / d * sds[t] / sds[s])
  })
  list(Y = scale(Y), beta_true = beta_true, d = d)
}


## ----------------------------------------------------------------
## 3. Metrics
## ----------------------------------------------------------------
compute_metrics <- function(est_adj, true_adj) {
  tv <- true_adj[upper.tri(true_adj)]; ev <- est_adj[upper.tri(est_adj)]
  TP <- sum(ev==1 & tv==1); FP <- sum(ev==1 & tv==0)
  TN <- sum(ev==0 & tv==0); FN <- sum(ev==0 & tv==1)
  TPR <- TP/max(TP+FN,1); FDP <- FP/max(TP+FP,1)
  F1  <- ifelse(2*TP+FP+FN > 0, 2*TP/(2*TP+FP+FN), 0)
  den <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC <- ifelse(den>0, (TP*TN - FP*FN)/den, 0)
  c(F1=F1, FDP=FDP, TPR=TPR, MCC=MCC, TP=TP, FP=FP, TN=TN, FN=FN)
}

## Area under the precision-recall and ROC curves from a symmetric score
## matrix; ties are handled by stepping over unique thresholds.
compute_auc <- function(score, true_adj) {
  tv <- true_adj[upper.tri(true_adj)]; sc <- score[upper.tri(score)]
  if (sum(tv) == 0 || sum(tv) == length(tv)) return(c(AUPRC=NA, AUROC=NA))
  thr <- sort(unique(sc), decreasing = TRUE)
  tp <- vapply(thr, function(h) sum(sc >= h & tv == 1), numeric(1))
  fp <- vapply(thr, function(h) sum(sc >= h & tv == 0), numeric(1))
  rec <- tp / sum(tv); prec <- tp / (tp + fp); fpr <- fp / sum(tv == 0)
  ap  <- sum(diff(c(0, rec)) * prec)
  roc <- sum(diff(c(0, fpr)) * (rec + c(0, rec[-length(rec)])) / 2)
  if (max(fpr) < 1) roc <- roc + (1 - max(fpr)) * (1 + max(rec)) / 2
  c(AUPRC = ap, AUROC = roc)
}

## Field recovery on the true directed edges: Pearson correlation with the
## truth and coverage of mean +/- 1.96 sd.
field_metrics <- function(mean_list, sd_list, beta_true, adj) {
  p <- ncol(adj); cors <- covs <- numeric(0)
  for (s in 1:p) {
    nb <- (1:p)[-s]
    for (pos in seq_along(nb)) {
      if (adj[s, nb[pos]] != 1) next
      bt <- beta_true[[s]][, pos]; bh <- mean_list[[s]][, pos]
      if (sd(bh) > 1e-10 && sd(bt) > 1e-10) cors <- c(cors, cor(bh, bt))
      sdj <- sd_list[[s]][, pos]
      covs <- c(covs, mean(abs(bt - bh) <= 1.96 * sdj))
    }
  }
  c(field_cor = if (length(cors)) mean(cors) else NA,
    coverage  = if (length(covs)) mean(covs) else NA)
}


## ----------------------------------------------------------------
## 4. Methods. Each returns list(adj, score, [fields]).
## ----------------------------------------------------------------

## 4a. Horseshoe family. basis = "hsgp" gives varying coefficients;
##     basis = "const" gives one constant coefficient per predictor.
run_horseshoe <- function(Y, coords, basis = "hsgp", grouped = TRUE) {
  n <- nrow(Y); p <- ncol(Y)
  if (basis == "hsgp") {
    hs <- build_hsgp_basis(coords, m = M_BASIS, nu = 1.5, rho = NULL)
    phi <- hs$phi; S_diag <- hs$S_diag
  } else {
    phi <- matrix(1, n, 1); S_diag <- 1
  }
  m2 <- ncol(phi)
  kap <- matrix(1, p, p)
  mean_list <- sd_list <- vector("list", p)
  for (s in 1:p) {
    fit <- gp_group_horseshoe_sampler(Y[, s], Y[, -s, drop = FALSE], phi,
                                      S_diag, nmc = NMC, burn = BURN,
                                      thin = THIN, grouped = grouped)
    summ <- summarize_node(fit, phi, S_diag, p - 1, m2, as.character((1:p)[-s]))
    kap[s, -s] <- summ$kappa_mean
    mean_list[[s]] <- summ$edge_maps; sd_list[[s]] <- summ$edge_maps_sd
  }
  kmax <- pmax(kap, t(kap)); diag(kmax) <- 1
  list(adj = (kmax < 0.5) * 1, score = 1 - kmax,
       mean_list = mean_list, sd_list = sd_list)
}

## 4b. Graphical lasso, BIC over 15 penalties. Score: the largest penalty
##     at which the edge is present.
run_glasso <- function(Y) {
  n <- nrow(Y); p <- ncol(Y); S <- cor(Y)
  lmax <- max(abs(S[upper.tri(S)]))
  lseq <- exp(seq(log(lmax), log(lmax * 0.02), length.out = 15))
  best_bic <- Inf; best_adj <- matrix(0, p, p); score <- matrix(0, p, p)
  for (lam in lseq) {
    Th <- glasso(S, rho = lam, penalize.diagonal = FALSE)$wi
    on <- abs(Th) > 1e-6; diag(on) <- FALSE; on <- on & t(on)
    score[on & score == 0] <- lam
    ll  <- as.numeric(determinant(Th, logarithm = TRUE)$modulus) - sum(S * Th)
    bic <- -n * ll + log(n) * sum(on[upper.tri(on)])
    if (bic < best_bic) { best_bic <- bic; best_adj <- on * 1 }
  }
  list(adj = best_adj, score = score)
}

## 4c. Nodewise lasso, node-specific five-fold CV, AND rule. Score: the
##     largest penalty at which the coefficient is nonzero (min over the
##     two directions).
run_nodewise_lasso <- function(Y) {
  p <- ncol(Y); adj <- matrix(0, p, p); ent <- matrix(0, p, p)
  for (s in 1:p) {
    cv <- cv.glmnet(Y[, -s], Y[, s], nfolds = 5, intercept = FALSE,
                    standardize = FALSE)
    B  <- as.matrix(cv$glmnet.fit$beta) != 0
    ent[s, -s] <- apply(B, 1, function(z) if (any(z)) max(cv$glmnet.fit$lambda[z]) else 0)
    adj[s, -s] <- as.numeric(coef(cv, s = "lambda.min"))[-1] != 0
  }
  list(adj = adj * t(adj), score = pmin(ent, t(ent)))
}

## 4d. Correlation threshold.
run_corr_threshold <- function(Y, threshold = 0.15) {
  R <- abs(cor(Y)); diag(R) <- 0
  list(adj = (R > threshold) * 1, score = R)
}

## 4e. Pairwise GP on partial residuals ("Exact GP"). For each ordered
##     pair, the response is adjusted for the other predictors by OLS and
##     an exact Matern-3/2 GP is fitted to r_j(l) = x_j(l) beta_j(l) + noise.
##     A directed pair is selected when the pointwise 95% interval excludes
##     zero at more than half of the locations.
run_exact_gp <- function(Y, coords) {
  n <- nrow(Y); p <- ncol(Y)
  hs <- build_hsgp_basis(coords, m = 2, nu = 1.5, rho = NULL)
  D  <- as.matrix(dist(hs$coords_norm))
  K  <- (1 + sqrt(3) * D / hs$rho) * exp(-sqrt(3) * D / hs$rho) + diag(1e-6, n)
  frac <- matrix(0, p, p)
  mean_list <- sd_list <- vector("list", p)
  for (s in 1:p) {
    y <- Y[, s]; X <- Y[, -s, drop = FALSE]; q <- p - 1
    sigma_sq <- max(var(lm.fit(X, y)$residuals), 1e-4)
    bm <- bs <- matrix(0, n, q)
    for (j in 1:q) {
      xj  <- X[, j]
      r_j <- y - X[, -j, drop = FALSE] %*% lm.fit(X[, -j, drop = FALSE], y)$coefficients
      xs  <- ifelse(abs(xj) < 1e-3, 1e-3 * ifelse(xj < 0, -1, 1), xj)
      Kn  <- K; diag(Kn) <- diag(Kn) + sigma_sq / xs^2
      R   <- chol(Kn)
      V   <- backsolve(R, K, transpose = TRUE)          # R^{-T} K
      bm[, j] <- crossprod(V, backsolve(R, r_j / xs, transpose = TRUE))
      bs[, j] <- sqrt(pmax(diag(K) - colSums(V^2), 1e-12))
      frac[s, (1:p)[-s][j]] <- mean(abs(bm[, j]) > 1.96 * bs[, j])
    }
    mean_list[[s]] <- bm; sd_list[[s]] <- bs
  }
  fmin <- pmin(frac, t(frac))
  list(adj = (fmin > 0.5) * 1, score = fmin,
       mean_list = mean_list, sd_list = sd_list)
}


## ----------------------------------------------------------------
## 5. One task = one dataset, all methods
## ----------------------------------------------------------------
run_task <- function(task) {
  out_file <- file.path(OUT_DIR, "tasks", sprintf("%s.rds", task$id))
  if (file.exists(out_file)) return(readRDS(out_file))
  set.seed(task$seed)
  coords <- matrix(runif(N * 2, 0, 10), N, 2)
  adj <- generate_scale_free_adj(task$p, task$pi)
  dat <- generate_dataset(adj, coords, task$scenario)
  Y <- dat$Y

  methods <- list(
    GP_GHS   = function() run_horseshoe(Y, coords, "hsgp", TRUE),
    HS_basis = function() run_horseshoe(Y, coords, "hsgp", FALSE),
    HS_const = function() run_horseshoe(Y, coords, "const", FALSE),
    Exact_GP = function() run_exact_gp(Y, coords),
    Glasso   = function() run_glasso(Y),
    Nodewise_Lasso = function() run_nodewise_lasso(Y),
    Corr_Threshold = function() run_corr_threshold(Y))

  ## The ungrouped ablation doubles the cost; it is run at p = 15 only.
  if (task$p > 15) methods$HS_basis <- NULL

  rows <- lapply(names(methods), function(nm) {
    t0 <- proc.time()["elapsed"]
    res <- methods[[nm]]()
    rt <- unname(proc.time()["elapsed"] - t0)
    fm <- if (!is.null(res$mean_list) && nm != "HS_const")
      field_metrics(res$mean_list, res$sd_list, dat$beta_true, adj)
    else c(field_cor = NA, coverage = NA)
    data.frame(id = task$id, p = task$p, scenario = task$scenario,
               pi = task$pi, rep = task$rep, method = nm,
               t(compute_metrics(res$adj, adj)), t(compute_auc(res$score, adj)),
               t(fm), runtime = rt, d = dat$d, n_edges = sum(adj) / 2)
  })
  out <- do.call(rbind, rows)
  saveRDS(out, out_file)
  out
}


## ----------------------------------------------------------------
## 6. Design and run
## ----------------------------------------------------------------
make_tasks <- function() {
  g <- rbind(
    expand.grid(p = 15, scenario = c("varying", "constant"),
                pi = c(0, 0.1, 0.3, 0.5, 0.7), rep = 1:20,
                stringsAsFactors = FALSE),
    expand.grid(p = 25, scenario = "varying",
                pi = c(0.1, 0.3, 0.5, 0.7), rep = 1:10,
                stringsAsFactors = FALSE))
  g$id   <- sprintf("p%d_%s_pi%02d_r%02d", g$p, g$scenario, round(10 * g$pi), g$rep)
  g$seed <- SEED_BASE + seq_len(nrow(g))
  g
}

## SIM_SKIP_RUN=1 lets another script source the functions above only.
if (!identical(Sys.getenv("SIM_SKIP_RUN", "0"), "1")) {

tasks <- make_tasks()
if (PILOT) tasks <- tasks[tasks$rep == 1 & tasks$pi %in% c(0.1, 0.5), ]
## longest tasks first so the cluster stays balanced
tasks <- tasks[order(-tasks$p), ]
task_list <- split(tasks, seq_len(nrow(tasks)))

cat(sprintf("Simulation v4 | %d tasks | %d cores | out: %s\n",
            length(task_list), N_CORES, OUT_DIR))
t0 <- proc.time()["elapsed"]
cl <- makeCluster(min(N_CORES, length(task_list)))
clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(MASS); library(glmnet); library(glasso) })
  source("gp_group_horseshoe_spatial_graph.R")
})
clusterExport(cl, setdiff(ls(globalenv()), c("cl", "task_list", "tasks")))
res <- clusterApplyLB(cl, task_list, function(tk)
  tryCatch(run_task(tk), error = function(e)
    data.frame(id = tk$id, error = conditionMessage(e))))
stopCluster(cl)
cat(sprintf("Done in %.1f min\n", (proc.time()["elapsed"] - t0) / 60))

ok <- vapply(res, function(r) is.null(r$error), logical(1))
if (any(!ok)) print(do.call(rbind, res[!ok]))
raw <- do.call(rbind, res[ok])
saveRDS(raw, file.path(OUT_DIR, if (PILOT) "pilot_raw.rds" else "simulation_raw_v4.rds"))

num <- c("F1","FDP","TPR","MCC","FP","AUPRC","AUROC","field_cor","coverage","runtime")
key <- interaction(raw$p, raw$scenario, raw$pi, raw$method, drop = TRUE)
summ <- do.call(rbind, lapply(split(raw, key), function(d) data.frame(
  p = d$p[1], scenario = d$scenario[1], pi = d$pi[1], method = d$method[1],
  n_rep = nrow(d),
  t(setNames(colMeans(d[, num], na.rm = TRUE), paste0(num, "_mean"))),
  t(setNames(apply(d[, num], 2, sd, na.rm = TRUE), paste0(num, "_sd"))))))
summ <- summ[order(summ$p, summ$scenario, summ$pi, summ$method), ]
rownames(summ) <- NULL
write.csv(summ, file.path(OUT_DIR, if (PILOT) "pilot_summary.csv" else "simulation_summary_v4.csv"),
          row.names = FALSE)
print(summ[, c("p","scenario","pi","method","F1_mean","FDP_mean","TPR_mean",
               "MCC_mean","AUPRC_mean","field_cor_mean","coverage_mean",
               "runtime_mean")], digits = 3)

}
