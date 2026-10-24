## ================================================================
## Spatial transcriptomics application, step 2: GP-GHS per zone
##
## For each pathologist zone of the slide prepared by st_prepare.R, the
## nodewise regression of gp_group_horseshoe_spatial_graph.R is fitted to
## every estimable gene: the spots of the zone are the locations, the
## log-normalised expression of the 100 selected genes, standardised
## within the zone, is the data. A gene is estimable in a zone when it is
## detected in at least ST_MIN_DET of the zone's spots; a gene that is not
## estimable is not a response and its column is set to zero.
##
## One task = one regression (zone, response gene). Tasks are saved one
## file each, so an interrupted run resumes where it stopped.
##
## Usage (R 4.4.1): Rscript st_fit.R                                  (primary, m = 2)
##                  ST_PERMUTE=1 ST_TAG=_permuted Rscript st_fit.R    (negative control)
##                  ST_SEED=777 ST_TAG=_seed2 Rscript st_fit.R        (second chain)
## ================================================================

suppressPackageStartupMessages(library(parallel))
source("gp_group_horseshoe_spatial_graph.R")

OUT_DIR  <- Sys.getenv("ST_OUT_DIR", "st_results")
M_BASIS  <- as.integer(Sys.getenv("ST_M_BASIS", "2"))
NMC      <- as.integer(Sys.getenv("ST_NMC", "3000"))
BURN     <- as.integer(Sys.getenv("ST_BURN", "1000"))
THIN     <- as.integer(Sys.getenv("ST_THIN", "5"))
MIN_DET  <- as.numeric(Sys.getenv("ST_MIN_DET", "0.10"))
N_CORES  <- as.integer(Sys.getenv("ST_N_CORES", "18"))
SEED     <- as.integer(Sys.getenv("ST_SEED", "20261004"))
## Negative control: ST_PERMUTE=1 permutes every gene independently over the
## spots of a zone, which keeps each marginal distribution and removes all
## dependence between genes. ST_TAG names the output folder of such a run.
PERMUTE  <- identical(Sys.getenv("ST_PERMUTE", "0"), "1")
FIT_DIR  <- file.path(OUT_DIR, sprintf("fit_m%d%s", M_BASIS, Sys.getenv("ST_TAG", "")))
dir.create(file.path(FIT_DIR, "tasks"), recursive = TRUE, showWarnings = FALSE)

prep  <- readRDS(file.path(OUT_DIR, "st_prepared.rds"))
genes <- prep$genes; p <- length(genes); K <- nrow(prep$zone_tab)

## per-zone data: standardised expression, basis, estimable genes
ZD <- lapply(seq_len(K), function(z) {
  idx <- which(prep$obs$zone == z)
  E <- prep$expr[idx, , drop = FALSE]
  if (PERMUTE) { set.seed(SEED + z); E <- apply(E, 2, sample) }
  est <- colMeans(E > 0) >= MIN_DET & apply(E, 2, sd) > 0
  Y <- scale(E); Y[, !est] <- 0; Y[!is.finite(Y)] <- 0
  coords <- as.matrix(prep$obs[idx, c("x", "y")])
  hs <- build_hsgp_basis(coords, m = M_BASIS, nu = 1.5, rho = NULL)
  list(idx = idx, Y = Y, est = est, coords = coords, phi = hs$phi, S_diag = hs$S_diag, rho = hs$rho)
})

run_task <- function(tk) {
  out_file <- file.path(FIT_DIR, "tasks", sprintf("z%d_g%03d.rds", tk$zone, tk$s))
  if (file.exists(out_file)) return(readRDS(out_file))
  d <- ZD[[tk$zone]]; s <- tk$s; m2 <- ncol(d$phi); q <- p - 1
  set.seed(SEED + 1000L * tk$zone + s)
  t0 <- proc.time()["elapsed"]
  fit <- gp_group_horseshoe_sampler(d$Y[, s], d$Y[, -s, drop = FALSE], d$phi, d$S_diag,
                                    nmc = NMC, burn = BURN, thin = THIN)
  theta_mean <- rowMeans(fit$theta_store)
  fld <- d$phi %*% matrix(theta_mean, m2, q)          # posterior mean coefficient functions, n x q
  fitted <- rowSums(d$Y[, -s, drop = FALSE] * fld)
  out <- list(zone = tk$zone, s = s,
              kappa_mean = rowMeans(fit$kappa_store), kappa_avg_mean = rowMeans(fit$kappa_avg_store),
              theta_mean = theta_mean, beta_mean = colMeans(fld), beta_rms = sqrt(colMeans(fld^2)),
              share_pos = colMeans(fld > 0),
              tau_mean = mean(sqrt(fit$tau_store)), sigma_mean = mean(sqrt(fit$sigma_store)),
              r2 = 1 - var(d$Y[, s] - fitted) / var(d$Y[, s]),
              runtime = unname(proc.time()["elapsed"] - t0))
  saveRDS(out, out_file)
  out
}

tasks <- do.call(rbind, lapply(seq_len(K), function(z)
  data.frame(zone = z, s = which(ZD[[z]]$est), n = length(ZD[[z]]$idx))))
tasks <- tasks[order(-tasks$n), ]                      # longest first
cat(sprintf("ST GP-GHS | slide %s | m = %d (%d coefficients per regression) | %d regressions | %d cores\n",
            prep$slide, M_BASIS, (p - 1) * M_BASIS^2, nrow(tasks), N_CORES))
for (z in seq_len(K)) cat(sprintf("  zone %d %-36s n = %3d | estimable genes %3d | rho = %.2f\n",
                                  z, prep$zone_tab$name[z], length(ZD[[z]]$idx), sum(ZD[[z]]$est), ZD[[z]]$rho))
t0 <- proc.time()["elapsed"]
cl <- makeCluster(min(N_CORES, nrow(tasks)))
invisible(clusterEvalQ(cl, source("gp_group_horseshoe_spatial_graph.R")))
clusterExport(cl, c("ZD", "FIT_DIR", "NMC", "BURN", "THIN", "SEED", "p", "run_task"))
res <- clusterApplyLB(cl, split(tasks, seq_len(nrow(tasks))), function(tk)
  tryCatch(run_task(tk), error = function(e) list(zone = tk$zone, s = tk$s, error = conditionMessage(e))))
stopCluster(cl)
bad <- vapply(res, function(r) !is.null(r$error), logical(1))
if (any(bad)) for (r in res[bad]) cat("ERROR zone", r$zone, "gene", r$s, ":", r$error, "\n")
res <- res[!bad]

## assemble per zone: directed p x p matrices (row = response, column = predictor)
fits <- lapply(seq_len(K), function(z) {
  blank <- matrix(NA_real_, p, p, dimnames = list(genes, genes))
  kap <- bmean <- brms <- spos <- blank
  tau <- sigma <- r2 <- setNames(rep(NA_real_, p), genes)
  theta <- vector("list", p)
  for (r in res[vapply(res, function(r) r$zone == z, logical(1))]) {
    s <- r$s; nb <- (1:p)[-s]
    kap[s, nb] <- r$kappa_mean; bmean[s, nb] <- r$beta_mean; brms[s, nb] <- r$beta_rms
    spos[s, nb] <- r$share_pos; tau[s] <- r$tau_mean; sigma[s] <- r$sigma_mean; r2[s] <- r$r2
    theta[[s]] <- r$theta_mean
  }
  est <- ZD[[z]]$est
  kap[!est, ] <- NA; kap[, !est] <- NA                  # pairs with a non-estimable gene
  kmax <- pmax(kap, t(kap))
  adj <- (kmax < 0.5) * 1; adj[is.na(adj)] <- 0; diag(adj) <- 0
  list(zone = z, name = prep$zone_tab$name[z], idx = ZD[[z]]$idx, est = est, coords = ZD[[z]]$coords,
       Y = ZD[[z]]$Y, phi = ZD[[z]]$phi, rho = ZD[[z]]$rho,
       kappa = kap, kappa_max = kmax, adj = adj, beta_mean = bmean, beta_rms = brms, share_pos = spos,
       tau = tau, sigma = sigma, r2 = r2, theta = theta)
})
saveRDS(list(fits = fits, genes = genes,
             settings = list(m_basis = M_BASIS, nmc = NMC, burn = BURN, thin = THIN,
                             min_det = MIN_DET, seed = SEED, permute = PERMUTE)),
        file.path(FIT_DIR, "st_fit.rds"))
rt <- vapply(res, function(r) r$runtime, numeric(1))
cat(sprintf("Done in %.1f min | seconds per regression (quartiles): %s\n",
            (proc.time()["elapsed"] - t0) / 60, paste(round(quantile(rt, c(.25, .5, .75)), 1), collapse = " ")))
for (f in fits) cat(sprintf("  zone %d %-36s edges %4d of %4d estimable pairs (%.3f) | median sigma %.2f | median R2 %.2f\n",
                            f$zone, f$name, sum(f$adj) / 2, choose(sum(f$est), 2),
                            sum(f$adj) / 2 / choose(sum(f$est), 2),
                            median(f$sigma, na.rm = TRUE), median(f$r2, na.rm = TRUE)))
