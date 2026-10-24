## ================================================================
## CRC Spatial Interaction Analysis
## Pipeline:
##   1. KDE-based cell type density features per image
##   2. GP-GHS spatial interaction graph per image
##   3. Differential edge testing: CLR (pathology=1) vs DII (pathology=2)
##
## Data: 140 images, 35 patients (4 images each), 15 retained cell types
##       pathology 1 = CLR (17 patients, 68 images)
##       pathology 2 = DII (18 patients, 72 images)
## ================================================================

source("gp_group_horseshoe_spatial_graph.R")
library(MASS)
library(parallel)

## ----------------------------------------------------------------
## 0. Config
## ----------------------------------------------------------------
DATA_PATH  <- "CRC_data_cleaned_corrected.csv"
OUT_DIR    <- Sys.getenv("CRC_OUT_DIR", "crc_results_v2")
N_CORES    <- as.integer(Sys.getenv(
  "CRC_N_CORES", max(1L, detectCores(logical = FALSE) - 1L)))

## HSGP + MCMC settings
M_BASIS    <- as.integer(Sys.getenv("CRC_M_BASIS", "5"))  # m^2 basis functions per field
NU         <- 1.5    # Matern 3/2
NMC        <- as.integer(Sys.getenv("CRC_NMC", "3000"))
BURN       <- as.integer(Sys.getenv("CRC_BURN", "1000"))
THIN       <- as.integer(Sys.getenv("CRC_THIN", "5"))

## KDE grid resolution (n_grid x n_grid per image)
N_GRID     <- as.integer(Sys.getenv("CRC_N_GRID", "30"))  # N_GRID^2 locations per image
## Local-abundance feature: "kde" (kernel density on the grid) or "count"
## (number of cells in each grid cell, square-root transformed).
FEATURE    <- Sys.getenv("CRC_FEATURE", "kde")
KDE_BW     <- "nrd"  # normal-reference bandwidth (MASS::bandwidth.nrd)
MIN_CELLS  <- as.integer(Sys.getenv("CRC_MIN_CELLS", "5"))
VAR_THR    <- as.numeric(Sys.getenv("CRC_VAR_THR", "1e-8"))
CHAIN_ID   <- as.integer(Sys.getenv("CRC_CHAIN_ID", "1"))
SEED_BASE  <- as.integer(Sys.getenv("CRC_SEED_BASE", "2024"))
DETREND    <- identical(Sys.getenv("CRC_DETREND", "0"), "1")

## Edge threshold for GP-GHS
KAPPA_THR  <- 0.5

## Differential testing FDR threshold
FDR_THR    <- 0.05

dir.create(OUT_DIR, showWarnings = FALSE)


## ================================================================
## STEP 1: KDE-BASED FEATURE MATRIX PER IMAGE
## ================================================================
## For each image, estimate the kernel density of each cell type
## over a regular grid covering the image domain.
## This converts point-pattern data (cell locations + labels)
## into a spatially referenced matrix Y (n_grid^2 x p):
##   Y[i, s] = estimated density of cell type s at grid point i
##
## Rationale: GP-GHS requires a continuous-valued matrix at observed
## spatial locations. Raw cell counts at individual spots are too
## sparse and noisy for regression. KDE smooths each cell type's
## spatial distribution onto a common grid, producing a dense
## matrix where each column is a smooth spatial abundance surface.
## The grid points serve as the spatial coordinates for GP-GHS.
## ================================================================

compute_kde_matrix <- function(cells, cell_types, n_grid = N_GRID,
                               bw = KDE_BW) {
  ## cells: data.frame with columns X, Y, type for ONE image
  ## Returns: list(Y = n_grid^2 x p matrix, coords = n_grid^2 x 2)

  p          <- length(cell_types)

  ## Define common grid over image bounding box
  x_range <- range(cells$X); y_range <- range(cells$Y)
  ## Expand slightly to avoid boundary effects
  x_pad   <- diff(x_range) * 0.05; y_pad <- diff(y_range) * 0.05
  x_seq   <- seq(x_range[1] - x_pad, x_range[2] + x_pad, length.out = n_grid)
  y_seq   <- seq(y_range[1] - y_pad, y_range[2] + y_pad, length.out = n_grid)
  grid    <- expand.grid(X = x_seq, Y = y_seq)

  Y <- matrix(0, nrow = nrow(grid), ncol = p,
               dimnames = list(NULL, cell_types))

  for (ct in cell_types) {
    pts <- cells[which(cells$type == ct), c("X", "Y")]
    pts <- pts[complete.cases(pts), ]
    if (nrow(pts) < MIN_CELLS) {
      ## Too few cells for KDE — fill with small constant
      Y[, ct] <- 1e-6
      next
    }
    ## 2D KDE via MASS::kde2d
    h <- tryCatch({
      hx <- bandwidth.nrd(pts$X); hy <- bandwidth.nrd(pts$Y)
      c(hx, hy)
    }, error = function(e) c(diff(x_range)/10, diff(y_range)/10))
    ## Guard: if sd~0 (all cells at same location) use domain fraction
    h[1] <- max(h[1], diff(x_range) * 0.03, na.rm = TRUE)
    h[2] <- max(h[2], diff(y_range) * 0.03, na.rm = TRUE)
    kde <- kde2d(pts$X, pts$Y,
                 h = h,
                 n = n_grid,
                 lims = c(x_range[1] - x_pad, x_range[2] + x_pad,
                           y_range[1] - y_pad, y_range[2] + y_pad))

    ## expand.grid varies x fastest, matching R's column-major vectorization
    ## of kde$z (rows=x, columns=y).
    Y[, ct] <- as.vector(kde$z)
  }

  ## Log-transform + scale: stabilizes variance, puts densities on
  ## comparable scale across cell types with very different abundances
  Y <- log1p(Y * 1e4)  # scale up before log to avoid log(tiny)
  Y <- scale(Y)

  list(Y = Y, coords = as.matrix(grid))
}


## Quadrat counts: the number of cells of each type in each cell of an
## n_grid x n_grid partition of the image bounding box, on the square-root
## scale. Unlike the KDE, neighbouring values share no cells, so the
## sampling noise of the point pattern is retained.
compute_count_matrix <- function(cells, cell_types, n_grid = N_GRID) {
  x_range <- range(cells$X); y_range <- range(cells$Y)
  xb <- seq(x_range[1], x_range[2], length.out = n_grid + 1)
  yb <- seq(y_range[1], y_range[2], length.out = n_grid + 1)
  xm <- (xb[-1] + xb[-(n_grid + 1)]) / 2; ym <- (yb[-1] + yb[-(n_grid + 1)]) / 2
  grid <- expand.grid(X = xm, Y = ym)
  Y <- matrix(0, nrow(grid), length(cell_types), dimnames = list(NULL, cell_types))
  for (ct in cell_types) {
    pts <- cells[which(cells$type == ct), c("X", "Y")]
    pts <- pts[complete.cases(pts), ]
    if (nrow(pts) < MIN_CELLS) next          # constant column, not estimable
    ix <- findInterval(pts$X, xb, all.inside = TRUE)
    iy <- findInterval(pts$Y, yb, all.inside = TRUE)
    Y[, ct] <- sqrt(as.vector(table(factor(ix, 1:n_grid), factor(iy, 1:n_grid))))
  }
  Y <- scale(Y)
  list(Y = Y, coords = as.matrix(grid))
}


## ================================================================
## STEP 2: GP-GHS PER IMAGE (parallelized over p nodes)
## ================================================================

## Persistent cluster for all images — created and closed in main pipeline
run_gp_ghs_image <- function(Y, coords, cl,
                              m     = M_BASIS,
                              nmc   = NMC,
                              burn  = BURN,
                              thin  = THIN,
                              kappa = KAPPA_THR) {

  p          <- ncol(Y)
  cell_types <- colnames(Y)

  ## Build HSGP basis once on master
  hsgp   <- build_hsgp_basis(coords, m = m, nu = NU, rho = NULL)
  phi    <- hsgp$phi
  S_diag <- hsgp$S_diag
  m2     <- hsgp$m2
  Y_sc   <- scale(Y)
  Y_sc[!is.finite(Y_sc)] <- 0   # guard: zero-variance cols -> NaN after scale
  if (DETREND) {
    ## Sensitivity analysis for an omitted response-specific smooth intercept:
    ## remove each field's ridge-regularized projection onto the HSGP basis.
    ridge <- crossprod(phi) + diag(1e-4, ncol(phi))
    trend_coef <- solve(ridge, crossprod(phi, Y_sc))
    Y_sc <- Y_sc - phi %*% trend_coef
    Y_sc <- scale(Y_sc)
    Y_sc[!is.finite(Y_sc)] <- 0
  }

  ## Package per-image data into a list and send directly to each worker.
  ## Using clusterCall avoids globalenv pollution across loop iterations:
  ## every worker always has the current image's exact data.
  img_data <- list(phi=phi, S_diag=S_diag, m2=m2, Y_sc=Y_sc,
                   p=p, cell_types=cell_types,
                   nmc=nmc, burn=burn, thin=thin, kappa=kappa)
  clusterCall(cl, function(d) {
    for (nm in names(d)) assign(nm, d[[nm]], envir=globalenv())
    NULL
  }, img_data)

  node_results <- parLapply(cl, 1:p, function(s) {
    y_s  <- as.vector(Y_sc[, s])
    X_s  <- Y_sc[, -s, drop = FALSE]
    q    <- ncol(X_s)

    ## Skip nodes with zero/near-zero response variance: constant KDE surface
    ## cannot be regressed on spatial predictors. Declare all edges absent.
    if (var(y_s) < VAR_THR || !all(is.finite(y_s))) {
      return(list(
        s           = s,
        nbr_idx     = (1:p)[-s],
        edge_active = rep(0L, q),
        kappa_mean  = rep(1.0, q),   # fully shrunk
        kappa_sd    = rep(0, q),
        edge_maps   = matrix(0, nrow(Y_sc), q),
        edge_maps_sd = matrix(0, nrow(Y_sc), q),
        field_l2_samples = matrix(0, q, floor((nmc - burn) / thin)),
        tau_mean = NA_real_, tau_sd = NA_real_,
        tau_q025 = NA_real_, tau_q975 = NA_real_,
        sigma_mean = NA_real_, tau_samples = numeric(0),
        sigma_samples = numeric(0), kappa_samples = matrix(numeric(0), q, 0)
      ))
    }

    fit  <- gp_group_horseshoe_sampler(y_s, X_s, phi, S_diag,
                                        nmc = nmc, burn = burn, thin = thin)
    summ <- summarize_node(fit, phi, S_diag, q, m2,
                            cell_types[-s], edge_threshold = kappa)
    field_l2_samples <- matrix(0, q, ncol(fit$theta_store))
    for (j in seq_len(q)) {
      idx <- ((j - 1L) * m2 + 1L):(j * m2)
      field_draws <- phi %*% fit$theta_store[idx, , drop = FALSE]
      field_l2_samples[j, ] <- sqrt(colMeans(field_draws^2))
    }

    list(
      s           = s,
      nbr_idx     = (1:p)[-s],
      edge_active = as.integer(summ$edge_summary$edge_active),
      kappa_mean  = summ$kappa_mean,
      kappa_sd    = summ$kappa_sd,
      edge_maps   = summ$edge_maps,
      edge_maps_sd = summ$edge_maps_sd,
      field_l2_samples = field_l2_samples,
      tau_mean = mean(sqrt(fit$tau_store)),
      tau_sd = sd(sqrt(fit$tau_store)),
      tau_q025 = unname(quantile(sqrt(fit$tau_store), 0.025)),
      tau_q975 = unname(quantile(sqrt(fit$tau_store), 0.975)),
      sigma_mean = mean(sqrt(fit$sigma_store)),
      tau_samples = sqrt(fit$tau_store),
      sigma_samples = sqrt(fit$sigma_store),
      kappa_samples = fit$kappa_store
    )
  })

  ## Assemble adjacency and kappa matrices
  adj   <- matrix(0, p, p, dimnames = list(cell_types, cell_types))
  kappa_mat <- matrix(NA, p, p, dimnames = list(cell_types, cell_types))

  for (r in node_results) {
    adj[r$s,   r$nbr_idx] <- r$edge_active
    kappa_mat[r$s, r$nbr_idx] <- r$kappa_mean
  }

  ## AND symmetrize adjacency
  adj_sym <- adj * t(adj); diag(adj_sym) <- 0

  ## Symmetrize kappa: take max (more conservative — edge requires BOTH nodes to agree)
  kappa_sym <- matrix(NA, p, p, dimnames = list(cell_types, cell_types))
  for (i in 1:p) for (j in 1:p) {
    if (!is.na(kappa_mat[i,j]) && !is.na(kappa_mat[j,i]))
      kappa_sym[i,j] <- max(kappa_mat[i,j], kappa_mat[j,i], na.rm = TRUE)
  }

  list(adj = adj_sym, kappa = kappa_sym,
       node_results = node_results, cell_types = cell_types)
}


## ================================================================
## STEP 3 (REVISED): DIFFERENTIAL EDGE TESTING
## ================================================================
## For each edge (s, s'), the LMM is:
##
##   kappa_{ss', i} = mu_{ss'} + beta_{ss'} * 1[group_i = DII]
##                    + b_{patient(i)} + epsilon_i
##
##   b_{patient} ~ N(0, sigma_b^2),  epsilon ~ N(0, sigma^2)
##
## beta_{ss'} is the primary inferential target: it is the
## group difference in kappa after accounting for within-patient
## correlation.  The LMM p-value and FDR are based on beta_hat.
## We also extract mu_hat (= estimated kappa for CLR) and
## mu_hat + beta_hat (= estimated kappa for DII) to use in
## summary figures, replacing the raw kappa means which do not
## account for patient clustering.
## ================================================================

run_differential_testing <- function(graph_list, meta) {
  
  imgs       <- names(graph_list)
  cell_types <- graph_list[[imgs[1]]]$cell_types
  p          <- length(cell_types)
  edges      <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  n_edges    <- nrow(edges)
  
  cat(sprintf("Testing %d edges across %d images...\n", n_edges, length(imgs)))
  
  results <- vector("list", n_edges)
  
  for (e in 1:n_edges) {
    i  <- edges[e, 1]; j <- edges[e, 2]
    ct_i <- cell_types[i]; ct_j <- cell_types[j]
    
    ## Collect per-image adj and kappa
    adj_vec <- sapply(imgs, function(img) {
      adj_mat <- graph_list[[img]]$adj
      cts     <- graph_list[[img]]$cell_types
      ii <- match(ct_i, cts); jj <- match(ct_j, cts)
      if (is.na(ii) || is.na(jj)) return(0L)
      adj_mat[ii, jj]
    })
    
    kappa_vec <- sapply(imgs, function(img) {
      kap_mat <- graph_list[[img]]$kappa
      cts     <- graph_list[[img]]$cell_types
      ii <- match(ct_i, cts); jj <- match(ct_j, cts)
      if (is.na(ii) || is.na(jj)) return(1.0)
      k <- kap_mat[ii, jj]
      if (is.na(k)) 1.0 else k
    })

    ## Directionally conservative spatial field magnitude.  Because all KDE
    ## columns are standardized, the two nodewise coefficient fields are on a
    ## common standardized scale.  Taking the smaller L2 magnitude requires
    ## both directional regressions to support a sizeable dependency.
    field_l2_vec <- sapply(imgs, function(img) {
      g <- graph_list[[img]]
      ii <- match(ct_i, g$cell_types); jj <- match(ct_j, g$cell_types)
      if (is.na(ii) || is.na(jj)) return(0)
      ni <- g$node_results[[ii]]; nj <- g$node_results[[jj]]
      si <- ni$field_l2_samples[match(jj, ni$nbr_idx), ]
      sj <- nj$field_l2_samples[match(ii, nj$nbr_idx), ]
      mean(pmin(si, sj))
    })
    
    df <- data.frame(
      img       = imgs,
      adj       = adj_vec,
      kappa     = kappa_vec,
      field_l2  = field_l2_vec,
      pathology = meta$pathology[match(imgs, meta$img)],
      patient   = meta$patient[match(imgs, meta$img)]
    )
    df$group <- factor(ifelse(df$pathology == 1, "CLR", "DII"),
                       levels = c("CLR", "DII"))   # CLR = reference
    
    ## (A) Fisher's exact test on edge prevalence
    tab <- table(df$adj, df$group)
    if (nrow(tab) < 2) {
      fisher_p <- 1.0
    } else {
      fisher_p <- fisher.test(tab)$p.value
    }
    prev_CLR <- mean(df$adj[df$group == "CLR"])
    prev_DII <- mean(df$adj[df$group == "DII"])
    
    ## (B) LMM on kappa: extract beta_hat and LMM-estimated marginals
    ## beta_hat > 0 means kappa is higher (more shrunk, weaker) in DII
    ## mu_hat     = estimated kappa for CLR (intercept)
    ## mu_hat_DII = mu_hat + beta_hat
    lmm_out <- tryCatch({
      if (requireNamespace("lme4", quietly = TRUE) &&
          requireNamespace("lmerTest", quietly = TRUE)) {
        
        suppressPackageStartupMessages({
          library(lme4)
          library(lmerTest)
        })
        eps <- 1e-4
        df$logit_kappa <- qlogis(pmin(pmax(df$kappa, eps), 1 - eps))
        fit_full <- lmer(logit_kappa ~ group + (1 | patient),
                         data = df, REML = FALSE)
        fit_null <- lmer(logit_kappa ~ 1 + (1 | patient),
                         data = df, REML = FALSE)
        lrt      <- anova(fit_null, fit_full)
        lrt_p    <- lrt[2, "Pr(>Chisq)"]
        
        ## Fixed effects: intercept = mu_CLR, groupDII = beta
        fe         <- fixef(fit_full)
        mu_hat     <- unname(fe["(Intercept)"])
        beta_hat   <- unname(fe["groupDII"])
        
        list(lmm_p    = lrt_p,
             beta_hat = beta_hat,
             mu_CLR   = plogis(mu_hat),
             mu_DII   = plogis(mu_hat + beta_hat))
        
      } else {
        ## Fallback: patient-mean aggregation then t-test
        ## No LMM-estimated marginals available; use patient-mean averages
        eps <- 1e-4
        df$logit_kappa <- qlogis(pmin(pmax(df$kappa, eps), 1 - eps))
        pt_means <- aggregate(logit_kappa ~ patient + group,
                              data = df, FUN = mean)
        tt       <- t.test(logit_kappa ~ group, data = pt_means)
        
        ## Approximate beta as difference in patient-level means
        eta_CLR <- mean(pt_means$logit_kappa[pt_means$group == "CLR"])
        eta_DII <- mean(pt_means$logit_kappa[pt_means$group == "DII"])
        
        list(lmm_p    = tt$p.value,
             beta_hat = eta_DII - eta_CLR,
             mu_CLR   = plogis(eta_CLR),
             mu_DII   = plogis(eta_DII))
      }
    }, error = function(e) {
      list(lmm_p    = NA_real_,
           beta_hat = NA_real_,
           mu_CLR   = NA_real_,
           mu_DII   = NA_real_)
    })

    ## Primary application analysis: magnitude of the estimated spatial field,
    ## rather than the local-scale-only kappa diagnostic.
    field_out <- tryCatch({
      df$log_field_l2 <- log(df$field_l2 + 1e-8)
      if (requireNamespace("lme4", quietly = TRUE) &&
          requireNamespace("lmerTest", quietly = TRUE)) {
        fit_field <- lmer(log_field_l2 ~ group + (1 | patient),
                          data = df, REML = FALSE)
        null_field <- lmer(log_field_l2 ~ 1 + (1 | patient),
                           data = df, REML = FALSE)
        fe <- fixef(fit_field)
        list(
          p = anova(null_field, fit_field)[2, "Pr(>Chisq)"],
          beta = unname(fe["groupDII"]),
          mean_CLR = exp(unname(fe["(Intercept)"])),
          mean_DII = exp(unname(fe["(Intercept)"] + fe["groupDII"]))
        )
      } else {
        pt <- aggregate(log_field_l2 ~ patient + group, df, mean)
        tt <- t.test(log_field_l2 ~ group, pt)
        a <- mean(pt$log_field_l2[pt$group == "CLR"])
        b <- mean(pt$log_field_l2[pt$group == "DII"])
        list(p = tt$p.value, beta = b - a,
             mean_CLR = exp(a), mean_DII = exp(b))
      }
    }, error = function(e) list(p = NA_real_, beta = NA_real_,
                                mean_CLR = NA_real_, mean_DII = NA_real_))
    
    ## Raw kappa means (diagnostic; NOT used for inference or plots)
    mean_kappa_CLR_raw <- mean(df$kappa[df$group == "CLR"], na.rm = TRUE)
    mean_kappa_DII_raw <- mean(df$kappa[df$group == "DII"], na.rm = TRUE)
    
    results[[e]] <- data.frame(
      ct_i              = ct_i,
      ct_j              = ct_j,
      edge              = paste(ct_i, "--", ct_j),
      prev_CLR          = round(prev_CLR, 3),
      prev_DII          = round(prev_DII, 3),
      delta_prev        = round(prev_CLR - prev_DII, 3),
      mu_CLR            = round(lmm_out$mu_CLR,   4),    # LMM: kappa for CLR
      mu_DII            = round(lmm_out$mu_DII,   4),    # LMM: kappa for DII
      beta_hat          = round(lmm_out$beta_hat, 4),    # LMM group difference
      ## raw means retained for diagnostics only
      kappa_raw_CLR     = round(mean_kappa_CLR_raw, 4),
      kappa_raw_DII     = round(mean_kappa_DII_raw, 4),
      field_l2_CLR      = round(field_out$mean_CLR, 4),
      field_l2_DII      = round(field_out$mean_DII, 4),
      field_beta_hat    = round(field_out$beta, 4),
      fisher_p          = fisher_p,
      lmm_p             = lmm_out$lmm_p,
      field_p           = field_out$p
    )
  }
  
  res <- do.call(rbind, results)
  
  ## FDR correction
  res$fisher_fdr <- p.adjust(res$fisher_p, method = "BH")
  res$lmm_fdr    <- p.adjust(res$lmm_p,   method = "BH")
  res$field_fdr  <- p.adjust(res$field_p, method = "BH")
  
  res$sig_fisher <- res$fisher_fdr < FDR_THR
  res$sig_lmm    <- res$lmm_fdr    < FDR_THR
  res$sig_field  <- res$field_fdr  < FDR_THR
  res$sig_both   <- res$sig_fisher & res$sig_lmm
  
  res[order(res$lmm_fdr), ]
}

## ================================================================
## STEP 4 (REVISED): SUMMARY PLOTS
## ================================================================
## Heatmaps now use LMM-estimated marginal means (mu_CLR, mu_DII)
## and beta_hat, rather than raw kappa averages.
## ================================================================

plot_network_summary <- function(graph_list, meta, diff_res, out_dir) {
  ## diff_res is now required: it contains mu_CLR, mu_DII, beta_hat
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible())
  library(ggplot2)
  
  imgs       <- names(graph_list)
  cell_types <- graph_list[[imgs[1]]]$cell_types
  p          <- length(cell_types)
  
  ## ---------- 4a: Edge prevalence heatmap (unchanged logic) ----------
  reindex_adj <- function(img) {
    cts_img <- graph_list[[img]]$cell_types
    idx     <- match(cell_types, cts_img)
    valid   <- which(!is.na(idx))
    m       <- matrix(0L, p, p)
    for (ii in valid)
      for (jj in valid)
        m[ii, jj] <- graph_list[[img]]$adj[idx[ii], idx[jj]]
    m
  }
  
  make_prev_mat <- function(group_imgs) {
    m <- matrix(0, p, p, dimnames = list(cell_types, cell_types))
    for (img in group_imgs) m <- m + reindex_adj(img)
    m / length(group_imgs)
  }
  
  clr_imgs <- intersect(meta$img[meta$pathology == 1], imgs)
  dii_imgs <- intersect(meta$img[meta$pathology == 2], imgs)
  
  prev_clr <- make_prev_mat(clr_imgs)
  prev_dii <- make_prev_mat(dii_imgs)
  
  mat_to_df <- function(m, group) {
    df <- expand.grid(ct_i = cell_types, ct_j = cell_types)
    df$prevalence <- as.vector(m)
    df$group <- group
    df[df$ct_i != df$ct_j, ]
  }
  
  df_prev <- rbind(mat_to_df(prev_clr, "CLR"), mat_to_df(prev_dii, "DII"))
  
  g1 <- ggplot(df_prev, aes(x = ct_j, y = ct_i, fill = prevalence)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_wrap(~group) +
    scale_fill_gradient2(low = "white", mid = "#fdbb84", high = "#d7301f",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Edge\nPrevalence") +
    labs(title = "Edge Prevalence by Pathology Group", x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y  = element_text(size = 7),
          strip.text   = element_text(face = "bold"),
          legend.position = "right")
  
  ## ---------- 4b: LMM-estimated kappa heatmaps ----------
  ## Reshape diff_res into p x p matrices for mu_CLR, mu_DII, beta_hat
  ## Only upper triangle is filled; symmetrize for display
  
  make_lmm_mat <- function(col_name) {
    m <- matrix(NA, p, p, dimnames = list(cell_types, cell_types))
    for (k in seq_len(nrow(diff_res))) {
      i <- match(diff_res$ct_i[k], cell_types)
      j <- match(diff_res$ct_j[k], cell_types)
      if (is.na(i) || is.na(j)) next
      val <- diff_res[[col_name]][k]
      m[i, j] <- val
      m[j, i] <- val   # symmetric display
    }
    m
  }
  
  mat_mu_CLR  <- make_lmm_mat("mu_CLR")
  mat_mu_DII  <- make_lmm_mat("mu_DII")
  mat_beta    <- make_lmm_mat("beta_hat")
  
  ## Fig: side-by-side mu_CLR and mu_DII (analogous to old raw-mean heatmaps)
  mat_to_df_lmm <- function(m, group) {
    df <- expand.grid(ct_i = cell_types, ct_j = cell_types)
    df$kappa_est <- as.vector(m)
    df$group <- group
    df[df$ct_i != df$ct_j, ]
  }
  
  df_mu <- rbind(mat_to_df_lmm(mat_mu_CLR, "CLR"),
                 mat_to_df_lmm(mat_mu_DII, "DII"))
  
  g2a <- ggplot(df_mu, aes(x = ct_j, y = ct_i, fill = kappa_est)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_wrap(~group) +
    scale_fill_gradient2(
      low = "#d7301f", mid = "white", high = "#2166ac",
      midpoint = 0.5, limits = c(0, 1),
      name = expression(hat(mu)[ss~"'"])
    ) +
    labs(
      title    = "LMM-Estimated Kappa by Pathology Group",
      subtitle = expression("Low "*hat(mu)*" (red) = strong interaction; high "*hat(mu)*" (blue) = shrunk"),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y   = element_text(size = 7),
          strip.text    = element_text(face = "bold"),
          legend.position = "right")
  
  ## Fig: beta_hat heatmap (estimated group difference, analogous to old delta_kappa)
  ## beta_hat > 0 => kappa higher (weaker interaction) in DII vs CLR
  ## beta_hat < 0 => kappa lower  (stronger interaction) in DII vs CLR
  dd_beta <- expand.grid(ct_i = cell_types, ct_j = cell_types)
  dd_beta$beta <- as.vector(mat_beta)
  dd_beta <- dd_beta[dd_beta$ct_i != dd_beta$ct_j, ]
  lim <- max(abs(dd_beta$beta), na.rm = TRUE)
  
  g2b <- ggplot(dd_beta, aes(x = ct_j, y = ct_i, fill = beta)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#d73027", mid = "white", high = "#2166ac",
      midpoint = 0, limits = c(-lim, lim),
      name = expression(hat(beta)[ss~"'"])
    ) +
    labs(
      title    = expression("LMM Group Effect "*hat(beta)*" (DII vs CLR)"),
      subtitle = expression("Positive (blue) = weaker in DII; negative (red) = stronger in DII"),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y   = element_text(size = 7))
  
  ## Save: edge prevalence (two panels)
  png(file.path(out_dir, "edge_prevalence_heatmap.png"),
      width = 1400, height = 600, res = 130)
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(g1, ncol = 1)
  } else print(g1)
  dev.off()
  
  ## Save: LMM-estimated kappa by group (two panels)
  png(file.path(out_dir, "kappa_lmm_heatmap.png"),
      width = 1400, height = 600, res = 130)
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(g2a, ncol = 1)
  } else print(g2a)
  dev.off()
  
  ## Save: beta_hat heatmap
  png(file.path(out_dir, "beta_hat_heatmap.png"),
      width = 700, height = 600, res = 130)
  print(g2b)
  dev.off()
  
  cat("Heatmaps saved.\n")
}

## ================================================================
## REVISED VOLCANO PLOT
## x-axis: beta_hat (LMM group effect) instead of raw delta_kappa
## ================================================================

plot_differential_results <- function(diff_res, out_dir, fdr_thr = FDR_THR) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible())
  library(ggplot2)
  
  diff_res$neg_log_fdr <- -log10(pmax(diff_res$lmm_fdr, 1e-10))
  diff_res$sig_color   <- ifelse(diff_res$sig_lmm, "Significant", "NS")
  
  g <- ggplot(diff_res, aes(x = beta_hat, y = neg_log_fdr, color = sig_color)) +
    geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed",
               color = "grey50", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    geom_point(size = 2.2, alpha = 0.8) +
    scale_color_manual(
      values = c("Significant" = "#d7301f", "NS" = "grey70"),
      name   = paste0("LMM FDR < ", fdr_thr)
    ) +
    labs(
      title    = "Differential Edge Strength: CLR vs DII",
      subtitle = expression(hat(beta) > 0~"= weaker in DII; "*hat(beta) < 0~"= stronger in DII"),
      x        = expression(hat(beta)[ss~"'"]~"(DII vs CLR, LMM)"),
      y        = expression(-log[10]*"(LMM FDR)")
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    library(ggrepel)
    g <- g + ggrepel::geom_text_repel(
      data         = diff_res[diff_res$sig_lmm, ],
      aes(label    = edge),
      size         = 2.8,
      max.overlaps = 30,
      box.padding  = 0.5,
      color        = "black",
      fontface     = "bold"
    )
  } else {
    g <- g + geom_text(
      data      = diff_res[diff_res$sig_lmm, ],
      aes(label = edge),
      size      = 2.5,
      vjust     = -0.8,
      color     = "black"
    )
  }
  
  png(file.path(out_dir, "differential_volcano.png"),
      width = 950, height = 750, res = 130)
  print(g)
  dev.off()
  cat("Volcano plot saved.\n")
}

## ================================================================
## MAIN PIPELINE
## ================================================================
cat("Loading data...\n")
dat <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

## Fixed analysis set used for every image. Sparse NK cells and unclassified
## records are excluded; absent retained types receive the documented floor.
CELL_TYPES <- c(
  "adipocytes", "B cells", "CD163+ macros", "CD4+ T cells",
  "CD68+ macros", "CD8+ T cells", "generic immune", "granulocytes",
  "memory CD4+ T", "plasma cells", "smooth muscle", "stroma", "Tregs",
  "tumor cells", "vasculature"
)
dat <- dat[dat$type %in% CELL_TYPES & !is.na(dat$type), ]

## Confirm pathology labels (1=CLR, 2=DII based on data structure)
dat$group <- ifelse(dat$pathology == 1, "CLR", "DII")

## Image-level metadata
meta <- unique(dat[, c("img", "pathology", "patient", "group")])
rownames(meta) <- NULL
imgs <- sort(unique(dat$img))
max_images <- as.integer(Sys.getenv("CRC_MAX_IMAGES", "0"))
if (max_images > 0L) imgs <- head(imgs, max_images)
meta <- meta[meta$img %in% imgs, , drop = FALSE]

cat(sprintf("Data: %d cells | %d images | %d patients\n",
            nrow(dat), length(imgs), length(unique(dat$patient))))
cat(sprintf("CLR: %d images | DII: %d images\n",
            sum(meta$pathology == 1), sum(meta$pathology == 2)))
cat(sprintf(paste0("Cell types: %d | chain: %d | seed base: %d | ",
                   "minimum cells: %d | detrend: %s\n\n"),
            length(unique(dat$type)), CHAIN_ID, SEED_BASE, MIN_CELLS, DETREND))

saveRDS(list(
  chain_id = CHAIN_ID, seed_base = SEED_BASE, min_cells = MIN_CELLS,
  var_thr = VAR_THR, detrend = DETREND, nmc = NMC, burn = BURN,
  thin = THIN, m_basis = M_BASIS, n_grid = N_GRID, feature = FEATURE,
  cell_types = CELL_TYPES
), file.path(OUT_DIR, "run_settings.rds"))


## ----------------------------------------------------------------
## 1. Build KDE feature matrices
## ----------------------------------------------------------------
cat("Computing KDE matrices for each image...\n")
kde_list <- vector("list", length(imgs))
names(kde_list) <- imgs

for (i in seq_along(imgs)) {
  if (i %% 10 == 0) cat(sprintf("  KDE: %d/%d\n", i, length(imgs)))
  cells <- dat[dat$img == imgs[i], c("X", "Y", "type")]
  kde_list[[imgs[i]]] <- if (FEATURE == "count")
    compute_count_matrix(cells, cell_types = CELL_TYPES, n_grid = N_GRID) else
    compute_kde_matrix(cells, cell_types = CELL_TYPES, n_grid = N_GRID)
}
cat("KDE computation complete.\n\n")


## ----------------------------------------------------------------
## 2. GP-GHS per image
## ----------------------------------------------------------------
cat(sprintf("Running GP-GHS on %d images (%d cores)...\n\n",
            length(imgs), N_CORES))

## Persistent cluster for all images
cl <- makeCluster(N_CORES)

## Export core GP-GHS functions once
clusterExport(cl, varlist = c(
  "gp_group_horseshoe_sampler", "summarize_node",
  "build_hsgp_basis", "matern_spectral_density",
  "M_BASIS", "NU", "NMC", "BURN", "THIN", "KAPPA_THR", "VAR_THR"
), envir = globalenv())
invisible(clusterEvalQ(cl, library(MASS)))

graph_list <- vector("list", length(imgs))
names(graph_list) <- imgs
checkpoint_path <- file.path(OUT_DIR, "graph_list_checkpoint.rds")
if (identical(Sys.getenv("CRC_RESUME", "0"), "1") &&
    file.exists(checkpoint_path)) {
  old <- readRDS(checkpoint_path)
  common <- intersect(names(old), names(graph_list))
  graph_list[common] <- old[common]
  cat(sprintf("Resuming from checkpoint with %d completed images.\n",
              sum(!sapply(graph_list, is.null))))
}

t_start <- proc.time()["elapsed"]
for (i in seq_along(imgs)) {
  img <- imgs[i]
  if (!is.null(graph_list[[img]])) {
    cat(sprintf("[%d/%d] Image: %s | already complete\n", i, length(imgs), img))
    next
  }
  clusterSetRNGStream(cl, iseed = SEED_BASE + 100000L * CHAIN_ID + i)
  cat(sprintf("[%d/%d] Image: %s | group: %s\n",
              i, length(imgs), img,
              meta$group[meta$img == img]))

  t_img <- proc.time()["elapsed"]
  graph_list[[img]] <- tryCatch(
    run_gp_ghs_image(
      Y      = kde_list[[img]]$Y,
      coords = kde_list[[img]]$coords,
      cl     = cl,
      m      = M_BASIS, nmc = NMC, burn = BURN, thin = THIN,
      kappa  = KAPPA_THR
    ),
    error = function(e) {
      message("  ERROR on image ", img, ": ", e$message)
      NULL
    }
  )
  dt <- proc.time()["elapsed"] - t_img
  n_edges_img <- if (!is.null(graph_list[[img]]))
    sum(graph_list[[img]]$adj[upper.tri(graph_list[[img]]$adj)])
  else NA
  cat(sprintf("  Done: %.1fs | Active edges: %s\n", dt,
              ifelse(is.na(n_edges_img), "ERROR", n_edges_img)))

  ## Checkpoint every 10 images
  if (i %% 10 == 0) {
    saveRDS(graph_list, checkpoint_path)
    cat(sprintf("  Checkpoint saved at image %d. Total elapsed: %.1fm\n",
                i, (proc.time()["elapsed"] - t_start)/60))
  }
}
stopCluster(cl)

## Remove failed images
failed <- sapply(graph_list, is.null)
if (any(failed)) {
  cat(sprintf("WARNING: %d images failed. Removing from analysis.\n",
              sum(failed)))
  graph_list <- graph_list[!failed]
  meta       <- meta[meta$img %in% names(graph_list), ]
}

saveRDS(graph_list, file.path(OUT_DIR, "graph_list_final.rds"))
cat(sprintf("\nGP-GHS complete. Total time: %.1f min\n\n",
            (proc.time()["elapsed"] - t_start)/60))

if (identical(Sys.getenv("CRC_FIT_ONLY", "0"), "1")) {
  tau_fit_only <- do.call(rbind, lapply(names(graph_list), function(img) {
    g <- graph_list[[img]]
    m <- meta[match(img, meta$img), ]
    do.call(rbind, lapply(g$node_results, function(x) data.frame(
      img = img, patient = m$patient, group = m$group,
      response = g$cell_types[x$s], tau_mean = x$tau_mean,
      tau_sd = x$tau_sd, tau_q025 = x$tau_q025, tau_q975 = x$tau_q975,
      sigma_mean = x$sigma_mean
    )))
  }))
  write.csv(tau_fit_only, file.path(OUT_DIR, "tau_by_image_node.csv"),
            row.names = FALSE)
  quit(save = "no", status = 0L)
}


## ----------------------------------------------------------------
## 3. Differential edge testing
## ----------------------------------------------------------------
cat("Running differential testing...\n")
diff_res <- run_differential_testing(graph_list, meta)

## Print significant edges
sig_edges <- diff_res[diff_res$sig_field | diff_res$sig_fisher |
                        diff_res$sig_lmm, ]
cat(sprintf("\nSignificant edges (Fisher FDR < %.2f OR LMM FDR < %.2f): %d\n",
            FDR_THR, FDR_THR, nrow(sig_edges)))

if (nrow(sig_edges) > 0) {
  cat("\nTop 20 by Fisher FDR:\n")
  print(head(diff_res[, c("edge","field_l2_CLR","field_l2_DII",
                           "field_beta_hat","field_fdr","lmm_fdr")], 20),
        row.names = FALSE, digits = 3)
}

tau_res <- do.call(rbind, lapply(names(graph_list), function(img) {
  g <- graph_list[[img]]
  m <- meta[match(img, meta$img), ]
  do.call(rbind, lapply(g$node_results, function(x) data.frame(
    img = img, patient = m$patient, group = m$group,
    response = g$cell_types[x$s], tau_mean = x$tau_mean,
    tau_sd = x$tau_sd, tau_q025 = x$tau_q025, tau_q975 = x$tau_q975,
    sigma_mean = x$sigma_mean
  )))
}))

write.csv(diff_res,   file.path(OUT_DIR, "differential_results.csv"),    row.names = FALSE)
write.csv(sig_edges,  file.path(OUT_DIR, "significant_edges.csv"),       row.names = FALSE)
write.csv(tau_res,    file.path(OUT_DIR, "tau_by_image_node.csv"),       row.names = FALSE)
cat("Results saved to", OUT_DIR, "\n")


## ----------------------------------------------------------------
## 4. Plots
## ----------------------------------------------------------------
cat("\nGenerating plots...\n")
plot_network_summary(graph_list, meta, diff_res, OUT_DIR)
plot_differential_results(diff_res, OUT_DIR)

cat("\nPipeline complete.\n")
