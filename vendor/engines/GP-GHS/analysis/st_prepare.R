## ================================================================
## Spatial transcriptomics application, step 1: data preparation
##
## One 10x Visium slide of colorectal cancer (Valdeolivas et al. 2024,
## npj Precision Oncology 8:10; h5ad from the cellxgene collection).
##   1. spots under tissue, quality filter, log-normalised expression
##   2. zones: the pathologist's annotation of each spot, as supplied
##      with the data (spots annotated "exclude" are dropped)
##   3. spatially variable genes: Moran's I on the spatial 6-nearest-
##      neighbour graph, Benjamini-Hochberg over genes
##   4. differential expression of the spatially variable genes across
##      the zones (Kruskal-Wallis), ranked by p-value; top 100
##
## Output (in ST_OUT_DIR): st_prepared.rds, gene_table.csv,
## zone_table.csv, prepare.log
##
## Usage (R 4.4.1): Rscript st_prepare.R
## ================================================================

suppressPackageStartupMessages({ library(hdf5r); library(Matrix); library(FNN) })

DATA_DIR <- Sys.getenv("ST_DATA_DIR", file.path(
  "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury",
  "Data for Analysis/CRC_CMS_Visium_cellxgene"))
SLIDE     <- Sys.getenv("ST_SLIDE", "S5_Rec_A121573_Rep2")
OUT_DIR   <- Sys.getenv("ST_OUT_DIR", "st_results")
N_GENES   <- 100L
MIN_UMI   <- 500      # spots with fewer counts are dropped
MIN_DET   <- 0.10     # a candidate gene is detected in at least this share of spots
K_SPACE   <- 6L       # neighbours in the spatial graph (Moran's I)
dir.create(OUT_DIR, showWarnings = FALSE)
log_con <- file(file.path(OUT_DIR, "prepare.log"), "w")
say <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = log_con) }

## display names of the pathologist's labels on this slide
ZONE_NAMES <- c(
  "epithelium&submucosa"        = "Epithelium and submucosa",
  "stroma_fibroblastic_IC high" = "Fibroblastic stroma, IC high",
  "tumor"                       = "Tumor",
  "tumor&stroma_IC med to high" = "Tumor and stroma, IC medium to high",
  "non neo epithelium"          = "Non-neoplastic epithelium",
  "submucosa"                   = "Submucosa")

## ---- 1. read the slide -------------------------------------------------------
f <- H5File$new(file.path(DATA_DIR, paste0(SLIDE, ".h5ad")), "r")
shape <- h5attr(f[["X"]], "shape")                       # spots x genes, CSR
counts <- sparseMatrix(i = f[["X/indices"]][] + 1L, p = f[["X/indptr"]][],
                       x = f[["X/data"]][], dims = c(shape[2], shape[1]))  # genes x spots
## categorical column; code -1 (spots off tissue) is missing
cat_col <- function(g) { cc <- g[["codes"]][]; lv <- g[["categories"]][]; lv[ifelse(cc < 0, NA, cc + 1L)] }
symbol  <- as.character(cat_col(f[["var/feature_name"]]))
ann_col <- grep("athologist", names(f[["obs"]]), value = TRUE)[1]   # name differs between slides
obs <- data.frame(
  barcode   = f[["obs/_index"]][],
  in_tissue = f[["obs/in_tissue"]][],
  array_row = f[["obs/array_row"]][], array_col = f[["obs/array_col"]][],
  annotation = as.character(cat_col(f[["obs"]][[ann_col]])),
  stringsAsFactors = FALSE)
xy <- t(f[["obsm/spatial"]][, ]); obs$x <- xy[, 1]; obs$y <- xy[, 2]
f$close_all()
rownames(counts) <- make.unique(symbol); colnames(counts) <- obs$barcode
say("slide: ", SLIDE, " | spots on the array: ", ncol(counts), " | genes: ", nrow(counts))

## ---- 2. spot filter and zones ---------------------------------------------------
umi <- colSums(counts)
keep <- obs$in_tissue == 1 & obs$annotation %in% names(ZONE_NAMES) & umi >= MIN_UMI
say("spots under tissue: ", sum(obs$in_tissue == 1),
    " | annotated 'exclude': ", sum(obs$in_tissue == 1 & obs$annotation %in% "exclude"),
    " | fewer than ", MIN_UMI, " counts: ", sum(obs$in_tissue == 1 & umi < MIN_UMI),
    " | kept: ", sum(keep))
stopifnot(all(obs$annotation[obs$in_tissue == 1] %in% c(names(ZONE_NAMES), "exclude")))
counts <- counts[, keep]; obs <- obs[keep, ]; umi <- umi[keep]
say("counts per spot (quartiles): ", paste(round(quantile(umi, c(.25, .5, .75))), collapse = " "),
    " | genes detected per spot (median): ", median(colSums(counts > 0)))
## zones ordered by size, largest first
ord <- names(sort(table(obs$annotation), decreasing = TRUE))
obs$zone <- match(obs$annotation, ord)
zone_tab <- data.frame(zone = seq_along(ord), annotation = ord, name = unname(ZONE_NAMES[ord]),
                       n_spots = as.numeric(table(obs$zone)),
                       median_umi = as.numeric(tapply(umi, obs$zone, median)))
nn_sp <- get.knn(as.matrix(obs[, c("x", "y")]), k = K_SPACE)$nn.index
same <- matrix(obs$zone[as.vector(nn_sp)] == rep(obs$zone, K_SPACE), ncol = K_SPACE)
zone_tab$same_zone_neighbours <- round(as.numeric(tapply(rowMeans(same), obs$zone, mean)), 3)
write.csv(zone_tab, file.path(OUT_DIR, "zone_table.csv"), row.names = FALSE)
say("zones (pathologist annotation):"); print(zone_tab, row.names = FALSE)
capture.output(print(zone_tab, row.names = FALSE), file = log_con)
say("share of spatial neighbours in the same zone: ", round(mean(same), 3),
    " | expected under random labels: ", round(sum((table(obs$zone) / nrow(obs))^2), 3))

## ---- 3. spatially variable genes: Moran's I ----------------------------------------
lib <- sweep(counts, 2, umi, "/") * 1e4
det <- rowMeans(counts > 0)
drop_tech <- grepl("^MT-|^RP[LS]|^MRP[LS]", rownames(counts))
cand <- det >= MIN_DET & !drop_tech
say("genes detected in >= ", 100 * MIN_DET, "% of spots: ", sum(det >= MIN_DET),
    " | mitochondrial or ribosomal among them (removed): ", sum(det >= MIN_DET & drop_tech),
    " | candidates: ", sum(cand))
E <- t(as.matrix(log1p(lib[cand, ])))                    # spots x genes, log1p(counts per 10,000)
n <- nrow(E)
W <- sparseMatrix(i = rep(seq_len(n), K_SPACE), j = as.vector(nn_sp), x = 1, dims = c(n, n))
W <- (W + t(W)) / 2                                       # symmetric weights
S0 <- sum(W); S1 <- 2 * sum(W^2); S2 <- sum((2 * rowSums(W))^2)
Z <- scale(E)
moran <- (n / S0) * colSums(Z * as.matrix(W %*% Z)) / colSums(Z^2)
EI <- -1 / (n - 1)
VI <- (n^2 * S1 - n * S2 + 3 * S0^2) / ((n^2 - 1) * S0^2) - EI^2   # normality assumption
moran_z <- (moran - EI) / sqrt(VI)
moran_p <- pnorm(moran_z, lower.tail = FALSE)
moran_fdr <- p.adjust(moran_p, "BH")
svg <- moran_fdr < 0.05
say("Moran's I (quartiles): ", paste(round(quantile(moran, c(.25, .5, .75)), 3), collapse = " "),
    " | spatially variable genes (BH < 0.05): ", sum(svg), " of ", length(svg))

## ---- 4. differential expression across zones: Kruskal-Wallis ---------------------
g <- obs$zone; ng <- as.numeric(table(g)); K <- length(ng)
kw <- apply(E[, svg, drop = FALSE], 2, function(x) {
  r <- rank(x); tt <- table(x)
  H <- 12 / (n * (n + 1)) * sum(tapply(r, g, sum)^2 / ng) - 3 * (n + 1)
  H / (1 - sum(tt^3 - tt) / (n^3 - n))
})
kw_logp <- pchisq(kw, K - 1, lower.tail = FALSE, log.p = TRUE) / log(10)
gene_tab <- data.frame(gene = colnames(E), detected = det[cand], moran_I = moran,
                       moran_z = moran_z, moran_fdr = moran_fdr, svg = svg,
                       kw_H = NA_real_, kw_log10p = NA_real_)
gene_tab$kw_H[svg] <- kw; gene_tab$kw_log10p[svg] <- kw_logp
gene_tab <- gene_tab[order(gene_tab$kw_log10p, -gene_tab$kw_H), ]
gene_tab$rank <- ifelse(is.na(gene_tab$kw_H), NA, seq_len(nrow(gene_tab)))
gene_tab$selected <- !is.na(gene_tab$rank) & gene_tab$rank <= N_GENES
## zone in which the mean expression is largest
top <- gene_tab$gene[gene_tab$selected]
mean_by <- sapply(seq_len(K), function(k) colMeans(E[g == k, top, drop = FALSE]))
colnames(mean_by) <- paste0("mean_zone", seq_len(K))
gene_tab$top_zone <- NA_integer_
gene_tab$top_zone[gene_tab$selected] <- max.col(mean_by, ties.method = "first")
write.csv(gene_tab, file.path(OUT_DIR, "gene_table.csv"), row.names = FALSE)
write.csv(data.frame(gene = top, round(mean_by, 3)), file.path(OUT_DIR, "selected_gene_zone_means.csv"),
          row.names = FALSE)
say("Kruskal-Wallis over ", sum(svg), " spatially variable genes, ", K - 1, " df | BH < 0.05: ",
    sum(p.adjust(exp(kw_logp * log(10)), "BH") < 0.05),
    " | H of gene ranked ", N_GENES, ": ", round(gene_tab$kw_H[N_GENES], 1),
    " (log10 p = ", round(gene_tab$kw_log10p[N_GENES], 1), ")")
say("selected genes: Moran's I range ", paste(round(range(gene_tab$moran_I[gene_tab$selected]), 3), collapse = " to "),
    " | detection range ", paste(round(range(gene_tab$detected[gene_tab$selected]), 2), collapse = " to "))
say("selected genes by zone of highest mean: ", paste(table(factor(gene_tab$top_zone, 1:K)), collapse = " "))
say("top ", N_GENES, " genes: ", paste(top, collapse = ", "))

saveRDS(list(slide = SLIDE, obs = obs, expr = E[, top], genes = top, gene_tab = gene_tab,
             zone_tab = zone_tab,
             settings = list(min_umi = MIN_UMI, min_det = MIN_DET, k_space = K_SPACE)),
        file.path(OUT_DIR, "st_prepared.rds"))
close(log_con)
