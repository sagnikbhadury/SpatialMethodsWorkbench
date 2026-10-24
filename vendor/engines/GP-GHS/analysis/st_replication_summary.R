## ================================================================
## Spatial transcriptomics: reproducibility of the gene networks across
## replicate sections and patients (all slides of Valdeolivas et al. 2024).
##
## A unit is one pathologist zone of one slide, with the network that
## st_fit.R estimated among the fixed panel of 100 genes. Every pair of
## units is compared on the gene pairs that are estimable in both.
##
## Output (st_results/replication): units.csv, unit_pairs.csv,
## category_summary.csv, replicated_edges.csv, summary.log
## Usage (R 4.4.1): Rscript st_replication_summary.R
## ================================================================

ROOT <- Sys.getenv("ST_ROOT", "st_results")
OUT  <- file.path(ROOT, "replication"); dir.create(OUT, showWarnings = FALSE)
sink(file.path(OUT, "summary.log"), split = TRUE)

harmonise_label <- function(x) {                       # as in st_prepare_panel.R
  x <- tolower(trimws(x)); x <- gsub("_", " ", x); x <- gsub("\\s+", " ", x)
  x <- sub("epitehlium", "epithelium", x)
  x <- sub("aggregregate|aggragate|aggreagate|aggreates", "aggregate", x)
  x[x == "tumor&stroma"] <- "tumor&stroma ic med to high"; x
}
index <- readRDS(file.path(ROOT, "st_prepared.rds"))
dirs <- c(setNames(ROOT, index$slide),
          sapply(list.dirs(file.path(ROOT, "slides"), recursive = FALSE), function(d) d, USE.NAMES = FALSE))
names(dirs)[-1] <- basename(dirs[-1])
units <- list()
for (sl in names(dirs)) {
  ff <- file.path(dirs[[sl]], "fit_m2", "st_fit.rds"); if (!file.exists(ff)) next
  fit <- readRDS(ff); prep <- readRDS(file.path(dirs[[sl]], "st_prepared.rds"))
  for (f in fit$fits) units[[length(units) + 1]] <- list(
    slide = sl, patient = sub("^(S[0-9]+)_.*", "\\1", sl), rep = sub(".*_(Rep[0-9])$", "\\1", sl),
    label = harmonise_label(prep$zone_tab$annotation[f$zone]), n = length(f$idx),
    est = f$est, adj = f$adj, kmax = f$kappa_max, bmean = f$beta_mean)
}
genes <- index$genes; p <- length(genes); up <- upper.tri(matrix(0, p, p)); U <- length(units)
ut <- do.call(rbind, lapply(units, function(u) data.frame(
  slide = u$slide, patient = u$patient, rep = u$rep, label = u$label, n_spots = u$n,
  n_genes = sum(u$est), n_pairs = choose(sum(u$est), 2), n_edges = sum(u$adj) / 2)))
ut$density <- ut$n_edges / ut$n_pairs
write.csv(ut, file.path(OUT, "units.csv"), row.names = FALSE)
cat("slides with a fit:", length(unique(ut$slide)), "| patients:", length(unique(ut$patient)),
    "| zones (units):", U, "| regressions:", sum(ut$n_genes), "\n")
cat("edges:", sum(ut$n_edges), "of", sum(ut$n_pairs), "pairs (", round(sum(ut$n_edges) / sum(ut$n_pairs), 4), ")\n\n")
print(ut[, c("slide", "label", "n_spots", "n_genes", "n_edges", "density")], digits = 3, row.names = FALSE)

## ---- all pairs of units ---------------------------------------------------------
rows <- list()
for (a in 1:(U - 1)) for (b in (a + 1):U) {
  A <- units[[a]]; B <- units[[b]]
  ok <- outer(A$est & B$est, A$est & B$est) & up
  x <- A$adj[ok] == 1; y <- B$adj[ok] == 1; N <- sum(ok)
  sa <- sign(A$bmean)[ok]; sb <- sign(B$bmean)[ok]
  same_pat <- A$patient == B$patient; same_slide <- A$slide == B$slide; same_lab <- A$label == B$label
  cat_ <- if (same_slide) "same slide, different zone" else
          if (same_pat && same_lab) "replicate sections, same zone" else
          if (same_pat) "replicate sections, different zone" else
          if (same_lab) "different patients, same zone" else "different patients, different zone"
  rows[[length(rows) + 1]] <- data.frame(
    a = a, b = b, slide_a = A$slide, slide_b = B$slide, label_a = A$label, label_b = B$label, category = cat_,
    n_pairs = N, edges_a = sum(x), edges_b = sum(y), both = sum(x & y),
    jaccard = sum(x & y) / max(sum(x | y), 1),
    overlap = sum(x & y) / max(min(sum(x), sum(y)), 1),              # share of the smaller network found in the other
    expected = sum(x) * sum(y) / N,                                   # edges in common if the two networks were unrelated
    fold = (sum(x & y) / N) / max((sum(x) / N) * (sum(y) / N), 1e-12),
    kappa_cor = suppressWarnings(cor(A$kmax[ok], B$kmax[ok], method = "spearman")),
    sign_agree = if (sum(x & y)) mean(sa[x & y] == sb[x & y]) else NA)
}
up_df <- do.call(rbind, rows)
write.csv(up_df, file.path(OUT, "unit_pairs.csv"), row.names = FALSE)
lev <- c("replicate sections, same zone", "replicate sections, different zone", "same slide, different zone",
         "different patients, same zone", "different patients, different zone")
cs <- do.call(rbind, lapply(lev, function(k) { d <- up_df[up_df$category == k, ]
  data.frame(category = k, n_comparisons = nrow(d), edges_in_common = sum(d$both), expected_by_chance = round(sum(d$expected), 1),
             fold_pooled = round(sum(d$both) / sum(d$expected), 1),
             jaccard_median = round(median(d$jaccard), 3), jaccard_q1 = round(quantile(d$jaccard, .25), 3),
             jaccard_q3 = round(quantile(d$jaccard, .75), 3),
             overlap_median = round(median(d$overlap), 3), overlap_pooled = round(sum(d$both) / sum(pmin(d$edges_a, d$edges_b)), 3),
             kappa_cor_median = round(median(d$kappa_cor, na.rm = TRUE), 3),
             sign_agree = round(weighted.mean(d$sign_agree, d$both, na.rm = TRUE), 3)) }))
write.csv(cs, file.path(OUT, "category_summary.csv"), row.names = FALSE)
cat("\nAgreement between the networks of two zones, by type of comparison:\n"); print(cs, row.names = FALSE)
cat("\nReplicate sections, same zone, one row per zone:\n")
rp <- up_df[up_df$category == lev[1], ]
print(rp[order(rp$slide_a), c("slide_a", "label_a", "edges_a", "edges_b", "both", "expected", "jaccard", "overlap", "kappa_cor", "sign_agree")],
      digits = 2, row.names = FALSE)

## ---- replicated edges: selected in both replicate sections of a patient -----------
rep_rows <- list()
for (k in seq_len(nrow(rp))) {
  A <- units[[rp$a[k]]]; B <- units[[rp$b[k]]]
  both <- which(A$adj == 1 & B$adj == 1 & up, arr.ind = TRUE)
  if (nrow(both)) rep_rows[[length(rep_rows) + 1]] <- data.frame(
    patient = A$patient, label = A$label, gene_i = genes[both[, 1]], gene_j = genes[both[, 2]],
    sign = ifelse(A$bmean[both] > 0 & t(A$bmean)[both] > 0, "positive",
                  ifelse(A$bmean[both] < 0 & t(A$bmean)[both] < 0, "negative", "mixed")))
}
re <- do.call(rbind, rep_rows); re$pair <- paste(re$gene_i, "--", re$gene_j)
write.csv(re, file.path(OUT, "replicated_edges.csv"), row.names = FALSE)
cat("\nreplicated edges (selected in both replicate sections of a patient, same zone):", nrow(re),
    "| distinct gene pairs:", length(unique(re$pair)), "\n")
cat("replicated edges per patient and zone:\n"); print(table(re$patient, re$label))
pp <- aggregate(patient ~ pair, re, function(z) length(unique(z))); names(pp)[2] <- "n_patients"
pz <- aggregate(label ~ pair, re, function(z) paste(sort(unique(z)), collapse = "; "))
pt <- aggregate(patient ~ pair, re, function(z) paste(sort(unique(z)), collapse = " "))
ps <- aggregate(sign ~ pair, re, function(z) paste(sort(unique(z)), collapse = "/"))
rec <- Reduce(function(x, y) merge(x, y, by = "pair"), list(pp, pt, pz, ps)); names(rec)[3] <- "patients"
rec <- rec[order(-rec$n_patients, rec$pair), ]
write.csv(rec, file.path(OUT, "replicated_edges_by_pair.csv"), row.names = FALSE)
cat("\ngene pairs replicated in k patients (k = 1..):", paste(table(rec$n_patients), collapse = " "), "\n")
cat("pairs replicated in at least three patients:\n")
print(rec[rec$n_patients >= 3, c("pair", "n_patients", "patients", "sign")], row.names = FALSE)

## zone-specific replicated edges: replicated in a zone of a patient and selected in no other zone of that patient
zs <- do.call(rbind, lapply(split(re, re$patient), function(d) {
  us <- which(sapply(units, function(u) u$patient == d$patient[1]))
  do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    gi <- match(d$gene_i[i], genes); gj <- match(d$gene_j[i], genes)
    other <- sum(sapply(units[us], function(u) u$label != d$label[i] && u$adj[gi, gj] == 1))
    data.frame(d[i, ], selected_in_other_zones = other) })) }))
write.csv(zs, file.path(OUT, "replicated_edges_zone_specificity.csv"), row.names = FALSE)
cat("\nreplicated edges that are selected in no other zone of the same patient:", sum(zs$selected_in_other_zones == 0),
    "of", nrow(zs), "\n")
print(table(patient = zs$patient, zone_specific = zs$selected_in_other_zones == 0))
saveRDS(list(units = units, ut = ut, up_df = up_df, re = re, rec = rec, zs = zs), file.path(OUT, "replication.rds"))
sink()
