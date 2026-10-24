## ================================================================
## Spatial transcriptomics, replication across slides: data preparation
##
## Prepares any slide of Valdeolivas et al. (2024) for st_fit.R with the
## gene panel fixed in advance: the 100 genes selected on the index slide
## by st_prepare.R. No gene selection is done on the new slide, so its
## networks are out-of-sample with respect to the choice of genes.
##
## Zones are the pathologist's labels of the slide. The spelling of the
## labels differs between slides; harmonise_label() maps the variants to
## one spelling so that the same zone can be matched across slides. Zones
## with fewer than ST_MIN_ZONE spots are not analysed.
##
## Output: <ST_ROOT>/slides/<slide>/st_prepared.rds, prepare.log
## Usage (R 4.4.1): ST_SLIDE=S5_Rec_A121573_Rep1 Rscript st_prepare_panel.R
## ================================================================

suppressPackageStartupMessages({ library(hdf5r); library(Matrix) })

DATA_DIR <- Sys.getenv("ST_DATA_DIR", file.path(
  "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury",
  "Data for Analysis/CRC_CMS_Visium_cellxgene"))
ROOT     <- Sys.getenv("ST_ROOT", "st_results")
SLIDE    <- Sys.getenv("ST_SLIDE")
MIN_UMI  <- 500
MIN_ZONE <- as.integer(Sys.getenv("ST_MIN_ZONE", "150"))
stopifnot(nzchar(SLIDE))
OUT_DIR <- file.path(ROOT, "slides", SLIDE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
log_con <- file(file.path(OUT_DIR, "prepare.log"), "w")
say <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = log_con) }

## One spelling per pathologist label (variants found in the 14 files).
harmonise_label <- function(x) {
  x <- tolower(trimws(x)); x <- gsub("_", " ", x); x <- gsub("\\s+", " ", x)
  x <- sub("epitehlium", "epithelium", x)
  x <- sub("aggregregate|aggragate|aggreagate|aggreates", "aggregate", x)
  x[x == "tumor&stroma"] <- "tumor&stroma ic med to high"     # slide S3 Rep2 omits the IC grade
  x
}

panel <- readRDS(file.path(ROOT, "st_prepared.rds"))$genes

f <- H5File$new(file.path(DATA_DIR, paste0(SLIDE, ".h5ad")), "r")
shape <- h5attr(f[["X"]], "shape")
counts <- sparseMatrix(i = f[["X/indices"]][] + 1L, p = f[["X/indptr"]][],
                       x = f[["X/data"]][], dims = c(shape[2], shape[1]))
cat_col <- function(g) { cc <- g[["codes"]][]; lv <- g[["categories"]][]; lv[ifelse(cc < 0, NA, cc + 1L)] }
symbol  <- as.character(cat_col(f[["var/feature_name"]]))
ann_col <- grep("athologist", names(f[["obs"]]), value = TRUE)[1]
obs <- data.frame(barcode = f[["obs/_index"]][], in_tissue = f[["obs/in_tissue"]][],
                  annotation_raw = as.character(cat_col(f[["obs"]][[ann_col]])),
                  stringsAsFactors = FALSE)
xy <- t(f[["obsm/spatial"]][, ]); obs$x <- xy[, 1]; obs$y <- xy[, 2]
f$close_all()
rownames(counts) <- make.unique(symbol)
obs$annotation <- harmonise_label(obs$annotation_raw)
umi <- colSums(counts)
ok <- obs$in_tissue == 1 & !is.na(obs$annotation) & obs$annotation != "exclude" & umi >= MIN_UMI
tab <- sort(table(obs$annotation[ok]), decreasing = TRUE)
zones <- names(tab)[tab >= MIN_ZONE]
say("slide: ", SLIDE, " | spots under tissue: ", sum(obs$in_tissue == 1), " | usable: ", sum(ok),
    " | labels: ", length(tab), " | zones with >= ", MIN_ZONE, " spots: ", length(zones))
for (z in names(tab)) say(sprintf("   %4d  %s%s", tab[[z]], z, if (tab[[z]] >= MIN_ZONE) "" else "   (not analysed)"))
keep <- ok & obs$annotation %in% zones
if (length(zones) == 0) { say("no zone is large enough; slide not analysed"); close(log_con); quit(save = "no") }
stopifnot(all(panel %in% rownames(counts)))
expr <- t(as.matrix(log1p(sweep(counts[panel, keep, drop = FALSE], 2, umi[keep], "/") * 1e4)))
obs <- obs[keep, ]; obs$zone <- match(obs$annotation, zones)
zone_tab <- data.frame(zone = seq_along(zones), annotation = zones, name = zones,
                       n_spots = as.numeric(table(obs$zone)))
say("analysed spots: ", nrow(obs), " | median counts per spot: ", round(median(umi[keep])))
saveRDS(list(slide = SLIDE, obs = obs, expr = expr, genes = panel, zone_tab = zone_tab,
             settings = list(min_umi = MIN_UMI, min_zone = MIN_ZONE, panel_from = "index slide")),
        file.path(OUT_DIR, "st_prepared.rds"))
close(log_con)
