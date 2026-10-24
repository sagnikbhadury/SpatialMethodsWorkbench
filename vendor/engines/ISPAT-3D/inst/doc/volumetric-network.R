## ----setup, message=FALSE-----------------------------------------------------
library(ISPAT3D)
set.seed(2026)
image <- ispat3d_example_image(n_per_section = 30, n_sections = 3,
                               bandwidth = 0.12, seed = 2026)
stopifnot(ncol(image$coords) == 3, nrow(image$Y) == 90)
table(image$zones, image$sections)

## ----map, fig.width=7, fig.height=2.7, fig.cap="Three registered synthetic sections. Color denotes annotated source cell type."----
old <- par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))
palette <- c(Tumor = "#C44E52", T_cell = "#4C72B0",
             Macrophage = "#55A868")
for (s in sort(unique(image$sections))) {
  at <- image$sections == s
  plot(image$coords[at, 1:2], xlim = c(0, 1), ylim = c(0, 1),
       xlab = "x", ylab = "y", main = paste("Section", s),
       pch = 19, cex = 0.6, col = palette[image$cell_type[at]])
}
par(old)

## ----sample-------------------------------------------------------------------
selected <- ispat3d_sample(image$coords, image$zones, image$sections,
                           budget = 36L, budget_kind = "per_zone",
                           seed = 2027L)
lengths(selected)
rows <- unlist(selected, use.names = FALSE)
stopifnot(length(rows) == 72L)

## ----fit3d--------------------------------------------------------------------
fit3d <- ispat3d_fit(
  Y = image$Y[rows, , drop = FALSE],
  coords = image$coords[rows, , drop = FALSE],
  zones = image$zones[rows],
  sections = image$sections[rows],
  rank = 2L,
  anchor_min = 12L, anchor_max = 24L,
  neighbors = 5L, gp_maxit = 3L, factor_maxit = 60L,
  threads = 1L, return_residuals = TRUE
)
fit3d$counts
table(fit3d$gp_log$status)

## ----inspect-network----------------------------------------------------------
names(fit3d$full)
round(fit3d$shared, 3)
round(fit3d$partial[["High"]], 3)
recomputed <- ispat3d_partial_correlation(fit3d$full[["High"]])
stopifnot(isTRUE(all.equal(recomputed, fit3d$partial[["High"]])))
ispat3d_edge_table(fit3d, zone = "High", threshold = 0.02)

## ----networks, fig.width=7, fig.height=3.5, fig.cap="Conditional cell-density association networks in two relative tumor-density zones."----
ispat3d_plot_zones(fit3d, threshold = 0.02, columns = 2,
                   label_cex = 0.8)

## ----real-input, eval=FALSE---------------------------------------------------
# # dat has X, Y, Z, section, zone, kde_Tumor, kde_T_cell, kde_Macrophage.
# coords <- as.matrix(dat[, c("X", "Y", "Z")])
# kde <- as.matrix(dat[, c("kde_Tumor", "kde_T_cell", "kde_Macrophage")])
# colnames(kde) <- c("Tumor", "T_cell", "Macrophage")
# Y <- log1p(1e9 * kde)
# selected <- ispat3d_sample(coords, dat$zone, dat$section,
#                            budget = 1000L, budget_kind = "per_zone")
# rows <- unlist(selected, use.names = FALSE)
# fit <- ispat3d_fit(Y[rows, , drop = FALSE], coords[rows, , drop = FALSE],
#                    dat$zone[rows], dat$section[rows], rank = 2L)
# ispat3d_plot_zones(fit, threshold = 0.05)

