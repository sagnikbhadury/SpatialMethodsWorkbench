## ----setup, message=FALSE-----------------------------------------------------
library(ISPAT3D)
image <- ispat3d_example_image(n_per_section = 30, n_sections = 3,
                               bandwidth = 0.12, seed = 2026)
table(image$zones, image$sections)

## ----sections, fig.width=7, fig.height=2.8, fig.cap="Synthetic planar cell locations, colored by relative tumor-density zone."----
old <- par(mfrow = c(1, 3), mar = c(3, 3, 2, 1))
zone_color <- c(Low = "#4C72B0", High = "#C44E52")
for (s in sort(unique(image$sections))) {
  at <- image$sections == s
  plot(image$coords[at, 1:2], xlim = c(0, 1), ylim = c(0, 1),
       xlab = "x", ylab = "y", main = paste("Section", s), pch = 19,
       cex = 0.6, col = zone_color[as.character(image$zones[at])])
}
par(old)

## ----roster-------------------------------------------------------------------
selected <- ispat3d_sample(image$coords, image$zones, image$sections,
                           budget = 36L, budget_kind = "per_zone",
                           seed = 2027L)
rows <- unlist(selected, use.names = FALSE)
table(image$zones[rows], image$sections[rows])
stopifnot(length(rows) == 72L)
Y <- image$Y[rows, , drop = FALSE]
xyz <- image$coords[rows, , drop = FALSE]
zone <- image$zones[rows]
section <- image$sections[rows]

## ----fit2d--------------------------------------------------------------------
common <- list(
  Y = Y, coords = xyz, zones = zone, sections = section,
  rank = 2L, anchor_min = 12L, anchor_max = 24L,
  neighbors = 5L, gp_maxit = 3L, factor_maxit = 60L,
  threads = 1L, return_residuals = TRUE
)
fit2d <- do.call(ispat3d_fit_2d, common)
table(fit2d$gp_log$status)
fit2d$counts

## ----network-values-----------------------------------------------------------
round(fit2d$full[["High"]], 3)
round(fit2d$partial[["High"]], 3)
stopifnot(isTRUE(all.equal(
  ispat3d_partial_correlation(fit2d$full[["High"]]),
  fit2d$partial[["High"]]
)))
ispat3d_edge_table(fit2d, zone = "High", threshold = 0.02)

## ----planar-networks, fig.width=7, fig.height=3.5, fig.cap="Matched planar conditional-density networks in the two simulated zones."----
ispat3d_plot_zones(fit2d, columns = 2, threshold = 0.02,
                   label_cex = 0.8)

## ----fit3d--------------------------------------------------------------------
fit3d <- do.call(ispat3d_fit, common)
stopifnot(identical(fit3d$counts, fit2d$counts))
delta_high <- fit3d$partial[["High"]] - fit2d$partial[["High"]]
round(delta_high, 3)

## ----comparison, fig.width=7, fig.height=3.5, fig.cap="The same selected cells analyzed with pooled 3D and section-wise planar GP adjustment."----
old <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
ispat3d_plot_network(fit3d, zone = "High", threshold = 0.02,
                     main = "Pooled 3D", label_cex = 0.8)
ispat3d_plot_network(fit2d, zone = "High", threshold = 0.02,
                     main = "Section-wise 2D", label_cex = 0.8)
par(old)

## ----real-input, eval=FALSE---------------------------------------------------
# # dat contains registered X,Y,Z, section, zone and cell-type KDE columns.
# coords <- as.matrix(dat[, c("X", "Y", "Z")])
# kde <- as.matrix(dat[, c("kde_Tumor", "kde_T_cell", "kde_Macrophage")])
# colnames(kde) <- c("Tumor", "T_cell", "Macrophage")
# Y <- log1p(1e9 * kde)
# selected <- ispat3d_sample(coords, dat$zone, dat$section,
#                            budget = 1000L, budget_kind = "per_zone")
# rows <- unlist(selected, use.names = FALSE)
# args <- list(Y = Y[rows, , drop = FALSE],
#              coords = coords[rows, , drop = FALSE],
#              zones = dat$zone[rows], sections = dat$section[rows],
#              rank = 2L)
# planar <- do.call(ispat3d_fit_2d, args)
# volumetric <- do.call(ispat3d_fit, args)
# ispat3d_plot_zones(planar, threshold = 0.05)

