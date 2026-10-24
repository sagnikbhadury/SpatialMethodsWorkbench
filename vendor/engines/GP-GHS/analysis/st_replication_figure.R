## Figure for the replication of the gene networks across sections and patients.
## (a) Agreement between the networks of two zones, by type of comparison.
## (b) The gene pairs selected in both replicate sections in at least two patients.
##
## Usage (R 4.4.1): Rscript st_replication_figure.R   (after st_replication_summary.R)

suppressPackageStartupMessages(library(igraph))
ROOT <- Sys.getenv("ST_ROOT", "st_results")
r <- readRDS(file.path(ROOT, "replication", "replication.rds"))
d <- r$up_df; rec <- r$rec
lev <- c("replicate sections, same zone", "replicate sections, different zone", "same slide, different zone",
         "different patients, same zone", "different patients, different zone")
lab <- c("Replicate\nsections,\nsame zone", "Replicate\nsections,\ndifferent zone", "Same\nsection,\ndifferent zone",
         "Different\npatients,\nsame zone", "Different\npatients,\ndifferent zone")
cols <- c("#2a78d6", "#1baf7a", "#eda100", "#eb6834", "#4a3aa7")
d$k <- match(d$category, lev)

cairo_pdf(file.path("manuscript", "Image_s", "st_replication.pdf"), width = 11, height = 5.4)
layout(matrix(1:2, 1), widths = c(1, 1.15))
## (a) share of the smaller network that is also selected in the other zone
par(mar = c(6.2, 4.4, 2.4, 0.8), mgp = c(2.6, 0.7, 0), las = 1)
plot(NA, xlim = c(0.5, 5.5), ylim = c(0, max(d$overlap) * 1.05), xaxt = "n", xlab = "",
     ylab = "Share of the smaller network also selected in the other", bty = "l")
set.seed(3)
for (k in 1:5) {
  v <- d$overlap[d$k == k]
  points(k + runif(length(v), -0.22, 0.22), v, pch = 16, cex = 0.55, col = adjustcolor(cols[k], 0.45))
  q <- quantile(v, c(.25, .5, .75))
  segments(k - 0.3, q[2], k + 0.3, q[2], lwd = 2.4); segments(k, q[1], k, q[3], lwd = 1.2)
  text(k, par("usr")[4], sprintf("n = %d", length(v)), cex = 0.75, pos = 1, xpd = NA)
}
ex <- sum(d$expected) / sum(pmin(d$edges_a, d$edges_b))
abline(h = ex, lty = 2, col = "grey40"); text(0.5, ex - 0.022, "expected if unrelated", pos = 4, cex = 0.75, col = "grey30", offset = 0.1)
axis(1, at = 1:5, labels = FALSE); text(1:5, par("usr")[3] - 0.04 * diff(par("usr")[3:4]), lab, srt = 0, xpd = NA, cex = 0.78, adj = c(0.5, 1))
title("(a) Agreement between the networks of two zones", cex.main = 1, adj = 0)

## (b) pairs replicated in at least two patients
e <- rec[rec$n_patients >= 2, ]
gi <- sub(" -- .*", "", e$pair); gj <- sub(".* -- ", "", e$pair)
g <- graph_from_data_frame(data.frame(gi, gj, n = e$n_patients, sign = e$sign), directed = FALSE)
set.seed(11); L <- layout_with_fr(g, niter = 3000)
## spread the crowded centre: keep each gene's angle, replace its radius by the rank of the radius
L <- sweep(L, 2, apply(L, 2, median)); rad <- sqrt(rank(sqrt(rowSums(L^2))) / nrow(L)); ang <- atan2(L[, 2], L[, 1])
L <- cbind(rad * cos(ang), rad * sin(ang))
ecol <- ifelse(E(g)$sign == "positive", "#e34948", ifelse(E(g)$sign == "negative", "#2a78d6", "grey50"))
par(mar = c(2.6, 0.3, 2.4, 0.3))
plot(g, layout = L, vertex.size = 4.5, vertex.color = "grey92", vertex.frame.color = "grey45",
     vertex.label.cex = 0.68, vertex.label.font = 2, vertex.label.color = "black", vertex.label.family = "sans",
     vertex.label.dist = 0.75, vertex.label.degree = -pi / 2,
     edge.width = 0.6 + 1.1 * (E(g)$n - 1), edge.color = adjustcolor(ecol, 0.8),
     edge.lty = ifelse(E(g)$sign == "positive", 1, ifelse(E(g)$sign == "negative", 5, 3)))
title(sprintf("(b) The %d gene pairs replicated in at least two patients", nrow(e)), cex.main = 1, adj = 0)
legend("bottom", horiz = TRUE, bty = "n", cex = 0.8, inset = -0.09, xpd = NA,
       legend = c("Positive", "Negative", "Sign differs", "2 patients", "6 patients"),
       col = c("#e34948", "#2a78d6", "grey50", "grey30", "grey30"), lty = c(1, 5, 3, 1, 1), lwd = c(1.7, 1.7, 1.7, 1.7, 6.1))
dev.off()
cat("pairs in panel b:", nrow(e), "| genes:", vcount(g), "| sign: positive", sum(e$sign == "positive"),
    "negative", sum(e$sign == "negative"), "other", sum(!e$sign %in% c("positive", "negative")), "\n")
cat("overlap expected if unrelated:", round(ex, 3), "\n")
## numbers quoted in the text
u <- r$ut; rp <- d[d$k == 1, ]
cat("replicate same zone: overlap by size of the smaller zone network; cor(overlap, log min edges):",
    round(cor(rp$overlap, log(pmin(rp$edges_a, rp$edges_b))), 2), "\n")
cat("zones:", nrow(u), "| spots:", sum(u$n_spots), "| edges per zone (quartiles):", quantile(u$n_edges, c(.25, .5, .75)), "\n")
cat("patients with at least one replicated zone:", length(unique(r$re$patient)), "| replicated zone pairs:", nrow(rp), "\n")
cat("share of a section's edges selected again in the replicate section, pooled over both directions:",
    round(2 * sum(rp$both) / sum(rp$edges_a + rp$edges_b), 3), "\n")
