## Workflow figure for the manuscript, drawn from real data and saved fits.
## Two inputs are shown: one CODEX image (cell types, quadrat counts) and one
## zone of a Visium slide (genes, log-normalised expression). Both feed the
## same model. Every panel is computed from the data or the saved fit.
## The imaging-only version of this figure is in archive_pre_20261003/.
##
## Usage (R 4.4.1): Rscript make_workflow_figure.R
##   (after crc_analysis.R with CRC_FEATURE=count, st_fit.R and st_network_summary.R)

suppressPackageStartupMessages({ library(ggplot2); library(gridExtra); library(grid) })

FIT_DIR <- Sys.getenv("CRC_OUT_DIR", "crc_results_counts")
IMG     <- Sys.getenv("WORKFLOW_IMG", "21_A")
TYPE_A  <- "CD163+ macros"; TYPE_B <- "CD8+ T cells"
ST_DIR  <- Sys.getenv("ST_OUT_DIR", "st_results")
ST_FIT  <- file.path(ST_DIR, "fit_m2")
ST_ZONE <- as.integer(Sys.getenv("WORKFLOW_ZONE", "1"))

base <- theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), plot.title = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 8), legend.key.size = unit(0.35, "cm"),
        axis.text = element_blank(), axis.ticks = element_blank(),
        strip.background = element_blank())
div_scale <- function(lim, name = "Coefficient")
  scale_fill_gradient2(low = "#2a78d6", mid = "#f4f4f2", high = "#e34948", midpoint = 0,
                       limits = c(-lim, lim), name = name, breaks = c(-1, 0, 1) * floor(10 * lim) / 10,
                       aesthetics = c("fill", "color"))
sign_cols <- c(Positive = "#e34948", Negative = "#2a78d6", Mixed = "grey55")
sign_lty  <- c(Positive = "solid", Negative = "longdash", Mixed = "dotted")
cat8 <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7")   # fixed categorical order

## ================= Input B (bottom row): multiplexed imaging ======================
set <- readRDS(file.path(FIT_DIR, "run_settings.rds")); NG <- set$n_grid
g   <- readRDS(file.path(FIT_DIR, "graph_list_final.rds"))[[IMG]]
CT  <- g$cell_types; p <- length(CT)
dat <- read.csv("CRC_data_cleaned_corrected.csv", stringsAsFactors = FALSE)
dat <- dat[dat$img == IMG & dat$type %in% CT, ]
pretty <- c("adipocytes" = "Adipocytes", "B cells" = "B cells",
            "CD163+ macros" = "CD163+ macrophages", "CD4+ T cells" = "CD4+ T cells",
            "CD68+ macros" = "CD68+ macrophages", "CD8+ T cells" = "CD8+ T cells",
            "generic immune" = "Generic immune", "granulocytes" = "Granulocytes",
            "memory CD4+ T" = "Memory CD4+ T", "plasma cells" = "Plasma cells",
            "smooth muscle" = "Smooth muscle", "stroma" = "Stroma", "Tregs" = "Tregs",
            "tumor cells" = "Tumor cells", "vasculature" = "Vasculature")

show <- c("CD163+ macros", "CD8+ T cells", "tumor cells", "stroma")
dat$shown <- ifelse(dat$type %in% show, pretty[dat$type], "Other types")
dat$shown <- factor(dat$shown, c(pretty[show], "Other types"))
a1 <- ggplot(dat[order(dat$shown == "Other types", decreasing = TRUE), ],
             aes(X, Y, color = shown, shape = shown)) +
  geom_point(size = 0.6, alpha = 0.85) +
  scale_color_manual(values = c(setNames(cat8[1:4], pretty[show]), "Other types" = "grey75"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18, 16), name = NULL) +
  coord_equal() + labs(x = NULL, y = NULL, title = "B1. Multiplexed image: phenotyped cells",
    subtitle = sprintf("One CODEX image: %s cells of %d types", format(nrow(dat), big.mark = ","),
                       length(unique(dat$type)))) +
  guides(color = guide_legend(override.aes = list(size = 2), nrow = 2)) +
  base + theme(legend.position = "bottom")

xb <- seq(min(dat$X), max(dat$X), length.out = NG + 1)
yb <- seq(min(dat$Y), max(dat$Y), length.out = NG + 1)
grid_a <- expand.grid(X = (xb[-1] + xb[-(NG + 1)]) / 2, Y = (yb[-1] + yb[-(NG + 1)]) / 2)
cnt <- function(ct) { d <- dat[dat$type == ct, ]
  as.vector(table(factor(findInterval(d$X, xb, all.inside = TRUE), 1:NG),
                  factor(findInterval(d$Y, yb, all.inside = TRUE), 1:NG))) }
qd <- rbind(data.frame(grid_a, type = unname(pretty[TYPE_A]), n = cnt(TYPE_A)),
            data.frame(grid_a, type = unname(pretty[TYPE_B]), n = cnt(TYPE_B)))
qd$type <- factor(qd$type, pretty[c(TYPE_A, TYPE_B)])
a2 <- ggplot(qd, aes(X, Y, fill = sqrt(n))) +
  geom_tile(width = diff(xb)[1], height = diff(yb)[1]) + facet_wrap(~ type, nrow = 1) +
  scale_fill_gradient(low = "#eef4fd", high = "#0d366b", name = expression(sqrt(count)),
                      breaks = scales::breaks_pretty(3)) +
  coord_equal() + labs(x = NULL, y = NULL, title = "B2. Variables: cell-type abundance",
    subtitle = sprintf("Counts in a %d x %d grid of quadrats (2 of %d cell types)", NG, NG, p)) +
  base + theme(legend.position = "bottom", legend.key.width = unit(0.8, "cm"))

i <- match(TYPE_A, CT); j <- match(TYPE_B, CT); ni <- g$node_results[[i]]
fld <- data.frame(grid_a, beta = ni$edge_maps[, match(j, ni$nbr_idx)])
a3 <- ggplot(fld, aes(X, Y, fill = beta)) +
  geom_tile(width = diff(xb)[1], height = diff(yb)[1]) + div_scale(max(abs(fld$beta))) +
  coord_equal() + labs(x = NULL, y = NULL, title = "B3. Coefficient function",
    subtitle = sprintf("%s given\n%s", pretty[TYPE_A], pretty[TYPE_B])) +
  base + theme(legend.position = "bottom", legend.key.width = unit(0.7, "cm"))

est <- !is.na(sapply(g$node_results, function(r) r$tau_mean))
ang <- pi / 2 - 2 * pi * (seq_len(p) - 0.5) / p
nodes <- data.frame(x = cos(ang), y = sin(ang), est = est, label = pretty[CT])
Ea <- which(upper.tri(g$adj) & g$adj == 1, arr.ind = TRUE)
edges_a <- do.call(rbind, lapply(seq_len(nrow(Ea)), function(e) {
  a <- Ea[e, 1]; b <- Ea[e, 2]; na <- g$node_results[[a]]; nb <- g$node_results[[b]]
  ma <- mean(na$edge_maps[, match(b, na$nbr_idx)]); mb <- mean(nb$edge_maps[, match(a, nb$nbr_idx)])
  data.frame(x = nodes$x[a], y = nodes$y[a], xend = nodes$x[b], yend = nodes$y[b],
             sign = if (ma > 0 && mb > 0) "Positive" else if (ma < 0 && mb < 0) "Negative" else "Mixed")
}))
edges_a$sign <- factor(edges_a$sign, names(sign_cols))
a4 <- ggplot() +
  geom_segment(data = edges_a, aes(x, y, xend = xend, yend = yend, color = sign, linetype = sign),
               linewidth = 0.6, show.legend = TRUE) +
  geom_point(data = nodes, aes(x, y, fill = est), shape = 21, size = 2.6, color = "grey20") +
  geom_text(data = nodes, aes(1.12 * x, 1.12 * y, label = label, hjust = ifelse(x > 0, 0, 1)), size = 2.3) +
  scale_color_manual(values = sign_cols, name = "Association", drop = FALSE) +
  scale_linetype_manual(values = sign_lty, name = "Association", drop = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "grey25", `FALSE` = "white"), guide = "none") +
  coord_equal(xlim = c(-2.4, 2.4), ylim = c(-1.15, 1.15)) +
  labs(title = "B4. Network of cell types",
       subtitle = sprintf("%d edges in this image; open nodes: fewer than %d cells", nrow(edges_a), set$min_cells)) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 8),
        legend.position = "bottom", legend.key.width = unit(0.8, "cm"))

## ================= Input A (top row): spatial transcriptomics ====================
prep <- readRDS(file.path(ST_DIR, "st_prepared.rds"))
sf   <- readRDS(file.path(ST_FIT, "st_fit.rds")); genes <- sf$genes; pg <- length(genes)
fo   <- readRDS(file.path(ST_FIT, "summary", "figure_objects.rds"))
z    <- sf$fits[[ST_ZONE]]; obs <- prep$obs
ex   <- fo$examples[2, ]; GENE_A <- ex$gene_i; GENE_B <- ex$gene_j
obs$zone_name <- factor(prep$zone_tab$name[obs$zone], prep$zone_tab$name)
b1 <- ggplot(obs, aes(x, -y, color = zone_name)) + geom_point(size = 0.45) +
  scale_color_manual(values = setNames(cat8[seq_len(nrow(prep$zone_tab))], prep$zone_tab$name), name = NULL) +
  coord_equal() +
  labs(x = NULL, y = NULL, title = "A1. Spatial transcriptomics: spots in zones",
       subtitle = sprintf("%s Visium spots, %d pathologist zones",
                          format(nrow(obs), big.mark = ","), nrow(prep$zone_tab))) +
  base + theme(legend.position = "none")

zo <- obs[z$idx, ]
xd <- rbind(data.frame(x = zo$x, y = zo$y, gene = GENE_A, v = prep$expr[z$idx, GENE_A]),
            data.frame(x = zo$x, y = zo$y, gene = GENE_B, v = prep$expr[z$idx, GENE_B]))
xd$gene <- factor(xd$gene, c(GENE_A, GENE_B))
b2 <- ggplot() + geom_point(data = obs[-z$idx, ], aes(x, -y), color = "grey90", size = 0.3) +
  geom_point(data = xd, aes(x, -y, color = v), size = 0.5) + facet_wrap(~ gene, nrow = 1) +
  scale_color_gradient(low = "#eef4fd", high = "#0d366b", name = "log expression",
                       breaks = scales::breaks_pretty(3)) +
  coord_equal() + labs(x = NULL, y = NULL, title = "A2. Variables: gene expression",
    subtitle = sprintf("Spots of the zone '%s' (2 of %d genes)", z$name, pg)) +
  base + theme(legend.position = "bottom", legend.key.width = unit(0.8, "cm"))

s <- match(GENE_A, genes); t <- match(GENE_B, genes); m2 <- ncol(z$phi)
bz <- as.vector(z$phi %*% matrix(z$theta[[s]], m2, pg - 1)[, match(t, (1:pg)[-s])])
b3 <- ggplot() + geom_point(data = obs[-z$idx, ], aes(x, -y), color = "grey90", size = 0.3) +
  geom_point(data = data.frame(zo, beta = bz), aes(x, -y, color = beta), size = 0.6) +
  div_scale(max(abs(bz))) + coord_equal() +
  labs(x = NULL, y = NULL, title = "A3. Coefficient function",
       subtitle = sprintf("%s given\n%s", GENE_A, GENE_B)) +
  base + theme(legend.position = "bottom", legend.key.width = unit(0.7, "cm"))

L <- fo$layout; Eb <- which(upper.tri(z$adj) & z$adj == 1, arr.ind = TRUE)
sg <- ifelse(z$beta_mean[Eb] > 0 & t(z$beta_mean)[Eb] > 0, "Positive",
             ifelse(z$beta_mean[Eb] < 0 & t(z$beta_mean)[Eb] < 0, "Negative", "Mixed"))
edges_b <- data.frame(x = L[Eb[, 1], 1], y = L[Eb[, 1], 2], xend = L[Eb[, 2], 1], yend = L[Eb[, 2], 2],
                      sign = factor(sg, names(sign_cols)))
deg <- rowSums(z$adj)
nodes_b <- data.frame(x = L[, 1], y = L[, 2], deg = deg, gene = genes)
lab_b <- nodes_b[order(-deg)[1:5], ]
b4 <- ggplot() +
  geom_segment(data = edges_b, aes(x, y, xend = xend, yend = yend, color = sign, linetype = sign),
               linewidth = 0.35, alpha = 0.8) +
  geom_point(data = nodes_b, aes(x, y, size = deg), shape = 21, fill = "grey25", color = "white", stroke = 0.2) +
  geom_label(data = lab_b, aes(x, y, label = gene), size = 2.1, fontface = "bold", vjust = -0.5,
             label.size = 0, label.padding = unit(0.08, "lines"), alpha = 0.8) +
  scale_size_continuous(range = c(0.5, 2.4), guide = "none") +
  scale_color_manual(values = sign_cols, name = "Association", drop = FALSE) +
  scale_linetype_manual(values = sign_lty, name = "Association", drop = FALSE) +
  coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.05, 1.12)) +
  labs(title = "A4. Network of genes",
       subtitle = sprintf("%d edges among %d genes in this zone", nrow(edges_b), pg)) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9.5, face = "bold"), plot.subtitle = element_text(size = 8),
        legend.position = "bottom", legend.key.width = unit(0.8, "cm"))

## ================= The model, common to both inputs ================================
eq <- c(
  "bold('Data: ') * n ~ 'locations' ~ l[1] * ',...,' * l[n] * ', ' ~ p ~ 'variables' ~ Y[1](l) * ',...,' * Y[p](l)",
  "bold('One regression per variable: ') * Y[s](l) == sum(beta[ss*minute](l)*Y[s*minute](l), s*minute != s, '') + epsilon[s](l)",
  "bold('Coefficient function, Hilbert space GP: ') * beta[ss*minute](l) == sum(theta[ss*minute*','*k]*phi[k](l), k == 1, m^2)",
  "bold('Group horseshoe, one local scale per function: ') * theta[ss*minute*','*k] %~% N(0*','~sigma[s]^2*tau[s]^2*lambda[ss*minute]^2*S(omega[k]^2))",
  "bold('Edge between ') * s * bold(' and ') * s*minute * bold(': ') * 'max' * group('{', list(bar(kappa)[s %->% s*minute], bar(kappa)[s*minute %->% s]), '}') < 1/2")
mdl <- ggplot() +
  annotate("text", x = 0.15, y = rev(seq_along(eq)), label = eq, parse = TRUE, hjust = 0, size = 3.1) +
  xlim(0, 10) + ylim(0.4, length(eq) + 0.6) +
  labs(title = "GP-GHS: the same model for both inputs (a variable is a gene in A and a cell type in B)") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9.5, face = "bold"),
        panel.border = element_rect(fill = NA, color = "grey40"), plot.margin = margin(4, 5.5, 4, 5.5))

w <- c(1, 1.55, 0.85, 1.45)
rowA <- arrangeGrob(a1, a2, a3, a4, nrow = 1, widths = w)
rowB <- arrangeGrob(b1, b2, b3, b4, nrow = 1, widths = w)
fig <- arrangeGrob(rowB, mdl, rowA, nrow = 3, heights = c(1, 0.42, 1))   # gene input on top
ggsave(file.path("manuscript", "Image_s", "gpghs_workflow.pdf"), fig, width = 13.5, height = 9.6,
       device = cairo_pdf)
cat("imaging: edges", nrow(edges_a), "| estimable types", sum(est), "| cells", nrow(dat),
    "\nST: zone", ST_ZONE, z$name, "| edges", nrow(edges_b), "| example", GENE_A, "given", GENE_B, "\n")
