## Figures and tables for simulation_study_v4.R.
## Reads the per-replicate results exactly as saved; no value is adjusted.
##
## Usage: Rscript simulation_figures_v4.R

suppressPackageStartupMessages(library(ggplot2))

OUT_DIR <- Sys.getenv("SIM_OUT_DIR", "simulation_v4")
FIG_DIR <- file.path("manuscript", "Image_s")
raw <- readRDS(file.path(OUT_DIR, "simulation_raw_v4.rds"))

method_levels <- c("GP_GHS", "HS_basis", "HS_const", "Exact_GP",
                   "Glasso", "Nodewise_Lasso", "Corr_Threshold")
method_labels <- c("GP-GHS", "HS-basis", "HS-const", "Pairwise GP",
                   "GLasso", "Nodewise lasso", "Correlation")
## fixed categorical order; shape is the secondary encoding
method_cols   <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                   "#e87ba4", "#008300", "#4a3aa7")
method_shapes <- c(16, 17, 15, 18, 4, 3, 8)
names(method_cols) <- names(method_shapes) <- method_labels
raw$method <- factor(raw$method, method_levels, method_labels)

metrics <- c(F1 = "F1", MCC = "MCC", TPR = "TPR", FDP = "FDP", AUPRC = "AUPRC")

long <- do.call(rbind, lapply(names(metrics), function(m) data.frame(
  raw[, c("p", "scenario", "pi", "rep", "method")],
  metric = metrics[[m]], value = raw[[m]])))
long <- long[long$pi > 0, ]
long$metric <- factor(long$metric, metrics)
agg <- aggregate(value ~ p + scenario + pi + method + metric, long,
                 function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
agg <- do.call(data.frame, agg)
names(agg)[names(agg) %in% c("value.mean", "value.se")] <- c("mean", "se")
agg$scenario <- factor(agg$scenario, c("varying", "constant"),
                       c("Varying coefficients", "Constant coefficients"))

make_fig <- function(d, file, height) {
  g <- ggplot(d, aes(pi, mean, color = method, shape = method, group = method)) +
    geom_line(aes(linewidth = method == "GP-GHS")) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.03,
                  linewidth = 0.35) +
    geom_point(size = 2.4) +
    scale_linewidth_manual(values = c(`TRUE` = 1.1, `FALSE` = 0.55),
                           guide = "none") +
    scale_color_manual(values = method_cols, drop = TRUE, name = NULL) +
    scale_shape_manual(values = method_shapes, drop = TRUE, name = NULL) +
    scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = expression("Graph density " * pi), y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
          legend.position = "bottom", strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
  g <- if (length(unique(d$scenario)) > 1)
    g + facet_grid(scenario ~ metric) else g + facet_wrap(~ metric, nrow = 1)
  g <- g + guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
  ggsave(file.path(FIG_DIR, file), g, width = 10.5, height = height, dpi = 300)
}

make_fig(agg[agg$p == 15, ], "simulation_v4_p15.png", 5.6)
make_fig(droplevels(agg[agg$p == 25, ]), "simulation_v4_p25.png", 3.5)

## ---- tables --------------------------------------------------------------
num <- c("F1", "MCC", "TPR", "FDP", "FP", "AUPRC", "AUROC", "field_cor",
         "coverage", "runtime")
tab <- aggregate(as.matrix(raw[, num]) ~ p + scenario + pi + method, raw,
                 mean, na.action = na.pass)
tab <- tab[order(tab$p, tab$scenario, tab$pi, tab$method), ]
write.csv(tab, file.path(OUT_DIR, "simulation_table_v4.csv"), row.names = FALSE)

## pooled over the four non-empty densities, for the main-text table
pooled <- aggregate(as.matrix(raw[raw$pi > 0, num]) ~ p + scenario + method,
                    raw[raw$pi > 0, ], mean, na.action = na.pass)
pooled <- pooled[order(pooled$p, pooled$scenario, pooled$method), ]
write.csv(pooled, file.path(OUT_DIR, "simulation_pooled_v4.csv"),
          row.names = FALSE)
print(format(pooled, digits = 3), row.names = FALSE)

cat("\nEmpty graph (pi = 0), mean number of selected pairs out of 105:\n")
print(aggregate(FP ~ scenario + method, raw[raw$pi == 0, ], mean), digits = 3)

cat("\nReplicates per cell:\n")
print(table(interaction(raw$p, raw$scenario), raw$pi)[, ] / 1)
