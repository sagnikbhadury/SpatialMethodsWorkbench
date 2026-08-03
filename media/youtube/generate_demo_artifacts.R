source("R/utils.R")
source("R/registry.R")
source("R/methods_core.R")
source("R/methods_advanced.R")
source("R/reporting.R")

output_dir <- file.path("media", "youtube", "captures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- demo_spatial_data(n = 180L, seed = 2026L)
mapping <- list(
  x = "x", y = "y", z = NULL, zone = "zone", id = "sample_id",
  outcome = "disease", exposure = NULL, mediator = NULL
)
features <- c("epithelial", "APC", "CTL", "Treg", "T_helper")
result <- run_analysis(
  "spatial_qc", data, mapping, features,
  list(feature = "epithelial", neighbors = 8L, permutations = 199L, seed = 1L)
)

ggplot2::ggsave(
  file.path(output_dir, "06_synthetic_result_figure.png"),
  result$plot + ggplot2::labs(caption = "Synthetic demonstration data — method training only") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 24, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 16),
      plot.caption = ggplot2::element_text(size = 13),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 13)
    ),
  width = 16, height = 9, dpi = 120, bg = "white"
)

utils::write.csv(result$table, file.path(output_dir, "synthetic_result_table.csv"), row.names = FALSE)
jsonlite::write_json(result_manifest(result), file.path(output_dir, "synthetic_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
bundle_path <- normalizePath(file.path("media", "youtube", "synthetic_reproducibility_bundle.zip"), mustWork = FALSE)
write_result_bundle(result, bundle_path)

cat(sprintf(
  "Moran's I: %.6f\nPermutation p-value: %.6f\n",
  as.numeric(result$summary[["moran_i"]]),
  as.numeric(result$summary[["p_value"]])
))
