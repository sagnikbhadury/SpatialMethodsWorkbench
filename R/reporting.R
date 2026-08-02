result_manifest <- function(result) {
  list(
    application = "Spatial Methods Workbench",
    application_version = "0.1.1",
    analysis_id = result$analysis_id,
    method = result$method,
    started_at = result$started_at,
    elapsed_seconds = result$elapsed_seconds,
    mapping = result$mapping,
    features = result$features,
    parameters = result$params,
    summary = as.list(result$summary),
    citations = result$citations,
    session = list(R = R.version.string, platform = R.version$platform)
  )
}

write_result_bundle <- function(result, destination) {
  bundle_dir <- tempfile("spatial-workbench-")
  dir.create(bundle_dir, recursive = TRUE)
  utils::write.csv(result$table %||% data.frame(), file.path(bundle_dir, "results.csv"), row.names = FALSE)
  jsonlite::write_json(result_manifest(result), file.path(bundle_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  writeLines(c(
    "Please cite the software and method references below in work that uses these results:",
    "", paste0("- ", result$citations)
  ), file.path(bundle_dir, "CITATION.txt"), useBytes = TRUE)
  writeLines(c(
    "RESULT INTERPRETATION IS NOT PROVIDED BY THIS APPLICATION",
    "",
    "For interpretation, study-specific guidance, or collaboration, contact:",
    "Sagnik Bhadury <bhadury@umich.edu>",
    "https://sagnikbhadury.github.io/work-with-me/"
  ), file.path(bundle_dir, "COLLABORATION.txt"), useBytes = TRUE)
  if (!is.null(result$plot)) ggplot2::ggsave(file.path(bundle_dir, "figure.png"), result$plot, width = 8, height = 6, dpi = 180)
  saveRDS(result, file.path(bundle_dir, "result_object.rds"))
  zip::zipr(destination, list.files(bundle_dir, full.names = TRUE), root = bundle_dir)
}
