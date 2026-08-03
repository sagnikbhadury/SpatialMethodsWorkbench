source(file.path("R", "utils.R"), local = FALSE)

output_dir <- file.path("inst", "extdata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(
  demo_spatial_data(n = 72, seed = 2026),
  file.path(output_dir, "spatial_demo.csv"), row.names = FALSE
)
utils::write.csv(
  demo_image_data(n = 48, features = 8, seed = 2026),
  file.path(output_dir, "image_demo.csv"), row.names = FALSE
)
utils::write.csv(
  demo_contour_data(contours = 6, points = 16, seed = 2026),
  file.path(output_dir, "contour_demo.csv"), row.names = FALSE
)
utils::write.csv(
  demo_spatial_3d_data(n = 72, volumes = 2, seed = 2026),
  file.path(output_dir, "spatial_3d_demo.csv"), row.names = FALSE
)
