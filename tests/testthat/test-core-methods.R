test_that("demo data validates and recommends compatible methods", {
  data <- demo_spatial_data(120)
  mapping <- coerce_mapping(data, list(x = "x", y = "y", z = NULL, zone = "zone",
                                       id = "sample_id", outcome = "disease",
                                       exposure = "epithelial", mediator = "APC"))
  features <- c("epithelial", "APC", "CTL", "Treg", "T_helper")
  validation <- validate_dataset(data, mapping)
  expect_true(validation$ok)
  rec <- recommend_methods(data, mapping, features)
  expect_true(rec$spatial_qc$compatible)
  expect_true(rec$conditional_network$compatible)
  expect_true(rec$spatial_ml$compatible)
  expect_false(rec$ispat3d$compatible)
})

test_that("spatial QC produces a permutation result", {
  data <- demo_spatial_data(90)
  mapping <- list(x = "x", y = "y")
  result <- run_spatial_qc(data, mapping, "epithelial",
                           list(feature = "epithelial", neighbors = 6, permutations = 19, seed = 1))
  expect_s3_class(result$plot, "ggplot")
  expect_equal(nrow(result$table), 4)
  expect_true(is.finite(as.numeric(result$summary["moran_i"])))
})

test_that("conditional networks are returned overall and by zone", {
  data <- demo_spatial_data(150)
  features <- c("epithelial", "APC", "CTL", "Treg", "T_helper")
  result <- run_conditional_network(data, list(zone = "zone"), features,
                                    list(threshold = .1, ridge = .05, by_zone = TRUE))
  expect_gte(length(result$matrices), 2)
  expect_true(all(c("source", "target", "weight") %in% names(result$table)))
})

test_that("prediction and mediation complete", {
  data <- demo_spatial_data(150)
  features <- c("epithelial", "APC", "CTL", "Treg", "T_helper")
  prediction <- run_spatial_ml(data, list(outcome = "disease", x = "x", y = "y", z = NULL),
                               features, list(seed = 3))
  expect_s3_class(prediction$plot, "ggplot")
  mediation <- run_mediation(data, list(exposure = "epithelial", mediator = "APC", outcome = "CTL"),
                             features, list(bootstrap = 20, seed = 3))
  expect_equal(nrow(mediation$table), 5)
})

test_that("shape PCA resamples multiple contours", {
  theta <- seq(0, 2 * pi, length.out = 25)[-25]
  shapes <- do.call(rbind, lapply(1:8, function(i) {
    data.frame(id = paste0("C", i), x = (1 + i / 20) * cos(theta), y = (1 - i / 30) * sin(theta))
  }))
  result <- run_shape_pca(shapes, list(id = "id", x = "x", y = "y"), character(), list(landmarks = 20))
  expect_equal(nrow(result$table), 8)
  expect_s3_class(result$plot, "ggplot")
})

test_that("result bundles contain reproducibility artifacts", {
  data <- demo_spatial_data(70)
  result <- run_spatial_qc(data, list(x = "x", y = "y"), "APC",
                           list(feature = "APC", neighbors = 5, permutations = 9, seed = 1))
  result$analysis_id <- "spatial_qc"; result$started_at <- "2026-08-02 UTC"; result$elapsed_seconds <- 0.1
  result$params <- list(); result$mapping <- list(x = "x", y = "y"); result$features <- "APC"
  result$citations <- method_citations("spatial_qc")
  path <- tempfile(fileext = ".zip")
  write_result_bundle(result, path)
  expect_true(file.exists(path))
  files <- utils::unzip(path, list = TRUE)$Name
  expect_true(all(c("results.csv", "manifest.json", "CITATION.txt", "COLLABORATION.txt", "figure.png", "result_object.rds") %in% files))
  expect_false(any(c("interpretation.md", "optional_llm_prompt.txt") %in% files))
})
