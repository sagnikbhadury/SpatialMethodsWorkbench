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

test_that("spatial AI workflows complete with held-out outputs", {
  data <- demo_spatial_data(160)
  features <- c("epithelial", "APC", "CTL", "Treg", "T_helper")
  clusters <- run_spatial_clustering(data, list(x = "x", y = "y"), features,
                                     list(clusters = 4, spatial_weight = .5, seed = 2))
  expect_equal(length(unique(clusters$table$cluster)), 4)
  expect_s3_class(clusters$plot, "ggplot")
  neural <- run_neural_prediction(data, list(outcome = "disease", x = "x", y = "y", z = NULL),
                                  features, list(hidden_units = 4, decay = .02, seed = 2))
  expect_true("accuracy" %in% names(neural$summary))
  expect_s3_class(neural$plot, "ggplot")
})

test_that("wide-image regression workflows complete", {
  set.seed(4)
  n <- 90; p <- 12
  input <- matrix(rnorm(n * p), n, p); output <- input[, 1:6] %*% matrix(rnorm(6 * p, sd = .2), 6, p) + matrix(rnorm(n * p), n, p)
  data <- data.frame(outcome = input[, 1] - input[, 2] + rnorm(n, sd = .3))
  data[paste0("pixel_", seq_len(p))] <- input
  data[paste0("input__", seq_len(p))] <- input
  data[paste0("output__", seq_len(p))] <- output
  scalar_features <- paste0("pixel_", seq_len(p))
  scalar <- run_scalar_image(data, list(outcome = "outcome"), scalar_features, list(image_ridge = 5, seed = 3))
  expect_equal(nrow(scalar$table), p)
  paired <- c(paste0("input__", seq_len(p)), paste0("output__", seq_len(p)))
  latent <- run_image_to_image(data, list(), paired, list(latent_factors = 4, seed = 3))
  expect_equal(nrow(latent$table), p)
  expect_true(is.finite(unname(latent$summary["overall_heldout_rmse"])))
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
  expect_match(workbench_citation(), "SpatialMethodsWorkbench", fixed = TRUE)
  expect_match(workbench_citation(), "Version 0.2.1", fixed = TRUE)
  expect_match(workbench_citation(), "10.5281/zenodo.21764196", fixed = TRUE)
})
