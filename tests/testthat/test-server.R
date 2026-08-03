test_that("Shiny server executes a compatible selected analysis", {
  app_path <- file.path("..", "..", "app.R")
  skip_if_not(file.exists(app_path), "The source-repository Shiny app is not included in the installed package check.")
  app_env <- new.env(parent = globalenv())
  old_user <- Sys.getenv("R_USER", unset = NA_character_)
  Sys.setenv(R_USER = tempdir())
  on.exit(if (is.na(old_user)) Sys.unsetenv("R_USER") else Sys.setenv(R_USER = old_user), add = TRUE)
  suppressWarnings(source(app_path, local = app_env))
  shiny::testServer(app_env$server, {
    session$setInputs(
      use_demo = TRUE,
      x_col = "x", y_col = "y", z_col = "", zone_col = "zone",
      id_col = "sample_id", outcome_col = "disease",
      exposure_col = "epithelial", mediator_col = "APC",
      feature_cols = c("epithelial", "APC", "CTL", "Treg", "T_helper"),
      citation_ack = TRUE,
      analysis = "spatial_qc", param_feature = "epithelial",
      param_neighbors = 6, param_permutations = 19, param_seed = 1
    )
    session$flushReact()
    expect_true(validation()$ok)
    expect_true(recommendations()$spatial_qc$compatible)
    session$setInputs(run = 1)
    session$flushReact()
    expect_null(result()$error)
    expect_s3_class(result()$plot, "ggplot")
  })
})
