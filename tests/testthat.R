library(testthat)
if (requireNamespace("SpatialMethodsWorkbench", quietly = TRUE)) {
  library(SpatialMethodsWorkbench)
  test_check("SpatialMethodsWorkbench")
} else {
  test_dir("tests/testthat")
}
