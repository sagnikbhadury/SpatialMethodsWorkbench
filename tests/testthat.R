library(testthat)
if (dir.exists("tests/testthat")) {
  test_dir("tests/testthat")
} else {
  library(SpatialMethodsWorkbench)
  test_check("SpatialMethodsWorkbench")
}
