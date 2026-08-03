#!/usr/bin/env Rscript

packages <- c(
  SpatialMethodsWorkbench = "SpatialMethodsWorkbench",
  ISPAT = "ISPAT",
  GPGHS = "GPGHS",
  ISPAT3D = "ISPAT3D",
  nnet = "nnet"
)

status <- data.frame(
  component = names(packages),
  package = unname(packages),
  installed = vapply(unname(packages), requireNamespace, logical(1), quietly = TRUE),
  row.names = NULL
)
print(status)

if (!requireNamespace("SpatialMethodsWorkbench", quietly = TRUE)) {
  stop("SpatialMethodsWorkbench is not installed. Run Rscript install-complete.R first.")
}

registry <- SpatialMethodsWorkbench::analysis_registry()
cat("\nRegistered Workbench paths:", length(registry), "\n")
print(data.frame(
  id = names(registry),
  label = vapply(registry, `[[`, character(1), "label"),
  engine = vapply(registry, `[[`, character(1), "engine"),
  row.names = NULL
))

if (!all(status$installed)) {
  stop("One or more required engines is unavailable; see the status table above.")
}
stopifnot(length(registry) == 12L)
cat("\nVerification passed: all 12 paths are registered and all direct engine packages are installed.\n")

