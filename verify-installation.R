#!/usr/bin/env Rscript

packages <- c(
  SpatialMethodsWorkbench = "SpatialMethodsWorkbench",
  ISPAT3D = "ISPAT3D",
  gpboost = "gpboost",
  GPGHS = "GPGHS",
  nnet = "nnet"
)
minimum <- c(ISPAT3D = "0.3.0", GPGHS = "1.0.0")

installed <- vapply(unname(packages), requireNamespace, logical(1), quietly = TRUE)
versions <- vapply(unname(packages), function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
}, character(1))
status <- data.frame(
  component = names(packages),
  package = unname(packages),
  installed = installed,
  version = versions,
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
for (p in names(minimum)) {
  if (utils::packageVersion(p) < minimum[[p]]) {
    stop(sprintf("%s %s is installed but version %s or newer is required.", p, utils::packageVersion(p), minimum[[p]]))
  }
}
stopifnot(length(registry) == 12L)
stopifnot(identical(registry$ispat$package, "ISPAT3D"), identical(registry$ispat3d$package, "ISPAT3D"))
cat("\nVerification passed: all 12 paths are registered; ISPAT3D (2D planar and 3D volumetric fits) and GPGHS are installed at compatible versions.\n")
