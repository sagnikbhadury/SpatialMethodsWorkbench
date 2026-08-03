#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
include_related <- "--include-related" %in% args
skip_vignettes <- "--skip-vignettes" %in% args

find_root <- function(start = getwd()) {
  candidate <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(candidate, "DESCRIPTION")) &&
        dir.exists(file.path(candidate, "vendor", "engines"))) return(candidate)
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  stop("Run this script from the extracted Spatial Methods Workbench release.")
}

root <- find_root()
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (.Platform$OS.type == "windows") {
  rtools_version <- paste0(
    R.version$major,
    strsplit(R.version$minor, "[.]")[[1]][1]
  )
  rtools_home <- c(
    Sys.getenv(paste0("RTOOLS", rtools_version, "_HOME")),
    paste0("C:/rtools", rtools_version)
  )
  rtools_home <- rtools_home[nzchar(rtools_home) & dir.exists(rtools_home)]
  if (length(rtools_home)) {
    tool_paths <- c(
      file.path(rtools_home[[1]], "usr", "bin"),
      file.path(rtools_home[[1]], "x86_64-w64-mingw32.static.posix", "bin")
    )
    Sys.setenv(PATH = paste(c(tool_paths, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  }
}

install_missing <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing CRAN dependencies: ", paste(missing, collapse = ", "))
    install.packages(
      missing,
      dependencies = c("Depends", "Imports", "LinkingTo")
    )
  }
  still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) {
    message(
      "Retrying unavailable or binary-incompatible packages from source: ",
      paste(still_missing, collapse = ", ")
    )
    install.packages(
      still_missing,
      dependencies = c("Depends", "Imports", "LinkingTo"),
      type = "source"
    )
    still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  }
  if (length(still_missing)) {
    stop(
      "Dependency installation failed for: ", paste(still_missing, collapse = ", "),
      ". Confirm that the R version and compiler toolchain match. On Windows, ",
      "install the Rtools version for this R release; rstan binaries built for a ",
      "newer R patch release may require upgrading R or compiling from source."
    )
  }
}

install_source <- function(relative_path, package_name) {
  path <- file.path(root, relative_path)
  if (!file.exists(file.path(path, "DESCRIPTION"))) {
    stop("Bundled source is missing DESCRIPTION: ", path)
  }
  message("Installing bundled ", package_name, " from ", relative_path)
  remotes::install_local(
    path,
    dependencies = FALSE,
    upgrade = "never",
    force = TRUE,
    build = TRUE,
    quiet = FALSE
  )
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop("Bundled package failed to install: ", package_name)
  }
}

runtime <- c(
  "shiny", "bslib", "DT", "ggplot2", "jsonlite", "markdown", "zip",
  "nnet", "knitr", "rmarkdown", "remotes", "testthat"
)
engine_dependencies <- c(
  "CVglasso", "doParallel", "foreach", "Matrix", "matrixcalc", "MASS",
  "mnormt", "rstan", "sparsepca"
)
related_dependencies <- c(
  "Rcpp", "RcppArmadillo", "BH", "BayesGPfit", "reshape2", "viridis",
  "RSpectra", "RANN", "imager", "ggpubr", "coda"
)

if (include_related) {
  message("Installing or refreshing related-workflow dependencies before compiled namespaces are loaded")
  install.packages(
    related_dependencies,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}
install_missing(runtime)
if (include_related) install_missing(related_dependencies)
install_missing(engine_dependencies)

install_source(file.path("vendor", "engines", "GP-GHS"), "GPGHS")
install_source(file.path("vendor", "engines", "ISPAT-3D"), "ISPAT3D")
install_source(file.path("vendor", "engines", "ISPAT"), "ISPAT")

if (include_related) {
  install_source(file.path("vendor", "reference-workflows", "STCAR"), "STCAR")
  if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Linux") {
    install_source(file.path("vendor", "reference-workflows", "SBLF"), "SBLF")
  } else {
    message("SBLF source is bundled but its upstream package supports Linux only; installation skipped on this platform.")
  }
}

message("Installing SpatialMethodsWorkbench package API")
remotes::install_local(
  root,
  dependencies = FALSE,
  build_vignettes = !skip_vignettes,
  upgrade = "never",
  force = TRUE
)

required <- c("SpatialMethodsWorkbench", "ISPAT", "GPGHS", "ISPAT3D")
status <- setNames(vapply(required, requireNamespace, logical(1), quietly = TRUE), required)
print(status)
stopifnot(all(status))

cat(
  "\nComplete Workbench installation verified.\n",
  "Launch from the extracted release directory with:\n",
  "  shiny::runApp(\".\", host = \"127.0.0.1\", port = 3838)\n",
  sep = ""
)
