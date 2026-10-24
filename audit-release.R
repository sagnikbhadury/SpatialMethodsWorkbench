#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required <- c(
  "app.R", "DESCRIPTION", "NAMESPACE", "README.md", "REFERENCES.md",
  "install-complete.R", "verify-installation.R", "vendor/manifest.json",
  "vignettes/app-walkthrough.Rmd", "vignettes/spatial-networks.Rmd",
  "vignettes/prediction-clustering-mediation.Rmd",
  "vignettes/imaging-shape-3d.Rmd", "vignettes/reproducibility-bundles.Rmd",
  "inst/extdata/spatial_demo.csv", "inst/extdata/spatial_3d_demo.csv",
  "inst/extdata/image_demo.csv", "inst/extdata/contour_demo.csv"
)

missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) stop("Release files missing: ", paste(missing, collapse = ", "))

release_version <- "0.4.0"
description <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]
manifest <- jsonlite::fromJSON(file.path(root, "vendor", "manifest.json"))
stopifnot(
  identical(unname(description[["Version"]]), release_version),
  identical(manifest$workbench_release, release_version),
  nrow(manifest$sources) == 8L,
  length(unique(manifest$sources$commit)) == 8L,
  all(nchar(manifest$sources$commit) == 40L)
)

engines <- list(
  GPGHS = list(path = "vendor/engines/GP-GHS", source = "GP-GHS", minimum = "1.0.0"),
  ISPAT3D = list(path = "vendor/engines/ISPAT-3D", source = "ISPAT-3D", minimum = "0.3.0")
)
for (package_name in names(engines)) {
  engine <- engines[[package_name]]
  fields <- read.dcf(file.path(root, engine$path, "DESCRIPTION"))[1, ]
  stopifnot(identical(unname(fields[["Package"]]), package_name))
  version <- unname(fields[["Version"]])
  if (package_version(version) < engine$minimum) {
    stop(sprintf("Bundled %s is version %s; the Workbench requires >= %s", package_name, version, engine$minimum))
  }
  recorded <- manifest$sources$version[manifest$sources$name == engine$source]
  if (!identical(recorded, version)) {
    stop(sprintf("vendor/manifest.json records %s version %s but the bundled DESCRIPTION says %s", engine$source, recorded, version))
  }
}
if ("ISPAT" %in% manifest$sources$name || dir.exists(file.path(root, "vendor", "engines", "ISPAT"))) {
  stop("The retired standalone ISPAT package must not be bundled; the 2D path uses ISPAT3D::ispat3d_fit_2d")
}

vendor_files <- list.files(file.path(root, "vendor"), recursive = TRUE, all.files = TRUE)
forbidden <- vendor_files[grepl("(^|/)[.]git(/|$)|[.](o|so|dll|a)$", vendor_files, ignore.case = TRUE)]
if (length(forbidden)) stop("Forbidden generated or Git files in vendor bundle: ", paste(forbidden, collapse = ", "))
startup_files <- vendor_files[basename(vendor_files) %in% c(".Rprofile", ".Renviron")]
if (length(startup_files)) stop("Project startup files are not allowed in the vendor bundle: ", paste(startup_files, collapse = ", "))

allowed_roots <- c("engines", "reference-workflows", "manifest.json", "README.md")
vendor_roots <- unique(sub("/.*$", "", vendor_files))
unexpected_roots <- setdiff(vendor_roots, allowed_roots)
if (length(unexpected_roots)) stop("Unexpected vendor content: ", paste(unexpected_roots, collapse = ", "))

info <- file.info(file.path(root, "vendor", vendor_files))
cat(
  "Release audit passed:\n",
  "  Workbench version: ", description[["Version"]], "\n",
  "  Registered public sources: ", nrow(manifest$sources), "\n",
  "  Bundled engines: ", paste(sprintf("%s %s", names(engines), vapply(engines, function(e) read.dcf(file.path(root, e$path, "DESCRIPTION"))[1, "Version"], character(1))), collapse = ", "), "\n",
  "  Vendor files: ", nrow(info), "\n",
  "  Vendor size (MiB): ", sprintf("%.2f", sum(info$size, na.rm = TRUE) / 1024^2), "\n",
  sep = ""
)
