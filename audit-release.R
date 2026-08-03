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

description <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]
manifest <- jsonlite::fromJSON(file.path(root, "vendor", "manifest.json"))
stopifnot(
  identical(unname(description[["Version"]]), "0.3.0"),
  identical(manifest$workbench_release, "0.3.0"),
  nrow(manifest$sources) == 9L,
  length(unique(manifest$sources$commit)) == 9L,
  all(nchar(manifest$sources$commit) == 40L)
)

engine_paths <- c(
  ISPAT = "vendor/engines/ISPAT",
  GPGHS = "vendor/engines/GP-GHS",
  ISPAT3D = "vendor/engines/ISPAT-3D"
)
for (package_name in names(engine_paths)) {
  fields <- read.dcf(file.path(root, engine_paths[[package_name]], "DESCRIPTION"))[1, ]
  stopifnot(identical(unname(fields[["Package"]]), package_name))
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
  "  Vendor files: ", nrow(info), "\n",
  "  Vendor size (MiB): ", sprintf("%.2f", sum(info$size, na.rm = TRUE) / 1024^2), "\n",
  sep = ""
)
