project_r <- normalizePath(file.path("..", "..", "R"), mustWork = TRUE)
invisible(lapply(list.files(project_r, pattern = "\\.R$", full.names = TRUE), source, local = FALSE))
