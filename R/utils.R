`%||%` <- function(x, y) if (is.null(x) || !length(x) || identical(x, "")) y else x

safe_names <- function(x) {
  make.unique(make.names(trimws(x), allow_ = TRUE))
}

read_uploaded_csv <- function(path) {
  out <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- safe_names(names(out))
  out
}

demo_spatial_data <- function(n = 180L, seed = 2026L) {
  set.seed(seed)
  x <- stats::runif(n, 0, 10)
  y <- stats::runif(n, 0, 10)
  zone <- cut(x + y, breaks = stats::quantile(x + y, probs = seq(0, 1, length.out = 4)),
              include.lowest = TRUE, labels = c("Low", "Intermediate", "High"))
  epithelial <- scale(sin(x / 2) + cos(y / 3) + stats::rnorm(n, sd = .45))[, 1]
  apc <- scale(.55 * epithelial + sin(y / 2) + stats::rnorm(n, sd = .7))[, 1]
  ctl <- scale(.6 * apc - .25 * epithelial + stats::rnorm(n, sd = .65))[, 1]
  treg <- scale(-.45 * ctl + .35 * epithelial + stats::rnorm(n, sd = .7))[, 1]
  helper <- scale(.45 * apc + .3 * treg + stats::rnorm(n, sd = .75))[, 1]
  disease <- factor(ifelse(.5 * epithelial - .45 * ctl + .25 * treg + stats::rnorm(n) > 0,
                           "Condition_B", "Condition_A"))
  data.frame(sample_id = sprintf("S%03d", seq_len(n)), x = x, y = y, zone = zone,
             disease = disease, epithelial = epithelial, APC = apc, CTL = ctl,
             Treg = treg, T_helper = helper, check.names = FALSE)
}

demo_spatial_3d_data <- function(n = 72L, volumes = 2L, seed = 2026L) {
  n <- max(40L, as.integer(n))
  volumes <- max(1L, min(as.integer(volumes), floor(n / 20L)))
  set.seed(seed)
  volume_id <- rep(sprintf("V%02d", seq_len(volumes)), length.out = n)
  x <- stats::runif(n, 0, 8)
  y <- stats::runif(n, 0, 8)
  z <- stats::runif(n, 0, 4) + (as.integer(factor(volume_id)) - 1L) * .25
  spatial_score <- x + y + 1.5 * z
  zone <- cut(spatial_score,
    breaks = stats::quantile(spatial_score, probs = seq(0, 1, length.out = 4)),
    include.lowest = TRUE, labels = c("Low", "Intermediate", "High")
  )
  epithelial <- scale(sin(x / 2) + cos(y / 2.5) + .2 * z + stats::rnorm(n, sd = .5))[, 1]
  apc <- scale(.5 * epithelial + sin(z) + stats::rnorm(n, sd = .7))[, 1]
  ctl <- scale(.55 * apc - .2 * epithelial + stats::rnorm(n, sd = .7))[, 1]
  treg <- scale(-.35 * ctl + .3 * epithelial + stats::rnorm(n, sd = .75))[, 1]
  helper <- scale(.4 * apc + .25 * treg + stats::rnorm(n, sd = .75))[, 1]
  data.frame(
    volume_id = volume_id, x = x, y = y, z = z, zone = zone,
    epithelial = epithelial, APC = apc, CTL = ctl, Treg = treg,
    T_helper = helper, check.names = FALSE
  )
}

demo_image_data <- function(n = 48L, features = 8L, seed = 2026L) {
  n <- max(36L, as.integer(n))
  features <- max(5L, min(as.integer(features), 40L))
  set.seed(seed)
  input <- matrix(stats::rnorm(n * features), nrow = n, ncol = features)
  loading <- matrix(stats::rnorm(features * features, sd = .18), features, features)
  output <- input %*% loading + matrix(stats::rnorm(n * features, sd = .65), n, features)
  outcome <- .8 * input[, 1] - .55 * input[, 2] + .3 * input[, 3] + stats::rnorm(n, sd = .5)
  out <- data.frame(subject_id = sprintf("I%03d", seq_len(n)), outcome = outcome, check.names = FALSE)
  out[paste0("pixel_", seq_len(features))] <- input
  out[paste0("input__", seq_len(features))] <- input
  out[paste0("output__", seq_len(features))] <- output
  out
}

demo_contour_data <- function(contours = 6L, points = 16L, seed = 2026L) {
  contours <- max(5L, as.integer(contours))
  points <- max(8L, as.integer(points))
  set.seed(seed)
  theta <- seq(0, 2 * pi, length.out = points + 1L)[-1L]
  do.call(rbind, lapply(seq_len(contours), function(i) {
    phase <- stats::runif(1, -.12, .12)
    radial <- 1 + .08 * sin(3 * theta + phase) + stats::rnorm(points, sd = .015)
    data.frame(
      contour_id = sprintf("C%02d", i),
      order = seq_len(points),
      x = radial * (1 + i / 35) * cos(theta),
      y = radial * (1 - i / 45) * sin(theta)
    )
  }))
}

numeric_columns <- function(data) names(data)[vapply(data, is.numeric, logical(1))]
categorical_columns <- function(data) names(data)[vapply(data, function(x) is.factor(x) || is.character(x), logical(1))]

coerce_mapping <- function(data, mapping) {
  keys <- c("x", "y", "z", "zone", "id", "outcome", "exposure", "mediator")
  mapping[keys] <- lapply(mapping[keys], function(value) {
    if (is.null(value) || !nzchar(value) || !value %in% names(data)) NULL else value
  })
  mapping
}

validate_dataset <- function(data, mapping = list()) {
  issues <- character()
  warnings <- character()
  if (!is.data.frame(data) || !nrow(data)) issues <- c(issues, "The dataset has no rows.")
  if (ncol(data) < 2L) issues <- c(issues, "At least two columns are required.")
  if (anyDuplicated(names(data))) issues <- c(issues, "Column names must be unique.")
  if (nrow(data) > 100000L) warnings <- c(warnings, "Large input: advanced Bayesian analyses may need a dedicated deployment.")
  missing_rate <- if (length(data)) mean(is.na(data)) else 0
  if (missing_rate > .2) warnings <- c(warnings, sprintf("%.1f%% of values are missing.", 100 * missing_rate))
  for (axis in c("x", "y", "z")) {
    column <- mapping[[axis]]
    if (!is.null(column) && !is.numeric(data[[column]])) issues <- c(issues, paste(toupper(axis), "coordinates must be numeric."))
  }
  list(ok = !length(issues), issues = issues, warnings = warnings,
       n = nrow(data), p = ncol(data), numeric = numeric_columns(data),
       categorical = categorical_columns(data), missing_rate = missing_rate)
}

sanitize_for_analysis <- function(data, columns) {
  columns <- unique(columns[!is.na(columns) & nzchar(columns) & columns %in% names(data)])
  out <- data[, columns, drop = FALSE]
  out[stats::complete.cases(out), , drop = FALSE]
}

download_text <- function(filename, text, mime = "text/plain") {
  shiny::downloadHandler(filename = function() filename,
    content = function(file) writeLines(text, file, useBytes = TRUE), contentType = mime)
}
