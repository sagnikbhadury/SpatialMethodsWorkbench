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

funding_acknowledgement <- function() {
  paste(
    "Research reported in this software was supported by the National Cancer Institute",
    "of the National Institutes of Health under Award Number R37CA214955-01A1.",
    "The content is solely the responsibility of the author and does not necessarily",
    "represent the official views of the National Institutes of Health."
  )
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

demo_complete_data <- function(contours = 8L, points = 16L, image_features = 6L,
                               seed = 2026L) {
  contours <- max(6L, as.integer(contours))
  points <- max(10L, as.integer(points))
  image_features <- max(5L, min(as.integer(image_features), 20L))
  set.seed(seed)

  theta <- rep(seq(0, 2 * pi, length.out = points + 1L)[-1L], contours)
  contour_index <- rep(seq_len(contours), each = points)
  point_order <- rep(seq_len(points), contours)
  contour_id <- sprintf("C%02d", contour_index)
  volume_index <- 1L + (contour_index > ceiling(contours / 2L))
  volume_id <- sprintf("V%02d", volume_index)
  center_x <- rep(seq(-3, 3, length.out = contours), each = points)
  center_y <- rep(c(-1.8, 1.8), length.out = contours)[contour_index]
  radial <- 1 + .07 * sin(3 * theta + contour_index / 7) +
    stats::rnorm(length(theta), sd = .015)
  x <- center_x + radial * (1 + contour_index / 45) * cos(theta)
  y <- center_y + radial * (1 - contour_index / 55) * sin(theta)
  z <- (volume_index - 1L) * 1.2 + point_order / points * .25

  exposure <- scale(sin(x / 2) + cos(y / 2.5) + .2 * z +
                      stats::rnorm(length(x), sd = .45))[, 1]
  mediator <- scale(.6 * exposure + sin(y / 2) +
                      stats::rnorm(length(x), sd = .65))[, 1]
  outcome <- scale(.45 * exposure + .55 * mediator - .15 * z +
                     stats::rnorm(length(x), sd = .6))[, 1]
  epithelial <- scale(exposure + stats::rnorm(length(x), sd = .35))[, 1]
  apc <- scale(mediator + stats::rnorm(length(x), sd = .4))[, 1]
  ctl <- scale(outcome - .25 * exposure + stats::rnorm(length(x), sd = .5))[, 1]
  treg <- scale(-.35 * ctl + .3 * epithelial + stats::rnorm(length(x), sd = .65))[, 1]
  helper <- scale(.4 * apc + .25 * treg + stats::rnorm(length(x), sd = .65))[, 1]
  disease <- ifelse(outcome + stats::rnorm(length(x), sd = .45) > 0,
                    "Condition_B", "Condition_A")
  spatial_score <- x + y + z
  zone <- cut(spatial_score,
    breaks = stats::quantile(spatial_score, probs = seq(0, 1, length.out = 4)),
    include.lowest = TRUE, labels = c("Low", "Intermediate", "High")
  )

  input <- matrix(stats::rnorm(length(x) * image_features),
                  nrow = length(x), ncol = image_features)
  input[, 1] <- scale(epithelial + stats::rnorm(length(x), sd = .35))[, 1]
  input[, 2] <- scale(apc + stats::rnorm(length(x), sd = .35))[, 1]
  loading <- matrix(stats::rnorm(image_features^2, sd = .2),
                    image_features, image_features)
  output <- input %*% loading +
    matrix(stats::rnorm(length(x) * image_features, sd = .6),
           nrow = length(x), ncol = image_features)

  result <- data.frame(
    volume_id = volume_id, contour_id = contour_id,
    sample_id = sprintf("S%03d", seq_along(x)), order = point_order,
    x = x, y = y, z = z, zone = zone, disease = disease,
    outcome = outcome, exposure = exposure, mediator = mediator,
    epithelial = epithelial, APC = apc, CTL = ctl, Treg = treg,
    T_helper = helper, check.names = FALSE
  )
  result[paste0("pixel_", seq_len(image_features))] <- input
  result[paste0("input__", seq_len(image_features))] <- input
  result[paste0("output__", seq_len(image_features))] <- output
  result
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
