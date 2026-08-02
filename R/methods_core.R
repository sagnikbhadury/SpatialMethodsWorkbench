knn_weights <- function(coords, k = 8L) {
  coords <- as.matrix(coords)
  d <- as.matrix(stats::dist(coords))
  diag(d) <- Inf
  k <- max(1L, min(as.integer(k), nrow(coords) - 1L))
  W <- matrix(0, nrow(coords), nrow(coords))
  for (i in seq_len(nrow(coords))) W[i, order(d[i, ])[seq_len(k)]] <- 1
  (W + t(W) > 0) * 1
}

moran_permutation <- function(values, coords, k = 8L, permutations = 199L, seed = 1L) {
  keep <- stats::complete.cases(values, coords)
  x <- as.numeric(values[keep]); xy <- as.matrix(coords[keep, , drop = FALSE])
  if (length(x) < 10L || stats::sd(x) == 0) stop("Moran's I requires at least 10 complete, non-constant observations.")
  if (length(x) > 4000L) { set.seed(seed); take <- sample(seq_along(x), 4000L); x <- x[take]; xy <- xy[take, , drop = FALSE] }
  W <- knn_weights(xy, k)
  z <- x - mean(x); scale_factor <- length(x) / sum(W)
  statistic <- scale_factor * sum(W * tcrossprod(z)) / sum(z^2)
  set.seed(seed)
  null <- replicate(as.integer(permutations), {
    zp <- sample(z); scale_factor * sum(W * tcrossprod(zp)) / sum(zp^2)
  })
  p <- (1 + sum(abs(null) >= abs(statistic))) / (length(null) + 1)
  list(statistic = statistic, p.value = p, null = null, n = length(x), k = k)
}

run_spatial_qc <- function(data, mapping, features, params) {
  feature <- params$feature %||% features[1]
  clean <- sanitize_for_analysis(data, c(mapping$x, mapping$y, feature))
  names(clean) <- c("x", "y", "value")
  test <- moran_permutation(clean$value, clean[, c("x", "y")],
                            k = params$neighbors %||% 8L,
                            permutations = params$permutations %||% 199L,
                            seed = params$seed %||% 1L)
  plot <- ggplot2::ggplot(clean, ggplot2::aes(x, y, color = value)) +
    ggplot2::geom_point(size = 2.2, alpha = .85) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_viridis_c(option = "C") +
    ggplot2::labs(title = paste("Spatial distribution of", feature), color = feature,
                  subtitle = sprintf("Moran's I = %.3f; permutation p = %.4f", test$statistic, test$p.value)) +
    ggplot2::theme_minimal(base_size = 13)
  list(method = "Spatial exploration & autocorrelation", plot = plot,
       table = data.frame(metric = c("Complete observations", "Neighbors", "Moran's I", "Permutation p-value"),
                          value = c(test$n, test$k, test$statistic, test$p.value)),
       summary = c(feature = feature, moran_i = test$statistic, p_value = test$p.value),
       notes = c("Moran's I measures global spatial autocorrelation under the selected k-nearest-neighbor graph.",
                 "The permutation p-value is descriptive unless the sampling design supports the exchangeability assumption."))
}

partial_correlation <- function(X, ridge = .05) {
  X <- as.matrix(X)
  X <- X[, apply(X, 2, stats::sd, na.rm = TRUE) > 0, drop = FALSE]
  S <- stats::cov(X, use = "pairwise.complete.obs")
  ridge_value <- ridge * mean(diag(S), na.rm = TRUE)
  precision <- solve(S + diag(ridge_value + 1e-8, ncol(S)))
  pcor <- -precision / sqrt(outer(diag(precision), diag(precision)))
  diag(pcor) <- 1
  pcor
}

matrix_edges <- function(mat, network = "Overall", threshold = 0) {
  idx <- which(upper.tri(mat) & abs(mat) >= threshold, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame(network = character(), source = character(), target = character(), weight = numeric()))
  data.frame(network = network, source = rownames(mat)[idx[, 1]], target = colnames(mat)[idx[, 2]],
             weight = mat[idx], sign = ifelse(mat[idx] >= 0, "Positive", "Negative"), row.names = NULL)
}

run_conditional_network <- function(data, mapping, features, params) {
  threshold <- params$threshold %||% .15
  ridge <- params$ridge %||% .05
  group_col <- mapping$zone
  groups <- if (!is.null(group_col) && isTRUE(params$by_zone)) unique(as.character(data[[group_col]])) else "Overall"
  matrices <- list(); edges <- list()
  for (group in groups) {
    subset <- if (identical(group, "Overall")) data else data[as.character(data[[group_col]]) == group, , drop = FALSE]
    X <- subset[, features, drop = FALSE]
    X <- X[stats::complete.cases(X), , drop = FALSE]
    if (nrow(X) <= length(features) + 2L) next
    mat <- partial_correlation(X, ridge)
    matrices[[group]] <- mat
    edges[[group]] <- matrix_edges(mat, group, threshold)
  }
  if (!length(matrices)) stop("No group has enough complete observations for the selected features.")
  edge_table <- do.call(rbind, edges)
  first <- matrices[[1]]
  heat <- expand.grid(row = rownames(first), column = colnames(first), stringsAsFactors = FALSE)
  heat$weight <- as.vector(first)
  plot <- ggplot2::ggplot(heat, ggplot2::aes(column, row, fill = weight)) +
    ggplot2::geom_tile(color = "white", linewidth = .3) +
    ggplot2::scale_fill_gradient2(low = "#d46245", mid = "#fffdf8", high = "#176b68", limits = c(-1, 1)) +
    ggplot2::coord_equal() + ggplot2::labs(title = paste(names(matrices)[1], "partial-correlation network"), x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  list(method = "Region-aware conditional network", plot = plot, table = edge_table,
       summary = c(networks = length(matrices), visible_edges = nrow(edge_table), threshold = threshold, ridge = ridge),
       matrices = matrices,
       notes = c("This fast exploratory workflow uses ridge-regularized precision matrices.",
                 "Use the Bayesian ISPat modules when their full input contract and server engine are available.",
                 "Conditional association does not establish physical contact, signaling, or causation."))
}

run_spatial_ml <- function(data, mapping, features, params) {
  outcome <- mapping$outcome
  predictors <- unique(c(features, mapping$x, mapping$y, mapping$z))
  clean <- sanitize_for_analysis(data, c(outcome, predictors))
  y <- clean[[outcome]]; X <- clean[, predictors, drop = FALSE]
  if (nrow(clean) < 30L) stop("Prediction requires at least 30 complete observations.")
  set.seed(params$seed %||% 1L)
  train_id <- sample(seq_len(nrow(clean)), floor(.75 * nrow(clean)))
  train <- clean[train_id, , drop = FALSE]; test <- clean[-train_id, , drop = FALSE]
  formula <- stats::reformulate(predictors, response = outcome)
  if (is.numeric(y)) {
    fit <- stats::lm(formula, data = train)
    prediction <- stats::predict(fit, newdata = test)
    observed <- test[[outcome]]
    metrics <- data.frame(metric = c("RMSE", "MAE", "R-squared"),
                          value = c(sqrt(mean((observed - prediction)^2)), mean(abs(observed - prediction)),
                                    stats::cor(observed, prediction)^2))
    table <- data.frame(observed = observed, predicted = prediction, residual = observed - prediction)
    plot <- ggplot2::ggplot(table, ggplot2::aes(observed, predicted)) + ggplot2::geom_point(alpha = .7, color = "#176b68") +
      ggplot2::geom_abline(linetype = 2, color = "#d46245") + ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(title = "Held-out predictions", subtitle = "25% test split")
  } else {
    y <- factor(y)
    if (nlevels(y) != 2L) stop("The current classification workflow supports a binary outcome.")
    train[[outcome]] <- factor(train[[outcome]], levels = levels(y)); test[[outcome]] <- factor(test[[outcome]], levels = levels(y))
    fit <- stats::glm(formula, data = train, family = stats::binomial())
    probability <- stats::predict(fit, newdata = test, type = "response")
    predicted <- factor(ifelse(probability >= .5, levels(y)[2], levels(y)[1]), levels = levels(y))
    observed <- test[[outcome]]
    accuracy <- mean(predicted == observed)
    metrics <- data.frame(metric = c("Accuracy", "Test observations"), value = c(accuracy, length(observed)))
    table <- data.frame(observed = observed, predicted = predicted, probability = probability)
    plot <- ggplot2::ggplot(table, ggplot2::aes(probability, fill = observed)) +
      ggplot2::geom_density(alpha = .45) + ggplot2::geom_vline(xintercept = .5, linetype = 2) +
      ggplot2::theme_minimal(base_size = 13) + ggplot2::labs(title = "Held-out class probabilities", fill = "Observed")
  }
  list(method = "Spatial machine-learning prediction", plot = plot, table = table, summary = stats::setNames(metrics$value, metrics$metric),
       model = fit, notes = c("This module uses a transparent generalized linear model with a held-out test split.",
                              "Spatial leakage can inflate performance; use spatially blocked validation for confirmatory studies."))
}

run_spatial_clustering <- function(data, mapping, features, params) {
  columns <- c(mapping$x, mapping$y, features)
  clean <- sanitize_for_analysis(data, columns)
  if (nrow(clean) < 20L) stop("Spatial clustering requires at least 20 complete observations.")
  k <- as.integer(params$clusters %||% 4L)
  if (k < 2L || k >= nrow(clean)) stop("Choose between 2 clusters and one fewer than the number of observations.")
  marker_matrix <- scale(as.matrix(clean[, features, drop = FALSE]))
  coord_matrix <- scale(as.matrix(clean[, c(mapping$x, mapping$y), drop = FALSE]))
  coord_weight <- as.numeric(params$spatial_weight %||% .5)
  design <- cbind(marker_matrix, coord_weight * coord_matrix)
  design[!is.finite(design)] <- 0
  set.seed(params$seed %||% 1L)
  fit <- stats::kmeans(design, centers = k, nstart = 25, iter.max = 100)
  table <- data.frame(x = clean[[mapping$x]], y = clean[[mapping$y]], cluster = factor(fit$cluster))
  plot <- ggplot2::ggplot(table, ggplot2::aes(x, y, color = cluster)) +
    ggplot2::geom_point(size = 2.2, alpha = .82) + ggplot2::coord_equal() +
    ggplot2::scale_color_brewer(palette = "Dark2") + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = "Spatially weighted phenotype clusters", color = "Cluster")
  list(method = "Spatially weighted phenotype clustering", plot = plot, table = table,
       summary = c(observations = nrow(clean), clusters = k, spatial_weight = coord_weight,
                   between_total_ratio = fit$betweenss / fit$totss), model = fit,
       notes = c("Coordinates and markers are standardized before clustering; the spatial weight controls their relative influence.",
                 "Clusters are exploratory partitions and are not validated biological cell states or clinical subtypes."))
}

standardize_train_test <- function(train, test) {
  center <- colMeans(train)
  scale_value <- apply(train, 2, stats::sd)
  scale_value[!is.finite(scale_value) | scale_value == 0] <- 1
  list(train = sweep(sweep(train, 2, center, "-"), 2, scale_value, "/"),
       test = sweep(sweep(test, 2, center, "-"), 2, scale_value, "/"),
       center = center, scale = scale_value)
}

run_neural_prediction <- function(data, mapping, features, params) {
  if (!requireNamespace("nnet", quietly = TRUE)) stop("The nnet engine is not installed on this deployment.")
  predictors <- unique(c(features, mapping$x, mapping$y, mapping$z))
  clean <- sanitize_for_analysis(data, c(mapping$outcome, predictors))
  if (nrow(clean) < 40L) stop("Neural prediction requires at least 40 complete observations.")
  if (length(predictors) > 250L) stop("Select no more than 250 predictors for this shallow-network workflow.")
  set.seed(params$seed %||% 1L)
  train_id <- sample(seq_len(nrow(clean)), floor(.75 * nrow(clean)))
  test_id <- setdiff(seq_len(nrow(clean)), train_id)
  scaled <- standardize_train_test(as.matrix(clean[train_id, predictors, drop = FALSE]),
                                   as.matrix(clean[test_id, predictors, drop = FALSE]))
  outcome <- clean[[mapping$outcome]]
  size <- as.integer(params$hidden_units %||% 5L)
  decay <- as.numeric(params$decay %||% .01)
  if (is.numeric(outcome)) {
    y_train <- outcome[train_id]; y_center <- mean(y_train); y_scale <- stats::sd(y_train)
    if (!is.finite(y_scale) || y_scale == 0) stop("The numeric outcome must vary in the training set.")
    fit <- nnet::nnet(scaled$train, (y_train - y_center) / y_scale, size = size, linout = TRUE,
                      decay = decay, maxit = 500, trace = FALSE, MaxNWts = 20000)
    prediction <- as.numeric(stats::predict(fit, scaled$test)) * y_scale + y_center
    observed <- outcome[test_id]
    table <- data.frame(observed = observed, predicted = prediction, residual = observed - prediction)
    metrics <- c(rmse = sqrt(mean((observed - prediction)^2)), mae = mean(abs(observed - prediction)),
                 r_squared = suppressWarnings(stats::cor(observed, prediction)^2))
    plot <- ggplot2::ggplot(table, ggplot2::aes(observed, predicted)) + ggplot2::geom_point(alpha = .7, color = "#176b68") +
      ggplot2::geom_abline(linetype = 2, color = "#d46245") + ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(title = "Held-out shallow-network predictions")
  } else {
    outcome <- factor(outcome)
    if (nlevels(outcome) != 2L) stop("Neural classification currently supports a binary outcome.")
    target <- as.numeric(outcome[train_id] == levels(outcome)[2])
    fit <- nnet::nnet(scaled$train, target, size = size, entropy = TRUE, decay = decay,
                      maxit = 500, trace = FALSE, MaxNWts = 20000)
    probability <- as.numeric(stats::predict(fit, scaled$test))
    predicted <- factor(ifelse(probability >= .5, levels(outcome)[2], levels(outcome)[1]), levels = levels(outcome))
    observed <- outcome[test_id]
    table <- data.frame(observed = observed, predicted = predicted, probability = probability)
    metrics <- c(accuracy = mean(observed == predicted), test_observations = length(test_id))
    plot <- ggplot2::ggplot(table, ggplot2::aes(probability, fill = observed)) + ggplot2::geom_density(alpha = .45) +
      ggplot2::geom_vline(xintercept = .5, linetype = 2) + ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(title = "Held-out shallow-network probabilities", fill = "Observed")
  }
  list(method = "Spatial shallow-neural prediction", plot = plot, table = table,
       summary = c(metrics, hidden_units = size, weight_decay = decay), model = fit,
       notes = c("This is a single-hidden-layer predictive workflow with standardized inputs and a held-out test set.",
                 "It is an applied spatial adaptation of public AI/deep-learning curricula, not a reproduction of their tutorial examples.",
                 "Use spatially blocked or external validation before confirmatory use."))
}

ridge_coefficients <- function(X, y, lambda) {
  if (ncol(X) <= nrow(X)) {
    solve(crossprod(X) + diag(lambda, ncol(X)), crossprod(X, y))
  } else {
    crossprod(X, solve(tcrossprod(X) + diag(lambda, nrow(X)), y))
  }
}

run_scalar_image <- function(data, mapping, features, params) {
  clean <- sanitize_for_analysis(data, c(mapping$outcome, features))
  outcome <- clean[[mapping$outcome]]
  if (!is.numeric(outcome)) stop("Wide-image scalar regression currently requires a numeric outcome.")
  if (nrow(clean) < 30L) stop("Wide-image regression requires at least 30 complete subjects.")
  if (length(features) > 5000L) stop("Select no more than 5,000 image features in the interactive workflow.")
  set.seed(params$seed %||% 1L)
  train_id <- sample(seq_len(nrow(clean)), floor(.75 * nrow(clean)))
  test_id <- setdiff(seq_len(nrow(clean)), train_id)
  scaled <- standardize_train_test(as.matrix(clean[train_id, features, drop = FALSE]),
                                   as.matrix(clean[test_id, features, drop = FALSE]))
  y_center <- mean(outcome[train_id]); y_train <- outcome[train_id] - y_center
  lambda <- as.numeric(params$image_ridge %||% 10)
  beta <- ridge_coefficients(scaled$train, y_train, lambda)
  prediction <- as.numeric(scaled$test %*% beta + y_center)
  observed <- outcome[test_id]
  coefficients <- data.frame(feature = features, coefficient = as.numeric(beta))
  coefficients <- coefficients[order(abs(coefficients$coefficient), decreasing = TRUE), , drop = FALSE]
  top <- head(coefficients, 30L); top$feature <- factor(top$feature, levels = rev(top$feature))
  plot <- ggplot2::ggplot(top, ggplot2::aes(feature, coefficient, fill = coefficient > 0)) +
    ggplot2::geom_col(show.legend = FALSE) + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("#d46245", "#176b68")) + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = "Largest standardized image-feature coefficients", x = NULL)
  list(method = "Wide-image scalar regression screen", plot = plot, table = coefficients,
       summary = c(subjects = nrow(clean), image_features = length(features), ridge = lambda,
                   heldout_rmse = sqrt(mean((observed - prediction)^2)), heldout_mae = mean(abs(observed - prediction))),
       model = list(coefficients = beta, feature_center = scaled$center, feature_scale = scaled$scale, outcome_center = y_center),
       notes = c("Each selected numeric column is treated as an aligned pixel, voxel, or image-derived feature.",
                 "This fast ridge screen is not the Bayesian SV-NN or ST-CAR posterior algorithm; their publications are included as related advanced methodology.",
                 "Image registration and common feature ordering must be established before upload."))
}

run_image_to_image <- function(data, mapping, features, params) {
  input_features <- features[startsWith(features, "input__")]
  output_features <- features[startsWith(features, "output__")]
  clean <- sanitize_for_analysis(data, c(input_features, output_features))
  if (nrow(clean) < 30L) stop("Image-to-image regression requires at least 30 complete subjects.")
  set.seed(params$seed %||% 1L)
  train_id <- sample(seq_len(nrow(clean)), floor(.75 * nrow(clean)))
  test_id <- setdiff(seq_len(nrow(clean)), train_id)
  X <- as.matrix(clean[, input_features, drop = FALSE]); Z <- as.matrix(clean[, output_features, drop = FALSE])
  x_scaled <- standardize_train_test(X[train_id, , drop = FALSE], X[test_id, , drop = FALSE])
  z_scaled <- standardize_train_test(Z[train_id, , drop = FALSE], Z[test_id, , drop = FALSE])
  components <- min(as.integer(params$latent_factors %||% 5L), nrow(x_scaled$train) - 2L,
                    ncol(x_scaled$train), ncol(z_scaled$train))
  x_pca <- stats::prcomp(x_scaled$train, center = FALSE, scale. = FALSE, rank. = components)
  z_pca <- stats::prcomp(z_scaled$train, center = FALSE, scale. = FALSE, rank. = components)
  x_train_scores <- x_pca$x[, seq_len(components), drop = FALSE]
  z_train_scores <- z_pca$x[, seq_len(components), drop = FALSE]
  coef <- qr.solve(cbind(1, x_train_scores), z_train_scores)
  x_test_scores <- x_scaled$test %*% x_pca$rotation[, seq_len(components), drop = FALSE]
  z_score_prediction <- cbind(1, x_test_scores) %*% coef
  z_scaled_prediction <- z_score_prediction %*% t(z_pca$rotation[, seq_len(components), drop = FALSE])
  prediction <- sweep(sweep(z_scaled_prediction, 2, z_scaled$scale, "*"), 2, z_scaled$center, "+")
  observed <- Z[test_id, , drop = FALSE]
  feature_rmse <- sqrt(colMeans((observed - prediction)^2))
  table <- data.frame(output_feature = output_features, heldout_rmse = feature_rmse)
  first <- data.frame(observed = observed[, 1], predicted = prediction[, 1])
  plot <- ggplot2::ggplot(first, ggplot2::aes(observed, predicted)) + ggplot2::geom_point(color = "#176b68", alpha = .75) +
    ggplot2::geom_abline(linetype = 2, color = "#d46245") + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = paste("Held-out prediction:", output_features[1]))
  list(method = "Latent image-to-image regression", plot = plot, table = table,
       summary = c(subjects = nrow(clean), input_features = length(input_features), output_features = length(output_features),
                   latent_factors = components, overall_heldout_rmse = sqrt(mean((observed - prediction)^2))),
       model = list(input_pca = x_pca, output_pca = z_pca, regression = coef),
       notes = c("Columns prefixed input__ and output__ are treated as aligned predictor- and outcome-image features.",
                 "This is a fast low-rank screening workflow, not the full Bayesian SBLF posterior engine.",
                 "All images must share registration, mask, resolution, and feature ordering before upload."))
}

run_mediation <- function(data, mapping, features, params) {
  columns <- c(mapping$exposure, mapping$mediator, mapping$outcome)
  clean <- sanitize_for_analysis(data, columns); names(clean) <- c("exposure", "mediator", "outcome")
  clean$exposure <- as.numeric(clean$exposure)
  if (!all(vapply(clean, is.numeric, logical(1)))) stop("The current mediation module requires numeric exposure, mediator, and outcome.")
  if (nrow(clean) < 30L) stop("Mediation requires at least 30 complete observations.")
  a_fit <- stats::lm(mediator ~ exposure, clean)
  b_fit <- stats::lm(outcome ~ exposure + mediator, clean)
  direct <- stats::coef(b_fit)["exposure"]; indirect <- stats::coef(a_fit)["exposure"] * stats::coef(b_fit)["mediator"]
  total <- stats::coef(stats::lm(outcome ~ exposure, clean))["exposure"]
  B <- as.integer(params$bootstrap %||% 300L); set.seed(params$seed %||% 1L)
  boot <- replicate(B, { d <- clean[sample(seq_len(nrow(clean)), replace = TRUE), ];
    stats::coef(stats::lm(mediator ~ exposure, d))["exposure"] * stats::coef(stats::lm(outcome ~ exposure + mediator, d))["mediator"] })
  ci <- stats::quantile(boot, c(.025, .975), na.rm = TRUE)
  table <- data.frame(effect = c("Indirect", "Direct", "Total", "Indirect CI lower", "Indirect CI upper"),
                      estimate = c(indirect, direct, total, ci[1], ci[2]))
  plot_data <- data.frame(indirect_effect = boot)
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(indirect_effect)) + ggplot2::geom_histogram(bins = 35, fill = "#176b68", color = "white") +
    ggplot2::geom_vline(xintercept = indirect, color = "#d46245", linewidth = 1) + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = "Bootstrap distribution of the indirect effect", subtitle = sprintf("Percentile 95%% CI: %.3f to %.3f", ci[1], ci[2]))
  list(method = "Mediation with sensitivity diagnostics", plot = plot, table = table,
       summary = c(indirect = indirect, direct = direct, total = total, ci_low = ci[1], ci_high = ci[2]),
       notes = c("These are model-based mediation estimates, not automatic causal effects.",
                 "Causal interpretation requires temporal ordering, consistency, positivity, and no unmeasured exposure–mediator, exposure–outcome, or mediator–outcome confounding."))
}

resample_contour <- function(x, y, landmarks = 40L) {
  x <- c(x, x[1]); y <- c(y, y[1]); segment <- sqrt(diff(x)^2 + diff(y)^2)
  distance <- c(0, cumsum(segment)); target <- seq(0, max(distance), length.out = landmarks + 1L)[-(landmarks + 1L)]
  xr <- stats::approx(distance, x, xout = target, ties = "ordered")$y
  yr <- stats::approx(distance, y, xout = target, ties = "ordered")$y
  xr <- xr - mean(xr); yr <- yr - mean(yr); scale <- sqrt(sum(xr^2 + yr^2))
  c(xr / scale, yr / scale)
}

run_shape_pca <- function(data, mapping, features, params) {
  clean <- sanitize_for_analysis(data, c(mapping$id, mapping$x, mapping$y))
  ids <- unique(clean[[mapping$id]])
  if (length(ids) < 5L) stop("Shape PCA requires at least five contours identified by the ID column.")
  landmarks <- as.integer(params$landmarks %||% 40L)
  shapes <- t(vapply(ids, function(id) { d <- clean[clean[[mapping$id]] == id, ];
    if (nrow(d) < 6L) stop("Every contour needs at least six ordered points.")
    resample_contour(d[[mapping$x]], d[[mapping$y]], landmarks) }, numeric(landmarks * 2L)))
  rownames(shapes) <- ids
  pca <- stats::prcomp(shapes, center = TRUE, scale. = FALSE)
  scores <- data.frame(id = ids, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  plot <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, label = id)) + ggplot2::geom_point(size = 2.5, color = "#176b68") +
    ggplot2::geom_text(check_overlap = TRUE, nudge_y = .02, size = 3) + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = "Aligned contour shape space", x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]))
  list(method = "Tumor contour / shape PCA", plot = plot, table = scores,
       summary = c(contours = length(ids), landmarks = landmarks, pc1_variance = variance[1], pc2_variance = variance[2]),
       model = pca, notes = c("Contours are centered, size-normalized, and resampled to equally spaced landmarks.",
                              "Input points must follow the contour boundary in order; rotation is not automatically standardized."))
}
