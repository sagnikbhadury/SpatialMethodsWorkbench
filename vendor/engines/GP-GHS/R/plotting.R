#' Plot one coefficient function over the tissue
#'
#' Shows the posterior mean coefficient of predictor `to` in the regression
#' of response `from` at every location, with a diverging scale centered at
#' zero. Uses ggplot2 when it is installed and base graphics otherwise.
#'
#' @param fit A `"gpghs"` object.
#' @param from,to Node names.
#' @param plot_sd Plot the posterior standard deviation instead of the mean.
#' @param point_size Point size.
#' @return A ggplot object, or `NULL` (base graphics) invisibly.
#' @export
plot_edge_map <- function(fit, from, to, plot_sd = FALSE, point_size = 1.5) {
  cm <- coefficient_map(fit, from, to)
  val <- if (plot_sd) cm$sd else cm$mean
  active <- fit$adj[from, to] == 1
  ttl <- paste0(to, " in the regression of ", from, if (active) " (selected)" else " (not selected)")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    df <- data.frame(x = cm$x, y = cm$y, value = val)
    g <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = value)) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::labs(title = ttl, x = NULL, y = NULL, color = if (plot_sd) "posterior SD" else "coefficient") +
      ggplot2::coord_equal() + ggplot2::theme_bw(base_size = 11)
    if (plot_sd) g <- g + ggplot2::scale_color_viridis_c()
    else g <- g + ggplot2::scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0)
    return(g)
  }
  pal <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
  r <- max(abs(val)); br <- 1 + round(50 * (val / ifelse(r > 0, r, 1) + 1))
  graphics::plot(cm$x, cm$y, col = pal[br], pch = 16, cex = 0.8, main = ttl, xlab = "", ylab = "", asp = 1)
  invisible(NULL)
}

#' Plot the estimated network
#'
#' Draws the selected edges with igraph: solid for positive edges, dashed
#' for negative and dotted for mixed (the convention of the article), with
#' node size increasing with degree.
#'
#' @param fit A `"gpghs"` object.
#' @param layout An igraph layout function (default
#'   `igraph::layout_with_fr`).
#' @param label_top Number of highest-degree nodes to label (`Inf` labels
#'   all).
#' @param ... Passed to `igraph::plot.igraph`.
#' @return The igraph object, invisibly.
#' @export
plot_network <- function(fit, layout = NULL, label_top = Inf, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) stop("plot_network needs the igraph package")
  ed <- gpghs_edges(fit)
  g <- igraph::graph_from_data_frame(if (nrow(ed)) ed[, c("from", "to", "sign", "score")] else
                                       data.frame(from = character(), to = character()),
                                     directed = FALSE,
                                     vertices = data.frame(name = rownames(fit$adj)[fit$estimable]))
  deg <- igraph::degree(g)
  lab <- names(deg); if (is.finite(label_top)) lab[rank(-deg, ties.method = "first") > label_top] <- ""
  lty <- if (nrow(ed)) c(positive = 1, negative = 2, mixed = 3)[igraph::E(g)$sign] else 1
  if (is.null(layout)) layout <- igraph::layout_with_fr
  igraph::plot.igraph(g, layout = layout, vertex.size = 3 + 2 * sqrt(deg), vertex.label = lab,
                      vertex.label.cex = 0.7, vertex.color = "grey80", vertex.frame.color = "grey40",
                      edge.lty = lty, edge.color = ifelse(lty == 1, "#b2182b", ifelse(lty == 2, "#2166ac", "grey50")),
                      edge.width = 1 + 2 * (if (nrow(ed)) igraph::E(g)$score else 1), ...)
  invisible(g)
}

#' Heatmap of pair scores with the selected edges outlined
#'
#' @param fit A `"gpghs"` object.
#' @return A ggplot object (ggplot2 required).
#' @export
plot_adj_heatmap <- function(fit) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("plot_adj_heatmap needs ggplot2")
  nodes <- rownames(fit$adj)
  df <- expand.grid(from = nodes, to = nodes, stringsAsFactors = FALSE)
  df$score <- as.vector(fit$score); df$selected <- as.vector(fit$adj) == 1
  df <- df[df$from != df$to, ]
  ggplot2::ggplot(df, ggplot2::aes(x = to, y = from, fill = score)) +
    ggplot2::geom_tile(color = "grey85") +
    ggplot2::geom_tile(data = df[df$selected, ], fill = NA, color = "black", linewidth = 0.8) +
    ggplot2::scale_fill_gradient(low = "white", high = "#67000d", na.value = "grey95", name = "1 - max kappa") +
    ggplot2::labs(x = NULL, y = NULL, title = "Pair scores; outlined cells are selected edges") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

if (getRversion() >= "2.15.1") utils::globalVariables(c("x", "y", "value", "to", "from", "score"))
