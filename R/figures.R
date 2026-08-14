utils::globalVariables(c("xend", "yend", "display_weight", "direction", "degree", "node"))

network_figure <- function(matrices, threshold = 0, title = "Estimated network",
                           top_n = 30L, max_networks = 6L, signed = TRUE,
                           weight_label = "Association strength") {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Network figures require the igraph package.")
  }
  if (!length(matrices)) stop("No network matrices are available to plot.")
  if (is.null(names(matrices))) names(matrices) <- paste("Network", seq_along(matrices))
  matrices <- matrices[seq_len(min(length(matrices), as.integer(max_networks)))]

  panels <- Map(function(mat, panel_name) {
    mat <- as.matrix(mat)
    node_names <- rownames(mat) %||% paste0("V", seq_len(nrow(mat)))
    rownames(mat) <- colnames(mat) <- node_names
    edges <- matrix_edges(mat, panel_name, threshold)
    edges <- edges[is.finite(edges$weight) & edges$weight != 0, , drop = FALSE]
    if (nrow(edges)) {
      edges <- edges[order(abs(edges$weight), decreasing = TRUE), , drop = FALSE]
      edges <- utils::head(edges, as.integer(top_n))
    }

    graph_edges <- edges[, c("source", "target"), drop = FALSE]
    graph <- igraph::graph_from_data_frame(graph_edges, directed = FALSE,
                                           vertices = data.frame(name = node_names))
    set.seed(2026)
    layout <- if (igraph::ecount(graph)) {
      igraph::layout_with_fr(graph, weights = pmax(abs(edges$weight), 1e-6))
    } else {
      igraph::layout_in_circle(graph)
    }
    layout <- igraph::norm_coords(layout, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
    nodes <- data.frame(network = panel_name, node = node_names,
                        x = layout[, 1], y = layout[, 2],
                        degree = igraph::degree(graph), stringsAsFactors = FALSE)
    if (!nrow(edges)) return(list(nodes = nodes, edges = data.frame()))

    from <- match(edges$source, nodes$node); to <- match(edges$target, nodes$node)
    edges$x <- nodes$x[from]; edges$y <- nodes$y[from]
    edges$xend <- nodes$x[to]; edges$yend <- nodes$y[to]
    edges$display_weight <- abs(edges$weight)
    edges$direction <- if (signed) ifelse(edges$weight >= 0, "Positive", "Negative") else "Selected"
    list(nodes = nodes, edges = edges)
  }, matrices, names(matrices))

  nodes <- do.call(rbind, lapply(panels, `[[`, "nodes"))
  edge_parts <- lapply(panels, `[[`, "edges")
  edge_parts <- edge_parts[vapply(edge_parts, nrow, integer(1)) > 0L]
  edges <- if (length(edge_parts)) do.call(rbind, edge_parts) else NULL

  plot <- ggplot2::ggplot()
  if (!is.null(edges)) {
    plot <- plot + ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                   linewidth = display_weight, color = direction),
      alpha = .76, lineend = "round"
    )
  }
  plot <- plot +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x, y, size = degree),
                        shape = 21, fill = "#f6f3ec", color = "#263a40", stroke = .7) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y, label = node),
                       size = 3.1, fontface = "bold", check_overlap = TRUE, vjust = -1.15) +
    ggplot2::facet_wrap(~network) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "Displayed degree") +
    ggplot2::scale_linewidth_continuous(range = c(.4, 2.7), name = weight_label) +
    ggplot2::scale_color_manual(values = c(Positive = "#d6604d", Negative = "#2166ac", Selected = "#176b68"),
                                name = if (signed) "Direction" else NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(color = "#52636a"),
      panel.spacing = grid::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(16, 30, 16, 30),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("Up to %d strongest edges per network meeting the configured threshold; node size is displayed degree.", as.integer(top_n))
    )
  plot
}

gpghs_coefficient_matrix <- function(fit) {
  adjacency <- as.matrix(fit$adj)
  coefficient_rms <- matrix(0, nrow(adjacency), ncol(adjacency), dimnames = dimnames(adjacency))
  selected <- which(upper.tri(adjacency) & adjacency != 0, arr.ind = TRUE)
  if (!nrow(selected)) return(coefficient_rms)
  directional_field <- function(source, target) {
    neighbors <- setdiff(seq_len(nrow(adjacency)), source)
    position <- match(target, neighbors)
    fields <- fit$edge_maps[[source]]
    if (is.null(fields) || is.na(position) || position > ncol(fields)) return(numeric())
    as.numeric(fields[, position])
  }
  rms <- function(field) {
    field <- field[is.finite(field)]
    if (!length(field)) return(NA_real_)
    sqrt(mean(field^2))
  }
  for (index in seq_len(nrow(selected))) {
    i <- selected[index, 1]; j <- selected[index, 2]
    directional_rms <- c(rms(directional_field(i, j)), rms(directional_field(j, i)))
    directional_rms <- directional_rms[is.finite(directional_rms)]
    value <- if (length(directional_rms)) min(directional_rms) else 0
    coefficient_rms[i, j] <- coefficient_rms[j, i] <- value
  }
  coefficient_rms
}
