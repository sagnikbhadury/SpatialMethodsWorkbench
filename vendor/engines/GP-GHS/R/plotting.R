#' Plot a spatially varying edge-strength map for one directed edge
#'
#' @param result Output of \code{\link{gp_group_horseshoe_graph}}.
#' @param node_s Source cell type name.
#' @param node_t Target (neighbor) cell type name.
#' @param plot_sd Currently unused; reserved for plotting posterior SD instead
#'   of the posterior mean.
#'
#' @return A \pkg{ggplot2} object if the \pkg{ggplot2} package is available,
#'   otherwise a base R plot is drawn and \code{NULL} is returned invisibly.
#' @export
plot_edge_map <- function(result, node_s, node_t, plot_sd = FALSE) {
  cell_types <- rownames(result$adj)
  coords     <- result$coords

  s_idx   <- which(cell_types == node_s)
  t_idx   <- which(cell_types == node_t)
  nbr_pos <- which((1:length(cell_types))[-s_idx] == t_idx)

  edge_strength <- result$edge_maps[[node_s]][, nbr_pos]
  is_active     <- result$adj[node_s, node_t] == 1

  df <- data.frame(x = coords[, 1], y = coords[, 2],
                   strength = edge_strength)

  title_str <- paste0("Edge: ", node_s, " -> ", node_t,
                      ifelse(is_active, "  [ACTIVE]", "  [shrunk]"))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = strength)) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::scale_color_gradient2(low = "steelblue", mid = "white",
                            high = "firebrick", midpoint = 0,
                            name = "Edge\nStrength") +
      ggplot2::labs(title = title_str, x = "X", y = "Y") +
      ggplot2::theme_bw(base_size = 12)
  } else {
    pal <- colorRampPalette(c("steelblue", "white", "firebrick"))(100)
    br  <- cut(edge_strength, breaks = 100, labels = FALSE)
    plot(coords, col = pal[br], pch = 16, cex = 0.8,
         main = title_str, xlab = "X", ylab = "Y")
    invisible(NULL)
  }
}

#' Plot all edge maps for one node as a facetted panel
#'
#' Requires the \pkg{ggplot2} and \pkg{tidyr} packages.
#'
#' @param result Output of \code{\link{gp_group_horseshoe_graph}}.
#' @param node_s Source cell type name.
#'
#' @return A \pkg{ggplot2} object.
#' @export
plot_all_edges <- function(result, node_s) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("ggplot2 and tidyr required for panel plot.")
  }

  cell_types <- rownames(result$adj)
  s_idx      <- which(cell_types == node_s)
  nbr_names  <- cell_types[-s_idx]
  coords     <- result$coords
  maps       <- result$edge_maps[[node_s]]   # n x q
  active     <- result$adj[node_s, -s_idx]

  df <- as.data.frame(maps)
  colnames(df) <- paste0(nbr_names, ifelse(active == 1, "*", ""))
  df$x <- coords[, 1]; df$y <- coords[, 2]

  df_long <- tidyr::pivot_longer(df, -c(x, y), names_to = "neighbor",
                          values_to = "strength")

  ggplot2::ggplot(df_long, ggplot2::aes(x = x, y = y, color = strength)) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::facet_wrap(~neighbor) +
    ggplot2::scale_color_gradient2(low = "steelblue", mid = "white",
                          high = "firebrick", midpoint = 0) +
    ggplot2::labs(title = paste("All edge maps for:", node_s),
         subtitle = "* = active edge (AND rule)") +
    ggplot2::theme_bw(base_size = 10)
}

#' Plot the global adjacency heatmap with shrinkage signal
#'
#' Requires the \pkg{ggplot2} and \pkg{reshape2} packages; falls back to a
#' base R \code{\link[graphics]{image}} plot of the raw adjacency if either
#' is unavailable.
#'
#' @param result Output of \code{\link{gp_group_horseshoe_graph}}.
#'
#' @return A \pkg{ggplot2} object, or (invisibly) \code{NULL} if falling back
#'   to the base R plot.
#' @export
plot_adj_heatmap <- function(result) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("reshape2", quietly = TRUE)) {
    image(result$adj, main = "Adjacency"); return(invisible())
  }

  signal <- 1 - result$kappa
  diag(signal) <- NA

  df     <- reshape2::melt(signal, na.rm = TRUE)
  adj_df <- reshape2::melt(result$adj, na.rm = FALSE)
  colnames(df)     <- c("From", "To", "signal")
  colnames(adj_df) <- c("From", "To", "active_val")

  adj_df <- adj_df[adj_df$From != adj_df$To, ]
  df <- merge(df, adj_df, by = c("From", "To"))
  df$active <- df$active_val == 1

  ggplot2::ggplot(df, ggplot2::aes(x = To, y = From, fill = signal)) +
    ggplot2::geom_tile(color = "gray85") +
    ggplot2::geom_tile(data = subset(df, active),
              fill = NA, color = "black", linewidth = 1.2) +
    ggplot2::scale_fill_gradient(low = "white", high = "darkred",
                        name = "1 - kappa\n(signal)") +
    ggplot2::labs(title = "Cell-Cell Interaction Graph",
         subtitle = "Black border = active edge | Color = signal strength") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot the HSGP spectral density diagnostic
#'
#' @param result Output of \code{\link{gp_group_horseshoe_graph}}.
#'
#' @return A \pkg{ggplot2} object if the \pkg{ggplot2} package is available,
#'   otherwise \code{NULL} invisibly.
#' @export
plot_spectral_density <- function(result) {
  df <- data.frame(
    basis_idx = 1:length(result$S_diag),
    S         = result$S_diag
  )
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggplot(df, ggplot2::aes(x = basis_idx, y = S)) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.8) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(title = paste0("Matern(nu=", result$hsgp_params$nu,
                          ") Spectral Densities"),
           subtitle = "High S = smooth (low-freq) basis allowed | Low S = rough basis suppressed",
           x = "Basis function index (ordered by frequency)",
           y = "S(omega_k) [log scale]") +
      ggplot2::theme_bw(base_size = 12)
  } else {
    invisible(NULL)
  }
}
