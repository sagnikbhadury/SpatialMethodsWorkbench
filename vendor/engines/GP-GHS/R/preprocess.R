#' Log-normalize spatial transcriptomics counts
#'
#' Scales the counts of every spot to `scale_to` total counts and takes
#' `log(1 + x)`, as in the gene application of the article (counts per
#' 10,000).
#'
#' @param counts A matrix of counts with genes in rows and spots in columns
#'   (dense or a `dgCMatrix`-like object that supports `colSums` and `%*%`
#'   coercion by `as.matrix`).
#' @param scale_to Library size after scaling.
#' @return A dense matrix with spots in rows and genes in columns.
#' @export
st_normalize <- function(counts, scale_to = 1e4) {
  counts <- as.matrix(counts)
  lib <- colSums(counts)
  if (any(lib == 0)) stop("some spots have zero counts; remove them first")
  t(log1p(sweep(counts, 2, lib, "/") * scale_to))
}

#' Spatial neighbor weights
#'
#' Symmetric binary weights of the graph that joins each location to its
#' `k` nearest neighbors.
#'
#' @param coords An `n x 2` matrix of locations.
#' @param k Number of nearest neighbors (6 for the hexagonal Visium array).
#' @return A dense `n x n` symmetric weight matrix.
#' @export
spatial_knn_weights <- function(coords, k = 6) {
  coords <- as.matrix(coords); n <- nrow(coords)
  D <- as.matrix(stats::dist(coords)); diag(D) <- Inf
  W <- matrix(0, n, n)
  for (i in seq_len(n)) W[i, order(D[i, ])[seq_len(k)]] <- 1
  (W + t(W)) / 2
}

#' Moran's I screen for spatially variable genes
#'
#' Computes Moran's I of every column of `E` on the neighbor graph `W`, a
#' one-sided p-value from the normal approximation, and Benjamini-Hochberg
#' adjusted values. This is the screen used in the article (6 nearest
#' neighbors, false discovery rate 0.05).
#'
#' @param E A matrix with locations in rows and genes in columns (for
#'   example the output of [st_normalize()]).
#' @param coords An `n x 2` matrix of locations, used when `W` is `NULL`.
#' @param W Optional weight matrix from [spatial_knn_weights()].
#' @param k Nearest neighbors used when `W` is `NULL`.
#' @return A data frame with `gene`, `moran_I`, `z`, `p` and `fdr`, in the
#'   column order of `E`.
#' @export
morans_i <- function(E, coords = NULL, W = NULL, k = 6) {
  E <- as.matrix(E); n <- nrow(E)
  if (is.null(W)) { if (is.null(coords)) stop("give coords or W"); W <- spatial_knn_weights(coords, k) }
  S0 <- sum(W); S1 <- 2 * sum(W^2); S2 <- sum((2 * rowSums(W))^2)
  Z <- scale(E); Z[!is.finite(Z)] <- 0
  I <- (n / S0) * colSums(Z * (W %*% Z)) / pmax(colSums(Z^2), 1e-12)
  EI <- -1 / (n - 1)
  VI <- (n^2 * S1 - n * S2 + 3 * S0^2) / ((n^2 - 1) * S0^2) - EI^2
  z <- (I - EI) / sqrt(VI)
  p <- stats::pnorm(z, lower.tail = FALSE)
  data.frame(gene = colnames(E), moran_I = as.numeric(I), z = as.numeric(z), p = p,
             fdr = stats::p.adjust(p, "BH"), stringsAsFactors = FALSE)
}

#' Select the genes for a network: spatially variable, then most different between zones
#'
#' Reproduces the three-step selection of the article: keep genes detected
#' at a minimum share of spots (and drop mitochondrial and ribosomal protein
#' genes by name), keep the spatially variable genes by Moran's I at a false
#' discovery rate, and rank them by a Kruskal-Wallis test across zones.
#' When `zone` is `NULL`, the genes are ranked by Moran's I instead.
#'
#' @param E Log-normalized expression, spots in rows and genes in columns.
#' @param coords An `n x 2` matrix of spot locations.
#' @param zone Optional vector of zone labels (one per spot).
#' @param n_genes Number of genes to keep.
#' @param min_detect Minimum share of spots at which a gene is nonzero.
#' @param fdr False discovery rate for the Moran's I screen.
#' @param drop_pattern Regular expression of gene names to drop before the
#'   screen (mitochondrial and ribosomal protein genes by default); `NULL`
#'   keeps all.
#' @param k Nearest neighbors for Moran's I.
#' @return A data frame with one row per candidate gene (`gene`, `detected`,
#'   `moran_I`, `moran_fdr`, `svg`, `kw_H`, `kw_log10p`, `rank`, `selected`),
#'   sorted by rank; the selected genes are the first `n_genes` rows.
#' @export
st_select_genes <- function(E, coords, zone = NULL, n_genes = 100, min_detect = 0.10,
                            fdr = 0.05, drop_pattern = "^MT-|^RP[LS]|^MRP[LS]", k = 6) {
  E <- as.matrix(E)
  det <- colMeans(E > 0)
  keep <- det >= min_detect
  if (!is.null(drop_pattern)) keep <- keep & !grepl(drop_pattern, colnames(E))
  E <- E[, keep, drop = FALSE]; det <- det[keep]
  mi <- morans_i(E, coords = coords, k = k)
  svg <- mi$fdr < fdr
  tab <- data.frame(gene = colnames(E), detected = det, moran_I = mi$moran_I, moran_fdr = mi$fdr,
                    svg = svg, kw_H = NA_real_, kw_log10p = NA_real_, stringsAsFactors = FALSE)
  if (!is.null(zone)) {
    g <- as.factor(zone); ng <- as.numeric(table(g)); K <- nlevels(g); n <- nrow(E)
    kw <- apply(E[, svg, drop = FALSE], 2, function(x) {
      r <- rank(x); tt <- table(x)
      H <- 12 / (n * (n + 1)) * sum(tapply(r, g, sum)^2 / ng) - 3 * (n + 1)
      H / (1 - sum(tt^3 - tt) / (n^3 - n))
    })
    tab$kw_H[svg] <- kw
    tab$kw_log10p[svg] <- stats::pchisq(kw, K - 1, lower.tail = FALSE, log.p = TRUE) / log(10)
    tab <- tab[order(!tab$svg, tab$kw_log10p, -tab$kw_H), ]
  } else {
    tab <- tab[order(!tab$svg, -tab$moran_I), ]
  }
  tab$rank <- ifelse(tab$svg, seq_len(nrow(tab)), NA)
  tab$selected <- !is.na(tab$rank) & tab$rank <= n_genes
  rownames(tab) <- NULL
  tab
}

#' Quadrat counts of cell types in a multiplexed image
#'
#' Turns a point pattern of phenotyped cells into the local abundance fields
#' used in the cell-type application of the article: the bounding box of the
#' cells is divided into an `n_grid x n_grid` grid of equal quadrats, the
#' cells of each type are counted in every quadrat, the counts are square-root
#' transformed and standardized within the image, and a type with fewer than
#' `min_cells` cells is marked as not estimable and set to zero.
#'
#' @param x,y Cell coordinates.
#' @param type Cell type of each cell (character or factor).
#' @param n_grid Quadrats per side (15 in the article).
#' @param min_cells Minimum number of cells of a type in the image for the
#'   type to be estimable (5).
#' @param types Optional vector of the types to include, in the order wanted
#'   for the columns; by default all types present.
#' @return A list with `Y` (the `n_grid^2 x p` standardized abundance matrix,
#'   square-root counts), `counts` (raw counts), `coords` (quadrat centers),
#'   `estimable` and `n_cells` (cells per type).
#' @examples
#' cells <- data.frame(x = runif(2000, 0, 1000), y = runif(2000, 0, 1000),
#'                     type = sample(c("T cell", "B cell", "tumor"), 2000, TRUE))
#' q <- quadrat_counts(cells$x, cells$y, cells$type, n_grid = 10)
#' dim(q$Y)
#' @export
quadrat_counts <- function(x, y, type, n_grid = 15, min_cells = 5, types = NULL) {
  type <- as.character(type)
  if (is.null(types)) types <- sort(unique(type))
  xb <- seq(min(x), max(x), length.out = n_grid + 1)
  yb <- seq(min(y), max(y), length.out = n_grid + 1)
  ix <- findInterval(x, xb, all.inside = TRUE); iy <- findInterval(y, yb, all.inside = TRUE)
  counts <- matrix(0, n_grid^2, length(types), dimnames = list(NULL, types))
  n_cells <- stats::setNames(integer(length(types)), types)
  for (ct in types) {
    sel <- type == ct; n_cells[ct] <- sum(sel)
    counts[, ct] <- as.vector(table(factor(ix[sel], 1:n_grid), factor(iy[sel], 1:n_grid)))
  }
  estimable <- n_cells >= min_cells
  Y <- scale(sqrt(counts)); Y[, !estimable] <- 0; Y[!is.finite(Y)] <- 0
  Y <- unclass(Y); attr(Y, "scaled:center") <- NULL; attr(Y, "scaled:scale") <- NULL
  cx <- (xb[-1] + xb[-(n_grid + 1)]) / 2; cy <- (yb[-1] + yb[-(n_grid + 1)]) / 2
  coords <- as.matrix(expand.grid(x = cx, y = cy))
  list(Y = Y, counts = counts, coords = coords, estimable = estimable, n_cells = n_cells)
}
