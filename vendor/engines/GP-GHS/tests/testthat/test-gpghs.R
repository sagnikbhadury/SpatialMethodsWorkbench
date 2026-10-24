test_that("the basis has the right dimensions and spectral weights decrease with frequency", {
  coords <- matrix(runif(200), 100, 2)
  hs <- build_hsgp_basis(coords, m = 3)
  expect_equal(dim(hs$phi), c(100, 9))
  expect_length(hs$S_diag, 9)
  expect_true(all(diff(hs$S_diag[order(hs$omega_sq)]) <= 0))
  expect_true(all(abs(hs$coords_norm) <= 1))
})

test_that("the sampler returns the requested number of draws and bounded shrinkage factors", {
  sim <- simulate_svggm(n = 120, p = 4, pi = 0.5, seed = 3)
  hs <- build_hsgp_basis(sim$coords, m = 2)
  fit <- gp_group_horseshoe_sampler(sim$Y[, 1], sim$Y[, -1], hs$phi, hs$S_diag,
                                    nmc = 300, burn = 100, thin = 2)
  expect_equal(ncol(fit$theta_store), 100)
  expect_equal(nrow(fit$theta_store), 3 * 4)
  expect_true(all(fit$kappa_store >= 0 & fit$kappa_store <= 1))
  expect_true(all(fit$sigma_store > 0))
})

test_that("gpghs_graph recovers a strong constant edge and selects nothing on an empty graph", {
  set.seed(10)
  n <- 200; coords <- matrix(runif(2 * n), n, 2)
  x1 <- rnorm(n); x2 <- 0.9 * x1 + rnorm(n, sd = 0.4); x3 <- rnorm(n); x4 <- rnorm(n)
  Y <- cbind(a = x1, b = x2, c = x3, d = x4)
  fit <- gpghs_graph(Y, coords, m = 2, nmc = 600, burn = 300, thin = 3, n_cores = 1,
                     seed = 1, verbose = FALSE)
  expect_s3_class(fit, "gpghs")
  expect_equal(fit$adj["a", "b"], 1)
  expect_equal(sum(fit$adj[c("c", "d"), ]), 0)
  ed <- gpghs_edges(fit)
  expect_true(any(ed$from == "a" & ed$to == "b" & ed$sign == "positive"))
  cm <- coefficient_map(fit, "a", "b")
  expect_equal(nrow(cm), n)
  expect_true(mean(cm$mean) > 0)
})

test_that("quadrat counts give one row per quadrat and respect the cell rule", {
  cells <- data.frame(x = runif(500, 0, 100), y = runif(500, 0, 100),
                      type = c(rep("A", 300), rep("B", 197), rep("C", 3)))
  q <- quadrat_counts(cells$x, cells$y, cells$type, n_grid = 8, min_cells = 5)
  expect_equal(dim(q$Y), c(64, 3))
  expect_equal(sum(q$counts), 500)
  expect_false(q$estimable["C"])
  expect_true(all(q$Y[, "C"] == 0))
  expect_equal(nrow(q$coords), 64)
})

test_that("Moran's I is large for a smooth field and near zero for noise", {
  set.seed(5)
  coords <- matrix(runif(600), 300, 2)
  smooth <- sin(4 * coords[, 1]) + cos(3 * coords[, 2]) + rnorm(300, sd = 0.2)
  noise <- rnorm(300)
  mi <- morans_i(cbind(smooth = smooth, noise = noise), coords = coords, k = 6)
  expect_true(mi$moran_I[1] > 0.5)
  expect_true(abs(mi$moran_I[2]) < 0.2)
  expect_true(mi$fdr[1] < 0.05)
})

test_that("network_overlap counts common edges", {
  sim <- simulate_svggm(n = 150, p = 5, pi = 0.4, scenario = "constant", seed = 7)
  f1 <- gpghs_graph(sim$Y, sim$coords, m = 2, nmc = 300, burn = 100, thin = 2, n_cores = 1,
                    seed = 1, verbose = FALSE)
  f2 <- gpghs_graph(sim$Y, sim$coords, m = 2, nmc = 300, burn = 100, thin = 2, n_cores = 1,
                    seed = 2, verbose = FALSE)
  ov <- network_overlap(f1, f2)
  expect_equal(ov$n_pairs, 10)
  expect_true(ov$common <= min(ov$edges1, ov$edges2))
  m <- graph_metrics(f1$adj, sim$adj)
  expect_true(all(c("TPR", "FDP", "F1", "MCC") %in% names(m)))
})
