## Posterior check of gp_group_horseshoe_sampler against quadrature.
##
## With theta and sigma^2 integrated out analytically, the posterior of
## (lambda_1, lambda_2, tau) is three-dimensional. Under u = atan(scale) the
## half-Cauchy prior is uniform on (0, pi/2), so a midpoint rule on the cube
## gives the posterior moments without any Markov chain. sigma^2 ~ IG(2, 2).
##
## Run from the project root:  Rscript tests/test_sampler_quadrature.R

source("gp_group_horseshoe_spatial_graph.R")
set.seed(11)

n <- 30; q <- 2; m2 <- 2
X_s <- matrix(rnorm(n * q), n, q)
phi <- matrix(rnorm(n * m2), n, m2)
S_diag <- c(1.5, 0.4)
a0 <- 2; b0 <- 2
Z <- cbind(X_s[, 1] * phi, X_s[, 2] * phi)
theta_true <- c(0.8, -0.5, 0, 0)          # block 1 active, block 2 null
y <- as.vector(Z %*% theta_true) + 0.7 * rnorm(n)

## log p(y | lambda, tau): y | scales, sigma^2 ~ N(0, sigma^2 (I + Z D Z')),
## sigma^2 ~ IG(a0, b0)
ZtZ <- crossprod(Z); Zty <- crossprod(Z, y); yty <- sum(y^2)
log_marg <- function(l1, l2, tau) {
  d <- tau^2 * rep(c(l1^2, l2^2), each = m2) * rep(S_diag, q)
  A <- ZtZ + diag(1 / d)
  R <- chol(A)
  quad <- yty - sum(backsolve(R, Zty, transpose = TRUE)^2)
  logdet <- 2 * sum(log(diag(R))) + sum(log(d))   # log|I + Z D Z'|
  -0.5 * logdet - (n / 2 + a0) * log(b0 + quad / 2)
}

G <- 70
u <- (seq_len(G) - 0.5) / G * (pi / 2)
s <- tan(u)
lw <- array(0, c(G, G, G))
for (i in 1:G) for (j in 1:G) for (k in 1:G)
  lw[i, j, k] <- log_marg(s[i], s[j], s[k])
w <- exp(lw - max(lw)); w <- w / sum(w)
quad_res <- c(atan_lam1 = sum(apply(w, 1, sum) * u),
              atan_lam2 = sum(apply(w, 2, sum) * u),
              atan_tau  = sum(apply(w, 3, sum) * u))

fit <- gp_group_horseshoe_sampler(y, X_s, phi, S_diag, nmc = 220000,
                                  burn = 20000, thin = 1,
                                  a_sigma = a0, b_sigma = b0)
gibbs <- rbind(atan(sqrt(fit$lambda_store)), atan(sqrt(fit$tau_store)))
nb <- 50; bs <- ncol(gibbs) %/% nb
se <- apply(gibbs, 1, function(x) sd(colMeans(matrix(x[1:(nb * bs)], bs, nb))) /
              sqrt(nb))
out <- data.frame(quadrature = quad_res, gibbs = rowMeans(gibbs), mcse = se)
out$z <- (out$gibbs - out$quadrature) / out$mcse
print(out, digits = 4)
cat("posterior mean kappa (block 1 active, block 2 null):",
    round(rowMeans(fit$kappa_store), 3), "\n")
