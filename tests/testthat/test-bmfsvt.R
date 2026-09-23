
# A tiny nested stack with known structure, reused across tests.
make_stack <- function(N = 30, P = 20, Tn = 4, r = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(rnorm(N * r), N, r)
  B <- matrix(rnorm(P * r), P, r)
  X0 <- 1 + tcrossprod(A, B) / sqrt(r)
  X <- vector("list", Tn)
  X[[1]] <- X0
  if (Tn > 1) {
    for (t in 2:Tn) {
      G <- tcrossprod(abs(matrix(rnorm(N * r), N, r)),
                      abs(matrix(rnorm(P * r), P, r)))
      X[[t]] <- X[[t - 1]] - G / mean(G)
    }
  }
  U <- matrix(runif(N * P), N, P)
  Y <- lapply(X, function(mat) 1 * (U < plogis(mat)))
  list(Y = Y, X = X, prob = lapply(X, plogis))
}


test_that("stack_forward and grad_backward invert the model parameterization", {
  sim <- make_stack()
  Y <- sim$Y
  mu <- vapply(Y, function(m) qlogis(mean(m)), numeric(1))
  nu <- diff(mu)
  Z <- lapply(seq_along(Y), function(t) matrix(rnorm(length(Y[[1]]), sd = 0.1),
                                               nrow(Y[[1]]), ncol(Y[[1]])))
  X <- DFBM:::stack_forward(Z, mu[1], nu)

  # X^(t) = mu_t + sum_{t' <= t} Z^(t')
  for (t in seq_along(Y)) {
    expect_equal(X[[t]], mu[t] + Reduce(`+`, Z[seq_len(t)]))
  }

  Psi <- DFBM:::grad_backward(X, Y)
  for (t in seq_along(Y)) {
    manual <- Reduce(`+`, lapply(seq(t, length(Y)),
                                 function(s) plogis(X[[s]]) - Y[[s]]))
    expect_equal(Psi[[t]], manual)
  }
})


test_that("grad_backward matches a finite difference of the stacked loss", {
  sim <- make_stack(N = 8, P = 6, Tn = 3, seed = 7)
  Y <- sim$Y
  mu <- vapply(Y, function(m) qlogis(mean(m)), numeric(1))
  nu <- diff(mu)
  set.seed(11)
  Z <- lapply(seq_along(Y), function(t) matrix(rnorm(48, sd = 0.3), 8, 6))

  f <- function(Zl) DFBM:::stack_ce(DFBM:::stack_forward(Zl, mu[1], nu), Y)
  Psi <- DFBM:::grad_backward(DFBM:::stack_forward(Z, mu[1], nu), Y)

  h <- 1e-6
  for (t in seq_along(Y)) {
    for (idx in c(1L, 17L, 40L)) {
      Zp <- Z; Zp[[t]][idx] <- Zp[[t]][idx] + h
      Zm <- Z; Zm[[t]][idx] <- Zm[[t]][idx] - h
      expect_equal((f(Zp) - f(Zm)) / (2 * h), Psi[[t]][idx], tolerance = 1e-5)
    }
  }
})


test_that("soft_svt reproduces an explicit soft thresholded SVD", {
  set.seed(3)
  M <- tcrossprod(matrix(rnorm(40 * 3), 40, 3), matrix(rnorm(25 * 3), 25, 3))
  thresh <- 2

  ref <- svd(M)
  keep <- ref$d > thresh
  expected <- ref$u[, keep, drop = FALSE] %*%
    ((ref$d[keep] - thresh) * t(ref$v[, keep, drop = FALSE]))

  for (method in c("full", "svds")) {
    got <- soft_svt(M, thresh, rank_guess = 1L, method = method)
    expect_equal(got$mat, expected, tolerance = 1e-8,
                 info = paste("method:", method))
    expect_equal(got$rank, sum(keep), info = paste("method:", method))
  }
})


test_that("lambda at or above lambda_max shrinks every block to zero", {
  sim <- make_stack(N = 40, P = 25, Tn = 3, seed = 5)
  lmax <- lambda_max_seq(sim$Y)
  fit <- bmfsvt(sim$Y, lambda = lmax * 1.001, max_iter = 30L, svd_method = "full")

  expect_true(all(vapply(fit$Z, function(m) max(abs(m)), numeric(1)) < 1e-8))
  # With Z = 0 the fitted probabilities collapse to the mask prevalences.
  for (t in seq_along(sim$Y)) {
    expect_equal(range(fit$prob[[t]]), rep(fit$pi[t], 2L), tolerance = 1e-8)
  }
})


test_that("T = 1 reduces to the single mask algorithm and beats the intercept", {
  sim <- make_stack(N = 200, P = 120, Tn = 1, seed = 9)
  fit <- bmfsvt(sim$Y, alpha = 0.5, max_iter = 200L, svd_method = "full")

  expect_length(fit$nu, 0L)
  expect_length(fit$b, 0L)
  baseline <- mean((fit$pi[1] - sim$prob[[1]])^2)
  expect_lt(mean((fit$prob[[1]] - sim$prob[[1]])^2), baseline)
})


test_that("clipping enforces monotone probabilities across thresholds", {
  sim <- make_stack(N = 60, P = 40, Tn = 5, seed = 13)
  fit <- bmfsvt(sim$Y, alpha = 0.05, max_iter = 150L, clip = TRUE,
                svd_method = "full")

  for (t in 2:length(fit$prob)) {
    expect_true(all(fit$prob[[t]] <= fit$prob[[t - 1]] + 1e-10),
                info = paste("threshold", t))
  }
  expect_length(fit$b, length(sim$Y) - 1L)
  expect_true(all(fit$b >= 0 & fit$b <= 1))
})


test_that("the fit recovers the truth better than the intercept only baseline", {
  sim <- make_stack(N = 120, P = 60, Tn = 5, seed = 21)
  fit <- bmfsvt(sim$Y, alpha = 0.3, max_iter = 300L, svd_method = "full")

  truth <- unlist(sim$prob)
  fitted <- unlist(fit$prob)
  baseline <- unlist(lapply(seq_along(sim$Y),
                            function(t) rep(fit$pi[t], length(sim$Y[[1]]))))

  expect_lt(sqrt(mean((fitted - truth)^2)), sqrt(mean((baseline - truth)^2)))
  expect_true(fit$converged)
})


test_that("weakening the shrinkage lowers the objective and raises the rank", {
  sim <- make_stack(N = 80, P = 40, Tn = 4, seed = 21)
  strong <- bmfsvt(sim$Y, alpha = 0.4, max_iter = 200L, svd_method = "full")
  weak <- bmfsvt(sim$Y, alpha = 0.05, max_iter = 200L, svd_method = "full")

  obj <- function(f) f$obj_trace[length(f$obj_trace)]
  # Both objectives include their own penalty, so compare the data fit only.
  expect_lt(DFBM:::stack_ce(weak$X, sim$Y),
            DFBM:::stack_ce(strong$X, sim$Y))
  expect_gt(sum(weak$ranks), sum(strong$ranks))
})


test_that("the objective decreases despite the clipping heuristic", {
  sim <- make_stack(N = 80, P = 40, Tn = 5, seed = 31)
  fit <- bmfsvt(sim$Y, alpha = 0.2, max_iter = 200L, clip = TRUE,
                svd_method = "full", restart = "gradient")

  trace <- fit$obj_trace
  expect_gt(length(trace), 5L)
  # Backtracking plus clipping can produce small upticks; require that the
  # trace is decreasing overall and that no single step increases much.
  expect_lt(trace[length(trace)], trace[1])
  expect_lt(max(diff(trace)) / abs(trace[1]), 1e-3)
})


test_that("non nested masks are rejected", {
  sim <- make_stack(N = 20, P = 10, Tn = 2, seed = 2)
  bad <- sim$Y
  bad[[2]][1, 1] <- 1
  bad[[1]][1, 1] <- 0
  expect_error(bmfsvt(bad), "not nested")
})


test_that("the objective trace is finite and the solver reports diagnostics", {
  sim <- make_stack(N = 50, P = 30, Tn = 4, seed = 17)
  fit <- bmfsvt(sim$Y, alpha = 0.1, max_iter = 100L, svd_method = "full")

  expect_true(all(is.finite(fit$obj_trace)))
  expect_length(fit$obj_trace, fit$n_iter)
  expect_gte(fit$n_svd, fit$n_iter)
  expect_true(is.numeric(fit$d) && all(fit$d >= 0))
})


test_that("entry_loss sums the per-threshold Brier loss for each entry", {
  sim <- make_stack(N = 12, P = 7, Tn = 3, seed = 41)
  X <- sim$X
  Y <- sim$Y
  idx <- c(3L, 17L, 60L)

  manual <- vapply(idx, function(i) {
    sum(vapply(seq_along(X), function(t) (plogis(X[[t]][i]) - Y[[t]][i])^2,
               numeric(1)))
  }, numeric(1))
  # One element per entry, in column major order.
  expect_length(DFBM:::entry_loss(X, Y), length(X[[1L]]))
  expect_equal(DFBM:::entry_loss(X, Y)[idx], manual)
})


test_that("svt_df matches the hard count and is smaller when weighted", {
  dvals <- list(c(5, 3, 1), c(2, 0.5), numeric(0))
  thresh <- c(1, 1, 1)
  N <- 40; P <- 25

  hard <- DFBM:::svt_df(dvals, thresh, N, P, weighted = FALSE)
  expect_equal(hard, c(3 * (N + P - 3), 2 * (N + P - 2), 0))

  wt <- DFBM:::svt_df(dvals, thresh, N, P, weighted = TRUE)
  expect_true(all(wt <= hard + 1e-8))

  # sum_{i<=r} (N + P - 2i + 1) == r(N + P - r), so with nothing shrunk away the
  # weighted count collapses onto the hard one exactly.
  wt0 <- DFBM:::svt_df(dvals, c(0, 0, 0), N, P, weighted = TRUE)
  expect_equal(wt0, hard)
  expect_gt(wt0[1], wt[1])

  # The threshold enters as d/(d + thresh), so a larger threshold means the
  # retained components were shrunk more and count for less.
  expect_lt(DFBM:::svt_df(dvals, c(10, 10, 10), N, P)[1], wt[1])
})


test_that("a per-alpha variance breaks the Cp surrogate but a fixed one does not", {
  # Synthetic path with the shape the real one has: as alpha weakens the Brier
  # falls, the degrees of freedom explode, and the plug-in variance collapses.
  # df is K x T, one column per slice.
  dfm <- cbind(c(0, 120, 500, 2000, 6000), c(0, 80, 300, 1000, 3000))
  path <- list(brier   = c(1000, 700, 500, 300, 100),
               df      = dfm,
               df_hard = dfm,
               v_tilde = c(0.21, 0.19, 0.12, 0.05, 0.01),
               v_t     = c(0.21, 0.09),
               Tn      = 2L,
               alpha   = c(1, 0.5, 0.25, 0.1, 0.01))

  fixed <- DFBM:::cp_surrogate(path)
  expect_equal(fixed$v_ref, path$v_t)
  expect_equal(fixed$cp, path$brier + 2 * as.vector(dfm %*% path$v_t))
  expect_lt(fixed$index, 5L)

  # A scalar v_ref is recycled, reproducing the superseded pooled form.
  pooled <- DFBM:::cp_surrogate(path, v_ref = 0.15)
  expect_equal(pooled$cp, path$brier + 2 * 0.15 * rowSums(dfm))

  # This is the failure the method note's earlier text would reintroduce.
  collapsing <- path$brier + 2 * path$v_tilde * rowSums(dfm)
  expect_equal(which.min(collapsing), 5L)
})


test_that("lambda_star_seq is reproducible, decreasing, and below lambda_max", {
  sim <- make_stack(N = 40, P = 25, Tn = 4, seed = 91)

  a <- lambda_star_seq(sim$Y, C = 5L, seed = 1L)
  b <- lambda_star_seq(sim$Y, C = 5L, seed = 1L)
  expect_equal(a, b)
  expect_false(isTRUE(all.equal(a, lambda_star_seq(sim$Y, C = 5L, seed = 2L))))

  expect_length(a, 4L)
  # Fewer residual terms are summed as t grows, so the floor falls.
  expect_true(all(diff(a) < 0))

  # The real data carries signal on top of the same noise, so its gradient norm
  # is strictly larger at every threshold. This is what makes the whole interval
  # [lambda*, lambda_max] the over-shrinking side of the useful range.
  expect_true(all(a < lambda_max_seq(sim$Y)))

  # It must not leave the global RNG stream disturbed.
  set.seed(99); before <- runif(1)
  set.seed(99); invisible(lambda_star_seq(sim$Y, C = 3L, seed = 5L))
  expect_equal(runif(1), before)
})


test_that("bootstrap masks stay nested and the covariance matches a direct sum", {
  sim <- make_stack(N = 40, P = 25, Tn = 3, seed = 51)
  path <- DFBM:::bmfsvt_path(sim$Y, alpha_grid = c(0.5, 0.2, 0.08),
                             max_iter = 40L, svd_method = "full")

  # The fitted probabilities are monotone across thresholds, so a single shared
  # uniform per entry yields a nested resample.
  set.seed(3)
  pd <- DFBM:::path_prob(path, 2L)
  U <- matrix(runif(40 * 25), 40, 25)
  Yb <- lapply(pd, function(p) 1 * (U <= p))
  expect_silent(DFBM:::as_mask_list(Yb, check_nested = TRUE))

  bt <- DFBM:::boot_covariance(path, index = 2L, B = 4L, window = 1L, seed = 7,
                               max_iter = 40L, svd_method = "full")
  expect_equal(bt$cand, 1:3)
  expect_length(bt$M, 3L)
  expect_equal(bt$M, path$brier[bt$cand] + 2 * bt$cov_total)

  # Recompute one candidate the slow, explicit way from the same seeds.
  set.seed(7)
  seeds <- sample.int(.Machine$integer.max, 4L)
  pis <- lapply(seeds, function(sd) {
    set.seed(sd)
    Ub <- matrix(runif(40 * 25), 40, 25)
    Yb <- lapply(pd, function(p) 1 * (Ub <= p))
    fb <- DFBM::bmfsvt(Yb, lambda = path$alpha[2] * path$lambda_star,
                       Z_init = path$Z[[2]], max_iter = 40L,
                       svd_method = "full")
    list(p = fb$prob, y = Yb)
  })
  direct <- 0
  for (t in seq_len(3)) {
    pmat <- vapply(pis, function(z) as.vector(z$p[[t]]), numeric(1000))
    ymat <- vapply(pis, function(z) as.vector(z$y[[t]]), numeric(1000))
    pbar <- rowMeans(pmat)
    direct <- direct + sum(rowSums((pmat - pbar) * (ymat - as.vector(pd[[t]])))) / 3
  }
  expect_equal(bt$cov_total[2], direct, tolerance = 1e-8)
})


test_that("boot_covariance gives the same answer on one core and several", {
  skip_on_os("windows")
  sim <- make_stack(N = 30, P = 20, Tn = 3, seed = 53)
  path <- DFBM:::bmfsvt_path(sim$Y, alpha_grid = c(0.4, 0.15),
                             max_iter = 30L, svd_method = "full")
  a <- DFBM:::boot_covariance(path, index = 1L, B = 4L, window = 1L, seed = 11,
                              ncores = 1L, max_iter = 30L, svd_method = "full")
  b <- DFBM:::boot_covariance(path, index = 1L, B = 4L, window = 1L, seed = 11,
                              ncores = 2L, max_iter = 30L, svd_method = "full")
  expect_equal(a$cov_total, b$cov_total)
})


test_that("tune.bmfsvt selects a grid alpha under both criteria", {
  sim <- make_stack(N = 60, P = 30, Tn = 3, seed = 23)
  grid <- c(0.5, 0.1, 0.02)

  cp <- tune.bmfsvt(sim$Y, alpha_grid = grid, criterion = "cp",
                    max_iter = 60L, svd_method = "full")
  expect_true(cp$alpha %in% grid)
  expect_equal(nrow(cp$path), 3L)
  expect_true(all(is.finite(cp$path$cp)))
  expect_true(all(is.na(cp$path$M)))          # no bootstrap was run
  expect_length(cp$fit$prob, 3L)
  expect_s3_class(cp$fit, "bmfsvt")
  expect_length(cp$fit$V, 3L)

  bo <- tune.bmfsvt(sim$Y, alpha_grid = grid, criterion = "boot", B = 3L,
                    window = 1L, seed = 2, max_iter = 60L, svd_method = "full")
  expect_true(bo$alpha %in% grid)
  expect_equal(bo$alpha_cp, cp$alpha)          # same surrogate, same anchor
  expect_true(any(is.finite(bo$path$M)))
})


test_that("bmfsvt_path reproduces individual fits and returns usable dvals", {
  sim <- make_stack(N = 50, P = 30, Tn = 3, seed = 57)
  grid <- c(0.4, 0.1)
  path <- DFBM:::bmfsvt_path(sim$Y, alpha_grid = grid, max_iter = 200L,
                             svd_method = "full")

  expect_equal(path$alpha, grid)
  expect_equal(dim(path$ranks), c(2L, 3L))
  # Probabilities rebuilt from the stored increments match a direct fit.
  direct <- bmfsvt(sim$Y, lambda = grid[1] * path$lambda_star, max_iter = 200L,
                   svd_method = "full")
  expect_equal(DFBM:::path_prob(path, 1L), direct$prob, tolerance = 1e-8)
  expect_length(direct$dvals, 3L)
  expect_equal(vapply(direct$dvals, length, integer(1)), direct$ranks)
})


test_that("soft_svt returns orthonormal right singular vectors for what it keeps", {
  set.seed(3)
  M <- tcrossprod(matrix(rnorm(40 * 3), 40, 3), matrix(rnorm(25 * 3), 25, 3))
  for (method in c("full", "svds")) {
    got <- soft_svt(M, 2, rank_guess = 1L, method = method)
    expect_equal(ncol(got$v), got$rank, info = paste("method:", method))
    expect_equal(crossprod(got$v), diag(got$rank), tolerance = 1e-8,
                 info = paste("method:", method))
  }
  expect_equal(dim(soft_svt(M, 1e6)$v), c(25L, 0L))
})


test_that("bmfsvt stores singular vectors spanning its unclipped blocks", {
  sim <- make_stack(N = 50, P = 30, Tn = 3, seed = 61)
  fit <- bmfsvt(sim$Y, alpha = 0.3, clip = FALSE, max_iter = 200L,
                svd_method = "full")

  expect_s3_class(fit, "bmfsvt")
  expect_equal(fit$offset, "scalar")
  expect_true(all(fit$ranks > 0L))
  for (t in 1:3) {
    expect_equal(ncol(fit$V[[t]]), fit$ranks[t])
    expect_length(fit$dvals[[t]], fit$ranks[t])
    expect_equal(fit$Z[[t]] %*% tcrossprod(fit$V[[t]]), fit$Z[[t]],
                 tolerance = 1e-8)
  }
  expect_error(bmfsvt(sim$Y, max_iter = 0L), "max_iter")
})


test_that("folding the training rows back in reproduces an unclipped fit", {
  sim <- make_stack(N = 60, P = 25, Tn = 3, seed = 63)
  for (offset in c("scalar", "column")) {
    fit <- bmfsvt(sim$Y, alpha = 0.3, clip = FALSE, tol = 1e-12,
                  max_iter = 5000L, svd_method = "full", offset = offset)
    expect_true(fit$converged, info = offset)
    pr <- predict(fit, sim$Y, tol = 1e-9)
    expect_true(all(pr$converged), info = offset)
    for (t in 1:3) {
      expect_lt(max(abs(pr$X[[t]] - fit$X[[t]])), 1e-4)
    }
  }
})


test_that("fold-in treats every row independently", {
  sim <- make_stack(N = 40, P = 20, Tn = 3, seed = 65)
  fit <- bmfsvt(sim$Y, alpha = 0.3, max_iter = 300L, svd_method = "full")
  all_rows <- predict(fit, sim$Y)

  sub <- predict(fit, lapply(sim$Y, function(y) y[1:5, , drop = FALSE]))
  one <- predict(fit, lapply(sim$Y, function(y) y[2, , drop = FALSE]))
  arr <- array(unlist(lapply(sim$Y, function(y) y[2, , drop = FALSE])),
               c(1L, 20L, 3L))
  one_arr <- predict(fit, arr)
  for (t in 1:3) {
    expect_lt(max(abs(sub$prob[[t]] - all_rows$prob[[t]][1:5, ])), 1e-6)
    expect_equal(dim(one$prob[[t]]), c(1L, 20L))
    expect_lt(max(abs(one$prob[[t]] - all_rows$prob[[t]][2, ])), 1e-6)
    expect_equal(one_arr$prob[[t]], one$prob[[t]])
  }
})


test_that("clipped fold-in is monotone and close to the in-sample fit", {
  sim <- make_stack(N = 60, P = 40, Tn = 5, seed = 13)
  fit <- bmfsvt(sim$Y, alpha = 0.2, max_iter = 300L, clip = TRUE,
                svd_method = "full")
  pr <- predict(fit, sim$Y)
  gap <- mean(abs(unlist(pr$prob) - unlist(fit$prob)))
  expect_lt(gap, 1e-2)

  # Masks the fit never saw.
  set.seed(4)
  U <- matrix(runif(10 * 40), 10, 40)
  pn <- predict(fit, lapply(c(0.2, 0.4, 0.6, 0.8, 0.9), function(q) 1 * (U > q)))
  for (res in list(pr, pn)) {
    for (t in 2:5) {
      expect_true(all(res$prob[[t]] <= res$prob[[t - 1]] + 1e-12))
    }
    expect_true(all(res$converged))
    expect_true(all(res$row_ce > 0))
  }
})


test_that("with every block at rank zero the fold-in returns the offsets", {
  sim <- make_stack(N = 40, P = 25, Tn = 3, seed = 5)
  fit <- bmfsvt(sim$Y, lambda = lambda_max_seq(sim$Y) * 1.001, max_iter = 30L,
                svd_method = "full")
  expect_true(all(vapply(fit$V, ncol, integer(1)) == 0L))

  pr <- predict(fit, sim$Y)
  expect_true(all(pr$converged))
  for (t in 1:3) {
    expect_equal(range(pr$prob[[t]]), rep(plogis(fit$mu[t]), 2L),
                 tolerance = 1e-12)
  }
  expect_error(predict(fit, sim$Y[1:2]), "masks")
  expect_error(predict(fit, lapply(sim$Y, function(y) y[, 1:10])), "columns")
})


test_that("column offsets calibrate a zero-inflated column", {
  sim <- make_stack(N = 80, P = 30, Tn = 4, seed = 71)
  set.seed(2)
  z <- runif(80) < 0.6
  # Zeroing the same rows in every mask keeps the stack nested.
  Y <- lapply(sim$Y, function(y) { y[z, 1] <- 0; y })
  lam <- 0.3 * lambda_star_seq(Y, C = 5L, seed = 1L)

  sc <- bmfsvt(Y, lambda = lam, max_iter = 300L, svd_method = "full")
  co <- bmfsvt(Y, lambda = lam, max_iter = 300L, svd_method = "full",
               offset = "column")
  expect_equal(dim(co$mu), c(4L, 30L))
  expect_equal(dim(co$nu), c(3L, 30L))
  expect_true(all(co$nu <= 0))

  err <- function(f) {
    sum(vapply(1:4, function(t) abs(mean(f$prob[[t]][, 1]) - mean(Y[[t]][, 1])),
               numeric(1)))
  }
  expect_lt(err(co), err(sc))

  default <- bmfsvt(Y, lambda = lam, max_iter = 300L, svd_method = "full")
  expect_identical(default$prob, sc$prob)
})
