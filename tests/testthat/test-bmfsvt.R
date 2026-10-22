
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


test_that("a very large lambda shrinks every block to zero", {
  sim <- make_stack(N = 40, P = 25, Tn = 3, seed = 5)
  fit <- bmfsvt(sim$Y, lambda = 1e6, max_iter = 30L, svd_method = "full")

  expect_true(all(vapply(fit$Z, function(m) max(abs(m)), numeric(1)) < 1e-8))
  # With Z = 0 the fitted probabilities collapse to the column offsets.
  for (t in seq_along(sim$Y)) {
    expect_equal(fit$prob[[t]], matrix(plogis(fit$mu[t, ]), 40, 25, byrow = TRUE),
                 tolerance = 1e-8)
  }
})


test_that("T = 1 reduces to the single mask algorithm and beats the intercept", {
  sim <- make_stack(N = 200, P = 120, Tn = 1, seed = 9)
  fit <- bmfsvt(sim$Y, alpha = 0.5, max_iter = 200L, svd_method = "full")

  expect_length(fit$nu, 0L)
  baseline <- mean((fit$pi[1] - sim$prob[[1]])^2)
  expect_lt(mean((fit$prob[[1]] - sim$prob[[1]])^2), baseline)
})


test_that("clipping enforces monotone probabilities across thresholds", {
  sim <- make_stack(N = 60, P = 40, Tn = 5, seed = 13)
  fit <- bmfsvt(sim$Y, alpha = 0.05, max_iter = 150L, svd_method = "full")

  for (t in 2:length(fit$prob)) {
    expect_true(all(fit$prob[[t]] <= fit$prob[[t - 1]] + 1e-10),
                info = paste("threshold", t))
  }
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
  fit <- bmfsvt(sim$Y, alpha = 0.2, max_iter = 200L, svd_method = "full",
                restart = "gradient")

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
  expect_true(all(unlist(fit$dvals) >= 0))
})


test_that("lambda_star_seq is reproducible and decreasing", {
  sim <- make_stack(N = 40, P = 25, Tn = 4, seed = 91)

  a <- lambda_star_seq(sim$Y, C = 5L, seed = 1L)
  b <- lambda_star_seq(sim$Y, C = 5L, seed = 1L)
  expect_equal(a, b)
  expect_false(isTRUE(all.equal(a, lambda_star_seq(sim$Y, C = 5L, seed = 2L))))

  expect_length(a, 4L)
  # Fewer residual terms are summed as t grows, so the floor falls.
  expect_true(all(diff(a) < 0))

  # It must not leave the global RNG stream disturbed.
  set.seed(99); before <- runif(1)
  set.seed(99); invisible(lambda_star_seq(sim$Y, C = 3L, seed = 5L))
  expect_equal(runif(1), before)
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


test_that("bmfsvt stores singular vectors paired with its singular values", {
  sim <- make_stack(N = 50, P = 30, Tn = 3, seed = 61)
  fit <- bmfsvt(sim$Y, alpha = 0.3, max_iter = 200L, svd_method = "full")

  expect_s3_class(fit, "bmfsvt")
  expect_true(all(fit$ranks > 0L))
  for (t in 1:3) {
    expect_equal(ncol(fit$V[[t]]), fit$ranks[t])
    expect_length(fit$dvals[[t]], fit$ranks[t])
  }
  # The first block is never clipped, so V spans its rows exactly.
  expect_equal(fit$Z[[1]] %*% tcrossprod(fit$V[[1]]), fit$Z[[1]],
               tolerance = 1e-8)
  expect_error(bmfsvt(sim$Y, max_iter = 0L), "max_iter")
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
  fit <- bmfsvt(sim$Y, alpha = 0.2, max_iter = 300L, svd_method = "full")
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
  fit <- bmfsvt(sim$Y, lambda = 1e6, max_iter = 30L, svd_method = "full")
  expect_true(all(vapply(fit$V, ncol, integer(1)) == 0L))

  pr <- predict(fit, sim$Y)
  expect_true(all(pr$converged))
  for (t in 1:3) {
    expect_equal(pr$prob[[t]], matrix(plogis(fit$mu[t, ]), 40, 25, byrow = TRUE),
                 tolerance = 1e-12)
  }
  expect_error(predict(fit, sim$Y[1:2]), "masks")
  expect_error(predict(fit, lapply(sim$Y, function(y) y[, 1:10])), "columns")
})


test_that("offsets are per column and non-increasing across thresholds", {
  sim <- make_stack(N = 80, P = 30, Tn = 4, seed = 71)
  set.seed(2)
  z <- runif(80) < 0.6
  # Zeroing the same rows in every mask keeps the stack nested.
  Y <- lapply(sim$Y, function(y) { y[z, 1] <- 0; y })
  lam <- 0.3 * lambda_star_seq(Y, C = 5L, seed = 1L)

  co <- bmfsvt(Y, lambda = lam, max_iter = 300L, svd_method = "full")
  expect_equal(dim(co$mu), c(4L, 30L))
  expect_equal(dim(co$nu), c(3L, 30L))
  expect_true(all(co$nu <= 0))
})
