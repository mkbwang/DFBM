# A nested stack with low-rank structure, as make_stack() in test-bmfsvt.R.
make_stack_cv <- function(N = 60, P = 20, Tn = 4, r = 2, seed = 1) {
  set.seed(seed)
  A <- matrix(rnorm(N * r), N, r)
  B <- matrix(rnorm(P * r), P, r)
  X <- vector("list", Tn)
  X[[1]] <- 1 + tcrossprod(A, B) / sqrt(r)
  if (Tn > 1) {
    for (t in 2:Tn) {
      G <- tcrossprod(abs(matrix(rnorm(N * r), N, r)),
                      abs(matrix(rnorm(P * r), P, r)))
      X[[t]] <- X[[t - 1]] - G / mean(G)
    }
  }
  U <- matrix(runif(N * P), N, P)
  lapply(X, function(mat) 1 * (U < plogis(mat)))
}


test_that("train_mask_generation balances every row and column", {
  N <- 103
  P <- 47
  W <- train_mask_generation(matrix(0, N, P), prop = 0.8, seed = 80)
  expect_equal(dim(W), c(N, P))
  expect_true(all(W %in% c(0, 1)))

  n_val <- round(N * 0.2)
  expect_true(all(colSums(W == 0) == n_val))
  row_val <- rowSums(W == 0)
  expect_lte(max(row_val) - min(row_val), 1)
  expect_equal(mean(W), 1 - n_val / N)

  # Far tighter than choosing the validation entries at random.
  set.seed(1)
  R <- matrix(rbinom(N * P, 1, 0.8), N, P)
  expect_lt(var(rowSums(W)), var(rowSums(R)) / 10)
  expect_lt(var(colSums(W)), var(colSums(R)) / 10)
})


test_that("train_mask_generation is reproducible and leaves the RNG alone", {
  X <- matrix(0, 50, 12)
  set.seed(5)
  before <- .Random.seed
  a <- train_mask_generation(X, seed = 9)
  expect_identical(.Random.seed, before)
  b <- train_mask_generation(X, seed = 9)
  expect_identical(a, b)
  expect_false(identical(a, train_mask_generation(X, seed = 10)))
  expect_error(train_mask_generation(X, prop = 1), "prop")
  expect_error(train_mask_generation(matrix(0, 3, 4), prop = 0.9), "one entry")
})


test_that("train_mask_generation scales to the HRS shape", {
  elapsed <- system.time(W <- train_mask_generation(matrix(0, 6212, 25),
                                                    seed = 1))[["elapsed"]]
  expect_lt(elapsed, 5)
  row_val <- rowSums(W == 0)
  expect_lte(max(row_val) - min(row_val), 1)
})


test_that("an all-ones training mask is the unmasked fit", {
  Y <- make_stack_cv()
  lam <- 0.5 * lambda_star_seq(Y, C = 5L, seed = 2)
  a <- bmfsvt(Y, lambda = lam, max_iter = 100L, svd_method = "full")
  b <- bmfsvt(Y, train = matrix(1, 60, 20), lambda = lam, max_iter = 100L,
              svd_method = "full")
  expect_identical(a$X, b$X)
  expect_identical(lambda_star_seq(Y, C = 5L, seed = 2),
                   lambda_star_seq(Y, C = 5L, seed = 2,
                                   train = matrix(1, 60, 20)))
})


test_that("held-out entries have no influence on the fit", {
  Y <- make_stack_cv()
  W <- train_mask_generation(Y[[1]], prop = 0.8, seed = 4)
  # Overwrite every held-out entry with a different nested pattern.
  Y2 <- lapply(Y, function(m) { m[W == 0] <- 0; m })
  for (offset in c("scalar", "column")) {
    # Both draw lambda_star from the global stream; reset it between them.
    set.seed(1)
    a <- bmfsvt(Y, train = W, alpha = 0.5, C = 5L, max_iter = 100L,
                offset = offset, svd_method = "full")
    set.seed(1)
    b <- bmfsvt(Y2, train = W, alpha = 0.5, C = 5L, max_iter = 100L,
                offset = offset, svd_method = "full")
    expect_identical(a$X, b$X)
    expect_identical(a$mu, b$mu)
    expect_equal(a$n_train, sum(W))
  }
})


test_that("the masked gradient matches a finite difference of the masked loss", {
  Y <- make_stack_cv(N = 8, P = 6, Tn = 3, seed = 7)
  W <- train_mask_generation(Y[[1]], prop = 0.7, seed = 3)
  mu <- vapply(Y, function(m) qlogis(mean(m)), numeric(1))
  nu <- diff(mu)
  set.seed(11)
  Z <- lapply(seq_along(Y), function(t) matrix(rnorm(48, sd = 0.3), 8, 6))

  f <- function(Zl) DFBM:::stack_ce(DFBM:::stack_forward(Zl, mu[1], nu), Y, W)
  Psi <- DFBM:::grad_backward(DFBM:::stack_forward(Z, mu[1], nu), Y, W)
  expect_true(all(Psi[[1]][W == 0] == 0))

  h <- 1e-6
  for (t in seq_along(Y)) {
    for (idx in c(1L, 17L, 40L)) {
      Zp <- Z; Zp[[t]][idx] <- Zp[[t]][idx] + h
      Zm <- Z; Zm[[t]][idx] <- Zm[[t]][idx] - h
      expect_equal((f(Zp) - f(Zm)) / (2 * h), Psi[[t]][idx], tolerance = 1e-5)
    }
  }
})


test_that("held-out entries are imputed and clipped to monotone probabilities", {
  Y <- make_stack_cv()
  W <- train_mask_generation(Y[[1]], prop = 0.8, seed = 4)
  set.seed(2)
  fit <- bmfsvt(Y, train = W, alpha = 0.3, C = 5L, max_iter = 150L)
  for (t in 2:length(Y)) {
    expect_true(all(fit$prob[[t]] <= fit$prob[[t - 1]] + 1e-12))
  }
  # The held-out probabilities come from the low-rank blocks, not the offsets.
  expect_gt(sd(fit$prob[[1]][W == 0]), 0.01)
})


test_that("the masked noise floor is below the full one", {
  Y <- make_stack_cv(N = 150, P = 60)
  W <- train_mask_generation(Y[[1]], prop = 0.8, seed = 4)
  full <- lambda_star_seq(Y, C = 10L, seed = 1)
  masked <- lambda_star_seq(Y, C = 10L, seed = 1, train = W)
  ratio <- masked / full
  expect_true(all(ratio < 1 & ratio > 0.8))
})


test_that("crps_discrete matches numerical integration", {
  set.seed(3)
  Tn <- 4
  for (rep in 1:20) {
    m <- sort(runif(Tn + 1, 0, 3))
    s <- sort(runif(Tn), decreasing = TRUE)
    y <- runif(1, -0.5, 3.5)
    Fx <- function(x) {
      # atoms at m, masses S_t - S_{t+1}
      cdf <- c(1 - s, 1)
      vapply(x, function(z) {
        k <- sum(m <= z)
        if (k == 0) 0 else cdf[k]
      }, numeric(1))
    }
    h <- function(x) (Fx(x) - (x >= y))^2
    knots <- sort(c(m, y))
    num <- sum(vapply(seq_len(length(knots) - 1), function(i) {
      integrate(h, knots[i], knots[i + 1])$value
    }, numeric(1)))
    got <- DFBM:::crps_discrete(as.list(s), matrix(m, Tn + 1, 1), y)
    expect_equal(got, num, tolerance = 1e-6)
  }

  # A degenerate distribution at atom k scores |m_k - y|: the first entry has
  # S = (1, 1, 0), all mass on m_2; the second S = (1, 0, 0), all on m_1.
  m <- c(0, 1, 2, 5)
  got <- DFBM:::crps_discrete(list(c(1, 1), c(1, 0), c(0, 0)),
                              matrix(m, 4, 2), c(0.3, 4))
  expect_equal(got, c(abs(2 - 0.3), abs(1 - 4)))
})


test_that("cv.bmfsvt selects a grid alpha and refits on every entry", {
  Y <- make_stack_cv(N = 80, P = 15, Tn = 3)
  grid <- c(1, 0.6, 0.3)
  cv <- cv.bmfsvt(Y, alpha_grid = grid, n_splits = 2L, C = 5L, seed = 1,
                  max_iter = 80L, svd_method = "full")
  expect_true(cv$alpha %in% grid)
  expect_equal(cv$cv$alpha, sort(grid, decreasing = TRUE))
  expect_equal(cv$alpha, cv$cv$alpha[which.min(cv$cv$rps)])
  expect_true(all(is.na(cv$cv$mse)))
  expect_equal(dim(cv$loss_matrix$rps), c(2L, 3L))
  expect_equal(dim(cv$lambda_star_splits), c(2L, 3L))
  expect_true(all(cv$lambda_star_splits < rep(cv$lambda_star, each = 2)))
  expect_equal(cv$fit$n_train, 80 * 15)
  expect_equal(cv$fit$lambda, cv$alpha * cv$lambda_star)

  # Reproducible, and the same whether split over cores or not.
  again <- cv.bmfsvt(Y, alpha_grid = grid, n_splits = 2L, C = 5L, seed = 1,
                     ncores = 2L, max_iter = 80L, svd_method = "full")
  expect_equal(again$cv, cv$cv)

  expect_error(cv.bmfsvt(Y, loss = "mse", n_splits = 1L), "needs `A`")
})


test_that("cv.bmfsvt scores values on the capped scale", {
  set.seed(1)
  n <- 80
  P <- 6
  score <- matrix(rnorm(n * 2), n, 2)
  load <- matrix(rnorm(P * 2), P, 2)
  A <- exp(tcrossprod(score, load) + rep(seq(2, -3, length.out = P), each = n) +
             matrix(rnorm(n * P, sd = 0.5), n, P))
  A <- A / rowSums(A)
  th <- choose_thresholds(A, levels = c(0.1, 0.4, 0.7, 0.9))
  Y <- make_masks(A, th$thresholds)

  grid <- c(1, 0.5)
  res <- lapply(c("rps", "mse", "crps"), function(l) {
    cv.bmfsvt(Y, alpha_grid = grid, loss = l, n_splits = 2L, A = A,
              thresholds = th, C = 5L, seed = 2, max_iter = 80L,
              svd_method = "full")
  })
  # The fits do not depend on the loss, so neither do the reported losses.
  expect_equal(res[[1]]$cv, res[[2]]$cv)
  for (i in 1:3) {
    l <- c("rps", "mse", "crps")[i]
    expect_false(anyNA(res[[i]]$cv[, c("rps", "mse", "crps")]))
    expect_equal(res[[i]]$alpha, grid[which.min(res[[i]]$cv[[l]])])
  }
  # Normalized losses are O(1), not dominated by the largest column.
  expect_lt(max(res[[1]]$cv$mse), 5)
})


test_that("interval_values with a training mask uses training rows only", {
  A <- matrix(c(1:10, 11:20), 10, 2)
  th <- matrix(c(3, 7, 13, 17), 2, 2)
  W <- matrix(1, 10, 2)
  W[1, 1] <- 0   # remove the value 1 from interval I0 of column 1
  M <- interval_values(A, th, cap = c(10, 20), train = W)
  expect_equal(M[[1, 1]], mean(2:3))
  expect_equal(M[, 2], interval_values(A, th, cap = c(10, 20))[, 2])
})
