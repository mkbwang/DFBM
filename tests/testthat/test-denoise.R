
# A small composition with low-rank structure and very different column scales.
make_composition <- function(n = 120, P = 8, seed = 1) {
  set.seed(seed)
  score <- matrix(rnorm(n * 2), n, 2)
  load <- matrix(rnorm(P * 2), P, 2)
  logits <- tcrossprod(score, load) +
    rep(seq(2, -3, length.out = P), each = n) +
    matrix(rnorm(n * P, sd = 0.5), n, P)
  A <- exp(logits)
  A <- A / rowSums(A)
  colnames(A) <- paste0("c", seq_len(P))
  A
}


test_that("thresholds on continuous columns are type 1 quantiles", {
  A <- make_composition()
  lev <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  th <- choose_thresholds(A, levels = lev)
  n <- nrow(A)

  expect_equal(dim(th$thresholds), c(5L, 8L))
  expect_equal(colnames(th$thresholds), colnames(A))
  for (j in seq_len(ncol(A))) {
    expect_equal(unname(th$thresholds[, j]),
                 unname(quantile(A[, j], lev, type = 1)))
  }
  # Order statistics with a strict inequality give exact, identical prevalences.
  Y <- make_masks(A, th)
  for (t in seq_along(lev)) {
    expect_equal(unname(colMeans(Y[[t]])),
                 rep(1 - ceiling(n * lev[t] - 1e-9) / n, ncol(A)))
  }
  expect_true(all(th$cap >= th$thresholds[5, ]))
})


test_that("a zero-heavy column gets a nonzero mask and positive thresholds", {
  A <- make_composition()
  set.seed(5)
  A[sample(nrow(A), 48), "c8"] <- 0
  lev <- seq(0.05, 0.95, by = 0.1)
  th <- choose_thresholds(A, levels = lev)
  d <- th$thresholds[, "c8"]

  expect_equal(unname(th$zero_frac["c8"]), 0.4)
  expect_equal(d[[1]], 1e-12)
  expect_true(all(d[-1] > 0))
  expect_false(is.unsorted(d, strictly = TRUE))
  expect_equal(unname(th$levels[, "c8"]),
               0.4 + 0.6 * (lev - lev[1]) / (1 - lev[1]))
  expect_equal(mean(make_masks(A, th)[[1]][, 8]), 0.6)
  # Columns with fewer zeros than the first level keep plain quantiles.
  expect_equal(unname(th$thresholds[, "c1"]),
               unname(quantile(A[, "c1"], lev, type = 1)))
})


test_that("residues at or below zero_tol count as zero in training and new data", {
  A <- make_composition()
  A[1:60, "c8"] <- 0
  A[61:70, "c8"] <- 1e-15
  th <- choose_thresholds(A)
  expect_equal(unname(th$zero_frac["c8"]), 70 / 120)
  expect_true(all(make_masks(A[61:65, ], th)[[1]][, 8] == 0))
})


test_that("discrete columns give nested masks and single rows stay matrices", {
  set.seed(9)
  A <- cbind(a = rpois(100, 2), b = runif(100))
  expect_warning(th <- choose_thresholds(A), "tied")
  Y <- make_masks(A, th)
  expect_silent(DFBM:::as_mask_list(Y, check_nested = TRUE))

  one <- make_masks(A[1, , drop = FALSE], th)
  expect_equal(dim(one[[1]]), c(1L, 2L))
})


test_that("masks are invariant to a monotone transform of the columns", {
  A <- make_composition()
  expect_equal(make_masks(A, choose_thresholds(A)),
               make_masks(sqrt(A), choose_thresholds(sqrt(A))))
})


test_that("interval_values summarizes capped values within each interval", {
  A <- matrix(c(0, 1, 2, 2.5, 4, 5, 6, 7, 8, 100), ncol = 1,
              dimnames = list(NULL, "x"))
  thr <- matrix(c(1, 4, 6), ncol = 1, dimnames = list(paste0("t", 1:3), "x"))

  # I0 = {0, 1}, I1 = {2, 2.5, 4}, I2 = {5, 6}, I3 = {7, 8, min(100, 9)}
  M <- interval_values(A, thr, cap = 9)
  expect_equal(unname(M[, 1]), c(0.5, 8.5 / 3, 5.5, 8))
  Mm <- interval_values(A, thr, cap = 9, summary = "median")
  expect_equal(unname(Mm[, 1]), c(0.5, 2.5, 5.5, 8))
  expect_equal(rownames(M), paste0("I", 0:3))
  expect_lte(max(M), 9)

  # Tied thresholds leave (4, 4] empty; it takes min(d_t, cap).
  thr2 <- matrix(c(1, 4, 4), ncol = 1, dimnames = list(paste0("t", 1:3), "x"))
  M2 <- interval_values(A, thr2, cap = 9)
  expect_equal(M2[3, 1], 4)
  expect_false(is.unsorted(M2[, 1]))
})


test_that("survival_expectation gives the binned value on masks and the interval sum", {
  A <- make_composition(n = 60, P = 5)
  th <- choose_thresholds(A, levels = c(0.2, 0.5, 0.8))
  M <- interval_values(A, th$thresholds, th$cap)
  Y <- make_masks(A, th)

  # On the observed masks each entry snaps to its own interval representative.
  bin <- Reduce(`+`, Y)
  expected <- matrix(M[cbind(as.vector(bin) + 1, rep(1:5, each = 60))], 60, 5)
  expect_equal(unname(DFBM:::survival_expectation(Y, M)), expected)

  # Abel form equals sum_t (S_t - S_{t+1}) m_t, bounded even when S is not monotone.
  set.seed(1)
  S <- replicate(3, matrix(runif(300), 60, 5), simplify = FALSE)
  Sx <- c(list(matrix(1, 60, 5)), S, list(matrix(0, 60, 5)))
  direct <- matrix(0, 60, 5)
  for (t in 0:3) {
    direct <- direct + (Sx[[t + 1]] - Sx[[t + 2]]) * rep(M[t + 1, ], each = 60)
  }
  got <- DFBM:::survival_expectation(S, M)
  expect_equal(unname(got), unname(direct))
  expect_true(all(got >= rep(M[1, ], each = 60) - 1e-12))
  expect_true(all(got <= rep(M[4, ], each = 60) + 1e-12))

  zeros <- replicate(3, matrix(0, 60, 5), simplify = FALSE)
  ones <- replicate(3, matrix(1, 60, 5), simplify = FALSE)
  expect_equal(unname(DFBM:::survival_expectation(zeros, M)),
               matrix(unname(rep(M[1, ], each = 60)), 60, 5))
  expect_equal(unname(DFBM:::survival_expectation(ones, M)),
               matrix(unname(rep(M[4, ], each = 60)), 60, 5))
})


test_that("choose_thresholds rejects invalid input", {
  A <- make_composition()
  bad <- A
  bad[1, 1] <- NA
  expect_error(choose_thresholds(bad), "missing")
  bad <- A
  bad[1, 1] <- -1
  expect_error(choose_thresholds(bad), "nonnegative")
  expect_error(choose_thresholds(A, levels = c(0.5, 0.2)), "increasing")
  expect_error(choose_thresholds(A, levels = c(0.2, 0.5), q_cap = 0.4), "q_cap")
})


test_that("dfbm fits with a fixed alpha and returns closed positive compositions", {
  A <- make_composition(n = 80, P = 6)
  A[1:40, 6] <- 0
  lev <- c(0.1, 0.4, 0.7, 0.9)
  obj <- dfbm(A, alpha = 0.5, levels = lev, seed = 3, C = 5L, max_iter = 200L,
              svd_method = "full")

  expect_s3_class(obj, "dfbm")
  expect_null(obj$cv)
  expect_equal(obj$alpha, 0.5)
  expect_equal(obj$fit$offset, "column")

  # The fixed mode is exactly one bmfsvt() call at alpha * lambda_star.
  Y <- make_masks(A, obj$thresholds)
  ref <- bmfsvt(Y, lambda = 0.5 * lambda_star_seq(Y, C = 5L, seed = 3),
                offset = "column", max_iter = 200L, svd_method = "full")
  expect_equal(obj$fit$prob, ref$prob)

  expect_equal(unname(rowSums(obj$denoised)), rep(1, 80))
  expect_true(all(obj$denoised > 0))
  expect_true(all(obj$denoised[1:40, 6] > 0))       # zeros filled in
  unclosed <- obj$denoised * obj$row_sum
  expect_true(all(unclosed <= rep(obj$cap, each = 80) + 1e-12))
  for (nm in c("row_sum", "binned_row_sum", "row_ce", "flag", "converged")) {
    expect_length(obj[[nm]], 80L)
  }
  expect_named(obj$insample_gap, c("mean", "max"))
})


test_that("dfbm cross-validates alpha on the grid when alpha is NULL", {
  A <- make_composition(n = 80, P = 6)
  obj <- dfbm(A, levels = c(0.1, 0.4, 0.7, 0.9), alpha_grid = c(1, 0.5),
              n_splits = 2L, seed = 3, C = 5L, max_iter = 100L,
              svd_method = "full")
  expect_true(obj$alpha %in% c(1, 0.5))
  expect_equal(nrow(obj$cv), 2L)
  # dfbm always supplies the values, so every loss is reported.
  expect_false(anyNA(obj$cv[, c("rps", "mse", "crps")]))
  expect_equal(obj$alpha, obj$cv$alpha[which.min(obj$cv$rps)])
})


test_that("predict.dfbm reproduces the training output and checks columns", {
  A <- make_composition(n = 80, P = 6)
  obj <- dfbm(A, alpha = 0.5, levels = c(0.1, 0.4, 0.7, 0.9), seed = 3,
              C = 5L, max_iter = 200L, svd_method = "full")

  again <- predict(obj, A)
  expect_equal(again$denoised, obj$denoised)
  sub <- predict(obj, A[1:3, ])
  expect_lt(max(abs(sub$denoised - obj$denoised[1:3, ])), 1e-6)
  expect_error(predict(obj, A[, 6:1]), "order")
  expect_error(predict(obj, A[, 1:5]), "columns")
})
