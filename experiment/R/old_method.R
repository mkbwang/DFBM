
# Wrappers putting the old and the new factorization on a common footing.
#
# Both return the SAME estimand -- the marginal survival S^(t) = P(value > d_t)
# -- so neither is scored on its own internal parameterization. The old method
# estimates conditional probabilities and multiplies them up; the new method
# estimates sigma(X^(t)) directly.
#
# The per-threshold loop below reproduces the one in the deleted R/dfbm.R
# (recoverable with `git show HEAD:R/dfbm.R`), including its row/column dropping
# fallback, so the comparison is against the method as it was actually used.


#' Fit the old conditional logistic collaborative filter to every threshold
#'
#' @param Y list of T nested binary masks
#' @param fix_K fixed rank; when NULL the rank is chosen by cross validation
#' @param max_K largest rank considered by cross validation
#' @param lambdas ridge penalties considered
#' @param ignore drop a column whose risk set falls to at most this proportion
#' @param ncores cores for the rank search
#'
#' @returns a list with `S` (marginal survival per threshold), `pi_cond`
#'   (conditional probabilities), `ranks`, `risk_set` sizes and `elapsed`.
fit_old <- function(Y, fix_K = NULL, max_K = 10L,
                    lambdas = c(0.01, 0.1, 1), ignore = 0, ncores = 1L) {
  Tn <- length(Y)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  full <- matrix(1, N, P)

  pi_cond <- vector("list", Tn)
  ranks <- integer(Tn)
  risk_set <- numeric(Tn)

  t0 <- proc.time()[["elapsed"]]
  for (t in seq_len(Tn)) {
    mask <- if (t == 1L) full else Y[[t - 1L]] * full
    obs <- Y[[t]]
    risk_set[t] <- sum(mask)

    pi_mat <- matrix(0, N, P)
    if (risk_set[t] == 0) {          # the risk set is exhausted
      pi_cond[[t]] <- pi_mat
      ranks[t] <- 0L
      next
    }

    col_means <- colMeans(mask)
    cols_retain <- if (any(col_means <= ignore)) which(col_means > ignore) else seq_len(P)
    mask_sub <- mask[, cols_retain, drop = FALSE]
    row_means <- rowMeans(mask_sub)
    rows_retain <- if (any(row_means == 0)) which(row_means > 0) else seq_len(N)
    mask_sub <- mask_sub[rows_retain, , drop = FALSE]
    obs_sub <- obs[rows_retain, cols_retain, drop = FALSE]

    if (length(rows_retain) == 0L || length(cols_retain) == 0L) {
      pi_cond[[t]] <- pi_mat
      ranks[t] <- 0L
      next
    }

    if (nrow(obs_sub) == 1L || ncol(obs_sub) == 1L) {
      # Degenerate shape: nothing left to factorize, fall back to a scalar.
      scalar <- sum(obs_sub) / max(sum(mask_sub), 1)
      est <- matrix(scalar, nrow(obs_sub), ncol(obs_sub))
      ranks[t] <- 0L
    } else if (is.null(fix_K)) {
      cvfit <- DFBM::cv.logisticcfR(X = obs_sub, Z = mask_sub, max_K = max_K,
                                    lambdas = lambdas, ncores = ncores)
      est <- cvfit$pi
      ranks[t] <- cvfit$selected_K
    } else {
      lfit <- DFBM::logisticcfR(X = obs_sub, Z = mask_sub, K = fix_K,
                                lambda = lambdas[1L])
      est <- lfit$pi
      ranks[t] <- fix_K
    }
    pi_mat[rows_retain, cols_retain] <- est
    pi_cond[[t]] <- pi_mat
  }
  elapsed <- proc.time()[["elapsed"]] - t0

  # Marginal survival is the running product of the conditional probabilities.
  S <- vector("list", Tn)
  running <- matrix(1, N, P)
  for (t in seq_len(Tn)) {
    running <- running * pi_cond[[t]]
    S[[t]] <- running
  }

  list(S = S, pi_cond = pi_cond, ranks = ranks, risk_set = risk_set,
       elapsed = elapsed, method = "old")
}


#' Fit the new joint Clip-SVT factorization
#'
#' @param Y list of T nested binary masks
#' @param alpha shrinkage relative to the noise floor
#'   [DFBM::lambda_star_seq()]; ignored when `tune = TRUE`
#' @param tune whether to select alpha with [DFBM::tune.bmfsvt()]
#' @param alpha_grid grid used when tuning
#' @param criterion selection criterion passed to [DFBM::tune.bmfsvt()]
#' @param ... further arguments passed to [DFBM::bmfsvt()]
#'
#' @returns a list with the same shape as [fit_old()], plus the solver
#'   diagnostics `b`, `d`, `n_iter`, `n_svd` and `obj_trace`.
fit_new <- function(Y, alpha = 0.55, tune = FALSE,
                    alpha_grid = exp(seq(log(1), log(0.2), length.out = 10L)),
                    criterion = "cp", ...) {
  t0 <- proc.time()[["elapsed"]]
  if (tune) {
    tuned <- DFBM::tune.bmfsvt(Y, alpha_grid = alpha_grid,
                               criterion = criterion, ...)
    fit <- tuned$fit
    alpha_used <- tuned$alpha
    path <- tuned$path
  } else {
    fit <- DFBM::bmfsvt(Y, alpha = alpha, ...)
    alpha_used <- alpha
    path <- NULL
  }
  elapsed <- proc.time()[["elapsed"]] - t0

  list(S = fit$prob, pi_cond = NULL, ranks = fit$ranks,
       risk_set = rep(length(Y[[1L]]), length(Y)),
       elapsed = elapsed, method = "new",
       alpha = alpha_used, path = path,
       b = fit$b, d = fit$d, n_iter = fit$n_iter, n_svd = fit$n_svd,
       n_backtrack = fit$n_backtrack, converged = fit$converged,
       obj_trace = fit$obj_trace, clip_frac_trace = fit$clip_frac_trace)
}


#' Intercept only reference: every entry gets the mask prevalence
#'
#' @param Y list of T nested binary masks
#' @returns an object shaped like [fit_old()]
#' @details
#' This is the floor any factorization must beat. Reporting it keeps "the new
#' method wins" honest: without it, two methods can be compared to each other
#' while both are worse than doing nothing.
fit_baseline <- function(Y) {
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  nobs <- N * P
  t0 <- proc.time()[["elapsed"]]
  S <- lapply(Y, function(mat) matrix(mean(mat), N, P))
  list(S = S, pi_cond = NULL, ranks = rep(0L, length(Y)),
       risk_set = rep(nobs, length(Y)),
       elapsed = proc.time()[["elapsed"]] - t0, method = "baseline")
}
