
#' Joint factorization of a nested sequence of binary masks
#'
#' Implements the Clip-SVT algorithm (Algorithm 2 of the method note): all
#' binary masks are factorized together under one objective, optimized by FISTA
#' with backtracking and block-wise singular value thresholding.
#'
#' @param Y a list of T binary matrices, or an N x P x T array. The masks must be
#'   nested and decreasing, i.e. \eqn{Y^{(0)}_{ij} \ge Y^{(1)}_{ij} \ge \cdots},
#'   which holds automatically when they come from thresholding one abundance
#'   matrix with an increasing threshold sequence.
#' @param Omega an N x P binary matrix marking observed entries, shared by every
#'   mask. `NULL` means all entries are observed.
#' @param lambda shrinkage parameter for each mask; a scalar is recycled. When
#'   `NULL` the value `alpha * lambda_max_seq(Y, Omega)` is used.
#' @param alpha shrinkage relative to the largest useful value, used only when
#'   `lambda` is `NULL`.
#' @param max_iter maximum number of outer FISTA iterations.
#' @param tol relative change in `Z` below which the algorithm stops.
#' @param gamma step size increment ratio for backtracking.
#' @param L_rewind factor by which the Lipschitz constants are optimistically
#'   relaxed at the start of each iteration. The method note rewinds by `gamma`,
#'   which keeps the step size adaptive but costs roughly one extra backtrack
#'   per iteration, and every backtrack repeats all T decompositions. Set to 1
#'   to make the constants non decreasing and trade adaptivity for speed.
#' @param delta small constant stabilizing the relative convergence criterion.
#' @param rank_init starting guess for the rank of each block.
#' @param rank_max largest rank retained for any block; `NULL` means no limit.
#' @param rank_step how many ranks to add when the partial SVD was too small.
#' @param clip whether to enforce the monotonicity restriction by clipping the
#'   SVT outcome. Setting this to `FALSE` solves the unconstrained stacked
#'   problem, which is useful as a reference for how hard the constraint binds.
#' @param L_init how to initialize the Lipschitz constants. `"stacked"` uses
#'   \eqn{\sum_{t' \ge t} \pi_{t'}(1-\pi_{t'})}, `"single"` reproduces the
#'   method note exactly with \eqn{\pi_t(1-\pi_t)}, `"worst"` uses the
#'   conservative bound \eqn{(T-t)/4}.
#' @param restart FISTA restart scheme, one of `"gradient"`, `"function"`,
#'   `"none"`.
#' @param svd_method one of `"auto"`, `"full"`, `"svds"`; see [soft_svt()].
#' @param track_objective how to record the objective trace. `"approx"` reuses
#'   the singular values already computed by the thresholding step, which is
#'   exact unless clipping actually altered entries; `"exact"` recomputes the
#'   nuclear norms of the clipped blocks at the cost of an extra SVD per block
#'   per iteration; `"none"` skips the trace.
#' @param Z_init optional list of T matrices used to warm start `Z`.
#' @param max_backtrack cap on backtracking steps within one outer iteration.
#' @param verbose whether to print per-iteration progress.
#'
#' @returns a list with components
#'   \describe{
#'     \item{X}{list of T natural parameter matrices}
#'     \item{prob}{list of T matrices, \eqn{\sigma(X^{(t)})}, the estimated
#'       marginal probability that an entry exceeds threshold t}
#'     \item{Z}{list of T increment matrices actually optimized}
#'     \item{mu, nu}{offsets and offset increments}
#'     \item{lambda, L}{shrinkage parameters and final Lipschitz constants}
#'     \item{ranks}{retained rank of each block}
#'     \item{b, d}{violation rate and maximum violation margin of an unclipped
#'       proximal step taken from the solution, for t = 1..T-1}
#'     \item{obj_trace, clip_frac_trace}{objective and clipped-entry fraction by
#'       iteration}
#'     \item{n_iter, n_svd, n_backtrack, converged}{solver diagnostics}
#'   }
#'
#' @details
#' The parameterization follows section 2.2. With
#' \eqn{X^{(t)} = X^{(0)} + \sum_{t' \le t} H^{(t')}} and \eqn{H^{(t)} \le 0},
#' writing \eqn{Z^{(0)} = X^{(0)} - \mu_0} and \eqn{Z^{(t)} = H^{(t)} - \nu_t}
#' gives \eqn{X^{(t)} = \mu_t + \sum_{t' \le t} Z^{(t')}} because the offsets
#' telescope. The gradient with respect to \eqn{Z^{(t)}} is therefore the
#' reverse cumulative sum \eqn{\sum_{t' \ge t}[\sigma(X^{(t')}) - Y^{(t')}]},
#' and the monotonicity restriction becomes \eqn{Z^{(t)} \le -\nu_t}.
#'
#' Unlike the conditional formulation used by [logisticcfR()], every observed
#' entry contributes to the likelihood of every mask. The set of entries
#' informing the fit does not shrink as thresholds grow.
#'
#' A single mask (`T = 1`) reduces this to Algorithm 1.
#'
#' @seealso [cv.bmfsvt()] for choosing `alpha`, [lambda_max_seq()] for the
#'   shrinkage scale.
#' @importFrom stats plogis qlogis
#' @export
bmfsvt <- function(Y, Omega = NULL, lambda = NULL, alpha = 0.1,
                   max_iter = 200L, tol = 1e-5, gamma = 1.1, L_rewind = gamma,
                   delta = 1e-6,
                   rank_init = 5L, rank_max = NULL, rank_step = 2L,
                   clip = TRUE,
                   L_init = c("stacked", "single", "worst"),
                   restart = c("gradient", "function", "none"),
                   svd_method = c("auto", "full", "svds"),
                   track_objective = c("approx", "exact", "none"),
                   Z_init = NULL, max_backtrack = 30L, verbose = FALSE) {

  L_init <- match.arg(L_init)
  restart <- match.arg(restart)
  svd_method <- match.arg(svd_method)
  track_objective <- match.arg(track_objective)

  Y <- as_mask_list(Y, check_nested = TRUE)
  Tn <- length(Y)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])

  if (!is.null(Omega)) {
    storage.mode(Omega) <- "double"
    if (!identical(dim(Omega), c(N, P))) {
      stop("`Omega` must have the same dimensions as the binary masks.")
    }
    if (any(Omega != 0 & Omega != 1)) stop("`Omega` must contain only 0 and 1.")
  }
  nobs <- if (is.null(Omega)) N * P else sum(Omega)
  if (nobs == 0) stop("`Omega` marks no observed entries.")

  # ---- offsets and Lipschitz constants (Algorithm 2, lines 1-8) -------------
  pis <- vapply(Y, function(mat) {
    if (is.null(Omega)) mean(mat) else sum(Omega * mat) / nobs
  }, numeric(1L))
  eps <- 1 / (2 * nobs)
  if (any(pis <= 0 | pis >= 1)) {
    warning("some masks are entirely 0 or entirely 1; their offsets are clamped.")
  }
  pis <- pmin(pmax(pis, eps), 1 - eps)
  mu <- stats::qlogis(pis)
  nu <- if (Tn > 1L) diff(mu) else numeric(0L)

  L <- switch(
    L_init,
    single  = pis * (1 - pis),
    stacked = rev(cumsum(rev(pis * (1 - pis)))),
    worst   = (Tn - seq.int(0L, Tn - 1L)) / 4
  )
  L <- pmax(L, 1e-8)

  # ---- shrinkage parameters -------------------------------------------------
  if (is.null(lambda)) {
    lambda <- alpha * lambda_max_seq(Y, Omega)
  } else if (length(lambda) == 1L) {
    lambda <- rep(lambda, Tn)
  }
  if (length(lambda) != Tn) {
    stop("`lambda` must have length 1 or length(Y).")
  }

  # ---- initialization (lines 9-17) -----------------------------------------
  zero <- matrix(0, N, P)
  Z <- if (is.null(Z_init)) replicate(Tn, zero, simplify = FALSE) else {
    if (length(Z_init) != Tn) stop("`Z_init` must have length(Y) elements.")
    lapply(Z_init, function(mat) {
      storage.mode(mat) <- "double"
      mat
    })
  }
  W <- Z
  X <- stack_forward(W, mu[1L], nu)
  ranks <- rep(as.integer(rank_init), Tn)
  retained <- rep(NA_integer_, Tn)
  s_curr <- 1

  obj_trace <- numeric(0L)
  clip_frac_trace <- numeric(0L)
  n_svd <- 0L
  n_backtrack <- 0L
  converged <- FALSE
  m <- 0L

  clip_level <- if (Tn > 1L) -nu else numeric(0L)

  # ---- main loop (lines 18-57) ---------------------------------------------
  while (m < max_iter) {
    L <- L / L_rewind                                 # optimistic decrease
    s_next <- (1 + sqrt(1 + 4 * s_curr^2)) / 2

    # Psi depends only on the extrapolation point W, so it is computed once
    # per outer iteration rather than once per backtrack.
    Psi <- grad_backward(X, Y, Omega)
    f_W <- stack_ce(X, Y, Omega)

    bt <- 0L
    repeat {
      Zp <- vector("list", Tn)
      dvals <- vector("list", Tn)
      new_ranks <- ranks
      retained <- integer(Tn)
      n_clipped <- 0L
      for (t in seq_len(Tn)) {
        svt <- soft_svt(W[[t]] - Psi[[t]] / L[t], thresh = lambda[t] / L[t],
                        rank_guess = ranks[t], rank_max = rank_max,
                        rank_step = rank_step, method = svd_method)
        n_svd <- n_svd + svt$n_svd
        new_ranks[t] <- svt$rank_next
        retained[t] <- svt$rank
        mat <- svt$mat
        dvals[[t]] <- svt$d
        if (clip && t > 1L) {
          viol <- mat > clip_level[t - 1L]
          n_clipped <- n_clipped + sum(viol)
          if (any(viol)) mat[viol] <- clip_level[t - 1L]
        }
        Zp[[t]] <- mat
      }

      Xp <- stack_forward(Zp, mu[1L], nu)
      f_new <- stack_ce(Xp, Y, Omega)

      # Q(Z' | W) = f(X_W) + <Z' - W, Psi> + sum_t L_t/2 ||Z'^(t) - W^(t)||_F^2
      quad <- f_W
      for (t in seq_len(Tn)) {
        diff_t <- Zp[[t]] - W[[t]]
        quad <- quad + sum(diff_t * Psi[[t]]) + L[t] / 2 * sum(diff_t^2)
      }
      Delta <- f_new - quad

      if (Delta > 0 && bt < max_backtrack) {
        L <- gamma * L
        bt <- bt + 1L
        n_backtrack <- n_backtrack + 1L
      } else {
        if (Delta > 0) {
          warning(sprintf(
            "backtracking hit `max_backtrack` (%d) at iteration %d.",
            max_backtrack, m + 1L))
        }
        break
      }
    }
    ranks <- new_ranks

    if (track_objective != "none") {
      obj_trace <- c(obj_trace, f_new +
                       nuclear_penalty(Zp, dvals, lambda, n_clipped,
                                       track_objective))
      clip_frac_trace <- c(clip_frac_trace, n_clipped / (N * P * max(Tn - 1L, 1L)))
    }

    # ---- convergence (lines 44-47) -----------------------------------------
    num <- 0
    den <- 0
    for (t in seq_len(Tn)) {
      num <- num + sum((Zp[[t]] - Z[[t]])^2)
      den <- den + sum(Z[[t]]^2)
    }
    m <- m + 1L
    if (verbose) {
      message(sprintf(
        "iter %3d  f=%.6g  rel.change=%.3e  backtracks=%d  ranks=%s",
        m, f_new, num / (den + delta), bt, paste(ranks, collapse = ",")))
    }
    if (num / (den + delta) < tol) {
      Z <- Zp
      X <- Xp
      converged <- TRUE
      break
    }

    # ---- FISTA extrapolation with optional restart (lines 48-55) -----------
    do_restart <- switch(
      restart,
      none = FALSE,
      gradient = {
        # O'Donoghue and Candes generalized gradient scheme: restart when the
        # step moves against the proximal direction. Cheap insurance, and the
        # clipping step breaks the exact proximal-gradient structure that plain
        # FISTA assumes.
        ip <- 0
        for (t in seq_len(Tn)) ip <- ip + sum((W[[t]] - Zp[[t]]) * (Zp[[t]] - Z[[t]]))
        ip > 0
      },
      `function` = length(obj_trace) > 1L &&
        obj_trace[length(obj_trace)] > obj_trace[length(obj_trace) - 1L]
    )

    if (do_restart) {
      s_next <- 1
      W <- Zp
    } else {
      fac <- (s_curr - 1) / s_next
      W <- if (fac == 0) Zp else
        lapply(seq_len(Tn), function(t) Zp[[t]] + fac * (Zp[[t]] - Z[[t]]))
    }
    s_curr <- s_next
    Z <- Zp
    X <- stack_forward(W, mu[1L], nu)
  }

  if (!converged) X <- stack_forward(Z, mu[1L], nu)

  # ---- constraint diagnostics (lines 58-63) --------------------------------
  diag_bd <- violation_diagnostics(Z, X, Y, Omega, lambda, L, nu, ranks,
                                   rank_max, rank_step, svd_method)

  list(X = X,
       prob = lapply(X, sigmoid),
       Z = Z,
       mu = mu, nu = nu, pi = pis,
       lambda = lambda, L = L,
       ranks = retained,
       rank_guess = ranks,
       b = diag_bd$b, d = diag_bd$d,
       obj_trace = obj_trace,
       clip_frac_trace = clip_frac_trace,
       n_iter = m, n_svd = n_svd, n_backtrack = n_backtrack,
       converged = converged,
       clip = clip)
}


#' Stacked cross entropy over all masks
#' @param X list of natural parameter matrices
#' @param Y list of binary matrices
#' @param Omega observed entry mask or NULL
#' @returns the summed negative log likelihood
#' @keywords internal
stack_ce <- function(X, Y, Omega = NULL) {
  total <- 0
  for (t in seq_along(X)) total <- total + logistic_ce(X[[t]], Y[[t]], Omega)
  total
}


#' Nuclear norm penalty of the current blocks
#'
#' @param Z list of blocks
#' @param dvals singular values returned by the thresholding step, before clipping
#' @param lambda shrinkage parameters
#' @param n_clipped how many entries clipping actually altered
#' @param mode either "approx" or "exact"
#' @returns the value of \eqn{\sum_t \lambda_t \|Z^{(t)}\|_*}
#' @details
#' When no entry was clipped the shrunk singular values are exactly those of
#' `Z`, so the approximation is not an approximation at all. Recomputing is only
#' necessary once clipping has bent the blocks away from the SVT output.
#' @keywords internal
nuclear_penalty <- function(Z, dvals, lambda, n_clipped, mode = "approx") {
  if (mode == "approx" || n_clipped == 0L) {
    return(sum(vapply(seq_along(Z), function(t) lambda[t] * sum(dvals[[t]]),
                      numeric(1L))))
  }
  sum(vapply(seq_along(Z), function(t) {
    lambda[t] * sum(svd(Z[[t]], nu = 0L, nv = 0L)$d)
  }, numeric(1L)))
}


#' Violation rate and margin of an unclipped step from the solution
#'
#' @param Z solution blocks
#' @param X solution natural parameters
#' @param Y list of binary masks
#' @param Omega observed entry mask or NULL
#' @param lambda shrinkage parameters
#' @param L Lipschitz constants
#' @param nu offset increments
#' @param ranks rank guesses
#' @param rank_max,rank_step,svd_method passed to [soft_svt()]
#' @returns a list with numeric vectors `b` and `d` of length T-1
#' @details
#' Implements lines 58 to 63: one proximal step is taken from the solution
#' *without* clipping, and `b` records the fraction of observed entries that
#' would leave the feasible set while `d` records the largest such excursion on
#' the logit scale.
#'
#' These measure how **active** the monotonicity restriction is, not how badly
#' the clipping heuristic fails. At any constrained optimum with active
#' constraints the unconstrained gradient step must point out of the feasible
#' set, since that is what the KKT multiplier balances. Evaluating this same
#' statistic at an exact solution of the constrained program (see
#' `experiment/03_constraint.R`) returns very nearly the same values as
#' clip-SVT does, so a large `b` is a statement about the data, not a defect.
#'
#' To judge the heuristic itself, compare against a solution obtained some other
#' way: the objective gap and the difference in fitted probabilities are the
#' quantities that speak to suboptimality.
#' @keywords internal
violation_diagnostics <- function(Z, X, Y, Omega, lambda, L, nu, ranks,
                                  rank_max, rank_step, svd_method) {
  Tn <- length(Z)
  if (Tn < 2L) return(list(b = numeric(0L), d = numeric(0L)))

  nobs <- if (is.null(Omega)) length(Y[[1L]]) else sum(Omega)
  Psi <- grad_backward(X, Y, Omega)
  b <- numeric(Tn - 1L)
  d <- numeric(Tn - 1L)
  for (t in seq.int(2L, Tn)) {
    svt <- soft_svt(Z[[t]] - Psi[[t]] / L[t], thresh = lambda[t] / L[t],
                    rank_guess = ranks[t], rank_max = rank_max,
                    rank_step = rank_step, method = svd_method)
    margin <- svt$mat + nu[t - 1L]
    if (!is.null(Omega)) margin[Omega == 0] <- -Inf
    b[t - 1L] <- sum(margin > 0) / nobs
    d[t - 1L] <- max(0, max(margin))
  }
  names(b) <- names(d) <- paste0("t", seq.int(1L, Tn - 1L))
  list(b = b, d = d)
}


#' Choose the shrinkage level for bmfsvt by held out cross entropy
#'
#' @param Y a list of T binary matrices, or an N x P x T array
#' @param Omega an N x P binary matrix of observed entries, or `NULL`
#' @param alpha_grid decreasing grid of shrinkage levels relative to
#'   `lambda_max_seq()`. Fitted in the given order with warm starts.
#' @param prop_holdout fraction of observed entries held out for validation
#' @param seed optional seed for the holdout split
#' @param selection `"1se"` takes the strongest shrinkage whose validation loss
#'   is within one standard error of the best, `"min"` takes the outright
#'   minimum. The default is `"1se"` for the reason given in Details.
#' @param refit whether to refit at the selected `alpha` using all of `Omega`
#' @param verbose whether to report the validation loss of each `alpha`
#' @param ... further arguments passed to [bmfsvt()]
#'
#' @returns a list with the selected fit in `fit`, the chosen `alpha`, and a
#'   data frame `path` of validation losses.
#'
#' @details
#' The validation criterion is the one from the method note: the cross entropy
#' of the held out entries summed over every mask. The grid is traversed from
#' strong to weak shrinkage with each fit warm starting the next, which is why
#' the whole path costs little more than a handful of individual fits.
#'
#' **Held out validation systematically under-shrinks here, and the default
#' `"1se"` is a partial remedy rather than a fix.** Denoising differs from
#' matrix completion: the output of interest is a fitted value at every entry,
#' including the observed ones, but a held out criterion can only score entries
#' the fit never saw. As `alpha` falls, the fit increasingly reproduces the
#' realized Bernoulli draws at the entries it was trained on -- the correlation
#' between the fitted error and the coin flip noise at those entries rises from
#' about 0.03 to above 0.95 across the grid -- while at held out entries it
#' stays near zero. The held out loss is therefore nearly blind to the
#' overfitting that most damages the denoised values, and its minimum sits at a
#' weaker shrinkage than the one minimizing error against the truth.
#'
#' Switching the criterion from cross entropy to the Brier score does not help,
#' since both are evaluated on the same protected entries. A criterion that
#' scores all entries with a degrees-of-freedom penalty would address it
#' properly.
#'
#' @export
cv.bmfsvt <- function(Y, Omega = NULL,
                      alpha_grid = 10^seq(0, -2, length.out = 10L),
                      prop_holdout = 0.1, seed = NULL,
                      selection = c("1se", "min"), refit = TRUE,
                      verbose = FALSE, ...) {
  selection <- match.arg(selection)

  Y <- as_mask_list(Y, check_nested = TRUE)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  if (is.null(Omega)) Omega <- matrix(1, N, P)
  storage.mode(Omega) <- "double"

  if (!is.null(seed)) set.seed(seed)
  split <- holdout_split(Omega, prop_holdout)
  Omega_train <- split$train
  Omega_val <- split$val
  if (sum(Omega_val) == 0) {
    stop("the holdout split produced no validation entries; ",
         "increase `prop_holdout` or supply a denser `Omega`.")
  }

  alpha_grid <- sort(alpha_grid, decreasing = TRUE)
  lmax <- lambda_max_seq(Y, Omega_train)
  val_idx <- which(Omega_val == 1)

  losses <- ses <- numeric(length(alpha_grid))
  Z_warm <- NULL
  for (i in seq_along(alpha_grid)) {
    fit <- bmfsvt(Y, Omega = Omega_train, lambda = alpha_grid[i] * lmax,
                  Z_init = Z_warm, ...)
    Z_warm <- fit$Z
    # Per-observation losses, so the spread of the criterion is available and
    # not just its mean.
    ce <- unlist(lapply(seq_along(Y), function(t) {
      x <- fit$X[[t]][val_idx]
      log1exp(x) - Y[[t]][val_idx] * x
    }))
    losses[i] <- mean(ce)
    ses[i] <- stats::sd(ce) / sqrt(length(ce))
    if (verbose) {
      message(sprintf("alpha=%.4g  validation cross entropy=%.6f (se %.6f)",
                      alpha_grid[i], losses[i], ses[i]))
    }
  }

  best <- which.min(losses)
  if (selection == "1se") {
    # alpha_grid is sorted decreasing, so the first index inside the band is the
    # strongest shrinkage that is statistically indistinguishable from the best.
    within <- which(losses <= losses[best] + ses[best])
    best <- within[1L]
  }
  alpha_best <- alpha_grid[best]

  final <- if (refit) {
    bmfsvt(Y, Omega = Omega, lambda = alpha_best * lambda_max_seq(Y, Omega), ...)
  } else {
    bmfsvt(Y, Omega = Omega_train, lambda = alpha_best * lmax, ...)
  }

  list(fit = final,
       alpha = alpha_best,
       alpha_min = alpha_grid[which.min(losses)],
       selection = selection,
       path = data.frame(alpha = alpha_grid, val_loss = losses, se = ses),
       Omega_train = Omega_train, Omega_val = Omega_val)
}


#' Split observed entries into training and validation sets
#'
#' @param Omega binary matrix of observed entries
#' @param prop fraction to hold out
#' @returns a list with `train` and `val` binary matrices
#' @details
#' Entries are held out uniformly at random, then any row or column left with no
#' training entry has one returned to it. Without that repair, a fully held out
#' row makes its factor unidentifiable and the validation loss for that row
#' measures nothing but the intercept.
#' @keywords internal
holdout_split <- function(Omega, prop = 0.1) {
  idx <- which(Omega == 1)
  n_val <- floor(prop * length(idx))
  train <- Omega
  if (n_val == 0L) {
    return(list(train = train, val = matrix(0, nrow(Omega), ncol(Omega))))
  }
  val_idx <- sample(idx, n_val)
  train[val_idx] <- 0

  repair <- function(train, val_idx, margin) {
    counts <- if (margin == 1L) rowSums(train) else colSums(train)
    empty <- which(counts == 0)
    for (k in empty) {
      pos <- arrayInd(val_idx, dim(train))[, margin]
      cand <- which(pos == k)
      if (length(cand) > 0L) {
        give_back <- val_idx[cand[1L]]
        train[give_back] <- 1
        val_idx <- setdiff(val_idx, give_back)
      }
    }
    list(train = train, val_idx = val_idx)
  }
  r <- repair(train, val_idx, 1L)
  r <- repair(r$train, r$val_idx, 2L)

  val <- matrix(0, nrow(Omega), ncol(Omega))
  val[r$val_idx] <- 1
  list(train = r$train, val = val)
}
