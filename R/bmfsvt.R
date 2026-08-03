
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
#'     \item{dvals}{retained singular values of each block after shrinkage, as
#'       produced by the final thresholding step and therefore measured before
#'       clipping, matching the `"approx"` convention of [nuclear_penalty()].
#'       Supplied so that [svt_df()] never has to recompute a decomposition.}
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
  dvals <- replicate(Tn, numeric(0L), simplify = FALSE)
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
       dvals = dvals,
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


#' Fit a whole shrinkage path with warm starts
#'
#' @param Y a list of T binary matrices, or an N x P x T array
#' @param Omega an N x P binary matrix of observed entries, or `NULL`
#' @param alpha_grid shrinkage levels relative to [lambda_max_seq()]; sorted
#'   decreasing internally so that each fit warm starts the next
#' @param ... further arguments passed to [bmfsvt()]
#'
#' @returns a list holding the grid, the per-alpha increment matrices `Z`, the
#'   shared offsets `mu` and `nu`, and one row of diagnostics per alpha:
#'   `brier` (summed over the observed entries and every mask), `v_tilde`,
#'   `df` and `df_hard`, `ranks`, `n_iter`, `n_svd`, `converged`, `elapsed`.
#'
#' @details
#' Every selection criterion in the package reads this one object, so no two
#' rules are ever compared on different fits.
#'
#' Only `Z` is retained per alpha, not the full fit: `mu` and `nu` depend on the
#' mask prevalences alone and so are shared by the whole path, which makes
#' `sigma(stack_forward(Z, mu, nu))` enough to recover any fitted probability on
#' demand. That is three stored matrices per alpha rather than nine. It also
#' removes the need for a refit once an alpha is chosen: the path was fit on all
#' of `Omega` to begin with, so the selected fit is already in hand.
#' @keywords internal
bmfsvt_path <- function(Y, Omega = NULL, alpha_grid, ...) {
  Y <- as_mask_list(Y, check_nested = TRUE)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  Tn <- length(Y)
  if (is.null(Omega)) Omega <- matrix(1, N, P)
  storage.mode(Omega) <- "double"

  alpha_grid <- sort(alpha_grid, decreasing = TRUE)
  K <- length(alpha_grid)
  lmax <- lambda_max_seq(Y, Omega)
  obs_idx <- which(Omega == 1)

  Zs <- vector("list", K)
  brier <- v_tilde <- df_w <- df_h <- elapsed <- numeric(K)
  n_iter <- n_svd <- integer(K)
  converged <- logical(K)
  ranks <- matrix(NA_integer_, K, Tn)
  mu <- nu <- NULL

  Z_warm <- NULL
  for (k in seq_len(K)) {
    t0 <- proc.time()[["elapsed"]]
    fit <- bmfsvt(Y, Omega = Omega, lambda = alpha_grid[k] * lmax,
                  Z_init = Z_warm, ...)
    elapsed[k] <- proc.time()[["elapsed"]] - t0
    Z_warm <- fit$Z
    Zs[[k]] <- fit$Z
    if (is.null(mu)) {
      mu <- fit$mu
      nu <- fit$nu
    }
    lam <- alpha_grid[k] * lmax
    brier[k] <- sum(entry_loss(fit$X, Y, obs_idx, loss = "brier"))
    v_tilde[k] <- mean_bernoulli_var(lapply(fit$prob, function(m) m[obs_idx]))
    # The proximal step subtracts lambda_t / L_t, not lambda_t, so that is what
    # has to be added back to recover the pre-shrinkage singular values.
    df_w[k] <- sum(svt_df(fit$dvals, lam / fit$L, N, P, weighted = TRUE))
    df_h[k] <- sum(svt_df(fit$dvals, lam / fit$L, N, P, weighted = FALSE))
    ranks[k, ] <- fit$ranks
    n_iter[k] <- fit$n_iter
    n_svd[k] <- fit$n_svd
    converged[k] <- fit$converged
  }

  list(Y = Y, Omega = Omega, obs_idx = obs_idx,
       alpha = alpha_grid, lambda_max = lmax, Z = Zs, mu = mu, nu = nu,
       N = N, P = P, Tn = Tn,
       brier = brier, v_tilde = v_tilde, df = df_w, df_hard = df_h,
       ranks = ranks, n_iter = n_iter, n_svd = n_svd,
       converged = converged, elapsed = elapsed)
}


#' Fitted probability stack at one point of a shrinkage path
#'
#' @param path an object from [bmfsvt_path()]
#' @param k index into `path$alpha`
#' @returns a list of T probability matrices
#' @keywords internal
path_prob <- function(path, k) {
  lapply(stack_forward(path$Z[[k]], path$mu[1L], path$nu), sigmoid)
}


#' Mallows Cp surrogate for the shrinkage path
#'
#' @param path an object from [bmfsvt_path()]
#' @param v_ref fixed variance estimate; `NULL` takes it from the strongest
#'   shrinkage on the grid
#' @param weighted whether to use the shrinkage weighted degrees of freedom
#' @returns a list with the criterion `cp`, the selected index `index`, and the
#'   `v_ref` actually used
#'
#' @details
#' \eqn{M'(\alpha) = \sum_{\Omega}\sum_t (\hat\pi^{(t)}_{ij}(\alpha) -
#' Y^{(t)}_{ij})^2 + 2 \tilde V \, \mathrm{df}(\alpha)}, used to narrow the grid
#' before the bootstrap stage rather than to make the final choice.
#'
#' Two departures from a plain reading of the method note are deliberate, and
#' both were confirmed on a 16 cell simulation grid against an oracle that
#' minimizes RMSE against the true survival probabilities.
#'
#' **`v_ref` is not re-estimated at each alpha.** Doing so breaks the surrogate:
#' as the shrinkage weakens the fit drives its probabilities toward 0 and 1, so
#' the plug-in variance collapses (0.174 to 0.018 across the K = 10 grid) and
#' the penalty vanishes exactly where the degrees of freedom explode. The Brier
#' term collapses in step, so the criterion is minimized at the weakest
#' shrinkage on the grid -- in **16 of 16** cells -- and the candidate window
#' handed to the bootstrap never contains a sensible alpha. Measured RMSE 0.316
#' against 0.141 for the oracle. In Mallows Cp the variance is a fixed estimate,
#' not one re-read off each candidate model. The default takes it from the
#' intercept-only end of the grid, where \eqn{Z = 0} and
#' \eqn{\hat\pi^{(t)} = \pi_t}; that is parameter free and an upper bound on the
#' mean Bernoulli variance, so it errs toward stronger shrinkage.
#'
#' **`weighted = TRUE` by default.** The hard count \eqn{r_t(N+P-r_t)} treats
#' every retained singular value as a whole free parameter and over-penalizes;
#' it selected one of the two strongest alphas on the grid in **16 of 16**
#' cells, RMSE 0.197. See [svt_df()].
#' @keywords internal
cp_surrogate <- function(path, v_ref = NULL, weighted = TRUE) {
  if (is.null(v_ref)) v_ref <- path$v_tilde[1L]
  df <- if (weighted) path$df else path$df_hard
  cp <- path$brier + 2 * v_ref * df
  list(cp = cp, index = which.min(cp), v_ref = v_ref)
}


#' Parametric bootstrap estimate of the Brier plus covariance criterion
#'
#' @param path an object from [bmfsvt_path()]
#' @param index anchor index on the grid, normally the Cp surrogate's choice
#' @param B number of bootstrap replicates
#' @param window candidate alphas are those within this many grid steps of
#'   `index`
#' @param ncores replicates to run in parallel; 1 runs serially
#' @param seed optional seed making the replicates reproducible
#' @param ... further arguments passed to [bmfsvt()]
#'
#' @returns a list with the candidate indices `cand`, the criterion `M`, the
#'   covariance total `cov_total`, and the selected `index`
#'
#' @details
#' Estimates \eqn{\mathrm{Cov}(\hat\pi^{(t)}_{ij}(\alpha), Y^{(t)}_{ij})}
#' directly, so unlike [cp_surrogate()] it needs neither a degrees of freedom
#' proxy nor a variance estimate.
#'
#' Two implementation points carry real weight.
#'
#' **One uniform per entry, shared by every threshold.** Drawing independently
#' per mask would break nesting and the resampled stack would be rejected by
#' [as_mask_list()]; sharing the uniform is also what happens when a real
#' abundance matrix is thresholded.
#'
#' **The covariance is accumulated as running sums, never as an array.** Writing
#' \eqn{\sum_b (\hat\pi_b - \bar{\hat\pi})(Y_b - c) = S_3 - S_1 S_2 / B} with
#' \eqn{S_1 = \sum_b \hat\pi_b}, \eqn{S_2 = \sum_b Y_b}, \eqn{S_3 = \sum_b
#' \hat\pi_b Y_b} needs three accumulators instead of the B by N by P by T array
#' the defining formula suggests, which at proteomics scale would be terabytes.
#' The identity also shows the centering constant \eqn{c} cancels, because
#' \eqn{\sum_b (\hat\pi_b - \bar{\hat\pi}) = 0}; subtracting
#' \eqn{\hat\pi(\alpha^\dagger)} as the note does is therefore harmless but has
#' no effect on the result.
#'
#' Candidates are processed one at a time and \eqn{S_2}, which needs no fitting,
#' is accumulated once outside the workers. Each worker therefore carries two
#' accumulators rather than \eqn{2|\mathrm{cand}| + 1}, which is what keeps the
#' memory bounded once it is multiplied by the number of cores.
#'
#' Replicates are warm started from the real-data fit at the same alpha rather
#' than from the previous replicate. That is both a better starting point and
#' what keeps the replicates independent, so they can be split across cores.
#' @importFrom stats runif
#' @keywords internal
boot_covariance <- function(path, index, B = 30L, window = 2L, ncores = 1L,
                            seed = NULL, ...) {
  if (B < 2L) stop("`B` must be at least 2 to form a covariance.")
  K <- length(path$alpha)
  cand <- seq.int(max(1L, index - window), min(K, index + window))
  N <- path$N
  P <- path$P
  Tn <- path$Tn
  Omega <- path$Omega

  # Resampling needs a monotone probability stack or the draws are not nested.
  # Clipping already guarantees it, but `clip = FALSE` reaches here through
  # `...`, so enforce it rather than rely on the caller.
  pd <- path_prob(path, index)
  if (Tn > 1L) {
    for (t in seq.int(2L, Tn)) pd[[t]] <- pmin(pd[[t]], pd[[t - 1L]])
  }

  if (!is.null(seed)) set.seed(seed)
  boot_seeds <- sample.int(.Machine$integer.max, B)

  zero_stack <- function() replicate(Tn, matrix(0, N, P), simplify = FALSE)
  draw <- function(sd) {
    set.seed(sd)
    U <- matrix(stats::runif(N * P), N, P)
    lapply(pd, function(p) 1 * (U <= p))
  }
  add_stack <- function(a, b) Map(`+`, a, b)

  # S2 needs no fitting, so it is accumulated once here rather than inside every
  # worker. Regenerating the draws from their seeds costs a `runif` and keeps
  # the B masks out of memory.
  S2 <- zero_stack()
  for (sd in boot_seeds) S2 <- add_stack(S2, draw(sd))

  chunks <- split(boot_seeds, rep(seq_len(max(1L, ncores)), length.out = B))
  parallel_ok <- ncores > 1L && .Platform$OS.type != "windows"

  # Candidates are handled one at a time. Holding all of them in flight would
  # make each worker carry 2|cand| + 1 stacks instead of two, which at
  # proteomics scale is the difference between hundreds of megabytes and tens of
  # gigabytes once it is multiplied by the number of cores.
  cov_total <- vapply(cand, function(k) {
    one_chunk <- function(seeds) {
      S1 <- zero_stack()
      S3 <- zero_stack()
      for (sd in seeds) {
        Yb <- draw(sd)
        fb <- bmfsvt(Yb, Omega = Omega,
                     lambda = path$alpha[k] * path$lambda_max,
                     Z_init = path$Z[[k]], ...)
        for (t in seq_len(Tn)) {
          S1[[t]] <- S1[[t]] + fb$prob[[t]]
          S3[[t]] <- S3[[t]] + fb$prob[[t]] * Yb[[t]]
        }
      }
      list(S1 = S1, S3 = S3)
    }
    parts <- if (parallel_ok) {
      parallel::mclapply(chunks, one_chunk, mc.cores = ncores)
    } else {
      lapply(chunks, one_chunk)
    }
    bad <- vapply(parts, inherits, logical(1L), "try-error")
    if (any(bad)) {
      stop("bootstrap replicate failed: ",
           conditionMessage(attr(parts[[which(bad)[1L]]], "condition")))
    }
    S1 <- Reduce(add_stack, lapply(parts, `[[`, "S1"))
    S3 <- Reduce(add_stack, lapply(parts, `[[`, "S3"))
    tot <- 0
    for (t in seq_len(Tn)) {
      contrib <- S3[[t]] - S1[[t]] * S2[[t]] / B
      tot <- tot + sum(contrib[path$obs_idx])
    }
    tot / (B - 1)
  }, numeric(1L))

  M <- path$brier[cand] + 2 * cov_total
  list(cand = cand, M = M, cov_total = cov_total, B = B,
       index = cand[which.min(M)])
}


#' Choose the shrinkage level for bmfsvt
#'
#' @param Y a list of T binary matrices, or an N x P x T array
#' @param Omega an N x P binary matrix of observed entries, or `NULL`
#' @param alpha_grid shrinkage levels relative to [lambda_max_seq()]. The
#'   default is the method note's grid, \eqn{10^{-2k/(K-1)}} with `K = 10`.
#' @param criterion `"boot"` for the full Brier plus covariance criterion,
#'   `"cp"` to stop at the Mallows Cp surrogate, `"holdout"` for the superseded
#'   held out cross entropy rule
#' @param B bootstrap replicates used by `"boot"`
#' @param window how many grid steps either side of the surrogate's choice the
#'   bootstrap examines
#' @param weighted_df whether the surrogate uses shrinkage weighted degrees of
#'   freedom; `FALSE` reproduces the method note's hard count
#' @param ncores bootstrap replicates to run in parallel
#' @param seed optional seed for the bootstrap draws, or the holdout split
#' @param prop_holdout fraction held out, used only by `criterion = "holdout"`
#' @param verbose whether to report the criterion at each alpha
#' @param ... further arguments passed to [bmfsvt()]
#'
#' @returns a list with the fit at the chosen shrinkage in `fit`, the chosen
#'   `alpha`, the per-alpha data frame `path`, and, where the criterion computed
#'   them, `alpha_cp` and the bootstrap diagnostics in `boot`.
#'
#' @details
#' Implements section 2.3 of the method note. The criterion is scored on **all**
#' observed entries rather than on a held out subset, which is the right
#' estimand for denoising: the output of interest is a fitted value at every
#' entry including the ones the fit saw, and a held out criterion can only ever
#' score entries it did not. Starting from
#' \eqn{E[(\hat\pi - Y^*)^2] = E[(\hat\pi - Y)^2] + 2\,\mathrm{Cov}(\hat\pi, Y)}
#' for an independent copy \eqn{Y^*},
#' \deqn{M(\alpha) = \sum_{(i,j)\in\Omega}\sum_t (\hat\pi^{(t)}_{ij}(\alpha) -
#'   Y^{(t)}_{ij})^2 + 2\sum_{(i,j)\in\Omega}\sum_t
#'   \mathrm{Cov}(\hat\pi^{(t)}_{ij}(\alpha), Y^{(t)}_{ij}).}
#' The covariance is not available in closed form because every fitted value
#' depends on every entry, so it is estimated in two stages: [cp_surrogate()]
#' narrows the grid, then [boot_covariance()] estimates the covariance by
#' parametric bootstrap on the surviving candidates.
#'
#' `criterion = "holdout"` is the rule this replaced -- held out cross entropy,
#' scored per (entry, threshold) pair, with the strongest shrinkage inside one
#' standard error of the best. It is retained so the stages can be compared on
#' the same data.
#'
#' **What the comparison measured, on a 16 cell grid.** Mean absolute distance
#' from the oracle alpha, in grid steps, with the resulting RMSE against the
#' truth and the wall clock of the selection:
#'
#' \tabular{lrrr}{
#'   stage \tab steps \tab RMSE \tab secs \cr
#'   oracle \tab 0 \tab 0.1413 \tab 26 \cr
#'   `"holdout"` \tab 0.06 \tab 0.1413 \tab 26 \cr
#'   `"cp"` \tab 0.94 \tab 0.1646 \tab 26 \cr
#'   `"boot"` \tab 1.00 \tab 0.1642 \tab 176
#' }
#'
#' Read that honestly: on this grid the superseded held out rule is the most
#' accurate of the three, hitting the oracle in 15 of 16 cells, and the
#' bootstrap costs about seven times as much as the surrogate for no gain in
#' RMSE over it. The held out rule's advantage is not robust for the reason
#' given above and its band width scales as \eqn{n^{-1/4}}, so it is calibrated
#' to a particular validation set size rather than to the problem; but that is
#' an argument from mechanism, and the number above is what was measured.
#'
#' `"boot"` does buy one thing the others do not: it is much the closest to the
#' *true rank*, 3.96 against 7.61 for `"holdout"` and 10.60 for `"cp"` in mean
#' absolute per-block error. Rank recovery and estimation accuracy pull in
#' opposite directions here -- the oracle alpha itself retains about 7 more
#' components per block than the truth has, because soft thresholding shrinks
#' what it keeps. Choose the criterion by which of the two you need.
#'
#' @seealso [bmfsvt()] for the fit itself, [lambda_max_seq()] for the shrinkage
#'   scale.
#' @export
tune.bmfsvt <- function(Y, Omega = NULL,
                        alpha_grid = 10^seq(0, -2, length.out = 10L),
                        criterion = c("boot", "cp", "holdout"),
                        B = 30L, window = 2L, weighted_df = TRUE,
                        ncores = 1L, seed = NULL, prop_holdout = 0.1,
                        verbose = FALSE, ...) {
  criterion <- match.arg(criterion)
  Y <- as_mask_list(Y, check_nested = TRUE)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  if (is.null(Omega)) Omega <- matrix(1, N, P)
  storage.mode(Omega) <- "double"

  if (criterion == "holdout") {
    return(tune_holdout(Y, Omega, alpha_grid, prop_holdout, seed, verbose, ...))
  }

  path <- bmfsvt_path(Y, Omega = Omega, alpha_grid = alpha_grid, ...)
  cp <- cp_surrogate(path, weighted = weighted_df)
  index <- cp$index
  boot <- NULL
  M <- rep(NA_real_, length(path$alpha))

  if (verbose) {
    message(sprintf("Cp surrogate picks alpha=%.4g (v_ref=%.4f)",
                    path$alpha[cp$index], cp$v_ref))
  }

  if (criterion == "boot") {
    boot <- boot_covariance(path, index = cp$index, B = B, window = window,
                            ncores = ncores, seed = seed, ...)
    index <- boot$index
    M[boot$cand] <- boot$M
    if (verbose) {
      message(sprintf("bootstrap picks alpha=%.4g", path$alpha[index]))
    }
  }

  # Warm started from the stored increments, so this converges almost at once
  # and buys a complete fit object including the constraint diagnostics.
  final <- bmfsvt(Y, Omega = Omega,
                  lambda = path$alpha[index] * path$lambda_max,
                  Z_init = path$Z[[index]], ...)

  list(fit = final,
       alpha = path$alpha[index],
       index = index,
       criterion = criterion,
       alpha_cp = path$alpha[cp$index],
       v_ref = cp$v_ref,
       path = data.frame(alpha = path$alpha, brier = path$brier,
                         v_tilde = path$v_tilde, df = path$df,
                         df_hard = path$df_hard, cp = cp$cp, M = M,
                         mean_rank = rowMeans(path$ranks),
                         n_iter = path$n_iter, n_svd = path$n_svd,
                         converged = path$converged, elapsed = path$elapsed),
       ranks = path$ranks,
       boot = boot)
}


#' The superseded held out cross entropy rule
#'
#' @param Y list of masks, already coerced
#' @param Omega observation mask, already coerced
#' @param alpha_grid shrinkage grid
#' @param prop_holdout fraction held out
#' @param seed optional seed for the split
#' @param verbose whether to report each alpha
#' @param ... passed to [bmfsvt()]
#' @returns the same shape as [tune.bmfsvt()]
#' @details
#' Kept so the method note's three stages can be compared against what came
#' before them, on identical data. Scores cross entropy per (entry, threshold)
#' pair on the held out entries and takes the strongest shrinkage within one
#' standard error of the best. The standard error treats the T terms belonging
#' to one entry as independent, which they are not, so the band it produces is
#' narrower than an honest one; that is part of what is being compared.
#' @keywords internal
tune_holdout <- function(Y, Omega, alpha_grid, prop_holdout = 0.1, seed = NULL,
                         verbose = FALSE, ...) {
  if (!is.null(seed)) set.seed(seed)
  split <- holdout_split(Omega, prop_holdout)
  if (sum(split$val) == 0) {
    stop("the holdout split produced no validation entries; ",
         "increase `prop_holdout` or supply a denser `Omega`.")
  }
  alpha_grid <- sort(alpha_grid, decreasing = TRUE)
  lmax <- lambda_max_seq(Y, split$train)
  val_idx <- which(split$val == 1)

  losses <- ses <- numeric(length(alpha_grid))
  Z_warm <- NULL
  for (i in seq_along(alpha_grid)) {
    fit <- bmfsvt(Y, Omega = split$train, lambda = alpha_grid[i] * lmax,
                  Z_init = Z_warm, ...)
    Z_warm <- fit$Z
    ce <- unlist(lapply(seq_along(Y), function(t) {
      x <- fit$X[[t]][val_idx]
      log1exp(x) - Y[[t]][val_idx] * x
    }))
    losses[i] <- mean(ce)
    ses[i] <- stats::sd(ce) / sqrt(length(ce))
    if (verbose) {
      message(sprintf("alpha=%.4g  held out cross entropy=%.6f (se %.6f)",
                      alpha_grid[i], losses[i], ses[i]))
    }
  }

  # `index_min` is the raw argmin, which is the quantity the under-shrinkage
  # complaint is actually about; `index` adds the one standard error band, which
  # pushes back toward stronger shrinkage. Both are returned because the two can
  # differ by a grid step or more and reporting only the second hides whether
  # the band is doing the work.
  index_min <- which.min(losses)
  within <- which(losses <= losses[index_min] + ses[index_min])
  index <- within[1L]                     # grid is decreasing: strongest inside

  final <- bmfsvt(Y, Omega = Omega,
                  lambda = alpha_grid[index] * lambda_max_seq(Y, Omega), ...)

  list(fit = final,
       alpha = alpha_grid[index],
       index = index,
       index_min = index_min,
       alpha_min = alpha_grid[index_min],
       criterion = "holdout",
       alpha_cp = NA_real_,
       v_ref = NA_real_,
       path = data.frame(alpha = alpha_grid, val_loss = losses, se = ses),
       ranks = NULL,
       boot = NULL)
}
