
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
#' @param train optional N x P 0/1 matrix, 1 marking the entries that enter the
#'   likelihood, shared by every mask. `NULL` (the default) uses every entry.
#'   Held-out entries are excluded from the gradient (line 21 of Algorithm 2),
#'   from \eqn{f} and \eqn{Q} in the backtracking test (line 38), from the
#'   offsets and from the default noise floor, but they are still clipped (line
#'   27) and receive fitted probabilities, imputed from the low-rank blocks. This
#'   is what [cv.bmfsvt()] scores.
#' @param lambda shrinkage parameter for each mask; a scalar is recycled. When
#'   `NULL` the value `alpha * lambda_star_seq(Y, C)` is used.
#' @param alpha shrinkage relative to the **noise floor** [lambda_star_seq()],
#'   used only when `lambda` is `NULL`. Section 2.2 of the method note puts the
#'   search range at \eqn{\alpha \in (0, 1]}, with `alpha = 1` thresholding
#'   exactly at the floor.
#'
#'   **This changed meaning.** It was formerly a multiplier on
#'   [lambda_max_seq()], where the useful values sat near 0.3-0.6; it is now a
#'   multiplier on [lambda_star_seq()], where the measured optimum is about
#'   **0.55** (argmin in 6 of 6 pilot cells, two shapes x three seeds). A value
#'   carried over from the old parameterization will be wrong, not merely
#'   suboptimal.
#' @param C pure noise replicates used by [lambda_star_seq()], used only when
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
#' @param offset `"scalar"` uses one unpenalized offset per mask,
#'   \eqn{\mu_t = \mathrm{logit}(\bar Y^{(t)})}, as in the method note.
#'   `"column"` uses one per mask **and column**,
#'   \eqn{\mu_{tj} = \mathrm{logit}(\bar Y^{(t)}_{\cdot j})}, so `mu` becomes a
#'   T x P matrix and `nu` a (T-1) x P matrix. Use `"column"` whenever column
#'   prevalences differ strongly: a column main effect is expensive for a
#'   nuclear-norm penalized block to absorb, and on the HRS composition data the
#'   scalar offsets fitted a column that is 49% nonzero at 0.81 on average. When
#'   the masks come from per-column quantile thresholds the column prevalences
#'   are identical except for zero-heavy columns, so the two settings then
#'   coincide everywhere else. [lambda_star_seq()] and the Cp variance stay
#'   pooled under either setting.
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
#'     \item{mu, nu}{offsets and offset increments: vectors of length T and
#'       T-1 under `offset = "scalar"`, T x P and (T-1) x P matrices under
#'       `"column"`}
#'     \item{offset}{which offset parameterization was used}
#'     \item{V}{list of T `P x ranks[t]` matrices of right singular vectors
#'       paired with `dvals`. Like `dvals` they describe the thresholded block
#'       before clipping. Used by [predict.bmfsvt()] to fold in new rows.}
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
#'     \item{n_train}{number of entries per mask that entered the likelihood}
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
#' The returned object has class `"bmfsvt"`.
#'
#' @seealso [cv.bmfsvt()] and [tune.bmfsvt()] for choosing `alpha`,
#'   [lambda_star_seq()] for the shrinkage scale, [predict.bmfsvt()] for new
#'   rows.
#' @importFrom stats plogis qlogis
#' @export
bmfsvt <- function(Y, train = NULL, lambda = NULL, alpha = 0.55, C = 20L,
                   max_iter = 200L, tol = 1e-5, gamma = 1.1, L_rewind = gamma,
                   delta = 1e-6,
                   rank_init = 5L, rank_max = NULL, rank_step = 2L,
                   clip = TRUE, offset = c("scalar", "column"),
                   L_init = c("stacked", "single", "worst"),
                   restart = c("gradient", "function", "none"),
                   svd_method = c("auto", "full", "svds"),
                   track_objective = c("approx", "exact", "none"),
                   Z_init = NULL, max_backtrack = 30L, verbose = FALSE) {

  offset <- match.arg(offset)
  L_init <- match.arg(L_init)
  restart <- match.arg(restart)
  svd_method <- match.arg(svd_method)
  track_objective <- match.arg(track_objective)
  # With no iteration the returned increments would carry no singular vectors.
  if (max_iter < 1L) stop("`max_iter` must be at least 1.")

  Y <- as_mask_list(Y, check_nested = TRUE)
  Tn <- length(Y)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  train <- check_train(train, N, P)
  nobs <- if (is.null(train)) N * P else sum(train)
  if (!is.null(train)) {
    n_row <- rowSums(train)
    n_col <- colSums(train)
    if (any(n_col == 0) && offset == "column") {
      stop("`train` leaves some column with no training entry, so its ",
           "column offset is undefined.")
    }
    if (any(n_row == 0) || any(n_col == 0)) {
      warning("`train` leaves ", sum(n_row == 0), " row(s) and ",
              sum(n_col == 0), " column(s) with no training entry; they are ",
              "fitted from the offsets alone.")
    }
  }

  # ---- offsets and Lipschitz constants (Algorithm 2, lines 1-8) -------------
  # Taken over the training entries only, or the held-out values would leak
  # into the fit through the offsets.
  pis <- mask_prevalence(Y, train)
  eps <- 1 / (2 * nobs)
  if (any(pis <= 0 | pis >= 1)) {
    warning("some masks are entirely 0 or entirely 1; their offsets are clamped.")
  }
  pis <- pmin(pmax(pis, eps), 1 - eps)
  if (offset == "scalar") {
    mu <- stats::qlogis(pis)
    nu <- if (Tn > 1L) diff(mu) else numeric(0L)
  } else {
    # Column prevalences are non-increasing in t because the masks are nested,
    # and clamping preserves the order, so nu <= 0 still holds entrywise.
    eps_col <- 1 / (2 * N)
    pic <- if (is.null(train)) {
      matrix(vapply(Y, colMeans, numeric(P)), P, Tn)
    } else {
      matrix(vapply(Y, function(mat) colSums(train * mat) / n_col,
                    numeric(P)), P, Tn)
    }
    pic <- pmin(pmax(pic, eps_col), 1 - eps_col)
    mu <- t(stats::qlogis(pic))
    dimnames(mu) <- NULL
    nu <- mu[-1L, , drop = FALSE] - mu[-Tn, , drop = FALSE]
  }

  L <- switch(
    L_init,
    single  = pis * (1 - pis),
    stacked = rev(cumsum(rev(pis * (1 - pis)))),
    worst   = (Tn - seq.int(0L, Tn - 1L)) / 4
  )
  L <- pmax(L, 1e-8)

  # ---- shrinkage parameters -------------------------------------------------
  # The anchor is the noise floor, not lambda_max: the whole interval
  # [lambda*, lambda_max] over-shrinks, so a fraction of lambda* is what the
  # search range needs to be. See lambda_star_seq().
  if (is.null(lambda)) {
    lambda <- alpha * lambda_star_seq(Y, C = C, train = train)
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
  mu0 <- offset_head(mu)
  X <- stack_forward(W, mu0, nu)
  ranks <- rep(as.integer(rank_init), Tn)
  retained <- rep(NA_integer_, Tn)
  dvals <- replicate(Tn, numeric(0L), simplify = FALSE)
  Vs <- replicate(Tn, matrix(0, P, 0L), simplify = FALSE)
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
    Psi <- grad_backward(X, Y, train)
    f_W <- stack_ce(X, Y, train)

    bt <- 0L
    repeat {
      Zp <- vector("list", Tn)
      dvals <- vector("list", Tn)
      # Reset alongside dvals so that V always pairs with the accepted step,
      # whichever way the loop exits.
      Vs <- vector("list", Tn)
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
        Vs[[t]] <- svt$v
        if (clip && t > 1L) {
          level <- if (is.matrix(clip_level)) {
            rep(clip_level[t - 1L, ], each = N)
          } else {
            clip_level[t - 1L]
          }
          viol <- mat > level
          n_clipped <- n_clipped + sum(viol)
          if (any(viol)) {
            mat[viol] <- if (length(level) > 1L) level[viol] else level
          }
        }
        Zp[[t]] <- mat
      }

      Xp <- stack_forward(Zp, mu0, nu)
      f_new <- stack_ce(Xp, Y, train)

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
    X <- stack_forward(W, mu0, nu)
  }

  if (!converged) X <- stack_forward(Z, mu0, nu)

  # ---- constraint diagnostics (lines 58-63) --------------------------------
  diag_bd <- violation_diagnostics(Z, X, Y, lambda, L, nu, ranks,
                                   rank_max, rank_step, svd_method,
                                   train = train)

  structure(
    list(X = X,
         prob = lapply(X, sigmoid),
         Z = Z,
         mu = mu, nu = nu, pi = pis, offset = offset,
         lambda = lambda, L = L,
         ranks = retained,
         rank_guess = ranks,
         dvals = dvals,
         V = Vs,
         b = diag_bd$b, d = diag_bd$d,
         obj_trace = obj_trace,
         clip_frac_trace = clip_frac_trace,
         n_iter = m, n_svd = n_svd, n_backtrack = n_backtrack,
         converged = converged,
         n_train = nobs,
         clip = clip),
    class = "bmfsvt")
}


#' Offset of the first mask
#'
#' @param mu offsets, a vector of length T or a T x P matrix
#' @returns `mu[1]`, or the first row of `mu` for column offsets
#' @keywords internal
offset_head <- function(mu) {
  if (is.matrix(mu)) mu[1L, ] else mu[1L]
}


#' Stacked cross entropy over all masks
#' @param X list of natural parameter matrices
#' @param Y list of binary matrices
#' @param train optional training mask, see [logistic_ce()]
#' @returns the summed negative log likelihood
#' @keywords internal
stack_ce <- function(X, Y, train = NULL) {
  total <- 0
  for (t in seq_along(X)) total <- total + logistic_ce(X[[t]], Y[[t]], train)
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
#' @param lambda shrinkage parameters
#' @param L Lipschitz constants
#' @param nu offset increments
#' @param ranks rank guesses
#' @param rank_max,rank_step,svd_method passed to [soft_svt()]
#' @param train optional training mask; the gradient is taken over it, as in
#'   the fit, while the rate `b` counts every entry since every entry is clipped
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
violation_diagnostics <- function(Z, X, Y, lambda, L, nu, ranks,
                                  rank_max, rank_step, svd_method,
                                  train = NULL) {
  Tn <- length(Z)
  if (Tn < 2L) return(list(b = numeric(0L), d = numeric(0L)))

  nobs <- length(Y[[1L]])
  Psi <- grad_backward(X, Y, train)
  b <- numeric(Tn - 1L)
  d <- numeric(Tn - 1L)
  for (t in seq.int(2L, Tn)) {
    svt <- soft_svt(Z[[t]] - Psi[[t]] / L[t], thresh = lambda[t] / L[t],
                    rank_guess = ranks[t], rank_max = rank_max,
                    rank_step = rank_step, method = svd_method)
    nu_t <- if (is.matrix(nu)) rep(nu[t - 1L, ], each = nrow(Z[[t]])) else
      nu[t - 1L]
    margin <- svt$mat + nu_t
    b[t - 1L] <- sum(margin > 0) / nobs
    d[t - 1L] <- max(0, max(margin))
  }
  names(b) <- names(d) <- paste0("t", seq.int(1L, Tn - 1L))
  list(b = b, d = d)
}


#' Balanced right factors of the fitted blocks
#'
#' @param V list of T `P x r_t` right singular vector matrices, as in `fit$V`
#' @param dvals list of T shrunken singular value vectors, as in `fit$dvals`
#' @returns list of T `P x r_t` matrices \eqn{B_t = V_t D_t^{1/2}}
#' @details
#' The square root is not a normalization choice. The nuclear norm has the
#' variational form \eqn{\|Z\|_* = \min_{Z = AB'} (\|A\|_F^2 + \|B\|_F^2)/2},
#' attained only at the balanced split \eqn{A = UD^{1/2}}, \eqn{B = VD^{1/2}}.
#' Only with that `B` is a plain ridge on a new row's scores that row's share of
#' the nuclear norm penalty. With \eqn{B = VD} the same ridge would shrink weak
#' components harder than strong ones; folding the training rows back in that
#' way was measured to miss the fit by 11.7 on the logit scale.
#'
#' Built with `rep(each = )` rather than `diag()`, which for a single component
#' returns an identity matrix of size `floor(sqrt(d))` instead.
#' @keywords internal
factor_blocks <- function(V, dvals) {
  lapply(seq_along(V), function(t) {
    V[[t]] * rep(sqrt(dvals[[t]]), each = nrow(V[[t]]))
  })
}


#' Ridge fold-in of new rows against fixed right factors
#'
#' @param Y list of T binary n x P masks for the new rows
#' @param B list of T `P x r_t` balanced right factors, see [factor_blocks()]
#' @param mu0 offset of the first mask, scalar or length P
#' @param nu offset increments, vector of length T-1 or (T-1) x P matrix
#' @param lambda the training shrinkage parameters, length T
#' @param tol a row stops once the Euclidean norm of its own gradient falls
#'   below this
#' @param max_iter maximum number of iterations
#' @returns a list with the scores `A` (list of T `n x r_t` matrices), the
#'   unclipped increments `Z` (list of T n x P matrices), and per row
#'   `grad_norm`, `n_iter` and `converged`
#' @details
#' Solves, for the new rows only,
#' \deqn{\min_{A_1..A_T} \sum_t \mathrm{CE}(X^{(t)}; Y^{(t)}) +
#'   \sum_t \frac{\lambda_t}{2}\|A_t\|_F^2, \qquad
#'   X^{(t)} = \mu_t + \sum_{t' \le t} A_{t'} B_{t'}^\top,}
#' whose gradient is \eqn{\Psi_t B_t + \lambda_t A_t} with \eqn{\Psi} from
#' [grad_backward()]. The problem separates by row and is strongly convex.
#'
#' **Step size.** The Hessian block for the pair (t, s) is
#' \eqn{B_t^\top \mathrm{diag}(\sum_{t' \ge \max(s,t)} w_{t'}) B_s} with
#' \eqn{w \le 1/4}, so its norm is at most
#' \eqn{(T - \max(s,t) + 1)/4 \cdot \sqrt{d_{t,1} d_{s,1}}}. Summing over s gives a
#' block Lipschitz constant \eqn{L_t} valid for every row at every iterate, so
#' no backtracking is needed; the step is \eqn{1/(L_t + \lambda_t)}.
#'
#' **Each row is solved independently.** Momentum and restarts are kept per
#' row, and a row is frozen at the iterate where its own gradient norm first
#' drops below `tol`. The answer for a row therefore does not depend on which
#' other rows share the batch; batching is only for speed.
#' @keywords internal
foldin_rows <- function(Y, B, mu0, nu, lambda, tol = 1e-6, max_iter = 2000L) {
  Tn <- length(Y)
  n <- nrow(Y[[1L]])
  ranks <- vapply(B, ncol, integer(1L))
  tt <- seq_len(Tn)

  # ||B_t||_2^2 is the largest shrunken singular value, since V is orthonormal.
  nb <- vapply(B, function(b) if (ncol(b) == 0L) 0 else sqrt(max(colSums(b^2))),
               numeric(1L))
  weight <- outer(tt, tt, function(s, t) (Tn - pmax(s, t) + 1) / 4)
  L <- pmax(rowSums(weight * outer(nb, nb)), 1e-8)
  step <- 1 / (L + lambda)

  A_out <- lapply(ranks, function(r) matrix(0, n, r))
  grad_norm <- rep(NA_real_, n)
  n_iter <- rep(as.integer(max_iter), n)
  converged <- logical(n)

  rows <- seq_len(n)
  Yw <- Y
  A <- lapply(ranks, function(r) matrix(0, n, r))
  W <- A
  s <- rep(1, n)
  gn <- rep(NA_real_, n)

  for (m in seq_len(max_iter)) {
    Zw <- lapply(tt, function(t) tcrossprod(W[[t]], B[[t]]))
    Psi <- grad_backward(stack_forward(Zw, mu0, nu), Yw)
    G <- lapply(tt, function(t) Psi[[t]] %*% B[[t]] + lambda[t] * W[[t]])
    gn <- sqrt(Reduce(`+`, lapply(G, function(g) rowSums(g^2))))

    done <- gn < tol
    if (any(done)) {
      idx <- rows[done]
      for (t in tt) {
        if (ranks[t] > 0L) A_out[[t]][idx, ] <- W[[t]][done, , drop = FALSE]
      }
      grad_norm[idx] <- gn[done]
      n_iter[idx] <- m - 1L
      converged[idx] <- TRUE
      keep <- !done
      rows <- rows[keep]
      if (length(rows) == 0L) break
      Yw <- lapply(Yw, function(y) y[keep, , drop = FALSE])
      A <- lapply(A, function(a) a[keep, , drop = FALSE])
      W <- lapply(W, function(w) w[keep, , drop = FALSE])
      G <- lapply(G, function(g) g[keep, , drop = FALSE])
      s <- s[keep]
      gn <- gn[keep]
    }

    A_new <- lapply(tt, function(t) W[[t]] - step[t] * G[[t]])
    # Gradient restart, per row: drop the momentum of any row whose step moved
    # against its proximal direction.
    ip <- Reduce(`+`, lapply(tt, function(t) {
      rowSums((W[[t]] - A_new[[t]]) * (A_new[[t]] - A[[t]]))
    }))
    s_next <- (1 + sqrt(1 + 4 * s^2)) / 2
    fac <- (s - 1) / s_next
    restart <- ip > 0
    fac[restart] <- 0
    s_next[restart] <- 1
    # `fac` has one element per row, so recycling it down the columns of an
    # n x r matrix scales each row by its own factor.
    W <- lapply(tt, function(t) A_new[[t]] + fac * (A_new[[t]] - A[[t]]))
    A <- A_new
    s <- s_next
  }

  if (length(rows) > 0L) {
    for (t in tt) {
      if (ranks[t] > 0L) A_out[[t]][rows, ] <- A[[t]]
    }
    grad_norm[rows] <- gn
    warning(sprintf(
      "fold-in did not reach `tol` for %d of %d rows within %d iterations.",
      length(rows), n, max_iter))
  }

  Z <- lapply(tt, function(t) tcrossprod(A_out[[t]], B[[t]]))
  list(A = A_out, Z = Z, grad_norm = grad_norm, n_iter = n_iter,
       converged = converged)
}


#' Denoise new rows with a fitted factorization
#'
#' @param object a fit from [bmfsvt()], e.g. `tune.bmfsvt(...)$fit`
#' @param newdata a list of T binary n x P matrices or an n x P x T array. The
#'   masks must be built with the **training** thresholds, so that mask t means
#'   the same event it meant when `object` was fitted.
#' @param tol per-row gradient norm at which a row stops
#' @param max_iter maximum number of fold-in iterations
#' @param ... unused
#'
#' @returns a list with components
#'   \describe{
#'     \item{prob}{list of T n x P matrices, the estimated
#'       \eqn{P(\text{value} > d_t)} for each new entry}
#'     \item{X, Z}{natural parameters and (clipped) increments}
#'     \item{A}{list of T `n x r_t` row score matrices}
#'     \item{row_ce}{mean cross entropy of each row's fitted probabilities
#'       against its own masks, a per-row misfit diagnostic}
#'     \item{grad_norm, n_iter, converged}{per-row solver diagnostics}
#'   }
#'
#' @details
#' Every training quantity is held fixed: the offsets `mu` and `nu`, the
#' shrinkage `lambda`, and the right factors \eqn{B_t = V_t D_t^{1/2}}. Only the
#' new rows' scores are solved for, by the ridge problem in [foldin_rows()].
#' Nothing is shared across new rows, so each is denoised on its own, the way a
#' principal component projection treats a new sample.
#'
#' **Why this is the fit's own map.** With \eqn{B} fixed the training problem
#' separates by row, and the stationarity condition of each row's ridge problem,
#' \eqn{\Psi_t V_t = -\lambda_t U_t}, is exactly the nuclear norm optimality
#' condition. Folding the training rows back into an unclipped fit therefore
#' reproduces `fit$X` up to solver tolerance. Under squared error loss the same
#' construction has the closed form \eqn{c_r = d_r/(d_r + \lambda)(V^\top x)_r}:
#' project onto `V`, then shrink each component. Logistic loss has no closed
#' form, hence the iterative solve.
#'
#' **Clipping is applied after solving, not during.** When `object` was fitted
#' with `clip = TRUE`, the increments are capped at \eqn{-\nu_t} once the ridge
#' problem is solved. Clipping inside the iterations would make each row's
#' problem non-convex. The in-sample fit clips inside its iterations, so on
#' training rows the two differ slightly (3.5e-5 mean absolute probability on
#' the HRS data).
#'
#' **What the fold-in cannot do.** A row can only be expressed through the
#' column patterns in `V`, which are shared across the training population. An
#' idiosyncratic profile carried by one sample is not represented; see the
#' flagging in [dfbm()].
#' @importFrom stats predict
#' @export
predict.bmfsvt <- function(object, newdata, tol = 1e-6, max_iter = 2000L, ...) {
  if (is.null(object$V)) {
    stop("`object` stores no singular vectors; refit it with the current bmfsvt().")
  }
  Y <- as_mask_list(newdata, check_nested = TRUE)
  Tn <- length(object$V)
  P <- nrow(object$V[[1L]])
  if (length(Y) != Tn) {
    stop(sprintf("`newdata` has %d masks but the fit has %d.", length(Y), Tn))
  }
  if (ncol(Y[[1L]]) != P) {
    stop(sprintf("`newdata` has %d columns but the fit has %d.",
                 ncol(Y[[1L]]), P))
  }

  mu0 <- offset_head(object$mu)
  B <- factor_blocks(object$V, object$dvals)
  sol <- foldin_rows(Y, B, mu0, object$nu, object$lambda,
                     tol = tol, max_iter = max_iter)
  Z <- if (isTRUE(object$clip)) clip_increments(sol$Z, object$nu) else sol$Z
  X <- stack_forward(Z, mu0, object$nu)

  ce <- Reduce(`+`, lapply(seq_len(Tn), function(t) {
    rowSums(log1exp(X[[t]]) - Y[[t]] * X[[t]])
  }))

  list(prob = lapply(X, sigmoid), X = X, Z = Z, A = sol$A,
       row_ce = ce / (Tn * P),
       grad_norm = sol$grad_norm, n_iter = sol$n_iter,
       converged = sol$converged)
}




#' Fit a whole shrinkage path with warm starts
#'
#' @param Y a list of T binary matrices, or an N x P x T array
#' @param alpha_grid shrinkage levels relative to [lambda_star_seq()]; sorted
#'   decreasing internally so that each fit warm starts the next
#' @param C pure noise replicates for [lambda_star_seq()]
#' @param seed optional integer making the [lambda_star_seq()] draw reproducible
#' @param ... further arguments passed to [bmfsvt()]
#'
#' @returns a list holding the grid, the per-alpha increment matrices `Z`, the
#'   shared offsets `mu` and `nu`, and one row of diagnostics per alpha:
#'   `brier` (summed over every entry and every mask), `v_tilde`, `df` and
#'   `df_hard` (each a K x T **matrix**, per slice), `ranks`, `n_iter`, `n_svd`,
#'   `converged`, `elapsed`. `v_t` holds the per-slice surrogate variance
#'   \eqn{\pi_t(1-\pi_t)}, which depends on the data alone and not on any fit.
#'
#' @details
#' Every selection criterion in the package reads this one object, so no two
#' rules are ever compared on different fits.
#'
#' Only `Z` is retained per alpha, not the full fit: `mu` and `nu` depend on the
#' mask prevalences alone and so are shared by the whole path, which makes
#' `sigma(stack_forward(Z, mu, nu))` enough to recover any fitted probability on
#' demand. That is three stored matrices per alpha rather than nine. It also
#' removes the need for a refit once an alpha is chosen: the path was fit on
#' every entry to begin with, so the selected fit is already in hand.
#'
#' `lambda_star_seq()` is evaluated **once** here, not once per candidate: the
#' whole grid rescales the same vector, and the Monte Carlo cost is `C` partial
#' decompositions per threshold.
#'
#' `df` is kept per slice rather than summed because the Cp penalty of section
#' 2.3 weights each slice by its own \eqn{\tilde V_t}; see [cp_surrogate()].
#' @keywords internal
bmfsvt_path <- function(Y, alpha_grid, C = 20L, seed = NULL, ...) {
  Y <- as_mask_list(Y, check_nested = TRUE)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  Tn <- length(Y)

  alpha_grid <- sort(alpha_grid, decreasing = TRUE)
  K <- length(alpha_grid)
  lstar <- lambda_star_seq(Y, C = C, seed = seed)

  # The surrogate variance of section 2.3, one per slice, from the marginal
  # prevalences. It must NOT be read off a fit: see cp_surrogate().
  pis <- vapply(Y, mean, numeric(1L))
  v_t <- pis * (1 - pis)

  Zs <- vector("list", K)
  brier <- v_tilde <- elapsed <- numeric(K)
  df_w <- df_h <- matrix(NA_real_, K, Tn)
  n_iter <- n_svd <- integer(K)
  converged <- logical(K)
  ranks <- matrix(NA_integer_, K, Tn)
  mu <- nu <- NULL

  Z_warm <- NULL
  for (k in seq_len(K)) {
    t0 <- proc.time()[["elapsed"]]
    fit <- bmfsvt(Y, lambda = alpha_grid[k] * lstar, Z_init = Z_warm, ...)
    elapsed[k] <- proc.time()[["elapsed"]] - t0
    Z_warm <- fit$Z
    Zs[[k]] <- fit$Z
    if (is.null(mu)) {
      mu <- fit$mu
      nu <- fit$nu
    }
    lam <- alpha_grid[k] * lstar
    brier[k] <- sum(entry_loss(fit$X, Y))
    # Retained only for the `cp_note_asis` control, which re-estimates the
    # variance at every alpha in order to demonstrate that doing so fails.
    v_tilde[k] <- mean_bernoulli_var(fit$prob)
    # The proximal step subtracts lambda_t / L_t, not lambda_t, so that is what
    # has to be added back to recover the pre-shrinkage singular values.
    df_w[k, ] <- svt_df(fit$dvals, lam / fit$L, N, P, weighted = TRUE)
    df_h[k, ] <- svt_df(fit$dvals, lam / fit$L, N, P, weighted = FALSE)
    ranks[k, ] <- fit$ranks
    n_iter[k] <- fit$n_iter
    n_svd[k] <- fit$n_svd
    converged[k] <- fit$converged
  }

  list(Y = Y,
       alpha = alpha_grid, lambda_star = lstar, Z = Zs, mu = mu, nu = nu,
       N = N, P = P, Tn = Tn,
       brier = brier, v_tilde = v_tilde, v_t = v_t,
       df = df_w, df_hard = df_h,
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
  lapply(stack_forward(path$Z[[k]], offset_head(path$mu), path$nu), sigmoid)
}


#' Mallows Cp surrogate for the shrinkage path
#'
#' @param path an object from [bmfsvt_path()]
#' @param v_ref fixed per-slice variance vector of length T; `NULL` uses
#'   `path$v_t`, i.e. \eqn{\pi_t(1-\pi_t)} from the marginal prevalences. A
#'   scalar is recycled, which reproduces the superseded pooled form.
#' @param weighted whether to use the shrinkage weighted degrees of freedom
#' @returns a list with the criterion `cp`, the selected index `index`, and the
#'   `v_ref` actually used
#'
#' @details
#' \eqn{M'(\alpha) = \sum_{ij}\sum_t (\hat\pi^{(t)}_{ij}(\alpha) -
#' Y^{(t)}_{ij})^2 + 2 \sum_t \tilde V_t \, \mathrm{df}_t(\alpha)}, section 2.3
#' of the method note.
#'
#' Three departures from a plain reading of an earlier draft are deliberate.
#'
#' **The penalty is a sum of products, not a product of sums.** Each slice is
#' weighted by its own \eqn{\tilde V_t = \pi_t(1-\pi_t)} rather than by a pooled
#' average times the total `df`. Pooling is a poor approximation here because
#' \eqn{\tilde V_t} and \eqn{\mathrm{df}_t} are *correlated* across `t`
#' (measured \eqn{+0.49} at the oracle alpha): \eqn{\tilde V_t} is hump shaped,
#' peaking where the prevalences cross 0.5, and `df` is front loaded under the
#' [lambda_star_seq()] anchor. **The per-slice form therefore makes the penalty
#' larger, not smaller** -- 9% at the oracle, 18% at weak shrinkage, since
#' Chebyshev's sum inequality runs the other way for like-ordered sequences.
#' Its value is conditioning, not level: the pooled criterion was flat across
#' half the grid (a 0.8% range over five candidates, its argmin winning by 0.1%,
#' i.e. noise), while the per-slice one has 10x to 20x more curvature and a
#' decisive minimum. **It is not a fix for over-shrinkage** and should not be
#' described as one.
#'
#' **`v_ref` is not re-estimated at each alpha.** Doing so breaks the surrogate:
#' as the shrinkage weakens the fit drives its probabilities toward 0 and 1, so
#' the plug-in variance collapses (0.174 to 0.018 across the K = 10 grid) and
#' the penalty vanishes exactly where the degrees of freedom explode. It then
#' selects the weakest alpha on the grid in **10 of 10** cells, retaining rank
#' 17.3 to 45.6 against a truth near 3.1. In Mallows Cp the variance is a fixed
#' estimate, not one re-read off each candidate model.
#'
#' Note that the variance comes from the *prevalences*, never from a fit on the
#' path. Reading it off the first candidate used to be equivalent, because the
#' old grid's first point (`alpha = 1` on [lambda_max_seq()]) collapsed every
#' block to zero. Under the [lambda_star_seq()] anchor `alpha = 1` still retains
#' rank, so that shortcut is now silently wrong.
#'
#' \eqn{\pi_t(1-\pi_t)} at the marginal prevalence is an upper bound on the true
#' mean Bernoulli variance, by concavity of \eqn{\pi(1-\pi)} and Jensen --
#' measured **1.40x** too large -- so the penalty is conservative and errs toward
#' stronger shrinkage, which is the direction the criterion is measurably biased.
#'
#' **Do not try to fix that bias by shrinking \eqn{\tilde V_t}, or by rescaling
#' `df`, or by swapping \eqn{L_t} for the average curvature. All three were
#' measured and none works.** Sweeping a scale `c` in
#' \eqn{\mathrm{Brier} + 2c\,\mathrm{penalty}} over `[0.05, 3]`, the selected
#' index runs 10, 9, 7, 6, 5, 3, 2, 1 on the square shape and 10, 7, 2, 1 on HRS:
#' the oracle is **never** selected, for any `c`. It is not a vertex of the lower
#' convex hull of the (penalty, Brier) points, and every constant reweighting --
#' the 1.40x Jensen factor included, which lands at `c ~ 0.71` and selects k = 9,
#' far *past* the oracle -- picks some hull vertex. Substituting the average
#' fitted curvature \eqn{\bar\nu_t} for \eqn{L_t} fails differently and worse: it
#' is 10x to 100x smaller than the backtracked \eqn{L_t}, so `df` collapses to
#' 5-15% of its value and the weakest alpha on the grid wins in every cell.
#'
#' The residual is a **shape** mismatch, not a level error: the divergence
#' \eqn{d_r/(d_r + \tau)} is borrowed from soft thresholding a Gaussian matrix
#' directly, whereas this estimator is a logistic MLE reached iteratively. For
#' k = 4 to win, the penalty would have to grow by less than 1256 from k = 3 to 4
#' and more than 1522 from k = 4 to 5; it grows by 1484 and 1793. An ~18% local
#' error in one interval, costing about 4.3% RMSE.
#'
#' **`weighted = TRUE` by default.** The hard count \eqn{r_t(N+P-r_t)} treats
#' every retained singular value as a whole free parameter and over-penalizes.
#' See [svt_df()].
#' @keywords internal
cp_surrogate <- function(path, v_ref = NULL, weighted = TRUE) {
  if (is.null(v_ref)) v_ref <- path$v_t
  if (length(v_ref) == 1L) v_ref <- rep(v_ref, path$Tn)
  df <- if (weighted) path$df else path$df_hard
  # df is K x T; weight each slice by its own variance, then sum over slices.
  cp <- path$brier + 2 * as.vector(df %*% v_ref)
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
        fb <- bmfsvt(Yb, lambda = path$alpha[k] * path$lambda_star,
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
      tot <- tot + sum(contrib)
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
#' @param alpha_grid shrinkage levels relative to [lambda_star_seq()]. The
#'   default is the method note's grid, \eqn{10^{-2k/(K-1)}} with `K = 10`.
#' @param criterion `"boot"` for the full Brier plus covariance criterion,
#'   `"cp"` to stop at the Mallows Cp surrogate
#' @param B bootstrap replicates used by `"boot"`
#' @param window how many grid steps either side of the surrogate's choice the
#'   bootstrap examines
#' @param weighted_df whether the surrogate uses shrinkage weighted degrees of
#'   freedom; `FALSE` reproduces the method note's hard count
#' @param ncores bootstrap replicates to run in parallel
#' @param seed optional seed for the bootstrap draws
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
#' \deqn{M(\alpha) = \sum_{ij}\sum_t (\hat\pi^{(t)}_{ij}(\alpha) -
#'   Y^{(t)}_{ij})^2 + 2\sum_{ij}\sum_t
#'   \mathrm{Cov}(\hat\pi^{(t)}_{ij}(\alpha), Y^{(t)}_{ij}).}
#' The covariance is not available in closed form because every fitted value
#' depends on every entry, so it is estimated in two stages: [cp_surrogate()]
#' narrows the grid, then [boot_covariance()] estimates the covariance by
#' parametric bootstrap on the surviving candidates.
#'
#' **Choosing between the two: use `"cp"`.** Measured under the
#' [lambda_star_seq()] anchor and the per-slice penalty (`T = 10`,
#' `pi_tail = 0.02`), as signed grid distance from the oracle alpha (the grid
#' point minimizing RMSE against the true survival probabilities), HRS 300 x 20
#' and square 150 x 100:
#'
#' \tabular{lrrr}{
#'   criterion \tab steps HRS / square \tab RMSE HRS / square \tab extra secs \cr
#'   oracle \tab 0 / 0 \tab 0.1273 / 0.1110 \tab -- \cr
#'   `"cp"` \tab -2 / -1 \tab 0.1444 / 0.1158 \tab 0 \cr
#'   `"boot"` \tab -3 / -2 \tab 0.1585 / 0.1250 \tab 28 / 73
#' }
#'
#' `"cp"` is closer to the oracle, lower in RMSE **and** free, so `"boot"` has
#' no remaining argument in its favour on this grid. It also no longer wins on
#' rank recovery: `"cp"` reaches 0.80 / 3.14 mean absolute per-block rank error
#' against the oracle's 2.08 / 6.34. Rank recovery and estimation accuracy still
#' pull in opposite directions -- the oracle alpha is the *worst* of the three at
#' recovering rank, because soft thresholding shrinks what it keeps -- so do not
#' tune by trying to hit the true rank.
#'
#' **Both criteria over-shrink, and that is not fixable by reweighting the
#' penalty.** Sweeping a scale `c` in \eqn{\mathrm{Brier} + 2c\,\mathrm{penalty}}
#' over `[0.05, 3]`, the selected index runs 10, 9, 7, 6, 5, 3, 2, 1 on the
#' square shape and 10, 7, 2, 1 on HRS: **the oracle is never selected for any
#' `c`.** It is not a vertex of the lower convex hull of the (penalty, Brier)
#' points, so no criterion of this form can reach it -- which rules out every
#' rescale of \eqn{\tilde V_t} and of `df` at once. The residual is a shape
#' mismatch in a divergence formula borrowed from the Gaussian sequence model,
#' not a level error. See [cp_surrogate()].
#'
#' @seealso [bmfsvt()] for the fit itself, [lambda_star_seq()] for the shrinkage
#'   scale.
#' @export
tune.bmfsvt <- function(Y,
                        alpha_grid = exp(seq(log(1), log(0.2), length.out = 10L)),
                        criterion = c("cp", "boot"),
                        B = 30L, window = 2L, weighted_df = TRUE,
                        C = 20L, ncores = 1L, seed = NULL, verbose = FALSE,
                        ...) {
  criterion <- match.arg(criterion)
  Y <- as_mask_list(Y, check_nested = TRUE)

  path <- bmfsvt_path(Y, alpha_grid = alpha_grid, C = C, seed = seed, ...)
  cp <- cp_surrogate(path, weighted = weighted_df)
  index <- cp$index
  boot <- NULL
  M <- rep(NA_real_, length(path$alpha))

  if (verbose) {
    message(sprintf("Cp surrogate picks alpha=%.4g", path$alpha[cp$index]))
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
  final <- bmfsvt(Y, lambda = path$alpha[index] * path$lambda_star,
                  Z_init = path$Z[[index]], ...)

  list(fit = final,
       alpha = path$alpha[index],
       index = index,
       criterion = criterion,
       alpha_cp = path$alpha[cp$index],
       v_ref = cp$v_ref,
       path = data.frame(alpha = path$alpha, brier = path$brier,
                         v_tilde = path$v_tilde,
                         df = rowSums(path$df),
                         df_hard = rowSums(path$df_hard), cp = cp$cp, M = M,
                         mean_rank = rowMeans(path$ranks),
                         n_iter = path$n_iter, n_svd = path$n_svd,
                         converged = path$converged, elapsed = path$elapsed),
       ranks = path$ranks,
       boot = boot)
}
