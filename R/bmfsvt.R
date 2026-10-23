
#' Joint factorization of a nested sequence of binary masks
#'
#' Implements the Clip-SVT algorithm: all
#' binary masks are factorized together under one objective, optimized by FISTA
#' with backtracking and block-wise singular value thresholding.
#'
#' @param Y a list of T binary matrices, or an N x P x T array. The masks must be
#'   nested and decreasing, i.e. \eqn{Y^{(0)}_{ij} \ge Y^{(1)}_{ij} \ge \cdots},
#'   which holds automatically when they come from thresholding one abundance
#'   matrix with an increasing threshold sequence.
#' @param train optional N x P 0/1 matrix, 1 marking the entries that enter the
#'   likelihood, shared by every mask. `NULL` (the default) uses every entry.
#'   Held-out entries are excluded from the gradient,
#'   from \eqn{f} and \eqn{Q} in the backtracking test, from the
#'   offsets and from the default noise floor, but they are still clipped
#'   and receive fitted probabilities, imputed from the low-rank blocks. This
#'   is what [cv.bmfsvt()] scores.
#' @param lambda shrinkage parameter for each mask; a scalar is recycled. When
#'   `NULL` the value `alpha * lambda_star_seq(Y, C)` is used.
#' @param alpha shrinkage relative to the **noise floor** [lambda_star_seq()],
#'   used only when `lambda` is `NULL`.
#'   The search range is \eqn{\alpha \in (0, 1]}, with `alpha = 1` thresholding
#'   exactly at the floor.
#'
#' @param C pure noise replicates for estimating noise floor used by [lambda_star_seq()],
#' used only when `lambda` is `NULL`.
#' @param max_iter maximum number of outer FISTA iterations.
#' @param tol relative change in `Z` below which the algorithm stops.
#' @param gamma step size increment ratio for backtracking.
#' @param delta small constant stabilizing the relative convergence criterion.
#' @param rank_init starting guess for the rank of each block.
#' @param rank_max largest rank retained for any block; `NULL` means no limit.
#' @param rank_step how many ranks to add when the partial SVD was too small.
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
#'     \item{mu, nu}{column offsets and their increments across thresholds,
#'       a T x P and a (T-1) x P matrix}
#'     \item{pi}{prevalence of each mask over the training entries, length T}
#'     \item{V}{list of T `P x ranks[t]` matrices of right singular vectors
#'       paired with `dvals`. Like `dvals` they describe the thresholded block
#'       before clipping. Used by [predict.bmfsvt()] to fold in new rows.}
#'     \item{lambda, L}{shrinkage parameters and final Lipschitz constants}
#'     \item{ranks}{retained rank of each block}
#'     \item{dvals}{retained singular values of each block after shrinkage, as
#'       produced by the final thresholding step and therefore measured before
#'       clipping, matching the `"approx"` convention of [nuclear_penalty()].}
#'     \item{obj_trace}{objective by iteration, empty when
#'       `track_objective = "none"`}
#'     \item{n_iter, converged}{solver diagnostics}
#'   }
#'
#' @details
#' The parameterization follows
#' \eqn{X^{(t)} = X^{(0)} + \sum_{t' \le t} H^{(t')}} and \eqn{H^{(t)} \le 0},
#' writing \eqn{Z^{(0)} = X^{(0)} - \mu_0} and \eqn{Z^{(t)} = H^{(t)} - \nu_t}
#' gives \eqn{X^{(t)} = \mu_t + \sum_{t' \le t} Z^{(t')}} because the offsets
#' telescope. The gradient with respect to \eqn{Z^{(t)}} is therefore the
#' reverse cumulative sum \eqn{\sum_{t' \ge t}[\sigma(X^{(t')}) - Y^{(t')}]},
#' and the monotonicity restriction becomes \eqn{Z^{(t)} \le -\nu_t}.
#'
#'
#' The returned object has class `"bmfsvt"`.
#'
#' @seealso [cv.bmfsvt()] for choosing `alpha`,
#'   [lambda_star_seq()] for the shrinkage scale, [predict.bmfsvt()] for new
#'   rows.
#' @importFrom stats plogis qlogis
#' @export
bmfsvt <- function(Y, train = NULL, lambda = NULL, alpha = 0.55, C = 20L,
                   max_iter = 200L, tol = 1e-5, gamma = 1.1,
                   delta = 1e-6,
                   rank_init = 5L, rank_max = NULL, rank_step = 2L,
                   restart = c("gradient", "function", "none"),
                   svd_method = c("auto", "full", "svds"),
                   track_objective = c("approx", "exact", "none"),
                   Z_init = NULL, max_backtrack = 30L, verbose = FALSE) {

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
    if (any(n_col == 0)) {
      stop("`train` leaves some column with no training entry, so its ",
           "column offset is undefined.")
    }
    if (any(n_row == 0)) {
      warning("`train` leaves ", sum(n_row == 0), " row(s) with no training ",
              "entry; they are fitted from the offsets alone.")
    }
  }

  # ---- offsets and Lipschitz constants -------------------------------------
  # Taken over the training entries only, or the held-out values would leak
  # into the fit through the offsets.
  pis <- mask_prevalence(Y, train)
  eps <- 1 / (2 * nobs)
  if (any(pis <= 0 | pis >= 1)) {
    warning("some masks are entirely 0 or entirely 1; their offsets are clamped.")
  }
  pis <- pmin(pmax(pis, eps), 1 - eps)

  # calculate offsets for each column separately
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


  L <- rev(cumsum(rev(pis * (1 - pis))))
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

  # ---- initialization ------------------------------------------------------
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
  ranks <- rep(as.integer(rank_init), Tn)
  retained <- rep(NA_integer_, Tn)
  dvals <- replicate(Tn, numeric(0L), simplify = FALSE)
  Vs <- replicate(Tn, matrix(0, P, 0L), simplify = FALSE)
  s_curr <- 1

  obj_trace <- numeric(0L)
  # clip_frac_trace <- numeric(0L)
  n_svd <- 0L
  n_backtrack <- 0L
  converged <- FALSE
  m <- 0L

  clip_level <- if (Tn > 1L) -nu else numeric(0L)

  # ---- main loop -----------------------------------------------------------
  while (m < max_iter) {
    s_next <- (1 + sqrt(1 + 4 * s_curr^2)) / 2

    # Psi depends only on the extrapolation point W, so it is computed once
    # per outer iteration rather than once per backtrack.
    gc_W <- stack_grad(W, mu0, nu, Y, train, ce = TRUE)
    Psi <- gc_W$Psi
    f_W <- gc_W$f

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
        svt <- soft_svt(prox_point_cpp(W[[t]], Psi[[t]], L[t]),
                        thresh = lambda[t] / L[t],
                        rank_guess = ranks[t], rank_max = rank_max,
                        rank_step = rank_step, method = svd_method)
        n_svd <- n_svd + svt$n_svd
        new_ranks[t] <- svt$rank_next
        retained[t] <- svt$rank
        mat <- svt$mat
        dvals[[t]] <- svt$d
        Vs[[t]] <- svt$v
        if (t > 1L) {
          cl <- clip_block(mat, if (is.matrix(clip_level)) {
            clip_level[t - 1L, ]
          } else {
            clip_level[t - 1L]
          })
          mat <- cl$mat
          n_clipped <- n_clipped + cl$n_clipped
        }
        Zp[[t]] <- mat
      }

      # Q(Z' | W) = f(X_W) + <Z' - W, Psi> + sum_t L_t/2 ||Z'^(t) - W^(t)||_F^2
      ev <- step_eval(Zp, W, Psi, L, mu0, nu, Y, train, f_W)
      f_new <- ev$f
      Delta <- f_new - ev$quad

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
      # clip_frac_trace <- c(clip_frac_trace, n_clipped / (N * P * max(Tn - 1L, 1L)))
    }

    # ---- convergence -----------------------------------------
    st <- fista_stats(Zp, Z, W)
    num <- st[["num"]]
    den <- st[["den"]]
    m <- m + 1L
    if (verbose) {
      message(sprintf(
        "iter %3d  f=%.6g  rel.change=%.3e  backtracks=%d  ranks=%s",
        m, f_new, num / (den + delta), bt, paste(ranks, collapse = ",")))
    }
    if (num / (den + delta) < tol) {
      Z <- Zp
      converged <- TRUE
      break
    }

    # ---- FISTA extrapolation with optional restart ------------------------
    do_restart <- switch(
      restart,
      none = FALSE,
      gradient = {
        # O'Donoghue and Candes generalized gradient scheme: restart when the
        # step moves against the proximal direction. Cheap insurance, and the
        # clipping step breaks the exact proximal-gradient structure that plain
        # FISTA assumes.
        st[["ip"]] > 0
      },
      `function` = length(obj_trace) > 1L &&
        obj_trace[length(obj_trace)] > obj_trace[length(obj_trace) - 1L]
    )

    if (do_restart) {
      s_next <- 1
      W <- Zp
    } else {
      fac <- (s_curr - 1) / s_next
      W <- if (fac == 0) Zp else extrapolate_cpp(Zp, Z, fac)
    }
    s_curr <- s_next
    Z <- Zp
  }

  X <- stack_forward(Z, mu0, nu)



  structure(
    list(X = X,
         prob = lapply(X, sigmoid),
         Z = Z,
         mu = mu, nu = nu, pi = pis,
         lambda = lambda, L = L,
         ranks = retained,
         dvals = dvals,
         V = Vs,
         obj_trace = obj_trace,
         n_iter = m,
         converged = converged),
    class = "bmfsvt")
}




#' Denoise new rows with a fitted factorization
#'
#' @param object a fit from [bmfsvt()], e.g. `cv.bmfsvt(...)$fit`
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
#'
#' **Clipping is applied after solving, not during.** The increments are capped
#' at \eqn{-\nu_t} once the ridge problem is solved, so the fitted
#' probabilities are non-increasing across thresholds. Clipping inside the
#' iterations would make each row's problem non-convex. The in-sample fit clips
#' inside its iterations, so on training rows the two differ slightly; [dfbm()]
#' reports the difference as `insample_gap`.
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
  Z <- clip_increments(sol$Z, object$nu)
  X <- stack_forward(Z, mu0, object$nu)

  ce <- Reduce(`+`, lapply(seq_len(Tn), function(t) {
    rowSums(log1exp(X[[t]]) - Y[[t]] * X[[t]])
  }))

  list(prob = lapply(X, sigmoid), X = X, Z = Z, A = sol$A,
       row_ce = ce / (Tn * P),
       grad_norm = sol$grad_norm, n_iter = sol$n_iter,
       converged = sol$converged)
}



