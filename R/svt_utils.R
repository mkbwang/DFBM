
# Helper kernels for the Clip-SVT binary matrix factorization (Algorithm 2).
#
# The three functions marked "hot kernel" below are the elementwise passes that
# dominate runtime whenever the SVD is cheap (e.g. the HRS shape, 9000 x 20).
# They are kept small and free of side effects so that they can be swapped for
# RcppArmadillo equivalents without touching bmfsvt() itself.


#' log one plus exponential
#'
#' @param x a vector
#' @returns log(1+exp(x))
#' @details
#' This function is used to avoid numerical issues when x is too large
#' @export
log1exp <- function(x){
  output <- x
  output[x<20] <- log(1+exp(x[x<20]))
  return(output)
}


#' Sigmoid with guard against overflow
#'
#' @param x a numeric vector or matrix
#' @returns 1/(1+exp(-x)), evaluated by [stats::plogis()]
#' @keywords internal
sigmoid <- function(x) {
  out <- x
  out[] <- stats::plogis(x)
  out
}


#' Coerce a stack of binary masks into a list of matrices
#'
#' @param Y either a list of T binary matrices or an N x P x T array
#' @param check_nested whether to verify that the masks are nested and decreasing
#' @returns a list of T numeric matrices
#' @details
#' The model in section 2.2 assumes the masks come from thresholding one
#' abundance matrix, hence \eqn{Y^{(0)} \ge Y^{(1)} \ge \cdots}. A stack that
#' violates nesting is not merely unusual, it invalidates the parameterization
#' \eqn{H^{(t)} \le 0}, so we stop rather than warn.
#' @keywords internal
as_mask_list <- function(Y, check_nested = TRUE) {
  if (is.array(Y) && length(dim(Y)) == 3L) {
    # Rebuild the matrix explicitly: `drop = TRUE` would turn a 1 x P x T array
    # (a single new sample) into plain vectors.
    Y <- lapply(seq_len(dim(Y)[3L]),
                function(t) matrix(Y[, , t], dim(Y)[1L], dim(Y)[2L]))
  }
  if (!is.list(Y)) {
    stop("`Y` must be a list of binary matrices or an N x P x T array.")
  }
  if (length(Y) == 0L) stop("`Y` is empty.")
  Y <- lapply(Y, function(mat) {
    if (!is.matrix(mat)) stop("every element of `Y` must be a matrix.")
    storage.mode(mat) <- "double"
    mat
  })
  dims <- vapply(Y, dim, integer(2L))
  if (any(dims[1L, ] != dims[1L, 1L]) || any(dims[2L, ] != dims[2L, 1L])) {
    stop("all binary masks must share the same dimensions.")
  }
  bad <- vapply(Y, function(mat) any(mat != 0 & mat != 1), logical(1L))
  if (any(bad)) stop("all binary masks must contain only 0 and 1.")

  if (check_nested && length(Y) > 1L) {
    for (t in seq.int(2L, length(Y))) {
      if (any(Y[[t]] > Y[[t - 1L]])) {
        stop(sprintf(
          paste0("mask %d is not nested inside mask %d: %d entries have ",
                 "Y[[%d]] > Y[[%d]]. The masks must come from thresholding a ",
                 "single matrix with an increasing threshold sequence."),
          t, t - 1L, sum(Y[[t]] > Y[[t - 1L]]), t, t - 1L))
      }
    }
  }
  Y
}


#' Validate a training mask
#'
#' @param train `NULL` or an N x P matrix of 0/1 (or logical) values, 1 marking
#'   the entries that enter the likelihood
#' @param N,P dimensions of the binary masks
#' @returns `NULL`, or `train` as a double matrix. A mask that selects every
#'   entry is returned as `NULL`, so the unmasked code path (and its cost) is
#'   used whenever nothing is held out.
#' @keywords internal
check_train <- function(train, N, P) {
  if (is.null(train)) return(NULL)
  if (!is.matrix(train) || !identical(dim(train), c(as.integer(N), as.integer(P)))) {
    stop("`train` must be an N x P matrix matching the binary masks.")
  }
  storage.mode(train) <- "double"
  if (anyNA(train) || any(train != 0 & train != 1)) {
    stop("`train` must contain only 0 and 1.")
  }
  if (all(train == 1)) return(NULL)
  if (!any(train == 1)) stop("`train` marks no training entries.")
  train
}


#' Prevalence of each mask over the training entries
#'
#' @param Y list of T binary matrices
#' @param train `NULL` or a validated training mask, see [check_train()]
#' @returns numeric vector of length T
#' @keywords internal
mask_prevalence <- function(Y, train = NULL) {
  if (is.null(train)) return(vapply(Y, mean, numeric(1L)))
  nobs <- sum(train)
  vapply(Y, function(mat) sum(train * mat) / nobs, numeric(1L))
}


#' Penalized cross entropy of one binary mask
#'
#' @param X real valued matrix of natural parameters
#' @param Y binary matrix
#' @param train optional 0/1 matrix of the same size marking the entries that
#'   enter the likelihood; `NULL` uses every entry
#' @returns the summed negative log likelihood over the training entries
#' @details
#' Hot kernel. Written as \eqn{-YX + \log(1+e^X)} so that [log1exp()] handles the
#' large-\eqn{X} branch; evaluating \code{log(sigmoid(X))} directly underflows.
#' @keywords internal
logistic_ce <- function(X, Y, train = NULL) {
  val <- log1exp(X) - Y * X
  if (is.null(train)) sum(val) else sum(train * val)
}


#' Per-entry Brier loss of a fitted mask stack, summed over thresholds
#'
#' @param X list of T natural parameter matrices
#' @param Y list of T binary matrices
#' @returns a numeric vector, one element per entry
#' @details
#' Hot kernel. One pass per mask with no allocation beyond the accumulator.
#'
#' The sum over \eqn{t} is taken *inside*, so the unit of observation of the
#' returned vector is the entry \eqn{(i,j)} rather than the pair
#' \eqn{(\text{entry}, t)}. That is what makes a standard error over the result
#' honest: the masks are nested and generated from a single latent value per
#' entry, so the T terms belonging to one entry are strongly dependent and
#' treating them as T independent observations understates the spread by roughly
#' \eqn{\sqrt{T}}.
#' @keywords internal
entry_loss <- function(X, Y) {
  out <- numeric(length(X[[1L]]))
  for (t in seq_along(X)) {
    out <- out + (stats::plogis(as.vector(X[[t]])) - as.vector(Y[[t]]))^2
  }
  out
}


#' Degrees of freedom of a soft thresholded singular value decomposition
#'
#' @param dvals list of T numeric vectors holding the retained singular values
#'   of each block *after* shrinkage, as returned by [soft_svt()]
#' @param thresh numeric vector of T thresholds **as actually applied by the
#'   proximal step**, that is \eqn{\lambda_t / L_t} and not \eqn{\lambda_t}
#' @param N,P matrix dimensions
#' @param weighted whether to weight each retained component by how little it
#'   was shrunk. `FALSE` reproduces the hard count \eqn{r_t(N+P-r_t)} of the
#'   method note and ignores `thresh`.
#' @returns a numeric vector of length T
#' @details
#' \eqn{\mathrm{df}_t = \sum_i s_i (N + P - 2i + 1)}, where
#' \eqn{(N+P-2i+1)} is the free parameter count of the i-th singular triplet, so
#' that \eqn{\sum_{i \le r}(N+P-2i+1) = r(N+P-r)} recovers the hard count
#' exactly when nothing is discounted.
#'
#' The hard count treats every retained singular value as a whole free
#' parameter. Soft thresholding shrinks what it retains, so the effective number
#' is smaller by \eqn{s_i = d_i/(d_i + \tau_t)}, with \eqn{d_i} the *shrunken*
#' value stored in `dvals` and \eqn{d_i + \tau_t} the value before thresholding.
#'
#' **`dvals` are the singular values before clipping, which is what this wants.**
#' Clipping is an elementwise \eqn{\min(\cdot, c)}, so its derivative with
#' respect to the data is 1 on untouched entries and 0 on capped ones: it can
#' only lower the sensitivity of the fit to `Y`, hence lower the degrees of
#' freedom. It also destroys low-rankness -- capping a rank 13 block was
#' measured to leave a matrix of rank 43 -- so feeding the clipped spectrum in
#' would count artifacts of the cap as free parameters and push `df` the wrong
#' way (3225 to 3340 on one fit). The exact adjustment is to scale by the
#' unclipped fraction, which was 0.9965 there and is not worth taking.
#'
#' **`thresh` must be \eqn{\lambda_t/L_t}.** That is what
#' `soft_svt(W - \Psi/L, thresh = lambda/L)` actually subtracts, so it is what
#' has to be added back to recover the original singular value. Passing
#' \eqn{\lambda_t} instead understates `df` -- the two differ by the Lipschitz
#' constant, measured around 2.7x on one block -- which weakens the Cp penalty
#' and biases selection toward too little shrinkage. Over 8 cells the correct
#' threshold gave a mean distance of 0.75 grid steps from the best available
#' alpha against 0.875, and RMSE 0.1505 against 0.1549. Both are still biased,
#' in opposite directions (-0.50 and +0.625 signed), which says the divergence
#' of an iterative estimator is not fully captured by the proximal step alone.
#' @keywords internal
svt_df <- function(dvals, thresh, N, P, weighted = TRUE) {
  vapply(seq_along(dvals), function(t) {
    d <- dvals[[t]]
    r <- length(d)
    if (r == 0L) return(0)
    if (!weighted) return(r * (N + P - r))
    s <- d / (d + thresh[t])
    sum(s * (N + P - 2 * seq_len(r) + 1))
  }, numeric(1L))
}


#' Mean Bernoulli variance of a fitted probability stack
#'
#' @param prob list of T probability matrices
#' @returns a single number, \eqn{\frac{1}{NPT}\sum \hat\pi(1-\hat\pi)}
#' @details
#' Named rather than inlined because the Mallows Cp surrogate is only correct
#' when this is evaluated once at a *fixed* reference fit and held constant
#' across the shrinkage grid. Re-estimating it from each candidate fit makes the
#' penalty vanish exactly where it is needed: as the shrinkage weakens the fit
#' drives the probabilities toward 0 and 1, so this quantity collapses (0.209 to
#' 0.019 across one K = 10 grid) while the degrees of freedom explode.
#' @keywords internal
mean_bernoulli_var <- function(prob) {
  mean(vapply(prob, function(m) mean(m * (1 - m)), numeric(1L)))
}


#' Accumulate the natural parameter stack from the increments
#'
#' @param Z list of T increment matrices
#' @param mu0 offset of the first mask: a scalar, or a vector of length P
#'   holding one offset per column
#' @param nu offset increments of masks 2..T: a numeric vector of length T-1, or
#'   a (T-1) x P matrix holding one increment per column
#' @returns list of T matrices \eqn{X^{(t)} = \mu_t + \sum_{t' \le t} Z^{(t')}}
#' @details
#' Hot kernel. Uses the telescoping identity \eqn{\mu_0 + \sum_{t'=1}^{t}\nu_{t'}
#' = \mu_t} so only a running sum is needed.
#'
#' Column offsets are broadcast with `rep(x, each = n)`. Plain recycling of a
#' length-P vector against an n x P matrix runs down the rows, which is wrong.
#' @keywords internal
stack_forward <- function(Z, mu0, nu) {
  Tn <- length(Z)
  X <- vector("list", Tn)
  if (is.matrix(nu) || length(mu0) > 1L) {
    n <- nrow(Z[[1L]])
    running <- Z[[1L]] + rep(mu0, each = n)
    X[[1L]] <- running
    if (Tn > 1L) {
      for (t in seq.int(2L, Tn)) {
        running <- running + Z[[t]] + rep(nu[t - 1L, ], each = n)
        X[[t]] <- running
      }
    }
    return(X)
  }
  running <- Z[[1L]] + mu0
  X[[1L]] <- running
  if (Tn > 1L) {
    for (t in seq.int(2L, Tn)) {
      running <- running + Z[[t]] + nu[t - 1L]
      X[[t]] <- running
    }
  }
  X
}


#' Enforce the monotonicity restriction on a stack of increments
#'
#' @param Z list of T increment matrices
#' @param nu offset increments, a vector of length T-1 or a (T-1) x P matrix, as
#'   in [stack_forward()]
#' @returns `Z` with \eqn{Z^{(t)} \leftarrow \min(Z^{(t)}, -\nu_t)} for t >= 2
#' @details
#' Pure elementwise pass. \eqn{Z^{(t)} \le -\nu_t} is exactly
#' \eqn{X^{(t)} \le X^{(t-1)}}, so the returned stack gives non-increasing
#' probabilities across thresholds.
#' @keywords internal
clip_increments <- function(Z, nu) {
  Tn <- length(Z)
  if (Tn < 2L) return(Z)
  n <- nrow(Z[[1L]])
  for (t in seq.int(2L, Tn)) {
    level <- if (is.matrix(nu)) rep(-nu[t - 1L, ], each = n) else -nu[t - 1L]
    Z[[t]] <- pmin(Z[[t]], level)
  }
  Z
}


#' Expected value implied by a stack of survival probabilities
#'
#' @param prob list of T n x P matrices, \eqn{S^{(t)}_{ij} = P(A_{ij} > d_{tj})}
#' @param M a (T+1) x P matrix of interval representatives: row 1 for
#'   \eqn{(-\infty, d_1]}, row t+1 for \eqn{(d_t, d_{t+1}]}, row T+1 for
#'   \eqn{(d_T, \infty)}
#' @returns an n x P matrix,
#'   \eqn{\hat A_{ij} = \sum_{t=0}^{T} (S_t - S_{t+1}) m_{tj}}
#' @details
#' Hot kernel. Evaluated in the Abel summation form
#' \eqn{\hat A_{ij} = m_{0j} + \sum_{t=1}^{T} S^{(t)}_{ij}(m_{tj} - m_{t-1,j})},
#' the discrete analogue of \eqn{E X = \int S(x)\,dx}, so only one pass per
#' mask is needed. Because the representatives are non-decreasing in t, the
#' result lies in \eqn{[m_{0j}, m_{Tj}]} for any probabilities in \eqn{[0, 1]},
#' monotone or not.
#'
#' Passing the binary masks themselves as `prob` returns each entry's own
#' interval representative, which is the binned control.
#' @keywords internal
survival_expectation <- function(prob, M) {
  Tn <- length(prob)
  n <- nrow(prob[[1L]])
  P <- ncol(prob[[1L]])
  if (!is.matrix(M) || nrow(M) != Tn + 1L || ncol(M) != P) {
    stop("`M` must be a (T+1) x P matrix matching `prob`.")
  }
  out <- matrix(rep(M[1L, ], each = n), n, P)
  for (t in seq_len(Tn)) {
    out <- out + prob[[t]] * rep(M[t + 1L, ] - M[t, ], each = n)
  }
  out
}


#' CRPS of the discrete distribution implied by a survival stack
#'
#' @param prob list of T probability vectors or matrices of a common length n,
#'   \eqn{S^{(t)} = P(A > d_t)}
#' @param M (T+1) x n matrix of atoms, column i holding \eqn{m_0 \le \cdots \le
#'   m_T} for entry i (for a column-wise representative matrix `Mcol` and entry
#'   columns `j`, pass `Mcol[, j]`)
#' @param y length n vector of observed values
#' @returns length n vector of CRPS values
#' @details
#' Hot kernel. The predictive distribution puts mass \eqn{S_t - S_{t+1}} (with
#' \eqn{S_0 = 1}, \eqn{S_{T+1} = 0}) on atom \eqn{m_t}, the distribution whose
#' mean is the denoised value of [survival_expectation()]. Its CDF equals
#' \eqn{1 - S_{k+1}} on \eqn{[m_k, m_{k+1})}, so
#' \deqn{\mathrm{CRPS} = \int (F(x) - 1\{x \ge y\})^2 dx =
#'   \sum_k \big[(1-S_{k+1})^2 b_k + S_{k+1}^2 (g_k - b_k)\big] +
#'   (m_0 - y)_+ + (y - m_T)_+,}
#' with gap \eqn{g_k = m_{k+1} - m_k} and \eqn{b_k = \min(\max(y - m_k, 0), g_k)}
#' the part of the gap lying below `y`. Exact and vectorized, in place of a
#' numerical integral per entry.
#'
#' Probabilities are not required to be monotone in t; the formula integrates
#' whatever step function they define.
#' @keywords internal
crps_discrete <- function(prob, M, y) {
  Tn <- length(prob)
  if (!is.matrix(M) || nrow(M) != Tn + 1L || ncol(M) != length(y)) {
    stop("`M` must be a (T+1) x n matrix matching `prob` and `y`.")
  }
  out <- pmax(M[1L, ] - y, 0) + pmax(y - M[Tn + 1L, ], 0)
  for (k in seq_len(Tn)) {
    s <- as.vector(prob[[k]])
    g <- M[k + 1L, ] - M[k, ]
    b <- pmin(pmax(y - M[k, ], 0), g)
    out <- out + (1 - s)^2 * b + s^2 * (g - b)
  }
  out
}


#' Reverse cumulative gradient of the stacked cross entropy
#'
#' @param X list of T natural parameter matrices
#' @param Y list of T binary matrices
#' @param train optional 0/1 matrix marking the training entries; `NULL` uses
#'   every entry
#' @returns list of T matrices \eqn{\Psi_t = \sum_{t' \ge t} (\sigma(X^{(t')}) - Y^{(t')})}
#' @details
#' Hot kernel. Because \eqn{X^{(t')}} depends on \eqn{Z^{(t)}} for every
#' \eqn{t' \ge t}, the gradient with respect to \eqn{Z^{(t)}} is a reverse
#' cumulative sum, computed here in a single backward pass.
#'
#' With `train`, held-out entries contribute a zero residual (line 21 of
#' Algorithm 2 restricted to the training entries). The proximal step then
#' leaves them at the current low-rank fit, so they are imputed as in softImpute.
#' @keywords internal
grad_backward <- function(X, Y, train = NULL) {
  Tn <- length(X)
  Psi <- vector("list", Tn)
  running <- NULL
  for (t in seq.int(Tn, 1L)) {
    resid <- sigmoid(X[[t]]) - Y[[t]]
    if (!is.null(train)) resid <- resid * train
    running <- if (is.null(running)) resid else running + resid
    Psi[[t]] <- running
  }
  Psi
}


#' Adaptive rank soft thresholded singular value decomposition
#'
#' @param M matrix to threshold
#' @param thresh soft thresholding level, i.e. lambda/L
#' @param rank_guess starting guess for the number of singular triplets
#' @param rank_max largest rank allowed
#' @param rank_step how many ranks to add when the guess was too small
#' @param method one of "auto", "full", "svds"
#' @param full_cutoff always use a full SVD when min(nrow, ncol) is at most this
#' @param full_frac use a full SVD when the requested rank is at least this
#'   fraction of min(nrow, ncol). Measured crossover is near 0.25 at 800 x 200
#'   and near 0.14 at 1000 x 1000, and being slightly conservative costs little
#'   because the two backends are within 15% of each other there.
#' @param max_partial give up searching for the rank after this many partial
#'   decompositions and compute the exact one instead
#' @returns a list with the thresholded matrix `mat`, the retained rank `rank`,
#'   the retained singular values `d`, their right singular vectors `v` (an
#'   `ncol(M) x rank` matrix, zero columns when nothing is retained), the rank
#'   guess to reuse next time `rank_next`, and the number of decompositions
#'   performed `n_svd`
#' @details
#' Follows the softImpute prescription referenced at the end of section 2.2:
#' compute a partial SVD of size `rank_guess`, and if its smallest singular
#' value still exceeds `thresh` the truncation discarded signal, so widen the
#' guess and refit. The returned `rank_next` warm starts the next outer
#' iteration.
#'
#' The `"auto"` rule picks a backend from the ratio of the requested rank to the
#' matrix dimension, not from the dimension alone. A partial SVD only pays off
#' while `k` is a small fraction of `min(nrow, ncol)`; once it approaches that
#' bound ARPACK is both slower and less reliable than a direct decomposition.
#' Choosing on dimension alone is a trap: a 800 x 200 block needing rank 10
#' would take a full 200 component decomposition and throw 190 of them away.
#' @importFrom RSpectra svds
#' @keywords internal
soft_svt <- function(M, thresh, rank_guess = 5L, rank_max = NULL,
                     rank_step = 2L, method = c("auto", "full", "svds"),
                     full_cutoff = 50L, full_frac = 0.25, max_partial = 3L) {
  method <- match.arg(method)
  nr <- nrow(M)
  nc <- ncol(M)
  mindim <- min(nr, nc)
  # `rank_max` caps how many components may be retained. It must not drive the
  # backend choice: its natural default is `mindim`, and testing `rank_max >=
  # mindim` would then send every matrix to the full decomposition.
  rank_cap <- if (is.null(rank_max)) mindim else
    max(1L, min(as.integer(rank_max), mindim))

  k <- max(1L, min(as.integer(rank_guess), rank_cap))
  # ARPACK needs k strictly below the smaller dimension, and loses its edge as k
  # approaches it. Thresholds set from end-to-end measurement, which disagreed
  # with the isolated kernel benchmark and is the one to trust: in isolation the
  # partial solver looked 2x to 16x faster everywhere, but inside a full fit the
  # gain is only about 1.2x at 800x200 and it is a net LOSS below roughly 50
  # columns, where a full LAPACK decomposition of a narrow matrix costs less
  # than ARPACK's per-call overhead. Hence a cutoff rather than "always
  # partial". The large wins are expected at proteomics scale (1000x1000 and
  # up); confirm there before tuning further.
  # Benchmarking on a flat-spectrum random matrix reverses the ordering
  # entirely, so it is the wrong thing to measure on.
  partial_worthwhile <- function(k) {
    k < mindim && (method == "svds" ||
                     (mindim > full_cutoff && k < full_frac * mindim))
  }
  use_full <- method == "full" || !partial_worthwhile(k)

  full_svd <- function() {
    sv <- svd(M)
    n_svd <<- n_svd + 1L
    sv
  }

  n_svd <- 0L
  if (use_full) {
    sv <- full_svd()
    d <- sv$d
  } else {
    attempt <- 0L
    repeat {
      # Fall back rather than fail: ARPACK can stall on a clustered spectrum,
      # and a slow fit beats an aborted one.
      sv <- tryCatch(RSpectra::svds(M, k = k), error = function(e) NULL)
      if (is.null(sv)) {
        sv <- full_svd()
        d <- sv$d
        break
      }
      n_svd <- n_svd + 1L
      attempt <- attempt + 1L
      d <- sv$d
      # If even the smallest computed singular value survives thresholding we
      # truncated too aggressively and are missing retained components.
      if (min(d) <= thresh || k >= rank_cap) break
      # Grow geometrically. Adding `rank_step` at a time needs ~50 partial
      # decompositions to walk from a guess of 5 up to a retained rank near 100,
      # which is far more expensive than the single full decomposition it was
      # meant to avoid.
      k <- min(max(k + rank_step, ceiling(k * 1.5)), rank_cap)
      # Give up on guessing after a few misses: repeatedly re-decomposing to
      # search for the rank costs more than computing it exactly once.
      if (attempt >= max_partial || !partial_worthwhile(k)) {
        sv <- full_svd()
        d <- sv$d
        break
      }
    }
  }
  keep <- which(d > thresh)
  # Warm start the next call one step above what we actually kept, so the rank
  # can grow and shrink with the iterates.
  rank_next <- max(1L, min(length(keep) + rank_step, rank_cap))

  if (length(keep) == 0L) {
    return(list(mat = matrix(0, nr, nc), rank = 0L, d = numeric(0),
                v = matrix(0, nc, 0L),
                rank_next = max(1L, min(rank_step, rank_cap)), n_svd = n_svd))
  }

  dshrunk <- d[keep] - thresh
  u <- sv$u[, keep, drop = FALSE]
  v <- sv$v[, keep, drop = FALSE]
  mat <- u %*% (dshrunk * t(v))

  # `v` is kept so that a fitted block can fold in new rows; see predict.bmfsvt().
  list(mat = mat, rank = length(keep), d = dshrunk, v = v,
       rank_next = rank_next, n_svd = n_svd)
}


#' Largest useful shrinkage parameter for each binary mask
#'
#' @param Y list of T binary matrices, or an N x P x T array
#' @returns a numeric vector of length T
#' @details
#' \eqn{\lambda_t^{max} = \| \sum_{t' \ge t} (\pi_{t'} - Y^{(t')}) \|_2}
#' is the spectral norm of the gradient at \eqn{Z = 0}, hence the smallest
#' shrinkage for which the estimate collapses to \eqn{\hat Z^{(t)} = 0}.
#'
#' **This no longer sets the search range.** It bounds the grid from above but
#' the whole interval \eqn{[\lambda^*_t, \lambda^{max}_t]} lies on the
#' over-shrinking side, so [lambda_star_seq()] is the anchor; see there. This is
#' kept because [violation_diagnostics()] and the reports still use it.
#' @seealso [lambda_star_seq()], which does set the range
#' @importFrom RSpectra svds
#' @export
lambda_max_seq <- function(Y) {
  Y <- as_mask_list(Y, check_nested = FALSE)
  Tn <- length(Y)
  pis <- vapply(Y, mean, numeric(1L))

  out <- numeric(Tn)
  running <- NULL
  for (t in seq.int(Tn, 1L)) {
    resid <- pis[t] - Y[[t]]
    running <- if (is.null(running)) resid else running + resid
    out[t] <- spectral_norm(running)
  }
  out
}


#' Noise floor shrinkage parameter for each binary mask
#'
#' @param Y list of T binary matrices, or an N x P x T array
#' @param C number of pure noise replicates to average
#' @param seed optional integer; makes the Monte Carlo draw reproducible
#' @param train optional 0/1 N x P matrix of training entries. The prevalences
#'   are then taken over the training entries and the noise residual is zeroed
#'   on the held-out ones, which is the noise floor of the gradient that
#'   [bmfsvt()] actually sees when fitted with the same `train`.
#' @returns a numeric vector of length T
#' @details
#' \eqn{\lambda^*_t = \frac{1}{C}\sum_c \| \sum_{t' \ge t}
#' (\pi_{t'} - Y^{(t')c*}) \|_2}, the spectral norm of the same gradient
#' evaluated on a stack carrying *no signal*. Section 2.2 of the method note.
#' The search range is \eqn{\lambda_t = \alpha \lambda^*_t} with
#' \eqn{\alpha \in (0, 1]}.
#'
#' **The noise stack must be nested, and one shared uniform per entry is what
#' makes it so.** Drawing \eqn{U_{ij} \sim U(0,1)} once and setting
#' \eqn{Y^{(t')c*}_{ij} = I[U_{ij} < \pi_{t'}]} for every \eqn{t'} gives
#' \eqn{Y^{(0)c*} \ge Y^{(1)c*} \ge \cdots} automatically. Drawing each mask
#' independently would break the nesting and understate the floor: the shared
#' uniform makes the residuals positively correlated across \eqn{t'}, so their
#' variances do not simply add. Measured, an independent-entry analytic floor
#' \eqn{\sqrt{\sum_{t' \ge t}\pi(1-\pi)}(\sqrt N + \sqrt P)} is 1.1x to 2.1x too
#' small, the gap widening toward the dense head where more terms are summed.
#'
#' `C = 20` is ample and does not need tuning. The spectral norm of a random
#' matrix has Tracy-Widom fluctuations of order \eqn{N^{-1/6}} (about 2% at
#' \eqn{N = 150}), so averaging 20 draws puts the Monte Carlo error near 0.5%,
#' well inside one step of the recommended grid (20%).
#'
#' **With `train`, the floor is that of the masked problem.** Removing a
#' fraction of the entries shrinks the noise gradient by roughly the square root
#' of the training fraction. Anchoring a cross-validation fit on the full-data
#' floor would shrink it about `1/prop` times harder than the refit at the same
#' `alpha`; anchoring on its own floor keeps `alpha` meaning the same fraction
#' of the problem's own noise level. See [cv.bmfsvt()].
#'
#' Cost is `C` partial decompositions per threshold, paid **once per data set**:
#' the whole \eqn{\alpha} grid rescales the same vector, so this must not be
#' called per candidate.
#' @seealso [lambda_max_seq()] for the upper bound
#' @importFrom stats runif
#' @export
lambda_star_seq <- function(Y, C = 20L, seed = NULL, train = NULL) {
  Y <- as_mask_list(Y, check_nested = FALSE)
  Tn <- length(Y)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  train <- check_train(train, N, P)
  pis <- mask_prevalence(Y, train)

  # This is called from inside bmfsvt() by default, so a `seed` must not leak
  # into the caller's RNG stream: isolate it and put the stream back on exit.
  # With `seed = NULL` the draws advance the stream normally, which is what a
  # caller who did not ask for reproducibility expects.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }

  acc <- numeric(Tn)
  for (cc in seq_len(C)) {
    U <- matrix(stats::runif(N * P), N, P)
    running <- NULL
    for (t in seq.int(Tn, 1L)) {
      # I[U < pi_t] is the pure noise mask; the same U drives every threshold.
      resid <- pis[t] - (U < pis[t])
      if (!is.null(train)) resid <- resid * train
      running <- if (is.null(running)) resid else running + resid
      acc[t] <- acc[t] + spectral_norm(running)
    }
  }
  acc / C
}


#' Largest singular value of a matrix
#'
#' @param M a matrix
#' @returns the spectral norm
#' @importFrom RSpectra svds
#' @keywords internal
spectral_norm <- function(M) {
  if (all(M == 0)) return(0)
  mindim <- min(dim(M))
  # Only the leading singular value is needed, so the partial solver wins as
  # soon as the matrix is bigger than trivial.
  if (mindim <= 50L) {
    max(svd(M, nu = 0L, nv = 0L)$d)
  } else {
    RSpectra::svds(M, k = 1L, nu = 0L, nv = 0L)$d[1L]
  }
}
