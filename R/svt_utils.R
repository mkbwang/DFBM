
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
    Y <- lapply(seq_len(dim(Y)[3L]), function(t) Y[, , t, drop = TRUE])
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


#' Penalized cross entropy of one binary mask
#'
#' @param X real valued matrix of natural parameters
#' @param Y binary matrix
#' @param Omega binary matrix of observed entries, or NULL when all are observed
#' @returns the summed negative log likelihood over the observed entries
#' @details
#' Hot kernel. Written as \eqn{-YX + \log(1+e^X)} so that [log1exp()] handles the
#' large-\eqn{X} branch; evaluating \code{log(sigmoid(X))} directly underflows.
#' @keywords internal
logistic_ce <- function(X, Y, Omega = NULL) {
  val <- log1exp(X) - Y * X
  if (is.null(Omega)) sum(val) else sum(Omega * val)
}


#' Per-entry loss of a fitted mask stack, summed over thresholds
#'
#' @param X list of T natural parameter matrices
#' @param Y list of T binary matrices
#' @param idx integer positions of the entries to score, in column major order,
#'   or NULL to score every entry
#' @param loss `"brier"` for squared error on the probability scale, `"ce"` for
#'   cross entropy
#' @returns a numeric vector, one element per scored entry
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
entry_loss <- function(X, Y, idx = NULL, loss = c("brier", "ce")) {
  loss <- match.arg(loss)
  n <- if (is.null(idx)) length(X[[1L]]) else length(idx)
  out <- numeric(n)
  for (t in seq_along(X)) {
    x <- if (is.null(idx)) as.vector(X[[t]]) else X[[t]][idx]
    y <- if (is.null(idx)) as.vector(Y[[t]]) else Y[[t]][idx]
    out <- out + if (loss == "brier") {
      (stats::plogis(x) - y)^2
    } else {
      log1exp(x) - y * x
    }
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
#' @returns a single number, \eqn{\frac{1}{|\Omega|T}\sum \hat\pi(1-\hat\pi)}
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
#' @param mu0 scalar offset of the first mask
#' @param nu numeric vector of length T-1 holding the offsets of masks 1..T-1
#' @returns list of T matrices \eqn{X^{(t)} = \mu_t + \sum_{t' \le t} Z^{(t')}}
#' @details
#' Hot kernel. Uses the telescoping identity \eqn{\mu_0 + \sum_{t'=1}^{t}\nu_{t'}
#' = \mu_t} so only a running sum is needed.
#' @keywords internal
stack_forward <- function(Z, mu0, nu) {
  Tn <- length(Z)
  X <- vector("list", Tn)
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


#' Reverse cumulative gradient of the stacked cross entropy
#'
#' @param X list of T natural parameter matrices
#' @param Y list of T binary matrices
#' @param Omega binary matrix of observed entries, or NULL when all are observed
#' @returns list of T matrices \eqn{\Psi_t = \sum_{t' \ge t} P_\Omega(\sigma(X^{(t')}) - Y^{(t')})}
#' @details
#' Hot kernel. Because \eqn{X^{(t')}} depends on \eqn{Z^{(t)}} for every
#' \eqn{t' \ge t}, the gradient with respect to \eqn{Z^{(t)}} is a reverse
#' cumulative sum, computed here in a single backward pass.
#' @keywords internal
grad_backward <- function(X, Y, Omega = NULL) {
  Tn <- length(X)
  Psi <- vector("list", Tn)
  running <- NULL
  for (t in seq.int(Tn, 1L)) {
    resid <- sigmoid(X[[t]]) - Y[[t]]
    if (!is.null(Omega)) resid <- resid * Omega
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
#'   the retained singular values `d`, the rank guess to reuse next time
#'   `rank_next`, and the number of decompositions performed `n_svd`
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
                rank_next = max(1L, min(rank_step, rank_cap)), n_svd = n_svd))
  }

  dshrunk <- d[keep] - thresh
  u <- sv$u[, keep, drop = FALSE]
  v <- sv$v[, keep, drop = FALSE]
  mat <- u %*% (dshrunk * t(v))

  list(mat = mat, rank = length(keep), d = dshrunk,
       rank_next = rank_next, n_svd = n_svd)
}


#' Largest useful shrinkage parameter for each binary mask
#'
#' @param Y list of T binary matrices, or an N x P x T array
#' @param Omega binary matrix of observed entries, or NULL when all are observed
#' @returns a numeric vector of length T
#' @details
#' \eqn{\lambda_t^{max} = \| \sum_{t' \ge t} P_\Omega(\pi_{t'} - Y^{(t')}) \|_2}
#' is the spectral norm of the gradient at \eqn{Z = 0}, hence the smallest
#' shrinkage for which the estimate collapses to \eqn{\hat Z^{(t)} = 0}.
#' @importFrom RSpectra svds
#' @export
lambda_max_seq <- function(Y, Omega = NULL) {
  Y <- as_mask_list(Y, check_nested = FALSE)
  Tn <- length(Y)
  nobs <- if (is.null(Omega)) length(Y[[1L]]) else sum(Omega)
  pis <- vapply(Y, function(mat) {
    if (is.null(Omega)) mean(mat) else sum(Omega * mat) / nobs
  }, numeric(1L))

  out <- numeric(Tn)
  running <- NULL
  for (t in seq.int(Tn, 1L)) {
    resid <- pis[t] - Y[[t]]
    if (!is.null(Omega)) resid <- resid * Omega
    running <- if (is.null(running)) resid else running + resid
    out[t] <- spectral_norm(running)
  }
  out
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
