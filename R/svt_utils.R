
# Helper kernels for the Clip-SVT binary matrix factorization.
#
# The functions marked "hot kernel" below are the elementwise passes that
# dominate runtime whenever the SVD is cheap (e.g. tall, narrow matrices).
# Their bodies are in C++ (src/svt_kernels.cpp), which reproduces the R
# expressions operation for operation; the R wrappers keep the signatures.


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
#' We need \eqn{Y^{(0)} \ge Y^{(1)} \ge \cdots}. A stack that
#' violates nesting is not merely unusual, it invalidates the parameterization
#' \eqn{H^{(t)} \le 0}, so we need to verify it.
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


#' Cross entropy of one binary mask
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
  logistic_ce_cpp(X, Y, train)
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
  if (!is.double(nu)) storage.mode(nu) <- "double"
  stack_forward_cpp(Z, as.double(mu0), nu)
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
  for (t in seq.int(2L, Tn)) {
    level <- if (is.matrix(nu)) -nu[t - 1L, ] else -nu[t - 1L]
    Z[[t]] <- clip_block(Z[[t]], level)$mat
  }
  Z
}


#' Cap one block at a per-column level
#'
#' @param mat a numeric matrix
#' @param level a scalar, or one level per column
#' @returns a list with `mat`, equal to `pmin(mat, level)` with `level`
#'   broadcast down the rows, and `n_clipped`, the number of entries changed
#' @details Hot kernel: the clipping step of [bmfsvt()] and [clip_increments()].
#' @keywords internal
clip_block <- function(mat, level) {
  if (!is.double(mat)) storage.mode(mat) <- "double"
  clip_block_cpp(mat, as.double(level))
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
#' monotone or not. Every entry whose mass sits in the top interval therefore
#' gets \eqn{m_{Tj}}, whatever its own value; see `beta` in [predict.dfbm()].
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


#' Probability the model leaves on each entry's observed interval
#'
#' @param prob list of T n x P matrices, \eqn{S^{(t)}_{ij} = P(A_{ij} > d_{tj})}
#' @param Y list of T nested n x P binary masks of the same entries
#' @returns an n x P matrix \eqn{\pi^{obs}_{ij} = S^{(k)}_{ij} - S^{(k+1)}_{ij}}
#'   with \eqn{k = \sum_t Y^{(t)}_{ij}} the observed interval, \eqn{S^{(0)} = 1}
#'   and \eqn{S^{(T+1)} = 0}
#' @details
#' Hot kernel. Because the masks are nested, \eqn{Y^{(t)} - Y^{(t+1)}} is the
#' indicator of the observed interval being t, so
#' \deqn{\pi^{obs} = (1 - Y^{(1)})(1 - S^{(1)}) +
#'   \sum_{t=1}^{T-1} (Y^{(t)} - Y^{(t+1)})(S^{(t)} - S^{(t+1)}) +
#'   Y^{(T)} S^{(T)}.}
#' Passing the masks as `prob` gives 1 everywhere. Negative values, possible
#' only for probabilities that are not monotone in t, are set to 0.
#' @keywords internal
observed_bin_prob <- function(prob, Y) {
  Tn <- length(prob)
  if (length(Y) != Tn) stop("`prob` and `Y` must have the same length.")
  out <- (1 - Y[[1L]]) * (1 - prob[[1L]])
  if (Tn > 1L) {
    for (t in seq_len(Tn - 1L)) {
      out <- out + (Y[[t]] - Y[[t + 1L]]) * (prob[[t]] - prob[[t + 1L]])
    }
  }
  out <- out + Y[[Tn]] * prob[[Tn]]
  pmax(out, 0)
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
#' With `train`, held-out entries contribute a zero residual. The proximal step
#' then leaves them at the current low-rank fit, so they are imputed as in
#' softImpute.
#' @keywords internal
grad_backward <- function(X, Y, train = NULL) {
  grad_backward_cpp(X, Y, train)
}


#' Gradient, and optionally cross entropy, at the stack built from increments
#'
#' @param Z list of T increment matrices
#' @param mu0,nu offsets, see [stack_forward()]
#' @param Y list of T binary matrices
#' @param train optional training mask, see [grad_backward()]
#' @param ce whether to also return the stacked cross entropy
#' @returns a list with `Psi`, equal to
#'   `grad_backward(stack_forward(Z, mu0, nu), Y, train)`, and `f`, equal to
#'   `stack_ce(stack_forward(Z, mu0, nu), Y, train)` when `ce` and `NA`
#'   otherwise
#' @details
#' Hot kernel. Each entry's \eqn{X^{(1..T)}} is built on the fly in the order
#' of operations of [stack_forward()], so the natural parameter stack is never
#' allocated. With `ce`, one exponential serves both the sigmoid and the cross
#' entropy, which moves the sigmoid by an ulp or two.
#' @keywords internal
stack_grad <- function(Z, mu0, nu, Y, train = NULL, ce = FALSE) {
  if (!is.double(nu)) storage.mode(nu) <- "double"
  stack_grad_cpp(Z, as.double(mu0), nu, Y, train, ce)
}


#' Objective and quadratic model at a candidate proximal step
#'
#' @param Zp list of T candidate increment matrices
#' @param W list of T extrapolation point matrices
#' @param Psi gradient at `W`, see [grad_backward()]
#' @param L Lipschitz constants, length T
#' @param mu0,nu offsets, see [stack_forward()]
#' @param Y list of T binary matrices
#' @param train optional training mask
#' @param f_W stacked cross entropy at `W`
#' @returns a list with `f = stack_ce(stack_forward(Zp, mu0, nu), Y, train)`
#'   and the backtracking bound
#'   \eqn{Q = f_W + \sum_t \langle Z'_t - W_t, \Psi_t\rangle +
#'   \frac{L_t}{2}\|Z'_t - W_t\|_F^2}
#' @details Hot kernel: one backtracking step of [bmfsvt()] in a single pass,
#' building each entry's natural parameters on the fly as [stack_grad()] does.
#' @keywords internal
step_eval <- function(Zp, W, Psi, L, mu0, nu, Y, train, f_W) {
  if (!is.double(nu)) storage.mode(nu) <- "double"
  step_eval_cpp(Zp, W, Psi, as.double(L), as.double(mu0), nu, Y, train, f_W)
}


#' Convergence and restart statistics of a FISTA step
#'
#' @param Zp,Z,W lists of T matrices: the new iterate, the previous iterate and
#'   the extrapolation point
#' @returns a named vector: `num` \eqn{\sum_t\|Z'_t - Z_t\|_F^2}, `den`
#'   \eqn{\sum_t\|Z_t\|_F^2} and `ip` \eqn{\sum_t\langle W_t - Z'_t, Z'_t -
#'   Z_t\rangle}
#' @keywords internal
fista_stats <- function(Zp, Z, W) {
  fista_stats_cpp(Zp, Z, W)
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
#' Follows the softImpute prescription: compute a partial SVD of size `rank_guess`, and if its smallest singular
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
  # approaches it.
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
    sv <- if (method == "auto" && mindim <= full_cutoff) gram_svd(M, thresh)
    if (is.null(sv)) sv <- full_svd() else n_svd <- n_svd + 1L
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
      # exit if the smallest singular value is smaller than threshold
      if (min(d) <= thresh || k >= rank_cap) break
      # Grow rank guess geometrically
      k <- min(max(k + rank_step, ceiling(k * 1.5)), rank_cap)
      # cap the number of guessing and go to full SVD if necessary
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
  v <- sv$v[, keep, drop = FALSE]
  mat <- if (isTRUE(sv$tall)) {
    # M V_k diag((d - thresh) / d) V_k' is the same matrix as U_k D_k V_k', and
    # for a narrow M costs one product with a ncol x ncol matrix.
    M %*% (v %*% ((dshrunk / d[keep]) * t(v)))
  } else {
    sv$u[, keep, drop = FALSE] %*% (dshrunk * t(v))
  }

  # `v` is kept so that a fitted block can fold in new rows; see predict.bmfsvt().
  list(mat = mat, rank = length(keep), d = dshrunk, v = v,
       rank_next = rank_next, n_svd = n_svd)
}



#' Singular value decomposition through the smaller Gram matrix
#'
#' @param M matrix to decompose
#' @param thresh soft thresholding level that will be applied, see [soft_svt()]
#' @param min_ratio give up, returning `NULL`, when `thresh` is below this
#'   fraction of the largest singular value
#' @returns `NULL`, or a list with all singular values `d` (decreasing), `v`
#'   and `u` holding at least the singular vectors of the components above
#'   `thresh` in their leading columns, and `tall`. When `tall` is `TRUE`, `u`
#'   is `NULL`: [soft_svt()] then rebuilds the block from `M` and `v` alone.
#' @details
#' Used by [soft_svt()] under `method = "auto"` when the smaller dimension is
#' at most `full_cutoff`. For the tall, narrow blocks of an abundance matrix the
#' symmetric eigendecomposition of the k x k Gram matrix, k = min(dim(M)), costs
#' a fraction of a full [svd()]. For a wide M the right vectors follow from one
#' matrix product, only for the components that are retained; a tall M needs
#' no left vectors at all.
#'
#' Squaring loses precision only in singular values far below the largest:
#' the absolute error of a soft thresholded component is about
#' \eqn{\epsilon d_1^2 / (2\,\mathrm{thresh})}. Below `min_ratio` that is no
#' longer negligible against \eqn{\epsilon d_1} and the caller falls back to
#' [svd()].
#' @keywords internal
gram_svd <- function(M, thresh, min_ratio = 1e-4) {
  tall <- nrow(M) >= ncol(M)
  e <- eigen(if (tall) crossprod(M) else tcrossprod(M), symmetric = TRUE)
  d <- sqrt(pmax(e$values, 0))
  if (d[1L] > 0 && thresh < min_ratio * d[1L]) return(NULL)
  if (tall) return(list(d = d, u = NULL, v = e$vectors, tall = TRUE))
  idx <- seq_len(sum(d > thresh))
  v <- crossprod(M, e$vectors[, idx, drop = FALSE]) /
    rep(d[idx], each = ncol(M))
  list(d = d, u = e$vectors, v = v, tall = FALSE)
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
#' evaluated on a stack carrying *no signal*.
#' The search range is \eqn{\lambda_t = \alpha \lambda^*_t} with
#' \eqn{\alpha \in (0, 1]}.
#'
#' **The noise stack must be nested, and one shared uniform per entry is what
#' makes it so.** Drawing \eqn{U_{ij} \sim U(0,1)} once and setting
#' \eqn{Y^{(t')c*}_{ij} = I[U_{ij} < \pi_{t'}]} for every \eqn{t'} gives
#' \eqn{Y^{(0)c*} \ge Y^{(1)c*} \ge \cdots} automatically.
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
