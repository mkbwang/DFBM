

# utility functions that include
# sanity checks (`check_abundance`, `check_columns`)
# choose thresholds(choose_thresholds) and generating masks (make_masks)


#' Check the within-interval weight
#'
#' @param beta a single number in \eqn{[0, 1]}
#' @returns `NULL`, invisibly; stops otherwise
#' @keywords internal
check_beta <- function(beta) {
  if (!is.numeric(beta) || length(beta) != 1L || is.na(beta) ||
      beta < 0 || beta > 1) {
    stop("`beta` must be a single number in [0, 1].")
  }
  invisible(NULL)
}

#' Validate an abundance matrix
#'
#' @param A a numeric matrix or data frame, samples in rows
#' @returns `A` as a double matrix
#' @keywords internal
check_abundance <- function(A) {
  if (is.data.frame(A)) A <- as.matrix(A)
  if (!is.matrix(A) || !is.numeric(A)) {
    stop("`A` must be a numeric matrix or a data frame of numeric columns.")
  }
  if (nrow(A) < 1L || ncol(A) < 1L) stop("`A` is empty.")
  if (anyNA(A)) stop("`A` contains missing values, which are not supported.")
  if (any(!is.finite(A))) stop("`A` contains infinite values.")
  if (any(A < 0)) stop("`A` must be nonnegative.")
  storage.mode(A) <- "double"
  A
}




#' Check that a matrix matches the columns thresholds were chosen on
#'
#' @param A abundance matrix
#' @param thresholds T x P threshold matrix
#' @returns `NULL`, invisibly; stops on a mismatch
#' @keywords internal
check_columns <- function(A, thresholds) {
  if (ncol(A) != ncol(thresholds)) {
    stop(sprintf("`A` has %d columns but the thresholds have %d.",
                 ncol(A), ncol(thresholds)))
  }
  ref <- colnames(thresholds)
  if (!is.null(colnames(A)) && !is.null(ref) && !identical(colnames(A), ref)) {
    stop("the column names of `A` differ from those the thresholds were ",
         "chosen on, or are in a different order.")
  }
  invisible(NULL)
}



#' Feature specific thresholds for the binary masks
#'
#' @param A training abundance matrix, samples in rows, nonnegative
#' @param levels base probability levels \eqn{0 < q_1 < \cdots < q_T < 1}, one
#'   per mask
#' @param zero_tol entries at or below this count as zero. The default absorbs
#'   floating point residues left by subtraction without touching real values.
#' @param q_cap probability level of the cap, above `max(levels)` and at most
#'   1. `NULL` (default) uses 0.99, or halfway between `max(levels)` and 1 when
#'   the levels reach 0.99, so that thresholds refined in the tail need no
#'   manual cap
#'
#' @returns a list with components
#'   \describe{
#'     \item{thresholds}{T x P matrix \eqn{d_{tj}}, column names from `A`}
#'     \item{cap}{length P vector \eqn{c_j}}
#'     \item{levels}{T x P matrix of the probability levels actually used}
#'     \item{cap_level}{length P vector of cap levels}
#'     \item{zero_frac}{length P training fraction of entries at or below
#'       `zero_tol`}
#'     \item{zero_tol}{the input}
#'   }
#'
#' @details
#' Column scales differ by orders of magnitude, so each column gets its own
#' thresholds from its own empirical quantiles. With
#' \eqn{z_j} the zero fraction of column j and \eqn{z^*_j = \max(z_j, q_1)},
#' \deqn{l_{tj} = z^*_j + (1 - z^*_j)\frac{q_t - q_1}{1 - q_1}, \qquad
#'   d_{tj} = \hat Q_j(l_{tj}),}
#' and \eqn{d_{1j}} is set to `zero_tol` when \eqn{z_j \ge q_1}.
#' A column with few zeros gets plain quantiles, \eqn{l_{tj} = q_t}.
#' A zero-heavy column gets a first mask meaning "nonzero" and spreads the remaining T-1
#' thresholds over its positive values. Without that adjustment such a column
#' would produce several identical masks, counting the zero/nonzero event
#' several times in the likelihood and pinning those blocks to the clipping
#' boundary.
#'
#'
#' The cap uses the same adjustment,
#' \eqn{c_j = \max(\hat Q_j(z^*_j + (1 - z^*_j)(q_{cap} - q_1)/(1 - q_1)),
#' d_{Tj})}; see [interval_values()].
#' @importFrom stats quantile
#' @export
choose_thresholds <- function(A, levels = seq(0.05, 0.95, by = 0.1),
                              zero_tol = 1e-12, q_cap = NULL) {
  A <- check_abundance(A)
  if (!is.numeric(levels) || length(levels) < 1L || anyNA(levels) ||
      any(levels <= 0 | levels >= 1)) {
    stop("`levels` must lie strictly between 0 and 1.")
  }
  if (is.unsorted(levels, strictly = TRUE)) {
    stop("`levels` must be strictly increasing.")
  }
  if (is.null(q_cap)) q_cap <- max(0.99, (1 + max(levels)) / 2)
  if (!is.numeric(q_cap) || length(q_cap) != 1L || q_cap <= max(levels) ||
      q_cap > 1) {
    stop("`q_cap` must be a single number in (max(levels), 1].")
  }

  P <- ncol(A)
  Tn <- length(levels)
  q1 <- levels[1L]
  zero_frac <- colMeans(A <= zero_tol) # count zero proportion with the zero tolerance
  zstar <- pmax(zero_frac, q1) # lowest starting quantile
  lev <- outer((levels - q1) / (1 - q1), 1 - zstar) + rep(zstar, each = Tn) # quantile levels, T*P
  cap_level <- zstar + (1 - zstar) * (q_cap - q1) / (1 - q1) # cap level, vector with length P

  thr <- matrix(NA_real_, Tn, P) # T threshold values for each of the P features
  cap <- numeric(P) # capping the largest value for all P features
  for (j in seq_len(P)) {
    q <- stats::quantile(A[, j], probs = c(lev[, j], cap_level[j]), type = 1,
                         names = FALSE)
    thr[, j] <- q[seq_len(Tn)]
    if (zero_frac[j] >= q1) thr[1L, j] <- zero_tol
    cap[j] <- max(q[Tn + 1L], thr[Tn, j])
  }

  tied <- which(apply(thr, 2L, function(d) any(diff(d) == 0)))
  if (length(tied) > 0L) {
    nm <- if (is.null(colnames(A))) as.character(tied) else colnames(A)[tied]
    warning("tied thresholds in column(s) ", paste(nm, collapse = ", "),
            ": these columns have too few distinct values for ", Tn,
            " masks, so some of their masks are identical.")
  }

  dimnames(thr) <- dimnames(lev) <- list(paste0("t", seq_len(Tn)), colnames(A))
  names(cap) <- names(cap_level) <- names(zero_frac) <- colnames(A)
  list(thresholds = thr, cap = cap, levels = lev, cap_level = cap_level,
       zero_frac = zero_frac, zero_tol = zero_tol)
}




#' Nested binary masks from an abundance matrix
#'
#' @param A abundance matrix, samples in rows
#' @param thresholds a T x P threshold matrix, or the list returned by
#'   [choose_thresholds()]. For new data these must be the **training**
#'   thresholds.
#' @returns a list of T n x P double matrices,
#'   \eqn{Y^{(t)}_{ij} = 1\{A_{ij} > d_{tj}\}}
#' @details
#' The thresholds are non-decreasing in t for every column, so the masks are
#' nested whatever data they are applied to.
#' @export
make_masks <- function(A, thresholds) {
  A <- check_abundance(A)
  if (is.list(thresholds)) thresholds <- thresholds$thresholds
  check_columns(A, thresholds)
  n <- nrow(A)
  P <- ncol(A)
  lapply(seq_len(nrow(thresholds)), function(t) {
    matrix(1 * (A > rep(thresholds[t, ], each = n)), n, P)
  })
}



#' Representative value of each interval between thresholds
#'
#' @param A training abundance matrix
#' @param thresholds T x P threshold matrix from [choose_thresholds()]
#' @param cap length P vector of caps from [choose_thresholds()]
#' @param summary `"mean"` for the mean of the capped values in each interval,
#'   `"median"` for their median
#' @param train optional n x P 0/1 matrix; each column is then summarized over
#'   its training entries only, so that held-out values never inform the
#'   representatives they are scored against (see [cv.bmfsvt()])
#' @returns a (T+1) x P matrix with rows `I0, ..., IT`: row 1 for
#'   \eqn{(-\infty, d_1]}, row t+1 for \eqn{(d_t, d_{t+1}]}, row T+1 for
#'   \eqn{(d_T, \infty)}
#' @details
#' \deqn{m_{tj} = \mathrm{summary}\{\min(A_{ij}, c_j) : A_{ij} \in I_{tj}\},}
#' computed on the training rows. An empty interval, which with type 1
#' quantiles happens only between tied thresholds, gets
#' \eqn{\min(d_{tj}, c_j)}.
#'
#' **Why the cap.** The top interval \eqn{(d_T, \infty)} holds the largest
#' training values, outliers included. Winsorizing at \eqn{c_j} before
#' summarizing keeps \eqn{m_{Tj} \le c_j}, so a single extreme value cannot move
#' the representative, and the denoised values built from these
#' representatives can never exceed \eqn{c_j}.
#'
#' @importFrom stats median
#' @export
interval_values <- function(A, thresholds, cap, summary = c("mean", "median"),
                            train = NULL) {
  summary <- match.arg(summary)
  A <- check_abundance(A)
  check_columns(A, thresholds)
  Tn <- nrow(thresholds)
  P <- ncol(A)
  if (length(cap) != P) stop("`cap` must have one element per column.")
  train <- check_train(train, nrow(A), P)
  f <- if (summary == "mean") mean else stats::median

  M <- matrix(NA_real_, Tn + 1L, P) # mean or median value in each of the T+1 interval
  for (j in seq_len(P)) {
    d <- thresholds[, j]
    a <- if (is.null(train)) A[, j] else A[train[, j] == 1, j]
    # bin k holds d_k < x <= d_{k+1}, i.e. exactly the entries whose masks
    # 1..k are one and k+1..T are zero.
    bin <- findInterval(a, d, left.open = TRUE) # open on the left and closed on the right
    xc <- pmin(a, cap[j])
    for (k in 0:Tn) {
      in_k <- bin == k
      M[k + 1L, j] <- if (any(in_k)) f(xc[in_k]) else
        min(d[max(k, 1L)], cap[j])
    }
  }
  if (any(diff(M) < -1e-12 * pmax(1, abs(M[-1L, , drop = FALSE])))) {
    stop("interval representatives are not non-decreasing; check that the ",
         "thresholds are non-decreasing in every column.")
  }
  dimnames(M) <- list(paste0("I", 0:Tn), colnames(thresholds))
  M
}



#' Close rows to unit sum
#'
#' @param A nonnegative matrix
#' @returns `A` with every row divided by its sum
#' @export
closure <- function(A) {
  rs <- rowSums(A)
  if (any(rs <= 0)) warning("some rows sum to zero and cannot be closed.")
  A / rs
}



