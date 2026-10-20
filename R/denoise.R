
# From an abundance matrix to denoised values: feature specific thresholds,
# nested binary masks, interval representatives, reconstruction from the fitted
# survival curves, and closure. dfbm() fits all of it on a training matrix;
# predict.dfbm() applies the fitted object to new rows without refitting, so a
# train/test comparison never lets test samples inform the denoiser.


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
#' @param q_cap probability level of the cap, above `max(levels)` and at most 1
#'
#' @returns a list with components
#'   \describe{
#'     \item{thresholds}{T x P matrix \eqn{d_{tj}}, column names from `A`}
#'     \item{cap}{length P vector \eqn{c_j}}
#'     \item{levels}{T x P matrix of the probability levels actually used}
#'     \item{cap_level}{length P vector of cap levels}
#'     \item{zero_frac}{length P training fraction of entries at or below
#'       `zero_tol`}
#'     \item{zero_tol, base_levels}{the inputs}
#'   }
#'
#' @details
#' Column scales differ by orders of magnitude, so each column gets its own
#' thresholds from its own empirical quantiles. With
#' \eqn{z_j} the zero fraction of column j and \eqn{z^*_j = \max(z_j, q_1)},
#' \deqn{l_{tj} = z^*_j + (1 - z^*_j)\frac{q_t - q_1}{1 - q_1}, \qquad
#'   d_{tj} = \hat Q_j(l_{tj}),}
#' and \eqn{d_{1j}} is set to `zero_tol` when \eqn{z_j \ge q_1}. A column with
#' few zeros therefore gets plain quantiles, \eqn{l_{tj} = q_t}. A zero-heavy
#' column gets a first mask meaning "nonzero" and spreads the remaining T-1
#' thresholds over its positive values. Without that adjustment such a column
#' would produce several identical masks, counting the zero/nonzero event
#' several times in the likelihood and pinning those blocks to the clipping
#' boundary.
#'
#' **Quantiles are type 1**, i.e. order statistics. With the strict inequality
#' of [make_masks()], the training prevalence of mask t is then exactly
#' \eqn{1 - \lceil n\,l_{tj} \rceil / n} in every continuous column, the masks
#' depend only on within-column ranks, and an interval between thresholds can
#' be empty only when two thresholds tie. Interpolating quantiles (type 7)
#' produce thresholds between zero and the smallest positive value that look
#' distinct but select the same entries.
#'
#' The cap uses the same adjustment,
#' \eqn{c_j = \max(\hat Q_j(z^*_j + (1 - z^*_j)(q_{cap} - q_1)/(1 - q_1)),
#' d_{Tj})}; see [interval_values()].
#' @importFrom stats quantile
#' @export
choose_thresholds <- function(A, levels = seq(0.05, 0.95, by = 0.1),
                              zero_tol = 1e-12, q_cap = 0.99) {
  A <- check_abundance(A)
  if (!is.numeric(levels) || length(levels) < 1L || anyNA(levels) ||
      any(levels <= 0 | levels >= 1)) {
    stop("`levels` must lie strictly between 0 and 1.")
  }
  if (is.unsorted(levels, strictly = TRUE)) {
    stop("`levels` must be strictly increasing.")
  }
  if (!is.numeric(q_cap) || length(q_cap) != 1L || q_cap <= max(levels) ||
      q_cap > 1) {
    stop("`q_cap` must be a single number in (max(levels), 1].")
  }

  P <- ncol(A)
  Tn <- length(levels)
  q1 <- levels[1L]
  zero_frac <- colMeans(A <= zero_tol)
  zstar <- pmax(zero_frac, q1)
  lev <- outer((levels - q1) / (1 - q1), 1 - zstar) + rep(zstar, each = Tn)
  cap_level <- zstar + (1 - zstar) * (q_cap - q1) / (1 - q1)

  thr <- matrix(NA_real_, Tn, P)
  cap <- numeric(P)
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
       zero_frac = zero_frac, zero_tol = zero_tol, base_levels = levels)
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
#' The representatives are non-decreasing in t within every column, which is
#' what makes [survival_expectation()] bounded and monotone in each survival
#' probability.
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

  M <- matrix(NA_real_, Tn + 1L, P)
  for (j in seq_len(P)) {
    d <- thresholds[, j]
    a <- if (is.null(train)) A[, j] else A[train[, j] == 1, j]
    # bin k holds d_k < x <= d_{k+1}, i.e. exactly the entries whose masks
    # 1..k are one and k+1..T are zero.
    bin <- findInterval(a, d, left.open = TRUE)
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


#' Denoise an abundance matrix by factorizing binary masks
#'
#' @param A training abundance matrix (or data frame), samples in rows,
#'   nonnegative, no missing values
#' @param alpha `NULL` to choose the shrinkage by entry-wise cross-validation,
#'   [cv.bmfsvt()], or a single number to skip tuning and fit once at
#'   \eqn{\lambda = \alpha \lambda^*}. A fixed `alpha` costs one fit instead of
#'   `n_splits` whole paths and is meant for a quick first look; the measured
#'   optimum in simulation is about 0.585.
#' @param levels,zero_tol,q_cap threshold settings, see [choose_thresholds()]
#' @param summary interval summary, see [interval_values()]
#' @param alpha_grid grid searched when `alpha` is `NULL`
#' @param loss validation loss selecting `alpha`, see [cv.bmfsvt()]. `"rps"`
#'   scores the survival probabilities and is free of column scale; `"mse"` and
#'   `"crps"` score the denoised values against the capped truth, each column
#'   divided by its own standard deviation. All three are reported in `cv`.
#' @param n_splits,train_prop number of train/validation splits and the
#'   training proportion of each row and column, see [cv.bmfsvt()]
#' @param offset offset parameterization passed to [bmfsvt()]. `"column"` is the
#'   default here because zero-heavy columns have prevalences far from the
#'   others.
#' @param close whether to divide each denoised row by its sum, for
#'   compositional data
#' @param flag_prob rows whose misfit `row_ce` exceeds this quantile of the
#'   training misfits are flagged
#' @param C pure noise replicates for [lambda_star_seq()]
#' @param seed optional seed making [lambda_star_seq()] and the
#'   cross-validation splits reproducible
#' @param tol per-row tolerance of the fold-in, see [predict.bmfsvt()]
#' @param ... further arguments passed to [bmfsvt()], or to [cv.bmfsvt()] when
#'   `alpha` is `NULL` (e.g. `ncores`, `rule`)
#'
#' @returns an object of class `"dfbm"`: the fitted settings (`thresholds`,
#'   `cap`, `levels`, `cap_level`, `zero_frac`, `zero_tol`, `base_levels`, the
#'   representatives `M`, `summary`, `close`, `tol`, `colnames`), the
#'   factorization `fit` with its `alpha` and, when tuned, `cv` (the per-alpha
#'   validation losses of [cv.bmfsvt()]), the
#'   misfit cutoff `row_ce_cut`, `insample_gap` (mean and max absolute
#'   difference between the fold-in and in-sample probabilities), and the
#'   training output of [predict.dfbm()]: `denoised`, `prob`, `row_sum`,
#'   `binned`, `binned_row_sum`, `row_ce`, `flag`, `converged`.
#'
#' @details
#' The steps are:
#' 1. per-column thresholds [choose_thresholds()] and nested masks
#'    [make_masks()];
#' 2. a joint factorization of the masks, [bmfsvt()] with `alpha` given or
#'    chosen by [cv.bmfsvt()], giving
#'    \eqn{S^{(t)}_{ij} = P(A_{ij} > d_{tj})};
#' 3. capped interval representatives [interval_values()];
#' 4. the denoised value
#'    \deqn{\hat A_{ij} = \sum_{t=0}^{T}(S^{(t)}_{ij} - S^{(t+1)}_{ij})\,m_{tj}
#'      = m_{0j} + \sum_{t=1}^{T} S^{(t)}_{ij}(m_{tj} - m_{t-1,j}),}
#'    with \eqn{S^{(0)} = 1} and \eqn{S^{(T+1)} = 0}, an estimate of
#'    \eqn{E\min(A_{ij}, c_j)}. It is positive wherever \eqn{m_{1j} > m_{0j}},
#'    so zeros are filled in, and it never exceeds \eqn{c_j}, so outliers are
#'    capped. An entry's own magnitude never enters, only its masks.
#' 5. optionally, closure to unit row sums.
#'
#' **The training output comes from the fold-in, not the in-sample fit.** Both
#' training and new rows go through [predict.dfbm()], so a model trained on the
#' denoised training data sees features built by exactly the map later applied
#' to test data. `insample_gap` reports how far this is from `fit$prob`.
#'
#' **Why the posterior median is not offered.** For a zero-heavy column most
#' zero entries have \eqn{S^{(1)} < 0.5}, so the median of their estimated
#' distribution is the bottom interval, i.e. zero again.
#'
#' **Binned control.** `binned` applies the same representatives to the
#' observed masks instead of the fitted probabilities, which snaps each entry
#' to the representative of its own interval. Comparing against it separates
#' the effect of the factorization from that of binning and capping.
#'
#' **Rows with a unique profile are erased, not merely capped, so they are
#' flagged.** Every row is expressed through column patterns shared by the
#' training population. A single row contributes at most
#' \eqn{\sqrt{P}(T - t + 1)} to the gradient of block t, which on the HRS data is
#' well below \eqn{\lambda_t}, so no pattern carried by one sample can enter the
#' fit. `row_ce` measures how badly a row's fitted probabilities match its own
#' masks and `flag` marks rows above the training `flag_prob` quantile; the
#' binned row sums (`binned_row_sum`) are a second signal, since the denoised
#' row sums stay near one even for an erased row.
#'
#' @seealso [predict.dfbm()] to apply the fit to new rows.
#' @export
dfbm <- function(A, alpha = NULL, levels = seq(0.05, 0.95, by = 0.1),
                 zero_tol = 1e-12, q_cap = 0.99,
                 summary = c("mean", "median"),
                 alpha_grid = exp(seq(log(1), log(0.2), length.out = 10L)),
                 loss = c("rps", "mse", "crps"), n_splits = 5L,
                 train_prop = 0.8,
                 offset = c("column", "scalar"), close = TRUE,
                 flag_prob = 0.999, C = 20L, seed = NULL, tol = 1e-6, ...) {
  summary <- match.arg(summary)
  loss <- match.arg(loss)
  offset <- match.arg(offset)
  A <- check_abundance(A)

  th <- choose_thresholds(A, levels = levels, zero_tol = zero_tol,
                          q_cap = q_cap)
  Y <- make_masks(A, th$thresholds)

  if (is.null(alpha)) {
    tuned <- cv.bmfsvt(Y, alpha_grid = alpha_grid, loss = loss,
                       n_splits = n_splits, train_prop = train_prop,
                       A = A, thresholds = th$thresholds, cap = th$cap,
                       summary = summary, C = C, seed = seed, offset = offset,
                       ...)
    fit <- tuned$fit
    alpha <- tuned$alpha
    cv <- tuned$cv
  } else {
    if (!is.numeric(alpha) || length(alpha) != 1L || !(alpha > 0)) {
      stop("`alpha` must be NULL or a single positive number.")
    }
    # lambda is formed here rather than inside bmfsvt(), which has no `seed`
    # and would draw the noise floor from the global RNG stream.
    lambda <- alpha * lambda_star_seq(Y, C = C, seed = seed)
    fit <- bmfsvt(Y, lambda = lambda, offset = offset, ...)
    cv <- NULL
  }

  obj <- structure(
    list(thresholds = th$thresholds, cap = th$cap, levels = th$levels,
         cap_level = th$cap_level, zero_frac = th$zero_frac,
         zero_tol = th$zero_tol, base_levels = th$base_levels,
         M = interval_values(A, th$thresholds, th$cap, summary = summary),
         summary = summary, close = close, tol = tol, colnames = colnames(A),
         fit = fit, alpha = alpha, cv = cv,
         flag_prob = flag_prob, row_ce_cut = NA_real_),
    class = "dfbm")

  out <- predict.dfbm(obj, A)
  obj$row_ce_cut <- stats::quantile(out$row_ce, flag_prob, names = FALSE)
  out$flag <- out$row_ce > obj$row_ce_cut

  gaps <- vapply(seq_along(out$prob), function(t) {
    dp <- abs(out$prob[[t]] - fit$prob[[t]])
    c(mean(dp), max(dp))
  }, numeric(2L))
  obj$insample_gap <- c(mean = mean(gaps[1L, ]), max = max(gaps[2L, ]))

  for (nm in names(out)) obj[[nm]] <- out[[nm]]
  obj
}


#' Denoise new rows with a fitted dfbm object
#'
#' @param object an object from [dfbm()]
#' @param newdata abundance matrix with the same columns, in the same order, as
#'   the training matrix
#' @param close whether to close the denoised rows to unit sum
#' @param tol per-row fold-in tolerance
#' @param max_iter maximum number of fold-in iterations
#' @param ... unused
#'
#' @returns a list with
#'   \describe{
#'     \item{denoised}{n x P denoised matrix, closed if `close`}
#'     \item{prob}{list of T n x P survival probability matrices}
#'     \item{row_sum}{row sums of the denoised values before closure}
#'     \item{binned, binned_row_sum}{the binned control and its row sums before
#'       closure, see [dfbm()]}
#'     \item{row_ce, flag}{per-row misfit and whether it exceeds the training
#'       cutoff}
#'     \item{converged}{per-row fold-in convergence}
#'   }
#'
#' @details
#' Everything is taken from the training fit: thresholds, interval
#' representatives, offsets, shrinkage and column factors. Each new row is
#' denoised on its own by [predict.bmfsvt()], so results for a row do not
#' depend on which other rows are passed with it.
#' @export
predict.dfbm <- function(object, newdata, close = object$close,
                         tol = object$tol, max_iter = 2000L, ...) {
  A <- check_abundance(newdata)
  Y <- make_masks(A, object$thresholds)
  pr <- predict.bmfsvt(object$fit, Y, tol = tol, max_iter = max_iter)

  den <- survival_expectation(pr$prob, object$M)
  bin <- survival_expectation(Y, object$M)
  row_sum <- rowSums(den)
  binned_row_sum <- rowSums(bin)
  if (close) {
    den <- closure(den)
    bin <- closure(bin)
  }
  dimnames(den) <- dimnames(bin) <- list(rownames(A), colnames(object$thresholds))

  list(denoised = den, prob = pr$prob,
       row_sum = row_sum,
       binned = bin, binned_row_sum = binned_row_sum,
       row_ce = pr$row_ce, flag = pr$row_ce > object$row_ce_cut,
       converged = pr$converged)
}
