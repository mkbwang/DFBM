
# From an abundance matrix to denoised values: feature specific thresholds,
# nested binary masks, interval representatives, reconstruction from the fitted
# survival curves, and closure. dfbm() fits all of it on a training matrix;
# predict.dfbm() applies the fitted object to new rows without refitting, so a
# train/test comparison never lets test samples inform the denoiser.


#' Denoise an abundance matrix by factorizing binary masks
#'
#' @param A training abundance matrix (or data frame), samples in rows,
#'   nonnegative, no missing values
#' @param alpha `NULL` to choose the shrinkage by entry-wise cross-validation,
#'   [cv.bmfsvt()], or a single number to skip tuning and fit once at
#'   \eqn{\lambda = \alpha \lambda^*}. A fixed `alpha` costs one fit instead of
#'   `n_splits` whole paths and is meant for a quick first look
#' @param levels,zero_tol,q_cap threshold settings, see [choose_thresholds()]
#' @param summary interval summary, see [interval_values()]
#' @param beta how far each entry's observed interval is represented by the
#'   entry's own capped value instead of the interval summary, in \eqn{[0, 1]};
#'   0 (default) is the plain survival expectation. Not tuned here, see
#'   [predict.dfbm()]. Stored in the object as the default for
#'   [predict.dfbm()]
#' @param alpha_grid grid searched when `alpha` is `NULL`
#' @param loss validation loss selecting `alpha`, see [cv.bmfsvt()]. `"rps"`
#'   scores the survival probabilities and is free of column scale; `"mse"` and
#'   `"crps"` score the denoised values against the capped truth, each column
#'   divided by its own standard deviation. All three are reported in `cv`.
#' @param n_splits,train_prop number of train/validation splits and the
#'   training proportion of each row and column, see [cv.bmfsvt()]
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
#'   `cap`, `levels`, `cap_level`, `zero_frac`, `zero_tol`, the
#'   representatives `M`, `summary`, `beta`, `close`, `tol`, `colnames`), the
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
#'    so zeros are filled in, and it lies in \eqn{[m_{0j}, m_{Tj}]}, so outliers
#'    are capped. An entry's own magnitude never enters, only its masks.
#' 5. with `beta > 0`, the within-interval correction of [predict.dfbm()];
#' 6. optionally, closure to unit row sums.
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
#' \eqn{\sqrt{P}(T - t + 1)} to the gradient of block t, which is typically
#' well below \eqn{\lambda_t}, so no pattern carried by one sample can enter the
#' fit. `row_ce` measures how badly a row's fitted probabilities match its own
#' masks and `flag` marks rows above the training `flag_prob` quantile; the
#' binned row sums (`binned_row_sum`) are a second signal, since the denoised
#' row sums stay near one even for an erased row.
#'
#' @seealso [predict.dfbm()] to apply the fit to new rows.
#' @export
dfbm <- function(A, alpha = NULL, levels = seq(0.05, 0.95, by = 0.1),
                 zero_tol = 1e-12, q_cap = NULL,
                 summary = c("mean", "median"), beta = 0,
                 alpha_grid = exp(seq(log(1), log(0.2), length.out = 10L)),
                 loss = c("rps", "mse", "crps"), n_splits = 5L,
                 train_prop = 0.8,
                 close = TRUE,
                 flag_prob = 0.999, C = 20L, seed = NULL, tol = 1e-6, ...) {
  summary <- match.arg(summary)
  loss <- match.arg(loss)
  A <- check_abundance(A)
  check_beta(beta)

  th <- choose_thresholds(A, levels = levels, zero_tol = zero_tol,
                          q_cap = q_cap)
  Y <- make_masks(A, th$thresholds)

  if (is.null(alpha)) {
    tuned <- cv.bmfsvt(Y, alpha_grid = alpha_grid, loss = loss,
                       n_splits = n_splits, train_prop = train_prop,
                       A = A, thresholds = th$thresholds, cap = th$cap,
                       summary = summary, C = C, seed = seed,
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
    fit <- bmfsvt(Y, lambda = lambda, ...)
    cv <- NULL
  }

  obj <- structure(
    list(thresholds = th$thresholds, cap = th$cap, levels = th$levels,
         cap_level = th$cap_level, zero_frac = th$zero_frac,
         zero_tol = th$zero_tol,
         M = interval_values(A, th$thresholds, th$cap, summary = summary),
         summary = summary, beta = beta, close = close, tol = tol,
         colnames = colnames(A),
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
#' @param beta weight in \eqn{[0, 1]} of the entry's own capped value within
#'   its observed interval; see Details. Changing it needs no refit.
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
#'
#' **Within-interval correction (`beta`).** The survival expectation
#' \eqn{\hat A = \sum_k \pi_k m_k}, with \eqn{\pi_k = S^{(k)} - S^{(k+1)}},
#' gives every entry whose mass sits in one interval the same value, so the
#' top interval \eqn{(d_T, \infty)}, which holds most of a skewed column's
#' spread, collapses to \eqn{m_T}. With \eqn{k} the entry's observed interval,
#' \eqn{\tilde x = \min(x, c_j)} and \eqn{\pi^{obs} = \pi_k}
#' ([observed_bin_prob()]), the representative of the observed interval alone
#' is moved toward the entry's own value:
#' \deqn{\hat A_\beta = \hat A + \beta\, \pi^{obs} (\tilde x - m_k).}
#' The thresholds fix the observed interval, the factorization decides how much
#' mass stays on it, and `beta` how much of that mass is placed at the entry's
#' own value. Properties:
#' - `beta = 0` is the plain expectation; with `beta = 1` and
#'   \eqn{\pi^{obs} = 1} the result is \eqn{\tilde x}.
#' - The result stays in \eqn{[\min(m_{0j}, \tilde x), c_j]}, so the cap still
#'   holds; it is \eqn{c_j} rather than \eqn{m_{Tj}} that now bounds it.
#' - Zeros in a zero-heavy column have \eqn{\tilde x = m_0 = 0}, so filled-in
#'   zeros are unchanged.
#' - Entries the model moves out of their observed interval keep little of
#'   their own value, in proportion to \eqn{\pi^{obs}}.
#'
#' **`beta` is not chosen by [cv.bmfsvt()], on purpose.** Its best value is set
#' by how much of the within-interval spread is noise, which a single noisy
#' copy of the data cannot reveal: a held-out entry's own value would be both
#' input and target, and any such criterion drifts to `beta = 1`. Choose it by
#' the downstream analysis, e.g. a grid such as `c(0, 0.5, 1)`, which costs one
#' fold-in per value and no refit.
#' @export
predict.dfbm <- function(object, newdata,
                         beta = if (is.null(object$beta)) 0 else object$beta,
                         close = object$close,
                         tol = object$tol, max_iter = 2000L, ...) {
  check_beta(beta)
  A <- check_abundance(newdata)
  Y <- make_masks(A, object$thresholds)
  pr <- predict.bmfsvt(object$fit, Y, tol = tol, max_iter = max_iter)

  den <- survival_expectation(pr$prob, object$M) # denoised value
  bin <- survival_expectation(Y, object$M) # represented by the bin summary
  if (beta > 0) {
    # Pull the observed interval's representative (bin) toward the entry's own
    # capped value, weighted by the mass the model leaves on that interval.
    xt <- pmin(A, rep(object$cap, each = nrow(A)))
    den <- den + beta * observed_bin_prob(pr$prob, Y) * (xt - bin)
  }
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


