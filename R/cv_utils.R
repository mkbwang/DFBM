
# Entry-wise cross-validation for the shrinkage level of bmfsvt(): balanced
# train/validation splits of the entries, validation losses on the held-out
# entries, and cv.bmfsvt() itself. Independent of tune.bmfsvt(), which selects
# alpha without holding anything out.


#' Evaluate an expression under a fixed seed without touching the caller's RNG
#'
#' @param seed `NULL` or an integer
#' @param expr the expression to evaluate
#' @returns the value of `expr`
#' @details
#' With `seed = NULL` the expression draws from the global stream as usual.
#' Otherwise the stream is put back afterwards, as in [lambda_star_seq()], so a
#' seeded call is reproducible and leaves the caller's own draws unchanged.
#' @keywords internal
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}


#' Balanced split of the entries of a matrix into training and validation
#'
#' @param X a matrix or data frame; only its dimensions are used
#' @param prop proportion of each column assigned to training
#' @param margin unused. Kept so that calls written for the earlier linear
#'   programming version still run: the construction below balances rows more
#'   tightly than any margin that version accepted.
#' @param seed optional integer; makes the split reproducible without
#'   disturbing the caller's random number stream
#' @returns an N x P matrix of 0/1, 1 marking training entries and 0 validation
#'   entries
#'
#' @details
#' Every column holds out exactly `round(N * (1 - prop))` entries, and the
#' number held out in each row differs across rows by **at most one**, so every
#' row and every column is split in (as nearly as integers allow) the same
#' proportion. Choosing validation entries uniformly at random over the whole
#' matrix gives no such guarantee: with a 20% hold-out and 25 columns, a
#' Binomial(25, 0.2) row count leaves about 1 row in 250 with 10 or more of its
#' 25 entries held out.
#'
#' The construction is greedy. Columns are visited in random order, and each
#' holds out the rows that have so far been held out least often, with ties
#' broken at random. If all row counts are `c` or `c + 1` before a column is
#' processed, they are all `c`, `c + 1` or `c + 2` after, and the smallest
#' value is gone whenever `c + 2` appears. So the spread never exceeds one.
#'
#' It replaces a 0/1 linear program whose dense constraint matrix has
#' \eqn{2(N+P) \times NP} entries, about 31 GB at the HRS training shape
#' (6212 x 25). This version needs O(NP) memory and O(NP log N) time.
#'
#' @importFrom stats runif
#' @export
train_mask_generation <- function(X, prop = 0.8, margin = 0.02, seed = NULL) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) stop("`X` must be a matrix or a data frame.")
  if (!is.numeric(prop) || length(prop) != 1L || !(prop > 0 && prop < 1)) {
    stop("`prop` must be a single number strictly between 0 and 1.")
  }
  N <- nrow(X)
  P <- ncol(X)
  n_val <- as.integer(round(N * (1 - prop)))
  if (n_val < 1L || n_val >= N) {
    stop(sprintf(paste0("`prop` = %g holds out %d of %d rows per column; ",
                        "at least one entry must fall on each side."),
                 prop, n_val, N))
  }

  with_seed(seed, {
    train <- matrix(1, N, P)
    held <- numeric(N)
    for (j in sample.int(P)) {
      # The uniform only breaks ties: it is below 1, so it never reorders rows
      # whose held-out counts differ.
      idx <- order(held + stats::runif(N))[seq_len(n_val)]
      train[idx, j] <- 0
      held[idx] <- held[idx] + 1
    }
    train
  })
}


#' Validation losses of one fit at the held-out entries
#'
#' @param prob list of T fitted probability matrices
#' @param Y list of T binary masks
#' @param val integer positions (column major) of the held-out entries
#' @param values `NULL`, or a list with the capped truth `y`, per-entry column
#'   scale `s` and per-entry representatives `M` ((T+1) x length(val)) of the
#'   held-out entries
#' @returns a named vector: mean `rps`, `mse` and `crps` over the held-out
#'   entries, the last two `NA` when `values` is `NULL`
#' @keywords internal
cv_losses <- function(prob, Y, val, values = NULL) {
  pv <- lapply(prob, function(p) p[val])
  rps <- 0
  for (t in seq_along(pv)) rps <- rps + (pv[[t]] - Y[[t]][val])^2
  out <- c(rps = mean(rps), mse = NA_real_, crps = NA_real_)
  if (!is.null(values)) {
    # One held-out entry per column of a 1 x n "matrix", so that each entry
    # carries its own column's representatives.
    den <- as.vector(survival_expectation(lapply(pv, function(p) matrix(p, 1L)),
                                          values$M))
    out[["mse"]] <- mean(((den - values$y) / values$s)^2)
    out[["crps"]] <- mean(crps_discrete(pv, values$M, values$y) / values$s)
  }
  out
}


#' Choose the shrinkage level for bmfsvt by entry-wise cross-validation
#'
#' @param Y a list of T binary matrices, or an N x P x T array
#' @param alpha_grid shrinkage levels relative to the noise floor
#'   [lambda_star_seq()], as in [tune.bmfsvt()]; sorted decreasing internally so
#'   that each fit warm starts the next
#' @param loss validation loss that selects `alpha`: `"rps"`, `"mse"` or
#'   `"crps"`. All that can be computed are reported whichever is chosen.
#' @param n_splits number of train/validation splits
#' @param train_prop proportion of each row and column used for training, see
#'   [train_mask_generation()]
#' @param train_masks optional list of N x P 0/1 training masks to use instead
#'   of generated ones; `n_splits` and `train_prop` are then ignored
#' @param A optional N x P abundance matrix the masks were built from. Needed,
#'   together with `thresholds` and `cap`, for `"mse"` and `"crps"`.
#' @param thresholds,cap the T x P thresholds and length P caps the masks were
#'   built with, as returned by [choose_thresholds()]
#' @param summary interval summary for the representatives, see
#'   [interval_values()]
#' @param rule `"min"` selects the alpha with the smallest mean validation loss;
#'   `"1se"` the strongest shrinkage within one standard error of it
#' @param C pure noise replicates for [lambda_star_seq()]
#' @param ncores splits to fit in parallel; 1 runs serially
#' @param seed optional seed making the splits, the noise floors and hence the
#'   whole selection reproducible
#' @param verbose whether to report progress
#' @param ... further arguments passed to [bmfsvt()], e.g. `offset`
#'
#' @returns a list with
#'   \describe{
#'     \item{fit}{the fit on **all** entries at the chosen alpha}
#'     \item{alpha, index}{the chosen alpha and its position in `cv`}
#'     \item{alpha_min, alpha_1se}{the choice under either rule}
#'     \item{loss, rule}{the settings used}
#'     \item{cv}{data frame, one row per alpha (strongest shrinkage first):
#'       mean and standard error over splits of each loss (`rps`, `rps_se`,
#'       `mse`, `mse_se`, `crps`, `crps_se`; the value-scale ones `NA` without
#'       `A`), `mean_rank`, `n_iter`, `converged`}
#'     \item{loss_matrix}{list of `n_splits x K` matrices, one per loss}
#'     \item{ranks}{`K x T` retained ranks averaged over splits}
#'     \item{lambda_star}{full-data noise floor used by the refit}
#'     \item{lambda_star_splits}{`n_splits x T` masked noise floors}
#'     \item{split_seeds, n_val}{per-split seeds and validation counts}
#'   }
#'
#' @details
#' For each split, [bmfsvt()] is fitted to the training entries of every mask
#' along the whole grid, and the held-out entries, whose probabilities are
#' imputed by the low-rank blocks, are scored. The losses are averaged over the
#' held-out entries of a split, then over splits. The chosen alpha is refitted on
#' every entry, since none of the split fits saw the whole matrix.
#'
#' **Each split has its own noise floor.** A split's `lambda` is
#' `alpha * lambda_star_seq(Y, train = W)`, the floor of the gradient restricted
#' to its training entries, while the refit uses `alpha` times the full-data
#' floor. Holding out a fraction of the entries shrinks the noise gradient by
#' about the square root of the training fraction. Anchoring the split fits on
#' the full-data floor would shrink them about `1/train_prop` times harder than
#' the refit at the same `alpha`, more than one grid step at 0.8. Some mismatch
#' remains, of the order of the held-out fraction, so a larger `train_prop`
#' transfers better at the cost of a noisier loss.
#'
#' **The losses.** With \eqn{S_t} the fitted survival probabilities of a
#' held-out entry and \eqn{Y_t} its masks:
#' \itemize{
#'   \item `"rps"`, the ranked probability score \eqn{\sum_t (S_t - Y_t)^2}. It is
#'     the CRPS on the probability scale of each column's quantiles. The
#'     thresholds are per-column quantiles, so the masks depend only on
#'     within-column ranks and every column counts equally whatever its units.
#'     It scores `S`, the estimand of the factorization, and needs no values.
#'   \item `"mse"`, \eqn{(\hat A - \min(A, c_j))^2 / s_j^2}, the denoised value
#'     of [survival_expectation()] against the **capped** truth.
#'   \item `"crps"`, the CRPS of the discrete distribution whose mean is
#'     \eqn{\hat A} ([crps_discrete()]) against the capped truth, divided by
#'     \eqn{s_j}.
#' }
#' The truth is capped because \eqn{\hat A} estimates \eqn{E\min(A, c_j)} and can
#' never exceed \eqn{c_j}. An uncapped outlier would add a large error that no
#' alpha can reduce. \eqn{s_j} is the standard deviation of the capped column
#' over all rows. It is fixed across splits and alphas, so it only weights the
#' columns against each other and keeps the large-valued ones from dominating.
#' The representatives \eqn{m_{tj}} are recomputed on each split's training
#' entries only, so a held-out value never informs its own prediction.
#'
#' The thresholds themselves are taken as given, i.e. chosen on all rows. They
#' fix only within-column ranks and are the same for every alpha, so they cannot
#' favour one alpha over another.
#'
#' **Held-out selection has been measured to under-shrink.** An earlier
#' held-out Brier criterion on the masks, which is `"rps"` here, landed on
#' average +0.75 grid steps from the oracle (+8.2% RMSE). That was measured
#' under the superseded `lambda_max` anchor with unbalanced splits, so it needs
#' re-measuring, but the mechanism is structural. The fit is scored only on
#' entries it did not see, whereas the output of interest includes the entries
#' it did. The `"1se"` rule was measured to cancel that bias only by
#' coincidence, so `"min"` is the default.
#'
#' **Parallelism.** Splits are independent and are forked with
#' `parallel::mclapply`. Their seeds are drawn in the parent, so serial and
#' parallel runs agree exactly. Pin the BLAS to one thread before forking; see
#' `experiment/run.sh`.
#'
#' @seealso [tune.bmfsvt()] for selection without holding entries out,
#'   [train_mask_generation()] for the splits.
#' @importFrom stats sd
#' @export
cv.bmfsvt <- function(Y,
                      alpha_grid = exp(seq(log(1), log(0.2), length.out = 10L)),
                      loss = c("rps", "mse", "crps"), n_splits = 5L,
                      train_prop = 0.8, train_masks = NULL,
                      A = NULL, thresholds = NULL, cap = NULL,
                      summary = c("mean", "median"), rule = c("min", "1se"),
                      C = 20L, ncores = 1L, seed = NULL, verbose = FALSE, ...) {
  loss <- match.arg(loss)
  summary <- match.arg(summary)
  rule <- match.arg(rule)
  Y <- as_mask_list(Y, check_nested = TRUE)
  N <- nrow(Y[[1L]])
  P <- ncol(Y[[1L]])
  Tn <- length(Y)

  # ---- values for the value-scale losses ------------------------------------
  if (is.list(thresholds)) {
    if (is.null(cap)) cap <- thresholds$cap
    thresholds <- thresholds$thresholds
  }
  has_values <- !is.null(A) && !is.null(thresholds) && !is.null(cap)
  if (loss != "rps" && !has_values) {
    stop(sprintf("`loss = \"%s\"` needs `A`, `thresholds` and `cap`.", loss))
  }
  if (has_values) {
    A <- check_abundance(A)
    check_columns(A, thresholds)
    if (nrow(A) != N || ncol(A) != P || nrow(thresholds) != Tn) {
      stop("`A` and `thresholds` must match the masks: N x P and T x P.")
    }
    if (length(cap) != P) stop("`cap` must have one element per column.")
    yc <- pmin(A, rep(cap, each = N))
    scale <- apply(yc, 2L, stats::sd)
    # A constant column has nothing to predict; any positive scale will do.
    scale[!(scale > 0)] <- 1
  }

  # ---- splits ----------------------------------------------------------------
  alpha_grid <- sort(unique(alpha_grid), decreasing = TRUE)
  K <- length(alpha_grid)
  if (!is.null(train_masks)) {
    if (is.matrix(train_masks)) train_masks <- list(train_masks)
    n_splits <- length(train_masks)
  }
  n_splits <- as.integer(n_splits)
  if (n_splits < 1L) stop("`n_splits` must be at least 1.")
  split_seeds <- with_seed(seed, sample.int(.Machine$integer.max, n_splits))

  one_split <- function(k) {
    sd_k <- split_seeds[k]
    W <- if (is.null(train_masks)) {
      train_mask_generation(Y[[1L]], prop = train_prop, seed = sd_k)
    } else {
      train_masks[[k]]
    }
    W <- check_train(W, N, P)
    if (is.null(W)) stop(sprintf("split %d holds out no entries.", k))
    val <- which(W == 0)
    lstar <- lambda_star_seq(Y, C = C, seed = sd_k, train = W)

    values <- NULL
    if (has_values) {
      Mk <- interval_values(A, thresholds, cap, summary = summary, train = W)
      jval <- (val - 1L) %/% N + 1L
      values <- list(y = yc[val], s = scale[jval],
                     M = Mk[, jval, drop = FALSE])
    }

    losses <- matrix(NA_real_, K, 3L,
                     dimnames = list(NULL, c("rps", "mse", "crps")))
    ranks <- matrix(NA_integer_, K, Tn)
    n_iter <- integer(K)
    converged <- logical(K)
    Z_warm <- NULL
    for (i in seq_len(K)) {
      fit <- bmfsvt(Y, train = W, lambda = alpha_grid[i] * lstar,
                    Z_init = Z_warm, ...)
      Z_warm <- fit$Z
      losses[i, ] <- cv_losses(fit$prob, Y, val, values)
      ranks[i, ] <- fit$ranks
      n_iter[i] <- fit$n_iter
      converged[i] <- fit$converged
    }
    if (verbose) {
      message(sprintf("split %d: %s argmin alpha=%.4g", k, loss,
                      alpha_grid[which.min(losses[, loss])]))
    }
    list(losses = losses, ranks = ranks, n_iter = n_iter,
         converged = converged, lstar = lstar, n_val = length(val))
  }

  parallel_ok <- ncores > 1L && n_splits > 1L && .Platform$OS.type != "windows"
  parts <- if (parallel_ok) {
    parallel::mclapply(seq_len(n_splits), one_split,
                       mc.cores = min(ncores, n_splits))
  } else {
    lapply(seq_len(n_splits), one_split)
  }
  bad <- vapply(parts, inherits, logical(1L), "try-error")
  if (any(bad)) {
    stop("cross-validation split failed: ",
         conditionMessage(attr(parts[[which(bad)[1L]]], "condition")))
  }

  # ---- aggregate ------------------------------------------------------------
  loss_matrix <- lapply(c(rps = "rps", mse = "mse", crps = "crps"), function(nm) {
    do.call(rbind, lapply(parts, function(p) p$losses[, nm]))
  })
  se <- function(m) {
    if (nrow(m) < 2L) rep(NA_real_, ncol(m)) else
      apply(m, 2L, stats::sd) / sqrt(nrow(m))
  }
  cv <- data.frame(
    alpha = alpha_grid,
    rps = colMeans(loss_matrix$rps), rps_se = se(loss_matrix$rps),
    mse = colMeans(loss_matrix$mse), mse_se = se(loss_matrix$mse),
    crps = colMeans(loss_matrix$crps), crps_se = se(loss_matrix$crps),
    mean_rank = rowMeans(Reduce(`+`, lapply(parts, `[[`, "ranks"))) / n_splits,
    n_iter = Reduce(`+`, lapply(parts, `[[`, "n_iter")) / n_splits,
    converged = Reduce(`&`, lapply(parts, `[[`, "converged")))

  crit <- cv[[loss]]
  crit_se <- cv[[paste0(loss, "_se")]]
  index_min <- which.min(crit)
  # The grid runs from strong to weak shrinkage, so the first index within the
  # band is the strongest shrinkage the band allows.
  index_1se <- if (is.na(crit_se[index_min])) index_min else
    min(which(crit <= crit[index_min] + crit_se[index_min]))
  index <- if (rule == "min") index_min else index_1se
  if (verbose) {
    message(sprintf("%s picks alpha=%.4g (min %.4g, 1se %.4g)", loss,
                    alpha_grid[index], alpha_grid[index_min],
                    alpha_grid[index_1se]))
  }

  # ---- refit on every entry ---------------------------------------------------
  lstar_full <- lambda_star_seq(Y, C = C, seed = seed)
  fit <- bmfsvt(Y, lambda = alpha_grid[index] * lstar_full, ...)

  list(fit = fit,
       alpha = alpha_grid[index], index = index,
       alpha_min = alpha_grid[index_min], alpha_1se = alpha_grid[index_1se],
       loss = loss, rule = rule,
       cv = cv,
       loss_matrix = loss_matrix,
       ranks = Reduce(`+`, lapply(parts, `[[`, "ranks")) / n_splits,
       lambda_star = lstar_full,
       lambda_star_splits = do.call(rbind, lapply(parts, `[[`, "lstar")),
       split_seeds = split_seeds,
       n_val = vapply(parts, `[[`, integer(1L), "n_val"))
}
