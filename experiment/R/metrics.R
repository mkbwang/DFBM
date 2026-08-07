
# Shared evaluation. Every fitting routine returns a list with `S`, the
# estimated marginal survival probability for each threshold, so all metrics
# live here rather than being duplicated per method.


#' Clamp probabilities away from the boundary before taking logs
#' @keywords internal
clamp01 <- function(p, eps = 1e-8) pmin(pmax(p, eps), 1 - eps)


#' Evaluate an estimated survival stack against the truth and held out data
#'
#' @param S_hat list of T estimated survival matrices
#' @param S_true list of T true survival matrices, or NULL
#' @param Y list of T observed binary masks
#' @param risk_set optional per-threshold risk set sizes from [fit_old()]
#'
#' @returns a data frame with one row per threshold plus a row for the overall
#'   aggregate (`t = NA`).
#'
#' @details
#' Reports both estimation error against the truth and out of sample predictive
#' accuracy, because the two can disagree: a method can be well calibrated on
#' average while badly misestimating individual entries.
eval_survival <- function(S_hat, S_true, Y, risk_set = NULL) {
  Tn <- length(S_hat)
  N <- nrow(S_hat[[1L]])
  P <- ncol(S_hat[[1L]])
  eval_idx <- seq_len(N * P)

  rows <- lapply(seq_len(Tn), function(t) {
    sh <- S_hat[[t]]
    yy <- Y[[t]][eval_idx]
    ph <- clamp01(sh[eval_idx])

    rmse <- mae <- NA_real_
    if (!is.null(S_true)) {
      err <- sh - S_true[[t]]
      rmse <- sqrt(mean(err^2))
      mae <- mean(abs(err))
    }

    data.frame(
      t = t,
      rmse = rmse,
      mae = mae,
      cross_entropy = mean(-(yy * log(ph) + (1 - yy) * log(1 - ph))),
      brier = mean((ph - yy)^2),
      auc = safe_auc(yy, ph),
      prevalence = mean(Y[[t]]),
      risk_set = if (is.null(risk_set)) NA_real_ else risk_set[t],
      mean_pred = mean(sh)
    )
  })
  per_t <- do.call(rbind, rows)

  overall <- data.frame(
    t = NA_integer_,
    rmse = if (is.null(S_true)) NA_real_ else
      sqrt(mean((unlist(S_hat) - unlist(S_true))^2)),
    mae = if (is.null(S_true)) NA_real_ else
      mean(abs(unlist(S_hat) - unlist(S_true))),
    cross_entropy = mean(per_t$cross_entropy),
    brier = mean(per_t$brier),
    auc = mean(per_t$auc, na.rm = TRUE),
    prevalence = NA_real_,
    risk_set = NA_real_,
    mean_pred = NA_real_
  )
  rbind(per_t, overall)
}


#' AUC that returns NA instead of failing on a degenerate held out set
#' @keywords internal
safe_auc <- function(y, p) {
  if (length(unique(y)) < 2L) return(NA_real_)
  suppressMessages(as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE))))
}


#' Fraction of entries whose estimated survival increases across thresholds
#'
#' @param S_hat list of T estimated survival matrices
#' @param tol numerical slack
#' @returns a single proportion
#' @details
#' Free for the old method, which multiplies probabilities in [0,1] and so is
#' monotone by construction. The new method has to earn it through clipping,
#' which is precisely why this is worth measuring.
monotonicity_violation <- function(S_hat, tol = 1e-10) {
  Tn <- length(S_hat)
  if (Tn < 2L) return(0)
  bad <- 0
  total <- 0
  for (t in seq.int(2L, Tn)) {
    bad <- bad + sum(S_hat[[t]] > S_hat[[t - 1L]] + tol)
    total <- total + length(S_hat[[t]])
  }
  bad / total
}


#' Calibration slope and intercept on held out entries
#'
#' @param S_hat list of estimated survival matrices
#' @param Y list of observed masks
#' @returns a data frame with one row per threshold
#' @details
#' Regresses the held out outcome on the estimated logit. A slope below one
#' indicates overconfident probabilities, which is the expected failure mode
#' when a rank is fit to a starved risk set.
calibration <- function(S_hat, Y) {
  eval_idx <- seq_along(S_hat[[1L]])
  do.call(rbind, lapply(seq_along(S_hat), function(t) {
    yy <- Y[[t]][eval_idx]
    lp <- stats::qlogis(clamp01(S_hat[[t]][eval_idx]))
    if (length(unique(yy)) < 2L || stats::sd(lp) < 1e-8) {
      return(data.frame(t = t, cal_intercept = NA_real_, cal_slope = NA_real_))
    }
    cf <- stats::coef(stats::glm(yy ~ lp, family = stats::binomial()))
    data.frame(t = t, cal_intercept = unname(cf[1L]), cal_slope = unname(cf[2L]))
  }))
}


# `error_vs_risk_set()` lived here. It stratified RMSE by the size of the
# conditional risk set the old method would have had, on the theory that the two
# engines should separate as that set starves. It was removed because it did not
# work: at small risk sets the true probabilities are near zero, so predicting
# zero everywhere is already near optimal and all methods collapse onto the
# intercept-only baseline (0.066 / 0.070 / 0.065 in the last run). The per
# threshold RMSE table in `eval_survival()` shows the real pattern instead --
# the engines tie at t = 1, where they see identical information, and diverge as
# conditioning compounds.


#' Error of the denoised expectation implied by an estimated survival stack
#'
#' @param S_hat list of estimated survival matrices
#' @param counts observed abundance matrix
#' @param thresholds threshold sequence used to build the masks
#' @param M_true true conditional mean matrix
#' @param cap values above this are left untouched, as in the original dfbm()
#'
#' @returns a one row data frame with the RMSE and correlation of the denoised
#'   expectation against the truth, plus the same for the raw counts.
#'
#' @details
#' Closes the loop back to what the survival probabilities are actually for.
#' Reuses the interval-mean reconstruction from the deleted `dfbm()`.
denoised_expectation_error <- function(S_hat, counts, thresholds, M_true,
                                       cap = NULL) {
  Tn <- length(thresholds)
  if (is.null(cap)) cap <- stats::quantile(counts, 0.99)

  interval_vals <- numeric(Tn)
  for (t in seq_len(Tn - 1L)) {
    sel <- counts > thresholds[t] & counts <= thresholds[t + 1L]
    interval_vals[t] <- if (any(sel)) mean(counts[sel]) else thresholds[t]
  }
  tail_sel <- counts > thresholds[Tn] & counts <= cap
  interval_vals[Tn] <- if (any(tail_sel)) mean(counts[tail_sel]) else thresholds[Tn]

  expected <- matrix(0, nrow(counts), ncol(counts))
  for (t in seq_len(Tn)) {
    s_next <- if (t == Tn) matrix(0, nrow(counts), ncol(counts)) else S_hat[[t + 1L]]
    expected <- expected + interval_vals[t] * (S_hat[[t]] - s_next)
  }
  expected[counts > cap] <- counts[counts > cap]

  data.frame(
    rmse_denoised = sqrt(mean((expected - M_true)^2)),
    rmse_raw = sqrt(mean((counts - M_true)^2)),
    cor_denoised = stats::cor(as.vector(expected), as.vector(M_true)),
    cor_raw = stats::cor(as.vector(counts), as.vector(M_true))
  )
}


#' Collect every metric for one fit into a tidy data frame
#'
#' @param fit an object returned by fit_old(), fit_new() or fit_baseline()
#' @param sim a simulation object
#' @param meta named list of design variables to attach to every row
#' @returns a data frame
summarize_fit <- function(fit, sim, meta = list()) {
  res <- eval_survival(fit$S, sim$S, sim$Y, fit$risk_set)
  cal <- calibration(fit$S, sim$Y)
  res <- merge(res, cal, by = "t", all.x = TRUE)

  res$method <- fit$method
  res$elapsed <- fit$elapsed
  res$mono_violation <- monotonicity_violation(fit$S)
  res$mean_rank <- mean(fit$ranks)
  res$n_iter <- if (is.null(fit$n_iter)) NA_integer_ else fit$n_iter
  res$n_svd <- if (is.null(fit$n_svd)) NA_integer_ else fit$n_svd

  for (nm in names(meta)) res[[nm]] <- meta[[nm]]
  res
}
