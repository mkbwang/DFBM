
# Choosing alpha, and building the tuned arms that get compared against the old
# engine.
#
# Two criteria: the Mallows Cp surrogate and the Brier plus bootstrap covariance
# rule of section 2.3. Both are scored against an `oracle`, the grid alpha that
# minimizes RMSE against the true survival probabilities -- a reference that only
# exists in simulation, which is the whole reason this lives here.
#
# The held out rule that used to be stage 1 is gone: selection no longer splits
# entries, so there is nothing left to compare it on.
#
# Everything runs off ONE call to DFBM:::bmfsvt_path, so no two stages are ever
# compared on different fits, and the package's own selection functions are
# called rather than reimplemented -- otherwise the experiment would validate a
# copy and leave the shipped code untested.


#' Fit the shrinkage path once and score every stage against the truth
#'
#' @param Y list of T nested binary masks
#' @param sim the simulation object, for `S` and `true_ranks`
#' @param alpha_grid decreasing shrinkage grid, as multiples of the noise floor
#'   [DFBM::lambda_star_seq()]
#' @param C pure noise replicates for the noise floor
#' @param do_boot whether to run the bootstrap covariance stage. Off by default:
#'   the bootstrap is set aside for now, it costs 3x to 4x the surrogate, and on
#'   the pilot grid it was the worse of the two on RMSE.
#' @param B bootstrap replicates for stage 3
#' @param window candidate half width for stage 3
#' @param ncores bootstrap replicates to run in parallel
#' @param seed seed for the bootstrap draws
#' @param stability whether to repeat stage 3 anchored one grid step either side
#'   of the surrogate's choice, the check the method note asks for. The
#'   resampling distribution is centred on the fit at the anchor, so the anchor
#'   can flatter itself; running all three shows whether the winner is a
#'   property of the criterion or of where it was started. Costs two extra
#'   bootstraps, so it is off by default.
#' @param ... further arguments passed to [DFBM::bmfsvt()]
#'
#' @returns a list with `per_alpha` (one row per alpha), `stages` (one row per
#'   stage), `ranks` (K x T retained ranks) and `boot`.
fit_selection_path <- function(Y, sim, alpha_grid, C = 20L, do_boot = FALSE,
                               B = 30L, window = 2L, ncores = 1L,
                               seed = NULL, stability = FALSE, ...) {
  t0 <- proc.time()[["elapsed"]]
  path <- DFBM:::bmfsvt_path(Y, alpha_grid = alpha_grid, C = C, seed = seed,
                             ...)
  t_path <- proc.time()[["elapsed"]] - t0
  cat("Original fits finished for all alphas\n")
  K <- length(path$alpha)
  Tn <- path$Tn

  # Truth based quantities, one pass over the path. RMSE is over EVERY entry,
  # which is the estimand: a denoised value is wanted everywhere, not only where
  # the fit was starved of data.
  rmse_all <- numeric(K)
  brier_t <- matrix(NA_real_, K, Tn)
  for (k in seq_len(K)) {
    prob <- DFBM:::path_prob(path, k)
    rmse_all[k] <- sqrt(mean((unlist(prob) - unlist(sim$S))^2))
    for (t in seq_len(Tn)) brier_t[k, t] <- mean((prob[[t]] - Y[[t]])^2)
  }

  # Stage 2, plus the two controls that show why the corrections are needed.
  cp_fixed <- DFBM:::cp_surrogate(path, weighted = TRUE)
  cp_hard <- DFBM:::cp_surrogate(path, weighted = FALSE)
  # The control: variance re-estimated at every alpha AND the hard df count.
  cp_asis <- path$brier + 2 * as.vector(path$df_hard %*% rep(1, Tn)) * path$v_tilde

  # Stage 3, skipped unless asked for.
  boot <- NULL
  t_boot <- 0
  M <- rep(NA_real_, K)
  if (do_boot) {
    t0 <- proc.time()[["elapsed"]]
    boot <- DFBM:::boot_covariance(path, index = cp_fixed$index, B = B,
                                   window = window, ncores = ncores,
                                   seed = seed, ...)
    t_boot <- proc.time()[["elapsed"]] - t0
    cat("Bootstraps finished for all alphas\n")
    M[boot$cand] <- boot$M
  }

  # The method note's stability check. The resample is drawn from the fit at the
  # anchor, so an over-smoothed anchor carries more Bernoulli variance than the
  # truth does; every candidate then chases more noise than it really would, the
  # weakly shrunk ones most of all, and they are over-penalized. Moving the
  # anchor therefore moves the answer, which is what this measures.
  stab <- NULL
  if (stability && do_boot) {
    anchors <- unique(pmax(1L, pmin(K, cp_fixed$index + c(-1L, 0L, 1L))))
    stab <- do.call(rbind, lapply(anchors, function(a) {
      bt <- if (a == cp_fixed$index) boot else
        DFBM:::boot_covariance(path, index = a, B = B, window = window,
                               ncores = ncores, seed = seed, ...)
      data.frame(anchor = a, anchor_alpha = path$alpha[a],
                 k = bt$cand, alpha = path$alpha[bt$cand], M = bt$M,
                 argmin = bt$index)
    }))
  }

  per_alpha <- data.frame(
    k = seq_len(K), alpha = path$alpha,
    rmse_all = rmse_all, brier = path$brier, v_tilde = path$v_tilde,
    df = rowSums(path$df), df_hard = rowSums(path$df_hard),
    cp = cp_fixed$cp, cp_hard = cp_hard$cp, cp_asis = cp_asis, M = M,
    mean_rank = rowMeans(path$ranks), n_iter = path$n_iter,
    n_svd = path$n_svd, converged = path$converged, elapsed = path$elapsed)

  # `elapsed` per stage is the cost of REACHING that alpha, which is what makes
  # the arms built below comparable against fit_old()'s self-timed CV. The Cp
  # surrogate needs only the path; the bootstrap needs the path plus its own
  # replicates, because the surrogate is what sets its anchor.
  picks <- list(
    oracle        = which.min(rmse_all),
    stage_2_cp    = cp_fixed$index,
    cp_note_asis  = which.min(cp_asis),
    cp_hard_df    = cp_hard$index)
  timing <- c(oracle = t_path, stage_2_cp = t_path,
              cp_note_asis = t_path, cp_hard_df = t_path)
  if (do_boot) {
    picks$stage_3_boot <- boot$index
    timing["stage_3_boot"] <- t_path + t_boot
  }

  k_or <- picks$oracle
  stages <- do.call(rbind, lapply(names(picks), function(nm) {
    k <- picks[[nm]]
    data.frame(stage = nm, k = k, alpha = path$alpha[k],
               k_minus_oracle = k - k_or,
               abs_k_minus_oracle = abs(k - k_or),
               rmse_all = rmse_all[k],
               excess_rmse = rmse_all[k] - rmse_all[k_or],
               mean_rank = mean(path$ranks[k, ]),
               mean_true_rank = mean(sim$true_ranks),
               rank_err = mean(path$ranks[k, ] - sim$true_ranks),
               abs_rank_err = mean(abs(path$ranks[k, ] - sim$true_ranks)),
               elapsed = unname(timing[nm]))
  }))

  # Per threshold rank recovery, long form, so the report can show which
  # thresholds a stage gets right rather than only the average.
  rank_long <- do.call(rbind, lapply(names(picks), function(nm) {
    k <- picks[[nm]]
    data.frame(stage = nm, t = seq_len(Tn), rank_fit = path$ranks[k, ],
               rank_true = sim$true_ranks)
  }))

  list(per_alpha = per_alpha, stages = stages, rank_long = rank_long,
       stability = stab,
       ranks = path$ranks, boot = boot, picks = picks, path = path,
       alpha_cp = path$alpha[cp_fixed$index],
       alpha_boot = if (do_boot) path$alpha[boot$index] else NA_real_,
       v_ref = cp_fixed$v_ref,
       t_path = t_path, t_boot = t_boot)
}


#' Turn one point of a fitted path into an arm the metrics code understands
#'
#' @param sel the object returned by [fit_selection_path()]
#' @param which `"cp"` or `"boot"`
#' @returns a list shaped like [fit_old()] / [fit_new()]
#' @details
#' Built from the path that was already fit rather than refitting, so `elapsed`
#' is exactly the cost of selecting that alpha and nothing is paid twice.
#'
#' This is the fix for a real defect: the tuned arm used to be constructed with
#' `fit_new(alpha = <already chosen>)`, whose timer wrapped only the final fit,
#' while `fit_old()` times its whole rank-and-lambda search. Comparing the two
#' flattered the new method by the entire cost of selection.
selection_arm <- function(sel, which = c("cp", "boot")) {
  which <- match.arg(which)
  k <- if (which == "cp") sel$picks$stage_2_cp else sel$picks$stage_3_boot
  path <- sel$path
  elapsed <- if (which == "cp") sel$t_path else sel$t_path + sel$t_boot
  nobs <- path$N * path$P
  list(S = DFBM:::path_prob(path, k),
       pi_cond = NULL,
       ranks = path$ranks[k, ],
       risk_set = rep(nobs, path$Tn),
       elapsed = elapsed,
       method = paste0("new_", which),
       alpha = path$alpha[k],
       n_iter = path$n_iter[k], n_svd = path$n_svd[k])
}


#' Attach design metadata to every row of a selection result
#'
#' @param sel the object returned by [fit_selection_path()]
#' @param meta named list of design variables
#' @returns a list of three data frames with `meta` bound on
tag_selection <- function(sel, meta) {
  bind <- function(df) {
    for (nm in names(meta)) df[[nm]] <- meta[[nm]]
    df
  }
  list(stages = bind(sel$stages),
       per_alpha = bind(sel$per_alpha),
       rank_long = bind(sel$rank_long),
       stability = if (is.null(sel$stability)) NULL else bind(sel$stability))
}
