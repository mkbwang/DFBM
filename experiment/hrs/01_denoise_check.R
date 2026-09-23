#!/usr/bin/env Rscript
#
# 01_denoise_check.R -- does the denoising pipeline behave on the real HRS data?
#
# A data-only check: no outcomes are needed. One 70/30 split of the merged leaf
# composition (8874 x 25). The denoiser is fit on the training rows only and
# applied to the test rows with predict(), exactly as the prediction study will
# do it. Checks, in order:
#
#   1. the fit: alpha, ranks, convergence, time, fold-in vs in-sample gap;
#   2. calibration: mean fitted P(> d_t) against mask prevalence, per column;
#   3. generalization: test Brier of the fold-in against the column offsets;
#   4. zero recovery against an EXTERNAL truth: pbas and peos come from a CBC
#      count rounded to 0.1 x10^9/L, so an observed zero means a count below
#      0.05. Denoised composition x row total should fall in (0, 0.05];
#   5. bounds and closure: positive, capped, rows sum to 1, CLR finite;
#   6. flags: are the known extreme rows (likely hematologic malignancy)
#      flagged, and what does denoising do to them.
#
# Usage, from the package root:
#   experiment/run.sh hrs/01_denoise_check.R          # fixed alpha = 0.5
#   experiment/run.sh hrs/01_denoise_check.R tuned    # adds the Cp-tuned fit
#
# The fixed-alpha run is about a minute; the tuned run fits a 10-point path and
# takes several minutes.

suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
run_tuned <- "tuned" %in% args
out_dir <- file.path("experiment", "hrs", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

comp <- as.matrix(read.csv("inst/extdata/flocyt_proj_composition_leaf_merged.csv"))
counts <- as.matrix(read.csv("inst/extdata/flocyt_proj_counts_leaf_merged.csv"))
colnames(comp) <- sub("_adj_pct$", "", colnames(comp))
total <- rowSums(counts)
n <- nrow(comp)
rownames(comp) <- seq_len(n)

known_extreme <- c(2830L, 8010L, 7147L, 388L, 3714L, 627L)

set.seed(2026)
train <- sort(sample(n, round(0.7 * n)))
test <- setdiff(seq_len(n), train)
cat(sprintf("rows: %d train, %d test; columns: %d\n\n",
            length(train), length(test), ncol(comp)))

rule <- function(title) cat("\n==", title, strrep("=", max(0, 70 - nchar(title))), "\n")

check <- function(label, ...) {
  rule(label)
  t0 <- proc.time()[["elapsed"]]
  obj <- dfbm(comp[train, ], seed = 1, ...)
  fit_secs <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  pt <- predict(obj, comp[test, ])
  pred_secs <- proc.time()[["elapsed"]] - t0

  # ---- 1. the fit ----------------------------------------------------------
  cat(sprintf("alpha %.4g | fit %.1f s | predict %.1f s | converged %s in %d iters\n",
              obj$alpha, fit_secs, pred_secs, obj$fit$converged, obj$fit$n_iter))
  cat("ranks:", obj$fit$ranks, "\n")
  cat(sprintf("fold-in vs in-sample |dp|: mean %.2e, max %.3f\n",
              obj$insample_gap[["mean"]], obj$insample_gap[["max"]]))
  cat(sprintf("fold-in converged: train %d/%d, test %d/%d\n",
              sum(obj$converged), length(train), sum(pt$converged), length(test)))
  if (!is.null(obj$tune_path)) print(obj$tune_path[, c("alpha", "brier", "df", "cp", "mean_rank")])

  # ---- 2. calibration ------------------------------------------------------
  calib <- function(prob, A) {
    Y <- make_masks(A, obj$thresholds)
    err <- vapply(seq_along(Y), function(t) colMeans(prob[[t]]) - colMeans(Y[[t]]),
                  numeric(ncol(A)))
    setNames(apply(abs(err), 1, max), colnames(A))
  }
  cal <- rbind(train = calib(obj$prob, comp[train, ]),
               test = calib(pt$prob, comp[test, ]))
  cat("\nmax over t of |mean fitted P(> d_t) - prevalence|, worst 6 columns (train) and pbas:\n")
  worst <- order(cal["train", ], decreasing = TRUE)[1:6]
  print(round(cal[, unique(c(worst, which(colnames(cal) == "pbas")))], 3))

  # ---- 3. generalization ---------------------------------------------------
  Yte <- make_masks(comp[test, ], obj$thresholds)
  brier <- function(prob) {
    mean(Reduce(`+`, lapply(seq_along(Yte), function(t) (prob[[t]] - Yte[[t]])^2)))
  }
  base_prob <- lapply(seq_along(Yte), function(t) {
    mu_t <- if (is.matrix(obj$fit$mu)) obj$fit$mu[t, ] else obj$fit$mu[t]
    matrix(rep(plogis(mu_t), each = length(test)), length(test), ncol(comp))
  })
  cat(sprintf("\ntest Brier per entry (summed over t): fold-in %.4f, offsets only %.4f (%.1f%% lower)\n",
              brier(pt$prob), brier(base_prob),
              100 * (1 - brier(pt$prob) / brier(base_prob))))

  # ---- 4. zero recovery against the rounded CBC counts ---------------------
  cat("\nimplied count (denoised composition x row total, x10^9/L), test rows:\n")
  for (col in c("pbas", "peos")) {
    obs <- counts[test, paste0(col, "_adj_count")]
    implied <- pt$denoised[, col] * total[test]
    grp <- cut(obs, c(-Inf, 0.025, 0.15, 0.25, Inf),
               labels = c("obs 0", "obs 0.1", "obs 0.2", "obs >= 0.3"))
    tab <- t(vapply(split(implied, grp), function(x) {
      c(n = length(x), quantile(x, c(0.1, 0.5, 0.9)))
    }, numeric(4)))
    cat(col, ":\n")
    print(signif(tab, 3))
    z <- obs == 0
    cat(sprintf("  observed zeros: %d; implied count in (0, 0.05]: %.1f%%, in (0, 0.10]: %.1f%%\n",
                sum(z), 100 * mean(implied[z] > 0 & implied[z] <= 0.05),
                100 * mean(implied[z] > 0 & implied[z] <= 0.10)))
    # Why: the fitted chance of clearing the first threshold, zero vs nonzero rows.
    s1 <- pt$prob[[1L]][, which(colnames(comp) == col)]
    cat(sprintf("  first threshold %.3g; mean fitted P(> d_1): zeros %.3f, nonzeros %.3f; train zero fraction %.3f\n",
                obj$thresholds[1L, col], mean(s1[z]), mean(s1[!z]),
                obj$zero_frac[[col]]))
  }

  # ---- 5. bounds and closure -----------------------------------------------
  unclosed <- pt$denoised * pt$row_sum
  cat(sprintf("\nall > 0: %s | unclosed <= cap: %s | max |row sum - 1|: %.1e | CLR finite: %s\n",
              all(pt$denoised > 0),
              all(unclosed <= rep(obj$cap, each = length(test)) * (1 + 1e-12)),
              max(abs(rowSums(pt$denoised) - 1)),
              all(is.finite(log(pt$denoised)))))
  cat("pre-closure row sums (test):", signif(quantile(pt$row_sum, c(0, .01, .5, .99, 1)), 4), "\n")
  cat("binned row sums (test):     ", signif(quantile(pt$binned_row_sum, c(0, .01, .5, .99, 1)), 4), "\n")

  # ---- 6. flags and the known extreme rows ---------------------------------
  cat(sprintf("\nrow_ce cutoff %.4f; flagged: train %d, test %d\n",
              obj$row_ce_cut, sum(obj$flag), sum(pt$flag)))
  cat("flagged test rows:", test[pt$flag], "\n")
  ce_all <- c(obj$row_ce, pt$row_ce)
  names(ce_all) <- c(train, test)
  bs_all <- c(obj$binned_row_sum, pt$binned_row_sum)
  names(bs_all) <- c(train, test)
  fl_all <- c(obj$flag, pt$flag)
  names(fl_all) <- c(train, test)
  den_all <- rbind(obj$denoised, pt$denoised)
  rownames(den_all) <- c(train, test)
  ex <- data.frame(row = known_extreme,
                   set = ifelse(known_extreme %in% train, "train", "test"),
                   row_ce = ce_all[as.character(known_extreme)],
                   ce_pctile = vapply(known_extreme, function(r) {
                     mean(ce_all <= ce_all[as.character(r)])
                   }, numeric(1)),
                   binned_row_sum = bs_all[as.character(known_extreme)],
                   flag = fl_all[as.character(known_extreme)])
  print(ex, row.names = FALSE, digits = 4)
  cat("\nlargest raw component of each known extreme row, raw vs denoised:\n")
  for (r in known_extreme) {
    j <- which.max(comp[r, ])
    j2 <- which.max(abs(comp[r, ] - den_all[as.character(r), ]))
    cat(sprintf("  row %d: %s raw %.4f -> %.4f | largest change %s raw %.4f -> %.4f\n",
                r, colnames(comp)[j], comp[r, j], den_all[as.character(r), j],
                colnames(comp)[j2], comp[r, j2], den_all[as.character(r), j2]))
  }

  summary <- list(label = label, alpha = obj$alpha, ranks = obj$fit$ranks,
                  fit_secs = fit_secs, pred_secs = pred_secs,
                  insample_gap = obj$insample_gap, calibration = cal,
                  brier = c(foldin = brier(pt$prob), offsets = brier(base_prob)),
                  extreme = ex, tune_path = obj$tune_path,
                  flagged_test = test[pt$flag])
  saveRDS(summary, file.path(out_dir, paste0("denoise_check_", label, ".rds")))
  invisible(summary)
}

check("fixed_alpha_0.5", alpha = 0.5)
if (run_tuned) check("tuned")
