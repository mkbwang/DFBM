#!/usr/bin/env Rscript
#
# 01_accuracy.R -- how accurately does each method recover P(value > d_t)?
#
# The primary design axis is `pi_tail`, the prevalence of the last mask. The old
# method conditions on survival to the previous threshold, so its effective
# sample size at threshold t is proportional to the prevalence at t-1; the new
# method uses every entry at every threshold. Driving pi_tail down is therefore
# the direct test of whether that structural difference matters.
#
# RUN THIS WITH experiment/run.sh, NOT bare Rscript, and not line by line in
# RStudio unless RhpcBLASctl is installed. Forked workers each run their own
# BLAS; with a threaded BLAS that is workers x cores threads and the whole
# machine thrashes. Measured 7.9x. Setting OPENBLAS_NUM_THREADS in .Renviron or
# with Sys.setenv() does NOT work -- OpenBLAS reads it before R starts. See
# check_blas_threads() in R/config.R.
#
# `alpha` is a multiple of the NOISE FLOOR, DFBM::lambda_star_seq(), not of
# lambda_max_seq(). The measured optimum is near 0.55. An alpha carried over
# from the old parameterization will be wrong, not merely suboptimal.
#
# Two questions are answered here: whether the new engine beats the old one on
# accuracy and on speed. The tuned arm is built from the fitted path
# (R/selection.R), so its `elapsed` includes the cost of selecting alpha --
# which is what makes it comparable against fit_old(), whose timer wraps its own
# rank-and-lambda CV.
#
# The bootstrap arm is off (cfg$selection$do_boot). Stages 02, 03 and 05 were
# retired; the suite is now 01_accuracy.R -> 04_render.R, with 06_serial_timing.R
# as an optional check. All entries are observed: there is no Omega anywhere.
#
# The count based end to end block that used to live at the foot of this script
# was removed. `simulate_nested_masks()` targets the survival probabilities
# directly, and those are the estimand; `simulate_from_counts()` and
# `denoised_expectation_error()` remain in experiment/R/ but nothing calls them.
#
# Usage:
#   EXP_CONFIG=quick experiment/run.sh 01_accuracy.R
#   EXP_CONFIG=full  experiment/run.sh 01_accuracy.R      # via slurm array

source("experiment/R/config.R")
source("experiment/R/grid_utils.R")
root <- setup_experiment(getwd())
cfg <- exp_config()
message("config: ", cfg$name)

acc <- cfg$accuracy
grid <- expand.grid(
  shape = names(acc$shapes),
  Tn = acc$Tn,
  rank = acc$rank,
  pi_tail = acc$pi_tail,
  rep = seq_len(acc$reps),
  stringsAsFactors = FALSE
)
grid <- slice_for_array(grid)
message(sprintf("%d design cells to run", nrow(grid)))

run_cell <- function(i) {
  g <- grid[i, ]
  dims <- acc$shapes[[g$shape]]
  seed <- 10000L * g$rep + as.integer(g$Tn) * 97L +
    round(1000 * g$pi_tail) + as.integer(g$rank)

  sim <- simulate_nested_masks(
    N = dims[["N"]], P = dims[["P"]], Tn = g$Tn, rank = g$rank,
    pi_head = cfg$pi_head, pi_tail = g$pi_tail, seed = seed)

  # No split: alpha selection scores every entry, and the headline metric is
  # squared error against the true S, which needs no held out set.
  meta <- list(shape = g$shape, N = dims[["N"]], P = dims[["P"]],
               Tn = g$Tn, rank = g$rank, pi_tail = g$pi_tail,
               rep = g$rep, config = cfg$name)

  fits <- list()
  fits$baseline <- fit_baseline(sim$Y)
  # Same core budget as the new method. Note the old engine gains nothing from
  # cores at these shapes (measured 1.00x and 0.78x), so this is fair rather
  # than generous -- see 06_serial_timing.R.
  fits$old <- try(fit_old(sim$Y, max_K = cfg$old_max_K,
                          lambdas = cfg$old_lambdas, ncores = cfg$ncores),
                  silent = TRUE)

  # One path serves the tuned arm and the stage table.
  sel <- try(fit_selection_path(sim$Y, sim,
                                alpha_grid = cfg$selection$alpha_grid,
                                C = cfg$selection$C,
                                do_boot = cfg$selection$do_boot,
                                B = cfg$selection$B,
                                window = cfg$selection$window,
                                stability = cfg$selection$stability,
                                ncores = cfg$ncores, seed = seed,
                                max_iter = cfg$max_iter), silent = TRUE)
  if (inherits(sel, "try-error")) {
    message(sprintf("  cell %d: selection failed -- %s", i,
                    conditionMessage(attr(sel, "condition"))))
    sel <- NULL
  } else {
    fits$new_cp <- selection_arm(sel, "cp")
    if (isTRUE(cfg$selection$do_boot)) fits$new_boot <- selection_arm(sel, "boot")
  }

  out <- list()
  for (nm in names(fits)) {
    f <- fits[[nm]]
    if (inherits(f, "try-error")) {
      message(sprintf("  cell %d: %s failed -- %s", i, nm, conditionMessage(attr(f, "condition"))))
      next
    }
    s <- summarize_fit(f, sim, meta = meta)
    s$fit_name <- nm
    s$alpha <- if (is.null(f$alpha)) NA_real_ else f$alpha
    out[[nm]] <- s
  }

  c(list(summary = do.call(rbind, out)),
    if (is.null(sel)) NULL else tag_selection(sel, meta))
}

t_start <- proc.time()[["elapsed"]]
results <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  message(sprintf("[%d/%d] %s T=%d rank=%d pi_tail=%.2f rep=%d",
                  i, nrow(grid), grid$shape[i], grid$Tn[i], grid$rank[i],
                  grid$pi_tail[i], grid$rep[i]))
  results[[i]] <- run_cell(i)
}
message(sprintf("total elapsed: %.1f min", (proc.time()[["elapsed"]] - t_start) / 60))

gather <- function(key) {
  parts <- lapply(results, `[[`, key)
  do.call(rbind, parts[!vapply(parts, is.null, logical(1L))])
}

saveRDS(list(summary = gather("summary"),
             selection = gather("stages"),
             sel_path = gather("per_alpha"),
             sel_ranks = gather("rank_long"),
             sel_stability = gather("stability"),
             config = cfg),
        out_path("accuracy", cfg$name))
message("done: ", out_path("accuracy", cfg$name))


#### load the result and plot Brier scores for comparison

library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 11))


method_labs <- c(baseline = "Intercept only", old = "Old",
                 new_cp = "New")
pal <- c("Intercept only" = "grey55", "Old" = "#C2453E",
         "New" = "#3A6FB0")
lab <- function(x) factor(method_labs[x], levels = names(pal))
acc <- readRDS("experiment/results/accuracy_quick.rds")

# overall accuracy
summary_df = acc$summary
summary_df <- summary_df %>% mutate(shape = case_match(shape,
                                                     "hrs" ~ "300*20",
                                                     "square" ~ "150*100",
                                                     .default = shape  # keeps other values unchanged
                                    ))



summary_df %>%
  filter(!is.na(t)) %>%
  mutate(method = lab(fit_name)) %>%
  group_by(method, t, shape) %>%
  summarise(sq_err = mean(rmse^2), .groups = "drop") %>%
  ggplot(aes(t, sq_err, colour = method)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  facet_wrap(~ shape, scales = "free_y") +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = 0:9)+
  labs(x = "Threshold",
       y = expression(mean((hat(pi) - pi[true])^2)), colour = NULL) +
  theme(legend.position = "bottom")




