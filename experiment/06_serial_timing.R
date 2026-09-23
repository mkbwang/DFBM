#!/usr/bin/env Rscript
#
# 06_serial_timing.R -- how much arithmetic does each engine actually do?
#
# RUN WITH experiment/run.sh, like every other stage.
#
# Stage 01 times both engines at the same `ncores`, which is the wall clock a
# user sees. That is the right headline number, but it conflates two things: how
# much work an algorithm does, and how well that work parallelizes.
#
# This stage removes that confound by running both at `ncores = 1`. Kept small:
# timing only, no accuracy, one rep per setting, because the point is a ratio
# and not a distribution.
#
# What it established, and why the stage 01 comparison is fair rather than
# generous: the OLD engine gains nothing from cores (measured 1.00x on 300x20
# and 0.78x -- i.e. slower -- on 150x100), because `cv.logisticcfR` can only
# spread `length(lambdas) = 3` fits per rank level and the fork does not repay
# itself. And `new_cp` never touches `ncores` at all: it reaches only
# `boot_covariance()`, and `bmfsvt_path()` is serial. So this stage's `speedup`
# column for `new_cp` is pure run-to-run noise (about 1.4x), which is the
# precision of any single-shot timing here -- prefer stage 01's replicate means.

source("experiment/R/config.R")
source("experiment/R/grid_utils.R")
root <- setup_experiment(getwd())
cfg <- exp_config()
message("config: ", cfg$name)

acc <- cfg$accuracy
grid <- expand.grid(shape = names(acc$shapes), ncores = c(1L, cfg$ncores),
                    stringsAsFactors = FALSE)

rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  dims <- acc$shapes[[g$shape]]
  Tn <- acc$Tn[1L]
  message(sprintf("[%d/%d] %s T=%d ncores=%d", i, nrow(grid), g$shape, Tn,
                  g$ncores))
  sim <- simulate_nested_masks(
    N = dims[["N"]], P = dims[["P"]], Tn = Tn, rank = acc$rank[1L],
    pi_head = cfg$pi_head, pi_tail = acc$pi_tail[1L], seed = 4242L)

  t0 <- proc.time()[["elapsed"]]
  fo <- try(fit_old(sim$Y, max_K = cfg$old_max_K,
                    lambdas = cfg$old_lambdas, ncores = g$ncores),
            silent = TRUE)
  t_old <- if (inherits(fo, "try-error")) NA_real_ else fo$elapsed

  sel <- try(fit_selection_path(sim$Y, sim,
                               alpha_grid = cfg$selection$alpha_grid,
                               C = cfg$selection$C, do_boot = FALSE,
                               stability = FALSE, ncores = g$ncores,
                               seed = 4242L, max_iter = cfg$max_iter),
             silent = TRUE)
  ok <- !inherits(sel, "try-error")

  rows[[i]] <- data.frame(
    shape = g$shape, N = dims[["N"]], P = dims[["P"]], Tn = Tn,
    ncores = g$ncores,
    old = t_old,
    new_cp = if (ok) sel$t_path else NA_real_,
    config = cfg$name)
  print(rows[[i]], row.names = FALSE, digits = 4)
}

out <- do.call(rbind, rows)
saveRDS(out, out_path("serial_timing", cfg$name))
message("done: ", out_path("serial_timing", cfg$name))
