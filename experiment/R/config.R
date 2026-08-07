
# One place to switch between a grid that runs on a laptop in minutes and the
# grid that goes to the cluster. Every driver script reads EXP_CONFIG from the
# environment, defaulting to "quick".

#' How many bootstrap replicates to fork, and a check on BLAS threading
#'
#' @returns an integer worker count
#' @details
#' Two settings interact and getting them wrong costs almost an order of
#' magnitude. Each forked worker runs its own BLAS, and a threaded BLAS defaults
#' to every core, so `n` workers ask for `n x cores` threads. Measured on one
#' bootstrap (150 x 100, T = 5, B = 8), seconds:
#'
#' \tabular{lrrrr}{
#'   BLAS threads \tab 1 worker \tab 4 \tab 8 \tab 14 \cr
#'   default (16) \tab 106.3 \tab 94.6 \tab 94.7 \tab 86.4 \cr
#'   1            \tab  40.7 \tab 15.3 \tab 10.9 \tab 11.1
#' }
#'
#' Single threaded BLAS with 8 workers is **7.9x** faster than 14 workers with a
#' threaded BLAS, and leaves half the machine free. Threading the BLAS hurts
#' even serially (106 vs 41), because the per block SVDs are small enough that
#' coordination costs more than it saves.
#'
#' Hence the default is capped at 8 rather than `detectCores() - 2`: past that
#' the curve is flat and the extra forks are pure memory. Override with
#' `EXP_NCORES`.
exp_ncores <- function() {
  env <- Sys.getenv("EXP_NCORES", "")
  if (nzchar(env)) return(max(1L, as.integer(env)))
  max(1L, min(8L, parallel::detectCores() - 2L))
}


#' Stop the BLAS from oversubscribing the machine once workers are forked
#'
#' @details
#' OpenBLAS fixes its thread count in its library constructor, which the dynamic
#' linker runs *before* R evaluates any startup file. Measured, on 3
#' multiplications of a 1500 x 1500 matrix:
#'
#' \tabular{lrl}{
#'   how the variable is set \tab time \tab single threaded? \cr
#'   not set \tab 0.26 s \tab no \cr
#'   exported before R starts \tab 1.38 s \tab yes \cr
#'   `Sys.setenv()` in-session \tab 0.26 s \tab no \cr
#'   `.Renviron` \tab 0.26 s \tab no \cr
#'   `RhpcBLASctl::blas_set_num_threads(1)` \tab 1.38 s \tab yes
#' }
#'
#' So neither `.Renviron` nor `Sys.setenv()` works -- only the environment at
#' process launch, or OpenBLAS's own runtime entry point via `RhpcBLASctl`,
#' which was checked to reach the same timing and to propagate into forked
#' workers. That is what makes an interactive RStudio session safe, since it
#' cannot be relaunched with the variable set.
check_blas_threads <- function() {
  if (any(nzchar(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"))))) {
    return(invisible(TRUE))
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    message("NOTE: BLAS threading was unset; RhpcBLASctl has pinned it to 1 ",
            "thread for this session.")
    return(invisible(TRUE))
  }
  message("NOTE: BLAS threading is unset, so every forked worker will try to ",
          "use every core (measured 7.9x on the bootstrap stage).\n",
          "      Either use experiment/run.sh, or install RhpcBLASctl so this ",
          "can be fixed in-session.\n",
          "      Setting it in .Renviron or with Sys.setenv() does NOT work: ",
          "OpenBLAS reads it before R starts.")
  invisible(FALSE)
}


exp_config <- function(which = Sys.getenv("EXP_CONFIG", "quick")) {
  which <- match.arg(which, c("quick", "full"))

  if (which == "quick") {
    list(
      name = "quick",
      # ---- 01_accuracy ----
      # Two settings, replicated, rather than a wide sweep. The binding cost
      # is cv.logisticcfR, which refits the old method once per rank per lambda
      # per threshold.
      accuracy = list(
        # Two design points only, chosen deliberately: the HRS shape and a
        # square one, at the sparsest tail, where the old method's conditional
        # risk set is most starved and the engines should separate. Replication
        # rather than breadth, because the seed is a function of the design so
        # every engine sees the SAME data within a rep -- the comparison is
        # paired, and 5 paired reps beat a wide unreplicated sweep.
        shapes = list(hrs = c(N = 300, P = 20), square = c(N = 150, P = 100)),
        Tn = 10L,
        rank = 3L,
        pi_tail = 0.02,
        reps = 5L
      ),
      # ---- alpha selection (inside 01_accuracy) ----
      # `alpha` is a multiple of the NOISE FLOOR `lambda_star_seq()`, not of
      # `lambda_max_seq()`. The measured optimum is near 0.55 (argmin in 6 of 6
      # pilot cells, two shapes x three seeds), so the grid is log spaced over
      # [0.2, 1] -- a step ratio of 1.20, bracketing 0.55 with headroom on the
      # over-shrinking side, which is the side the Cp criterion errs toward.
      # `alpha = 1` means "threshold exactly at the noise floor".
      #
      # `boot` is off: the bootstrap is set aside for now and the accuracy
      # comparison runs baseline / old / new_cp. `B`, `window` and `stability`
      # are kept so turning it back on needs no edit here.
      selection = list(
        alpha_grid = exp(seq(log(1), log(0.2), length.out = 10L)),
        C = 20L,
        do_boot = FALSE,
        B = 10L,
        window = 2L,
        stability = FALSE
      ),
      pi_head = 0.75,
      max_iter = 200L,
      # The full config's tuning, not a cheaper one. A speed or accuracy win
      # against a deliberately under-tuned baseline is not defensible, and this
      # is the setting that gives the old engine its best shot.
      old_max_K = 10L,
      old_lambdas = c(0.01, 0.1, 1),
      ncores = exp_ncores()
    )
  } else {
    list(
      name = "full",
      accuracy = list(
        shapes = list(hrs = c(N = 2000, P = 20),
                      square = c(N = 500, P = 500),
                      proteomics = c(N = 1000, P = 1000)),
        Tn = c(5L, 10L, 20L),
        rank = c(3L, 6L),
        pi_tail = c(0.20, 0.05, 0.01),
        reps = 20L
      ),
      selection = list(
        alpha_grid = exp(seq(log(1), log(0.2), length.out = 10L)),
        C = 20L,
        do_boot = FALSE,
        B = 30L,
        window = 2L,
        stability = TRUE
      ),
      pi_head = 0.75,
      max_iter = 400L,
      old_max_K = 10L,
      old_lambdas = c(0.01, 0.1, 1),
      ncores = exp_ncores()
    )
  }
}


#' Load the package and every experiment helper
setup_experiment <- function(root = NULL) {
  if (is.null(root)) {
    root <- tryCatch(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
                     error = function(e) getwd())
  }
  suppressMessages(devtools::load_all(root, quiet = TRUE))
  for (f in c("simulate.R", "old_method.R", "metrics.R", "selection.R")) {
    source(file.path(root, "experiment", "R", f))
  }
  check_blas_threads()
  invisible(root)
}
