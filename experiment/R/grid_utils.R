
# Utilities shared by the driver scripts so that the same script runs both
# locally over the whole grid and on a cluster as one shard of a job array.


#' Keep only the rows of a design grid belonging to this array task
#'
#' @param grid a data frame of design cells
#' @param n_tasks total number of array tasks; read from the environment when NULL
#' @param task_id this task's index, 1 based; read from the environment when NULL
#' @returns the subset of `grid` assigned to this task
#' @details
#' Rows are dealt out round robin rather than in contiguous blocks so that slow
#' cells (large N, many thresholds) spread evenly across tasks instead of
#' piling into one.
slice_for_array <- function(grid, n_tasks = NULL, task_id = NULL) {
  if (is.null(task_id)) {
    raw <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
    task_id <- if (nzchar(raw)) as.integer(raw) else NA_integer_
  }
  if (is.null(n_tasks)) {
    raw <- Sys.getenv("SLURM_ARRAY_TASK_COUNT", "")
    n_tasks <- if (nzchar(raw)) as.integer(raw) else NA_integer_
  }
  if (is.na(task_id) || is.na(n_tasks) || n_tasks <= 1L) return(grid)

  keep <- which((seq_len(nrow(grid)) - 1L) %% n_tasks == (task_id - 1L))
  message(sprintf("array task %d/%d: %d of %d cells",
                  task_id, n_tasks, length(keep), nrow(grid)))
  grid[keep, , drop = FALSE]
}


#' Path for a result file, tagged by config and array task
#'
#' @param stem short name of the result set
#' @param config_name "quick" or "full"
#' @param dir output directory
#' @returns a file path
out_path <- function(stem, config_name, dir = "experiment/results") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  task <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  suffix <- if (nzchar(task)) paste0("_task", task) else ""
  file.path(dir, sprintf("%s_%s%s.rds", stem, config_name, suffix))
}


#' Wall clock time of an expression, returning both value and seconds
#' @keywords internal
timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, elapsed = proc.time()[["elapsed"]] - t0)
}
