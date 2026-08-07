
# Generators for nested binary mask stacks with known ground truth.
#
# Two properties have to hold simultaneously and neither is automatic:
#   1. every X^(t) is (approximately) low rank, so that a factorization method
#      has something to find;
#   2. the entrywise probabilities decrease across thresholds, so that the
#      masks are nested and the H^(t) <= 0 parameterization is well posed.
#
# The obvious construction X^(t) = Theta - delta_t satisfies both but is
# degenerate: it makes Z^(t) identically zero for t >= 1, so any method that
# shrinks hard looks perfect. We instead give each step its own entrywise
# nonnegative low rank gap, and expose `gap_hetero` to interpolate between the
# degenerate case (0) and fully heterogeneous gaps (1).


#' Numerical rank of each block in the model's own parameterization
#'
#' @param X list of T true natural parameter matrices
#' @param k_max largest rank that will be looked for; the generator's blocks are
#'   low rank by construction, so a partial decomposition suffices and a full
#'   one would cost T SVDs of the whole matrix per simulated cell
#' @param tol singular values below this fraction of the largest are treated as
#'   zero
#' @returns an integer vector of length T
#' @details
#' `bmfsvt()` penalizes the increments \eqn{Z^{(t)}}, not the levels
#' \eqn{X^{(t)}}, so "the true rank of mask t" only means something once the
#' same reparameterization is applied to the truth: the first block is the
#' centered signal and the rest are successive differences. For the default
#' generator this comes out as `rank + 1` throughout -- the low rank gap plus
#' the rank one offset that the model absorbs into \eqn{\mu_0} and \eqn{\nu_t}.
#' @keywords internal
block_ranks <- function(X, k_max = 12L, tol = 1e-8) {
  blocks <- c(list(X[[1L]] - mean(X[[1L]])),
              lapply(seq_along(X)[-1L], function(t) X[[t]] - X[[t - 1L]]))
  mindim <- min(dim(blocks[[1L]]))
  k <- min(as.integer(k_max), mindim)
  vapply(blocks, function(M) {
    d <- if (k < mindim) {
      tryCatch(RSpectra::svds(M, k = k, nu = 0L, nv = 0L)$d,
               error = function(e) svd(M, nu = 0L, nv = 0L)$d)
    } else {
      svd(M, nu = 0L, nv = 0L)$d
    }
    sum(d > tol * max(d))
  }, integer(1L))
}


#' Entrywise nonnegative low rank matrix with unit mean
#' @keywords internal
nonneg_lowrank <- function(N, P, rank) {
  U <- abs(matrix(rnorm(N * rank), N, rank))
  V <- abs(matrix(rnorm(P * rank), P, rank))
  G <- tcrossprod(U, V)
  G / mean(G)
}


#' Simulate a nested stack of binary masks with known probabilities
#'
#' @param N,P matrix dimensions
#' @param Tn number of thresholds
#' @param rank rank of the latent signal
#' @param pi_head target prevalence of the first mask
#' @param pi_tail target prevalence of the last mask
#' @param signal_sd standard deviation of the latent signal on the logit scale
#' @param gap_hetero in [0, 1]; 0 gives identical gaps for every entry (the
#'   degenerate, maximally constraint-binding case), 1 gives fully
#'   heterogeneous low rank gaps
#' @param seed optional RNG seed
#'
#' @returns a list with the mask list `Y`, the true natural parameters `X`, the
#'   true survival probabilities `S`,
#'   the shared uniforms `U`, and the realized prevalences.
#'
#' @details
#' The single most important detail is that all T Bernoulli draws are coupled
#' through **one** uniform per entry. Drawing independently for each threshold
#' would break nesting and make the conditional likelihood used by the old
#' method undefined; coupling is also exactly what happens when a real
#' abundance matrix is thresholded.
simulate_nested_masks <- function(N, P, Tn, rank = 3L,
                                  pi_head = 0.7, pi_tail = 0.05,
                                  signal_sd = 1.5, gap_hetero = 1,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(Tn >= 1L, pi_tail <= pi_head, gap_hetero >= 0, gap_hetero <= 1)

  # Latent low rank signal, centered and scaled to a controlled logit spread.
  A <- matrix(rnorm(N * rank), N, rank)
  B <- matrix(rnorm(P * rank), P, rank)
  Theta <- tcrossprod(A, B) / sqrt(rank)
  Theta <- (Theta - mean(Theta)) / stats::sd(as.vector(Theta)) * signal_sd

  # Solve for the intercept so the head prevalence is exactly on target rather
  # than approximately, which keeps the design axis interpretable.
  mu0 <- stats::uniroot(
    function(m) mean(stats::plogis(m + Theta)) - pi_head,
    interval = c(-30, 30))$root
  X <- vector("list", Tn)
  X[[1L]] <- mu0 + Theta

  if (Tn > 1L) {
    # Cumulative nonnegative gaps, each step mixing a constant with a
    # heterogeneous low rank component.
    cum_gap <- vector("list", Tn - 1L)
    running <- matrix(0, N, P)
    for (t in seq_len(Tn - 1L)) {
      G <- (1 - gap_hetero) + gap_hetero * nonneg_lowrank(N, P, rank)
      running <- running + G
      cum_gap[[t]] <- running
    }
    # One global scale so that the tail prevalence lands on target.
    scale_c <- stats::uniroot(
      function(cc) mean(stats::plogis(X[[1L]] - cc * cum_gap[[Tn - 1L]])) - pi_tail,
      interval = c(1e-8, 200))$root
    for (t in seq_len(Tn - 1L)) X[[t + 1L]] <- X[[1L]] - scale_c * cum_gap[[t]]
  }

  S <- lapply(X, stats::plogis)

  # One uniform per entry, shared by every threshold: this is what makes the
  # masks nested.
  U <- matrix(runif(N * P), N, P)
  Y <- lapply(S, function(p) 1 * (U < p))

  list(Y = Y, X = X, S = S, U = U,
       true_ranks = block_ranks(X),
       prevalence = vapply(Y, mean, numeric(1L)),
       mu0 = mu0, params = list(N = N, P = P, Tn = Tn, rank = rank,
                                pi_head = pi_head, pi_tail = pi_tail,
                                signal_sd = signal_sd,
                                gap_hetero = gap_hetero))
}


#' Simulate an abundance matrix and threshold it the way dfbm does
#'
#' @param N,P matrix dimensions
#' @param Tn number of thresholds
#' @param rank rank of the latent log mean
#' @param size negative binomial dispersion
#' @param zi_prob probability that an entry is technically zeroed out
#' @param seed optional RNG seed
#'
#' @returns a list adding `counts`, `thresholds` and the true conditional mean
#'   `M` to the fields returned by [simulate_nested_masks()].
#'
#' @details
#' Ground truth here is the entry level survival function of the negative
#' binomial, which is exact, so this generator supports the end to end
#' denoised-expectation metric as well as the probability metrics. Thresholds
#' follow the quantile pooling rule from the original `dfbm()`.
simulate_from_counts <- function(N, P, Tn, rank = 3L, size = 2,
                                 zi_prob = 0.1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  A <- matrix(rnorm(N * rank), N, rank)
  B <- matrix(rnorm(P * rank), P, rank)
  log_mu <- tcrossprod(A, B) / sqrt(rank)
  log_mu <- (log_mu - mean(log_mu)) / stats::sd(as.vector(log_mu))
  mu_mat <- exp(1.5 + log_mu)

  counts <- matrix(rnbinom(N * P, size = size, mu = as.vector(mu_mat)), N, P)
  dropout <- matrix(rbinom(N * P, 1, zi_prob), N, P)
  observed <- counts * (1 - dropout)

  # Thresholds: pooled column deciles, thinned so that each step drops a
  # reasonable share of the surviving entries (the dfbm() heuristic).
  pool <- sort(unique(as.vector(apply(observed, 2, stats::quantile,
                                      probs = seq(0.1, 0.9, 0.1)))))
  pool <- pool[pool < stats::quantile(observed, 0.99)]
  if (length(pool) < Tn) pool <- sort(unique(c(pool, seq_len(Tn))))
  thresholds <- unique(stats::quantile(pool, probs = seq(0, 1, length.out = Tn),
                                       type = 1))
  Tn <- length(thresholds)

  Y <- lapply(thresholds, function(d) 1 * (observed > d))
  # True survival accounts for the zero inflation: P(count > d) = (1-zi) * NB tail.
  S <- lapply(thresholds, function(d) {
    (1 - zi_prob) * stats::pnbinom(d, size = size, mu = mu_mat, lower.tail = FALSE)
  })

  list(Y = Y, S = S, X = lapply(S, stats::qlogis),
       counts = observed, truth_counts = counts,
       M = (1 - zi_prob) * mu_mat, thresholds = thresholds,
       prevalence = vapply(Y, mean, numeric(1L)),
       params = list(N = N, P = P, Tn = Tn, rank = rank, size = size,
                     zi_prob = zi_prob))
}
