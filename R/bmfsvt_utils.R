

# utility functions for binary matrix factorization with singular value thresholding algorithm


#' Offset of the first mask
#'
#' @param mu offsets, a vector of length T or a T x P matrix
#' @returns `mu[1]`, or the first row of `mu` for column offsets
#' @keywords internal
offset_head <- function(mu) {
  if (is.matrix(mu)) mu[1L, ] else mu[1L]
}


#' Stacked cross entropy over all masks
#' @param X list of natural parameter matrices
#' @param Y list of binary matrices
#' @param train optional training mask, see [logistic_ce()]
#' @returns the summed negative log likelihood
#' @keywords internal
stack_ce <- function(X, Y, train = NULL) {
  total <- 0
  for (t in seq_along(X)) total <- total + logistic_ce(X[[t]], Y[[t]], train)
  total
}


#' Nuclear norm penalty of the current blocks
#'
#' @param Z list of blocks
#' @param dvals singular values returned by the thresholding step, before clipping
#' @param lambda shrinkage parameters
#' @param n_clipped how many entries clipping actually altered
#' @param mode either "approx" or "exact"
#' @returns the value of \eqn{\sum_t \lambda_t \|Z^{(t)}\|_*}
#' @details
#' When no entry was clipped the shrunk singular values are exactly those of
#' `Z`, so the approximation is not an approximation at all. Recomputing is only
#' necessary once clipping has bent the blocks away from the SVT output.
#' @keywords internal
nuclear_penalty <- function(Z, dvals, lambda, n_clipped, mode = "approx") {
  if (mode == "approx" || n_clipped == 0L) {
    return(sum(vapply(seq_along(Z), function(t) lambda[t] * sum(dvals[[t]]),
                      numeric(1L))))
  }
  sum(vapply(seq_along(Z), function(t) {
    lambda[t] * sum(svd(Z[[t]], nu = 0L, nv = 0L)$d)
  }, numeric(1L)))
}



#' Balanced right factors of the fitted blocks
#'
#' @param V list of T `P x r_t` right singular vector matrices, as in `fit$V`
#' @param dvals list of T shrunken singular value vectors, as in `fit$dvals`
#' @returns list of T `P x r_t` matrices \eqn{B_t = V_t D_t^{1/2}}
#' @details
#' The square root is not a normalization choice. The nuclear norm has the
#' variational form \eqn{\|Z\|_* = \min_{Z = AB'} (\|A\|_F^2 + \|B\|_F^2)/2},
#' attained only at the balanced split \eqn{A = UD^{1/2}}, \eqn{B = VD^{1/2}}.
#' Only with that `B` is a plain ridge on a new row's scores that row's share of
#' the nuclear norm penalty. With \eqn{B = VD} the same ridge would shrink weak
#' components harder than strong ones; folding the training rows back in that
#' way was measured to miss the fit by 11.7 on the logit scale.
#'
#' Built with `rep(each = )` rather than `diag()`, which for a single component
#' returns an identity matrix of size `floor(sqrt(d))` instead.
#' @keywords internal
factor_blocks <- function(V, dvals) {
  lapply(seq_along(V), function(t) {
    V[[t]] * rep(sqrt(dvals[[t]]), each = nrow(V[[t]]))
  })
}



#' Ridge fold-in of new rows against fixed right factors
#'
#' @param Y list of T binary n x P masks for the new rows
#' @param B list of T `P x r_t` balanced right factors, see [factor_blocks()]
#' @param mu0 offset of the first mask, scalar or length P
#' @param nu offset increments, vector of length T-1 or (T-1) x P matrix
#' @param lambda the training shrinkage parameters, length T
#' @param tol a row stops once the Euclidean norm of its own gradient falls
#'   below this
#' @param max_iter maximum number of iterations
#' @returns a list with the scores `A` (list of T `n x r_t` matrices), the
#'   unclipped increments `Z` (list of T n x P matrices), and per row
#'   `grad_norm`, `n_iter` and `converged`
#' @details
#' Solves, for the new rows only,
#' \deqn{\min_{A_1..A_T} \sum_t \mathrm{CE}(X^{(t)}; Y^{(t)}) +
#'   \sum_t \frac{\lambda_t}{2}\|A_t\|_F^2, \qquad
#'   X^{(t)} = \mu_t + \sum_{t' \le t} A_{t'} B_{t'}^\top,}
#' whose gradient is \eqn{\Psi_t B_t + \lambda_t A_t} with \eqn{\Psi} from
#' [grad_backward()]. The problem separates by row and is strongly convex.
#'
#' **Step size.** The Hessian block for the pair (t, s) is
#' \eqn{B_t^\top \mathrm{diag}(\sum_{t' \ge \max(s,t)} w_{t'}) B_s} with
#' \eqn{w \le 1/4}, so its norm is at most
#' \eqn{(T - \max(s,t) + 1)/4 \cdot \sqrt{d_{t,1} d_{s,1}}}. Summing over s gives a
#' block Lipschitz constant \eqn{L_t} valid for every row at every iterate, so
#' no backtracking is needed; the step is \eqn{1/(L_t + \lambda_t)}.
#'
#' **Each row is solved independently.** Momentum and restarts are kept per
#' row, and a row is frozen at the iterate where its own gradient norm first
#' drops below `tol`. The answer for a row therefore does not depend on which
#' other rows share the batch; batching is only for speed.
#' @keywords internal
foldin_rows <- function(Y, B, mu0, nu, lambda, tol = 1e-6, max_iter = 2000L) {
  Tn <- length(Y)
  n <- nrow(Y[[1L]])
  ranks <- vapply(B, ncol, integer(1L))
  tt <- seq_len(Tn)

  # / ||B_t||_2^2 is the largest shrunken singular value, since V is orthonormal.
  nb <- vapply(B, function(b) if (ncol(b) == 0L) 0 else sqrt(max(colSums(b^2))),
               numeric(1L))
  weight <- outer(tt, tt, function(s, t) (Tn - pmax(s, t) + 1) / 4)
  L <- pmax(rowSums(weight * outer(nb, nb)), 1e-8)
  step <- 1 / (L + lambda)

  A_out <- lapply(ranks, function(r) matrix(0, n, r))
  grad_norm <- rep(NA_real_, n)
  n_iter <- rep(as.integer(max_iter), n)
  converged <- logical(n)

  rows <- seq_len(n)
  Yw <- Y
  A <- lapply(ranks, function(r) matrix(0, n, r))
  W <- A
  s <- rep(1, n)
  gn <- rep(NA_real_, n)

  for (m in seq_len(max_iter)) {
    Zw <- lapply(tt, function(t) tcrossprod(W[[t]], B[[t]]))
    Psi <- stack_grad(Zw, mu0, nu, Yw)$Psi
    G <- lapply(tt, function(t) Psi[[t]] %*% B[[t]] + lambda[t] * W[[t]])
    gn <- sqrt(Reduce(`+`, lapply(G, function(g) rowSums(g^2))))

    done <- gn < tol
    if (any(done)) {
      idx <- rows[done]
      for (t in tt) {
        if (ranks[t] > 0L) A_out[[t]][idx, ] <- W[[t]][done, , drop = FALSE]
      }
      grad_norm[idx] <- gn[done]
      n_iter[idx] <- m - 1L
      converged[idx] <- TRUE
      keep <- !done
      rows <- rows[keep]
      if (length(rows) == 0L) break
      Yw <- lapply(Yw, function(y) y[keep, , drop = FALSE])
      A <- lapply(A, function(a) a[keep, , drop = FALSE])
      W <- lapply(W, function(w) w[keep, , drop = FALSE])
      G <- lapply(G, function(g) g[keep, , drop = FALSE])
      s <- s[keep]
      gn <- gn[keep]
    }

    A_new <- lapply(tt, function(t) W[[t]] - step[t] * G[[t]])
    # Gradient restart, per row: drop the momentum of any row whose step moved
    # against its proximal direction.
    ip <- Reduce(`+`, lapply(tt, function(t) {
      rowSums((W[[t]] - A_new[[t]]) * (A_new[[t]] - A[[t]]))
    }))
    s_next <- (1 + sqrt(1 + 4 * s^2)) / 2
    fac <- (s - 1) / s_next
    restart <- ip > 0
    fac[restart] <- 0
    s_next[restart] <- 1
    # `fac` has one element per row, so recycling it down the columns of an
    # n x r matrix scales each row by its own factor.
    W <- lapply(tt, function(t) A_new[[t]] + fac * (A_new[[t]] - A[[t]]))
    A <- A_new
    s <- s_next
  }

  if (length(rows) > 0L) {
    for (t in tt) {
      if (ranks[t] > 0L) A_out[[t]][rows, ] <- A[[t]]
    }
    grad_norm[rows] <- gn
    warning(sprintf(
      "fold-in did not reach `tol` for %d of %d rows within %d iterations.",
      length(rows), n, max_iter))
  }

  Z <- lapply(tt, function(t) tcrossprod(A_out[[t]], B[[t]]))
  list(A = A_out, Z = Z, grad_norm = grad_norm, n_iter = n_iter,
       converged = converged)
}

