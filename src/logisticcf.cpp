#include <RcppArmadillo.h>
#include "logisticcf.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


double mse_loss(const arma::mat &estim, const arma::mat &target, const arma::umat &mask){
  arma::mat diff = estim - target;
  arma::mat diff2 = pow(diff, 2);
  unsigned n_obs = accu(mask);
  double loss = accu(diff2.elem(find(mask))) / n_obs;
  return loss;
}


struct cfoutput cf(const arma::mat &S, const arma::umat &Z, arma::mat &init_A, arma::mat &init_B,
        unsigned K, double lambda, unsigned max_iter, double tol, int cap){

  bool express = all(vectorise(Z) == 1); // if no missing entries, we have a faster implementation
  unsigned N = S.n_rows;
  unsigned P = S.n_cols;
  arma::mat A, B;
  if (init_A.is_zero()){ // no A or B from previous iterations
    A = randn(K, N) / 5;
    B = randn(K, P) / 5;
    B.row(0) = sum(S % Z, 0) / sum(Z, 0);
  } else{
    A = init_A;
    B = init_B;
  }
  unsigned iter = 0;
  arma::mat S_lowrank = A.t() * B;
  double loss = mse_loss(S_lowrank, S, Z);
  vec loss_trace = zeros<vec>(max_iter);
  arma::mat new_A, new_B;
  while (iter < max_iter){
    if (express){ // all entries are observed
      new_A = solve(B * B.t() + lambda * eye(K, K), B * S.t(),
                    solve_opts::likely_sympd);
      new_B = solve(new_A * new_A.t() + lambda * eye(K, K), new_A * S,
                    solve_opts::likely_sympd);
    } else {
      new_A = ALS_update(S, Z, A, B, K, lambda, true);
      new_B = ALS_update(S, Z, new_A, B, K, lambda, false);
      // for (unsigned i = 0; i < N; i++) { // update A
      //   uvec subrow = {i};
      //   uvec subcol = find(Z.row(i) == 1);
      //   arma::mat B_filter = B.cols(subcol);
      //   vec S_selected = vectorise(S.row(i));
      //   vec S_filter = S_selected.elem(subcol);
      //   new_A.col(i) = solve(B_filter * B_filter.t() + lambda * eye(K, K),
      //             B_filter * S_filter,
      //             solve_opts::likely_sympd);
      // }
      // for (unsigned j = 0; j < P; j++) { // update B
      //   uvec subrow = find(Z.col(j) == 1);
      //   uvec subcol = {j};
      //   arma::mat A_filter = new_A.cols(subrow);
      //   vec S_selected = vectorise(S.col(j));
      //   vec S_filter = S_selected.elem(subrow);
      //   new_B.col(j) = solve(A_filter * A_filter.t() + lambda * eye(K, K),
      //             A_filter * S_filter,
      //             solve_opts::likely_sympd);
      // }
    }
    arma::mat new_S_lowrank = new_A.t() * new_B;
    A = new_A;
    B = new_B;
    double new_loss = mse_loss(new_S_lowrank, S, Z);
    loss_trace(iter) = new_loss;
    iter++;

    if (loss - new_loss < tol){// check convergence
      loss = new_loss;
      S_lowrank = new_S_lowrank;
      break;
    } else{
      loss = new_loss;
      S_lowrank = new_S_lowrank;
    }
  }

  S_lowrank.elem(find(S_lowrank > cap)).fill(cap);
  S_lowrank.elem(find(S_lowrank < -cap)).fill(-cap);

  struct cfoutput result{A, B, S_lowrank, loss_trace(span(0, iter-1))};

  return result;

}


double logit_loss(const arma::mat &estim, const arma::umat &target, const arma::umat &mask){
  arma::mat loglik = target % log(estim) + (1 - target) % log(1 - estim);
  unsigned n_obs = accu(mask);
  double loss = -accu(loglik.elem(find(mask))) / n_obs;
  return loss;
}


// [[Rcpp::export]]
List logisticcf(const arma::umat &X, const arma::umat &Z,
                unsigned K, double lambda, unsigned max_iter, double tol, int cap){

  arma::mat W = 2.0*conv_to<arma::mat>::from(X) - 1; // convert 0-1 to -1, 1
  arma::mat init_A = {0};
  arma::mat init_B = {0};
  struct cfoutput init_fit = cf(4*W, Z, init_A, init_B, K, lambda, max_iter, tol, cap); // initial fit
  arma::mat A = init_fit.A;
  arma::mat B = init_fit.B;
  arma::mat S = init_fit.representation;
  arma::mat pi = 1/(1+exp(-S));

  double loss = logit_loss(pi, X, Z);
  vec loss_trace = zeros<vec>(max_iter);

  unsigned iter = 0;
  while (iter < max_iter){
    arma::mat S_star = S + 4*W/(1+exp(W%S)); // MM algorithm
    struct cfoutput new_fit = cf(S_star, Z, A, B, K, lambda, max_iter, tol, cap);
    A = new_fit.A;
    B = new_fit.B;
    S = new_fit.representation;
    pi = 1/(1+exp(-S));
    double new_loss = logit_loss(pi, X, Z);
    loss_trace(iter) = new_loss;
    iter++;
    if (loss - new_loss < tol){
      loss = new_loss;
      break;
    } else{
      loss = new_loss;
    }
  }

  Rcpp::List result= Rcpp::List::create(Rcpp::Named("A") = A,
                                        Rcpp::Named("B") = B,
                                        Rcpp::Named("pi") = pi,
                                        Rcpp::Named("loss_trace") = loss_trace(span(0, iter - 1)));

  return result;
}

