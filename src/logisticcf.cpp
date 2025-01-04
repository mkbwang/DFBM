#include <RcppArmadillo.h>
#include "logisticcf.h"

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
List logisticcf(const arma::mat &X, const arma::mat &Z,
                unsigned K, double lambda, unsigned max_iter, double tol){

  unsigned N = X.n_rows;
  unsigned P = X.n_cols;
  arma::mat W = 2*X- 1;
  arma::mat XZ = X % Z;
  arma::mat WZ = W % Z;

  // initialize A and B
  arma::mat A = 0.1 * arma::randn(K, N);
  arma::vec rowmeanW = sum(WZ, 1) / sum(Z, 1);
  A.row(0) = rowmeanW.t();
  arma::mat B = 0.1 * arma::randn(K, P);
  B.row(0) = sum(WZ, 0) / sum(Z, 0);

  arma::mat S = A.t() * B;
  arma::mat pi = 1.0 / (1.0 + exp(-S));

  arma::mat A_penalty = mat(A.n_rows, A.n_cols, fill::value(lambda));
  A_penalty.row(0).fill(0); // first row does not have penalty
  arma::mat B_penalty = mat(B.n_rows, B.n_cols, fill::value(lambda));
  B_penalty.row(0).fill(0); // first row does not have penalty

  double loss = penalized_loss(X, Z, A, B, A_penalty, B_penalty);
  arma::vec loss_trace(max_iter+1, fill::zeros);
  loss_trace(0) = loss;
  unsigned iter = 0;

  while (iter < max_iter){

    // update A
    arma::mat g_star = -B * XZ.t()  + B * (Z.t() % pi.t()) + 2 * (A_penalty % A);
    arma::mat f_star = 0.25 * (B % B) * Z.t() + 2 * A_penalty;
    A = A - g_star / (2*f_star);

    // update B
    S = A.t() * B;
    pi = 1.0 / (1.0 + exp(-S));
    arma::mat g_dagger = -A * XZ + A * (Z % pi) + 2 * (B_penalty % B);
    arma::mat f_dagger = 0.25 * (A % A) * Z + 2 * B_penalty;
    B = B - g_dagger / (2*f_dagger);

    S = A.t() * B;
    pi = 1.0 / (1.0 + exp(-S));
    double new_loss = penalized_loss(X, Z, A, B, A_penalty, B_penalty);

    iter++;
    loss_trace(iter) = new_loss;
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
                                        Rcpp::Named("loss_trace") = loss_trace(span(0, iter)));

  return result;
}

