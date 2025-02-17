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
  arma::mat X_t = X.t();
  arma::mat Z_t = Z.t();
  uvec nonzero_z_index = find(Z > 0);
  arma::mat XZ = mat(N, P, fill::zeros);
  XZ.elem(nonzero_z_index) = X.elem(nonzero_z_index);
  uvec nonzero_xz_index = find(XZ > 0);
  arma::mat XZ_t = XZ.t();
  arma::mat WZ = mat(N, P, fill::zeros);
  WZ.elem(nonzero_z_index) = W.elem(nonzero_z_index);

  // initialize A and B
  arma::mat A = 0.1 * arma::randn(K, N);
  arma::vec rowmeanW = sum(WZ, 1) / sum(Z, 1);
  A.row(0) = rowmeanW.t();
  arma::mat A2 = square(A);
  arma::mat B = 0.1 * arma::randn(K, P);
  B.row(0) = sum(WZ, 0) / sum(Z, 0);
  arma::mat B2 = square(B);

  arma::mat S = A.t() * B;
  arma::mat pi = 1.0 / (1.0 + exp(-S));
  arma::mat pi_Z = mat(N, P, fill::zeros);
  pi_Z.elem(nonzero_z_index) = pi.elem(nonzero_z_index);

  arma::vec lambda_vec = lambda * arma::ones(K);
  lambda_vec(0) = 0; // first row does not have penalty
  // arma::mat A_penalty = mat(A.n_rows, A.n_cols, fill::value(lambda));
  // A_penalty.row(0).fill(0); // first row does not have penalty
  // arma::mat B_penalty = mat(B.n_rows, B.n_cols, fill::value(lambda));
  // B_penalty.row(0).fill(0); // first row does not have penalty

  double loss = penalized_loss(S, nonzero_xz_index, nonzero_z_index,
                               lambda_vec, A2, B2);

  arma::vec loss_trace(max_iter+1, fill::zeros);
  loss_trace(0) = loss;
  unsigned iter = 0;

  while (iter < max_iter){

    // update A
    arma::mat g_star = B * (pi_Z.t() - XZ_t);
    g_star += 2 * (A.each_col() % lambda_vec);
    arma::mat f_star = 0.25 * B2 * Z_t;
    f_star.each_col() += 2*lambda_vec;
    A = A - g_star / (2*f_star);
    A2 = square(A);

    // update B
    S = A.t() * B;
    pi = 1.0 / (1.0 + exp(-S));
    pi_Z.elem(nonzero_z_index) = pi.elem(nonzero_z_index);
    arma::mat g_dagger = A * (pi_Z - XZ);
    g_dagger += 2*(B.each_col() % lambda_vec);
    arma::mat f_dagger = 0.25 * A2 * Z;
    f_dagger.each_col() += 2*lambda_vec;
    B = B - g_dagger / (2*f_dagger);
    B2 = square(B);

    S = A.t() * B;
    pi = 1.0 / (1.0 + exp(-S));
    pi_Z.elem(nonzero_z_index) = pi.elem(nonzero_z_index);
    double new_loss = penalized_loss(S, nonzero_xz_index, nonzero_z_index,
                                     lambda_vec, A2, B2);

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

