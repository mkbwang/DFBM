#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace std;
using namespace Rcpp;

#ifndef LOGISTICCF_LOGISTICCF_H
#define LOGISTICCF_LOGISTICCF_H


arma::mat log1exp(const arma::mat &input){

  arma::mat output = input;
  uvec mask = find(input < 20);
  output.elem(mask) = log(1 + exp(input.elem(mask)));
  return output;

}

// ridge logloss
double penalized_loss(const arma::mat &S, const arma::uvec &nonzero_xz_index,
                      const arma::uvec & nonzero_z_index, const arma::vec &lambda_vec,
                         const arma::mat& A2, const arma::mat& B2){

  arma::mat log1S = log1exp(S);
  double logloss =  -accu(S.elem(nonzero_xz_index)) + accu(log1S.elem(nonzero_z_index));
  double ridge = accu(A2.each_col() % lambda_vec) + accu(B2.each_col() % lambda_vec);
  // arma::mat S = A.t() * B;
  // double logloss = accu(-Z % (X % S - log1exp(S)));
  // double ridge = accu(A_penalty % A % A) + accu(B_penalty % B % B);
  double mean_loss = (logloss + ridge) / nonzero_z_index.n_elem;
  return mean_loss;

}


Rcpp::List logisticcf(const arma::mat &X, const arma::mat &Z,
                int K, double lambda, unsigned max_iter=1000, double tol=1e-5);


#endif
