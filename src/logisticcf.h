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
double penalized_loss(const arma::mat &X, const arma::mat &Z,
                         const arma::mat &A, const arma::mat &B,
                         const arma::mat& A_penalty, const arma::mat& B_penalty){

  arma::mat S = A.t() * B;
  double logloss = accu(-Z % (X % S - log1exp(S)));
  double ridge = accu(A_penalty % A % A) + accu(B_penalty % B % B);
  double mean_loss = (logloss + ridge) / accu(Z);
  return mean_loss;

}


Rcpp::List logisticcf(const arma::mat &X, const arma::mat &Z,
                int K, double lambda, unsigned max_iter=1000, double tol=1e-5);


#endif
