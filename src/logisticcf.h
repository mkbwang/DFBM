#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

#ifndef LOGISTICCF_LOGISTICCF_H
#define LOGISTICCF_LOGISTICCF_H


struct cfoutput{
  arma::mat A;
  arma::mat B;
  arma::mat representation; // product of A and B
  arma::vec loss_trace; // trace of loss

  cfoutput(const arma::mat &result_A, const arma::mat &result_B,
           const arma::mat &result_representation, const arma::vec &losses):
    A(result_A), B(result_B), representation(result_representation), loss_trace(losses){}
};


double mse_loss(const arma::mat &estim, const arma::mat &target, const arma::umat &mask);
double logit_loss(const arma::mat &estim, const arma::umat &target, const arma::umat &mask);

struct cfoutput cf(const arma::mat &S, const arma::umat &Z, arma::mat &init_A, arma::mat &init_B,
          unsigned K, double lambda, unsigned max_iter=1000, double tol=1e-5, int cap=10);


Rcpp::List logisticcf(const arma::umat &X, const arma::umat &Z,
                int K, double lambda, unsigned max_iter=1000, double tol=1e-5, int cap=10);


#endif
