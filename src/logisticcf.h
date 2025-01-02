#include <RcppArmadillo.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
using namespace arma;
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

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


struct ALS_a: public RcppParallel::Worker{
  const arma::mat &S; // target
  const arma::umat &Z; // mask
  const arma::mat &B;
  const double lambda;
  const unsigned K;

  arma::mat &A;

  ALS_a (const mat &target, const umat &mask, mat &init_A, const mat &init_B,
       const double lambda,  const unsigned K):
    S(target), Z(mask), A(init_A), B(init_B), lambda(lambda), K(K){}

  void operator()(size_t begin, size_t end){
      for(size_t i=begin; i<end; i++){
        uvec subcol = find(Z.row(i) == 1);
        mat B_filter = B.cols(subcol);
        vec S_selected = vectorise(S.row(i));
        vec S_filter = S_selected.elem(subcol);
        A.col(i) = solve(B_filter * B_filter.t() + lambda * eye(K, K),
                  B_filter * S_filter,
                  solve_opts::likely_sympd);
      }
  }
};

arma::mat ALS_updateA(const arma::mat &S, const arma::umat &Z, arma::mat &A, const arma::mat &B,
                     unsigned K, double lambda){
    arma::mat updated_A = A;
    size_t N = S.n_rows;
    ALS_a als_obj(S, Z, updated_A, B, lambda, K);
    parallelFor(0, N, als_obj, 10);
    return updated_A;
}



struct ALS_b: public RcppParallel::Worker{
  const arma::mat &S; // target
  const arma::umat &Z; // mask
  const arma::mat &A;
  const double lambda;
  const unsigned K;

  arma::mat &B;

  ALS_b (const mat &target, const umat &mask,const mat &init_A, mat &init_B,
         const double lambda,  const unsigned K):
    S(target), Z(mask), A(init_A), B(init_B), lambda(lambda), K(K){}

  void operator()(size_t begin, size_t end){
    for (size_t j = begin; j < end; j++) { // update B
      uvec subrow = find(Z.col(j) == 1);
      arma::mat A_filter = A.cols(subrow);
      vec S_selected = vectorise(S.col(j));
      vec S_filter = S_selected.elem(subrow);
      B.col(j) = solve(A_filter * A_filter.t() + lambda * eye(K, K),
                A_filter * S_filter,
                solve_opts::likely_sympd);
    }
  }
};


arma::mat ALS_updateB(const arma::mat &S, const arma::umat &Z, const arma::mat &A, arma::mat &B,
                      unsigned K, double lambda){
  arma::mat updated_B = B;
  size_t P = S.n_cols;
  ALS_b als_obj(S, Z, A, updated_B, lambda, K);
  parallelFor(0, P, als_obj, 10);
  return updated_B;
}

double mse_loss(const arma::mat &estim, const arma::mat &target, const arma::umat &mask);
double logit_loss(const arma::mat &estim, const arma::umat &target, const arma::umat &mask);

struct cfoutput cf(const arma::mat &S, const arma::umat &Z, const arma::mat &init_A, const arma::mat &init_B,
          unsigned K, double lambda, unsigned max_iter=1000, double tol=1e-5, int cap=10);


Rcpp::List logisticcf(const arma::umat &X, const arma::umat &Z,
                int K, double lambda, unsigned max_iter=1000, double tol=1e-5, int cap=10);


#endif
