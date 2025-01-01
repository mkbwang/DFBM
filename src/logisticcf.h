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


struct ALS: public RcppParallel::Worker{
  const arma::mat &S; // target
  const arma::umat &Z; // mask
  arma::mat &A;
  arma::mat &B;
  const double lambda;
  const unsigned K;
  bool update_A;

  ALS (const mat &target, const umat &mask, mat &init_A, mat &init_B,
       const double lambda,  const unsigned K, bool update_A):
    S(target), Z(mask), A(init_A), B(init_B), lambda(lambda), K(K), update_A(update_A){}

  void operator()(size_t begin, size_t end){
    if(update_A){
      for(unsigned i=begin; i<end; i++){
        uvec subrow = {i};
        uvec subcol = find(Z.row(i) == 1);
        mat B_filter = B.cols(subcol);
        vec S_selected = vectorise(S.row(i));
        vec S_filter = S_selected.elem(subcol);
        A.col(i) = solve(B_filter * B_filter.t() + lambda * eye(K, K),
                  B_filter * S_filter,
                  solve_opts::likely_sympd);
      }
    }else{
      for(unsigned j=begin; j<end; j++){
        uvec subrow = find(Z.col(j) == 1);
        uvec subcol = {j};
        mat A_filter = A.cols(subrow);
        vec S_selected = vectorise(S.col(j));
        vec S_filter = S_selected.elem(subrow);
        B.col(j) = solve(A_filter * A_filter.t() + lambda * eye(K, K),
                  A_filter * S_filter,
                  solve_opts::likely_sympd);
      }
    }
  }
};

arma::mat ALS_update(const arma::mat &S, const arma::umat &Z, arma::mat &A, arma::mat &B,
                     unsigned K, double lambda, bool update_A){
  if(update_A){
    arma::mat updated_A = A;
    size_t N = S.n_rows;
    ALS als_obj(S, Z, updated_A, B, lambda, K, true);
    parallelFor(0, N, als_obj, 20);
    return updated_A;
  } else{
    arma::mat updated_B = B;
    size_t P = S.n_cols;
    ALS als_obj(S, Z, A, updated_B, lambda, K, false);
    parallelFor(0, P, als_obj, 20);
    return updated_B;
  }

}

double mse_loss(const arma::mat &estim, const arma::mat &target, const arma::umat &mask);
double logit_loss(const arma::mat &estim, const arma::umat &target, const arma::umat &mask);

struct cfoutput cf(const arma::mat &S, const arma::umat &Z, arma::mat &init_A, arma::mat &init_B,
          unsigned K, double lambda, unsigned max_iter=1000, double tol=1e-5, int cap=10);


Rcpp::List logisticcf(const arma::umat &X, const arma::umat &Z,
                int K, double lambda, unsigned max_iter=1000, double tol=1e-5, int cap=10);


#endif
