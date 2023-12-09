#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Soft thresholding function
  //'
  //' @param x A numeric value.
  //' @param a A numeric value for threshold.
  //' @return The soft-thresholded value.
  //'
  //' @examples
  //' soft_thresholding(1.5, 1)
  //' soft_thresholding(-1.5, 1)
  //' soft_thresholding(0.5, 1)
double soft_thresholding(double x, double a){
  if (x > a) return (x - a);
  if (x < -a) return (x + a);
  return 0.0;
}

//' Coordinate Gradient Descent
  //'
  //' @param X A matrix of predictors of dimension n x p.
  //' @param b A vector of responses of length p.
  //' @param lambda A numeric value for regularization.
  //' @param tol A numeric value for tolerance.
  //' @param maxit An integer for maximum number of iterations.
  //' @return A vector of coefficients of length p.
  //'
  //' @examples
  //' \donotrun {
  //' set.seed(2023)
  //' n <- 100
  //' p <- 5
  //' b <- rep(0, p)
  //' b[1] <- 1
  //' X <- matrix(rnorm(n * p), nrow = n, ncol = p) 
  //' lambda <- 0.2
  //' cg_cd_result = cg_cd(X, b, lambda)
  //' print(cg_cd_result)
  //' }
// [[Rcpp::export]]
arma::vec cg_cd(arma::mat X, arma::vec b, double lambda, double tol=1e-4, int maxit=1000){
  int p = X.n_cols;
  arma::vec u(p);
  u.zeros();
  
  arma::mat A = (1.0/X.n_rows)*X.t()*X + 0.01*arma::eye(p,p);
  arma::vec u_old = u;
  
  for (int j = 0; j < maxit; j++){
    for (int k = 0; k < p; k++){
      double r_k = b[k] - dot(A.row(k), u) + A(k,k)*u[k];
      double z_k = soft_thresholding(r_k, lambda) / A(k,k);
      u[k] = z_k;
    }
    
    if (norm(u_old - u, "fro") < tol) { break; }
    u_old = u;
  } 
  
  return u;
}

//' Coordinate Gradient Descent with Elastic Net Regularization
 //'
 //' @param X A matrix of predictors of dimension n x p.
 //' @param b A vector of responses of length p.
 //' @param lambda A numeric value for regularization.
 //' @param alpha A numeric value for elastic net mixing parameter.
 //' @param tol A numeric value for tolerance.
 //' @param maxit An integer for maximum number of iterations.
 //' @return A vector of coefficients of length p.
 //'
 //' @examples
 //' \donotrun {
 //' set.seed(2023)
 //' n <- 100
 //' p <- 5
 //' b <- rep(0, p)
 //' b[1] <- 1
 //' X <- matrix(rnorm(n * p), nrow = n, ncol = p) 
 //' lambda <- 0.2
 //' cg_cd_result = cg_cd_alpha(X, b, lambda)
 //' print(cg_cd_result)
 //' }
// [[Rcpp::export]]
arma::vec cg_cd_alpha(arma::mat X, arma::vec b, double lambda, double tol=1e-4, int maxit=1000, double alpha=0){
  int p = X.n_cols;
  arma::vec u(p);
  u.zeros();
  
  arma::mat A = (1.0/X.n_rows)*X.t()*X + 0.01*arma::eye(p,p);
  arma::vec u_old = u;
  
  for (int j = 0; j < maxit; j++){
    for (int k = 0; k < p; k++){
      double r_k = b[k] - dot(A.row(k), u) + A(k,k)*u[k];
      double z_k = soft_thresholding(r_k, lambda*(1-alpha)) / (A(k,k) + lambda*alpha);
      u[k] = z_k;
    }
    
    if (norm(u_old - u, "fro") < tol) { break; }
    u_old = u;
  } 
  
  return u;
}
