#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Soft thresholding function
double soft_thresholding(double x, double a){
  if (x > a) return (x - a);
  if (x < -a) return (x + a);
  return 0.0;
}

// [[Rcpp::export]]
arma::vec cg_cd(arma::mat X, arma::vec b, double lambda, double tol=1e-4, int maxiter=1000){
  int p = X.n_cols;
  arma::vec u(p);
  u.zeros();
  
  arma::mat A = (1.0/X.n_rows)*X.t()*X + 0.01*arma::eye(p,p);
  arma::vec u_old = u;
  
  for (int j = 0; j < maxiter; j++){
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

// [[Rcpp::export]]
arma::vec cg_cd_alpha(arma::mat X, arma::vec b, double lambda, double alpha=0, double tol=1e-4, int maxiter=1000){
  int p = X.n_cols;
  arma::vec u(p);
  u.zeros();
  
  arma::mat A = (1.0/X.n_rows)*X.t()*X + 0.01*arma::eye(p,p);
  arma::vec u_old = u;
  
  for (int j = 0; j < maxiter; j++){
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
