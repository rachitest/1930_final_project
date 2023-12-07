#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::colvec soft_threshold(arma::colvec x, double lambda)
{
  int p = x.n_rows;
  return (arma::max(x - lambda, arma::zeros<arma::colvec>(p)) - arma::max(-x - lambda, arma::zeros<arma::colvec>(p)));
}

// [[Rcpp::export]]
arma::colvec cg_admm(arma::mat X, arma::colvec b, double lambda, double tol, int maxit) {
  int m = X.n_rows;
  int n = X.n_cols;
  
//initial variables
    arma::mat u = arma::zeros(n, 1);
    arma::mat z = arma::zeros(n, 1);
    arma::mat v = arma::zeros(n, 1);
    double rho = 0.9;
    
//identity matrix
    arma::mat I = arma::eye(n,n);
      
//pre-compute A and A^Tb
    arma::mat A = (1/m) * X.t() * X + .01 * I;
    arma::mat Atb = A.t() * b;
        
//ADMM iterations
        for (int k = 1; k <= maxit; k++) {
//update u
          u = arma::inv(A + rho * I) * (Atb + rho * (z - v));
          
//update z (soft-thresholding for L1 penalty)
          z = soft_threshold(u + v, lambda / rho);
            
//update v
            v = v + u - z;
            
//check for convergence
            if (arma::norm(u - z) < tol) {
              break;
            }
        }  
//calculate the objective value
        //double objective_value = 0.5 * u.t() * ((1/n) * X.t() %*% X + 0.01 * I + arma::eye(p)) * u - u.t() * b + lambda * sum(abs(u))
          return u;
}