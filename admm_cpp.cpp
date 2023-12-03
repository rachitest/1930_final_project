#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::colvec soft_threshold(arma::colvec x, double tau)
{
  int p = x.n_rows;
  return (arma::sign(x) % arma::max(arma::zeros<arma::colvec>(p), arma::abs(x) -
          tau));
}

// [[Rcpp::export]]
arma::colvec cg_admm_cpp(arma::mat X, arma::colvec b, double lambda, double tol, int maxit) {
  int n = X.n_rows;
  int p = X.n_cols;
  double rho = .9;
  
//initial variables
    arma::mat u = arma::zeros(p, 1);
    arma::mat z = arma::zeros(p, 1);
    arma::mat w = arma::zeros(p, 1);
    
    arma::mat b_matrix = arma::repmat(b, n/p,1);
    
//identity matrix
    arma::mat I = arma::eye(p,p);
      
//pre-compute matrices for efficiency
    arma::mat XtX_over_n_plus_I_inv = arma::inv((1/n) * X.t() * X + 0.01 * I + arma::eye(p,p));
        
//ADMM loop
        for (int iter = 1; iter <= maxit; iter++) {
//print(iter)
//update u
          u = XtX_over_n_plus_I_inv * (X.t() * b_matrix + rho * (z - w));
          
//update z (soft-thresholding for L1 penalty)
          z = soft_threshold(u + w, lambda / rho);
            
//update w
            w = w + u - z;
            
//check for convergence
            if (norm(u - z) < tol) {
              break;
            }
        }  
//calculate the objective value
        //double objective_value = 0.5 * u.t() * ((1/n) * X.t() %*% X + 0.01 * I + arma::eye(p)) * u - u.t() * b + lambda * sum(abs(u))
          
          return(u);
}