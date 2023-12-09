#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Soft thresholding function
 //'
 //' @param x A column vector.
 //' @param lambda A numeric value for threshold.
 //' @return The soft-thresholded value.
 //'
 //' @examples
 //' soft_threshold(1.5, 1)
 //' soft_threshold(-1.5, 1)
 //' soft_threshold(0.5, 1)
 arma::colvec soft_threshold(arma::colvec x, double lambda)
 {
   int p = x.n_rows;
   return (arma::max(x - lambda, arma::zeros<arma::colvec>(p)) - arma::max(-x - lambda, arma::zeros<arma::colvec>(p)));
 }

//' Alternating Direction Method of Multipliers
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
 //' cg_admm_result = cg_admm(X, b, lambda, tol=1e-4, maxit=1000)
 //' print(cg_admm_result)
 //' }
// [[Rcpp::export]]
arma::colvec cg_admm(arma::mat X, arma::colvec b, double lambda, double tol=1e-4, int maxit=1000) {
  int m = X.n_rows;
  int n = X.n_cols;
  
//initial variables
    arma::mat u = arma::zeros(n, 1);
    arma::mat z = arma::zeros(n, 1);
    arma::mat v = arma::zeros(n, 1);
    double rho = .01;
    
//identity matrix
    arma::mat I = arma::eye(n,n);
      
//pre-compute A and A^Tb
    arma::mat A = (1.0/m) * X.t() * X + .01 * I;
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
 
 //' Alternating Direction Method of Multipliers with Elastic Net Regularization
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
 //' cg_admm_result = cg_admm_alpha(X, b, lambda)
 //' print(cg_admm_result)
 //' }
 // [[Rcpp::export]]
 arma::colvec cg_admm_alpha(arma::mat X, arma::colvec b, double lambda, double tol=1e-4, int maxit=1000, double alpha=0) {
   int m = X.n_rows;
   int n = X.n_cols;
   
   //initial variables
   arma::mat u = arma::zeros(n, 1);
   arma::mat z = arma::zeros(n, 1);
   arma::mat v = arma::zeros(n, 1);
   double rho = .01;
   
   //identity matrix
   arma::mat I = arma::eye(n,n);
   
   //pre-compute A and A^Tb
   arma::mat A = (1.0/m) * X.t() * X + .01 * I;
   arma::mat Atb = A.t() * b;
   
   //ADMM iterations
   for (int k = 1; k <= maxit; k++) {
     //update u
     u = arma::inv(A + rho * I) * (b - lambda + rho*z);
     
     //update z (soft-thresholding for L1 penalty)
     //z = soft_threshold(rho*u, lambda*alpha) / (rho + lambda*alpha);
     //z = soft_threshold(u+v,(1-alpha)*lambda/rho) / (rho + lambda*alpha);
     z = soft_threshold(lambda + rho*u, lambda*(1-alpha)) / (rho + lambda*alpha);
     
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