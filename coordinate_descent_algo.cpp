#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Soft Thresholding function
double soft_thresholding_cpp(double x, double a) {
  if (x > a)
    return x - a;
  else if (x < -a)
    return x + a;
  else
    return 0;
}

// Coordinate Gradient Descent function
// [[Rcpp::export]]
arma::vec cg_cd(arma::mat X, arma::vec b, double lambda, double tol, int maxiter){
  // Get the number of columns of X
  int p = X.n_cols;
  // Get the number of rows of X
  int n = X.n_rows;
  
  // Initialize a vector u of zeros with the same size as b
  arma::vec u = arma::zeros<arma::vec>(b.n_elem);
  // Compute the matrix A
  arma::mat A = (1.0/n) * X.t() * X + arma::eye<arma::mat>(p,p) * 0.01;
  
  // Start the main loop
  for (int j = 0; j < maxiter; j++){
    // Loop over each element of u
    for (int k = 0; k < u.n_elem; k++){
      // Apply the soft thresholding operator to b[k]
      double soft_threshold = soft_thresholding_cpp(b[k], b.n_elem * lambda);
      // Compute the inverse of A[k,k]
      double inv_tran_a = 1.0 / A(k,k);
      // Update u[k]
      u[k] = inv_tran_a * soft_threshold;
    }
    // Check the stopping criterion
    if (arma::norm(u - u, "fro") < tol) { break; }
  }
  
  // Return the final u
  return u;
}
