######### Functions for 1930 Project
## Created by: Caroline Schaefer
## Last Update: 30NOV2023

## functions
# soft_thresholding
# cg_cd
# cg_admm



soft_thresholding <- function(x,a){
  if (x>a) return (x-a)
  if (x< -a) return (x+a)
  return(0)
}

soft_threshold <- function(x, lambda) {
  pmax(0, x - lambda) - pmax(0, -x - lambda)
}

################# cg_cd #################
cg_cd <- function(X,b,lambda,tol=1e-4,maxiter=1000,quiet=FALSE){
  u <- rep(0, length(b))
  u <- as.matrix(u); X <- as.matrix(X)
  ulist <- list(length=(maxiter+1))
  ulist[[1]] <- u
  
  A<-(1/n)*t(X)%*%X+diag(0.01,p)
  
  for (j in 1:maxiter){
    for (k in 1:length(u)){
      soft_threshold<-as.matrix(rep(soft_thresholding(b[k],length(b)*lambda),length(b)))
      inv_tran_a<-1/(t(A[,k]))
      u[k] <- inv_tran_a%*%soft_threshold 
    }
    ulist[[(j+1)]] <- u
    
    if (norm(ulist[[j]] - u,"F") < tol) { break }
  } 
  
  return(list(u=u)) 
}     


################# cg_admm #################
cg_admm <- function(X, b, lambda, tol = 1e-4, maxit = 1000) {
  n <- nrow(X)
  p <- ncol(X)
  
  #initial variables
  u <- matrix(0, nrow = p, ncol = 1)
  z <- matrix(0, nrow = p, ncol = 1)
  w <- matrix(0, nrow = p, ncol = 1)
  
  b_matrix <- matrix(b, nrow = n, ncol=1, byrow = F)
  
  #identity matrix
  I <- diag(p)
  
  #pre-compute matrices for efficiency
  XtX_over_n_plus_I_inv <- solve((1/n) * t(X) %*% X + 0.01 * I + diag(p))
  
  #ADMM loop
  for (iter in 1:maxit) {
    #print(iter)
    #update u
    u <- XtX_over_n_plus_I_inv %*% (t(X) %*% b_matrix + rho * (z - w))
    
    #update z (soft-thresholding for L1 penalty)
    z <- soft_threshold(u + w, lambda / rho)
    
    #update w
    w <- w + u - z
    
    #check for convergence
    if (norm(u - z) < tol) {
      break
    }
  }
  
  #calculate the objective value
  objective_value <- 0.5 * t(u) %*% ((1/n) * t(X) %*% X + 0.01 * I + diag(p)) %*% u - t(u) %*% b + lambda * sum(abs(u))
  
  return(u = u)
}