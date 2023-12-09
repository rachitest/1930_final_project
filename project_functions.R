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

################# cg_cd #################
cg_cd_r <- function(X,b,lambda,tol=1e-4,maxiter=1000,quiet=FALSE){
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
#soft-thresholding function
soft_threshold <- function(x, lambda) {
  pmax(x - lambda, 0) - pmax(-x - lambda, 0)
}

#function to solve regression using ADMM
admm_r <- function(X, b, lambda, rho, max_iter=1000, tolerance=1e-4) {
  rho = 4
  n <- ncol(X)
  m <- nrow(X)
  
  #initial variables
  u <- matrix(0, n, 1)
  z <- matrix(0, n, 1)
  v <- matrix(0, n, 1)
  
  #precompute A and A^Tb
  A <- (1/m) * t(X) %*% X + 0.01 * diag(n)
  Atb <- t(A) %*% b
  
  #ADMM iterations
  for (k in 1:max_iter) {
    #u-update
    u = solve(A + rho * diag(n)) %*% (Atb + rho * (z - v))
    
    #z-update
    z <- soft_threshold(u + v, lambda / rho)
    
    #v-update
    v = v + u - z
    
    #check for convergence
    if (sum(abs(u - z)) < tolerance) {
      break
    }
  }
  
  return(list(u = u, z = z))
}