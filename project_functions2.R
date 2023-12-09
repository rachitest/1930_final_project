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
  
  A<-(1/n)*t(X)%*%X+diag(0.01,p)
  
  #initial variables
  u <- matrix(0, nrow = p, ncol = 1)
  z <- matrix(0, nrow = p, ncol = 1)
  w <- matrix(0, nrow = p, ncol = 1)
  
  b_matrix <- matrix(b, nrow = n, ncol=1, byrow = F)
  
  #identity matrix
  I <- diag(p)
  
  #pre-compute matrices for efficiency
  #XtX_over_n_plus_I_inv <- solve((1/n) * t(X) %*% X + 0.01 * I + diag(p))
  
  #ADMM loop
  for (iter in 1:maxit) {
    for (k in 1:length(u)){
    #print(iter)
    #update u
    soft_threshold<-as.matrix(rep(soft_thresholding(b[k],length(b)*lambda),length(b)))
    inv_tran_a<-1/(t(A[,k]))
    u[k] <- inv_tran_a%*%soft_threshold
    #u <- XtX_over_n_plus_I_inv %*% (t(X) %*% b_matrix + rho * (z - w))
    
    #update z (soft-thresholding for L1 penalty)
    #z <- soft_threshold(u + w, lambda / rho)
    z <- soft_thresholding(u[k] + w[k], lambda / rho)
    
    #update w
    w <- w + u - z
    
    #check for convergence
    if (norm(u - z) < tol) {
      break
      }
    }
  }
  
  #calculate the objective value
  objective_value <- 0.5 * t(u) %*% ((1/n) * t(X) %*% X + 0.01 * I + diag(p)) %*% u - t(u) %*% b + lambda * sum(abs(u))
  
  return(list(u = u, objective_value = objective_value))
}



#set parameters for simulation
rho <- 0.9
n <- 200
p <- 4

#create covariance matrix
cov_matrix <- matrix(0, nrow = p, ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    cov_matrix[i, j] <- rho^abs(i - j)
  }
}

#simulate data
set.seed(100)  # for reproducibility
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cov_matrix)

#set b
b <- c(1, rep(0, p - 1))

#lambda value
lambda_value <- 0.01

#ADMM  and CD algorithms and compare their results
result_admm <- cg_admm(X, b, lambda_value)
re <- cg_cd(X, b, lambda_value)
