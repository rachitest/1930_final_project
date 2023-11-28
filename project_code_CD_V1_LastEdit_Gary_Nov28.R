##Edit by Gary C Nov28, 4:40PM

soft_thresholding <- function(x,a){
  if (x>a) return (x-a)
  if (x< -a) return (x+a)
  return(0)
} 


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

#Simulation for 4 parametes and 20 samples 
set.seed(100)
n <- 20 #Sample Size

p=4 # number of parameter 

#generate variance covariance matrix pxp
cov_m<-matrix(nrow=p,ncol=p) 
for (i in 1:p){
  for (j in 1:p){
    cov_m[i, j] <- 0.9^(abs(i-j))
  }
}

library(MASS)
X <- mvrnorm(n=5, mu = rep(0, p), Sigma = cov_m) #Generate  nxp of X matrix 
b<-c(1,rep(0,p-1)) 

# find u for lambda = 0.01
re <- cg_cd(X, b, 0.01)
print(re)

