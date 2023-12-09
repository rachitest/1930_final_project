library(MASS)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("/Users/cschaef3/Documents/GitHub/1930_final_project/coordinate_descent_algo_alt.cpp")
sourceCpp("/Users/cschaef3/Documents/GitHub/1930_final_project/admm_cpp_v3.cpp")
source("/Users/cschaef3/Documents/GitHub/1930_final_project/project_functions.R")

obj <- function(S, b, u, lambda = 0.2) {
  ucol <- matrix(u, ncol = 1)
  0.5*t(ucol) %*% S %*% ucol - t(ucol) %*% matrix(b, ncol = 1) + lambda * sum(abs(ucol))
}

###### Scenario 1
# Test 1
set.seed(2023)
n <- 100
p <- 5
b <- rep(0, p)
b[1] <- 1
X <- matrix(rnorm(n * p), nrow = n, ncol = p) 
S <- t(X) %*% X / n + diag(0.01, p)
lambda <- 0.2

cd <- cg_cd(X, b, lambda)
admm <- cg_admm(X, b, lambda)
admm_res <- admm_r(X, b, lambda) 
cd_elastic <- cg_cd_alpha(X, b, lambda)
admm_elastic <- cg_admm_alpha(X, b, lambda)

print(list(cd_cpp=cd,admm_cpp=admm,admm_r =admm_res$u,cd_elastic=cd_elastic, admm_elastic=admm_elastic))

# Test 2
set.seed(2023)
n <- 100
p <- 5
b <- rep(0, p)
b[1] <- 1
X <- matrix(rnorm(n * p), nrow = n, ncol = p) 
S <- t(X) %*% X / n + diag(0.01, p)
lambda <- 0.4

cd <- cg_cd(X, b, lambda)
admm <- cg_admm(X, b, lambda)
elastic <- cg_cd_alpha(X, b, lambda)

print(list(cd=cd,admm=admm,elastic=elastic))

###### Scenario 2
# Test 1
set.seed(2023)
n <- 100
p <- 5
b <- rep(0, p)
b[1] <- 1
Sigma <- diag(0.2, p) + 0.8
X <- matrix(rnorm(n * p), nrow = n, ncol = p)%*% chol(Sigma)
S <- t(X) %*% X / n + diag(0.01, p)
lambda <- 0.2

cd <- cg_cd(X, b, lambda)
admm <- cg_admm(X, b, lambda)
elastic <- cg_cd_alpha(X, b, lambda)

print(list(cd=cd,admm=admm,elastic=elastic))