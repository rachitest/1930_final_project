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

cd1 <- cg_cd(X, b, lambda)
admm1 <- cg_admm(X, b, lambda)
cd_elastic1 <- cg_cd_alpha(X, b, lambda, alpha=.05)
admm_elastic1 <- cg_admm_alpha(X, b, lambda, alpha=.05)

print(data.frame(cd_cpp=cd1,admm_cpp=admm1,cd_elastic=cd_elastic1, admm_elastic=admm_elastic1))


# Test 2
set.seed(2023)
n <- 100
p <- 5
b <- rep(0, p)
b[1] <- 1
X <- matrix(rnorm(n * p), nrow = n, ncol = p) 
S <- t(X) %*% X / n + diag(0.01, p)
lambda <- 0.4

cd2 <- cg_cd(X, b, lambda)
admm2 <- cg_admm(X, b, lambda)
cd_elastic2 <- cg_cd_alpha(X, b, lambda, alpha=.05)
admm_elastic2 <- cg_admm_alpha(X, b, lambda, alpha=.05)

print(data.frame(cd_cpp=cd2,admm_cpp=admm2,cd_elastic=cd_elastic2, admm_elastic=admm_elastic2))

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

cd3 <- cg_cd(X, b, lambda)
admm3 <- cg_admm(X, b, lambda)
cd_elastic3 <- cg_cd_alpha(X, b, lambda, alpha=.05)
admm_elastic3 <- cg_admm_alpha(X, b, lambda, alpha=.05)

print(data.frame(cd_cpp=cd3,admm_cpp=admm3,cd_elastic=cd_elastic3, admm_elastic=admm_elastic3))

# Test 2
set.seed(2023)
n <- 100
p <- 5
b <- rep(0, p)
b[1] <- 1
Sigma <- diag(0.2, p) + 0.8
X <- matrix(rnorm(n * p), nrow = n, ncol = p)%*% chol(Sigma)
S <- t(X) %*% X / n + diag(0.01, p)
lambda <- 0.4

cd4 <- cg_cd(X, b, lambda)
admm4 <- cg_admm(X, b, lambda)
cd_elastic4 <- cg_cd_alpha(X, b, lambda, alpha=.05)
admm_elastic4 <- cg_admm_alpha(X, b, lambda, alpha=.05)

print(data.frame(cd_cpp=cd4,admm_cpp=admm4,cd_elastic=cd_elastic4, admm_elastic=admm_elastic4))