\name{ADMM function}
\alias{cg_admm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ADMM algorithm
}
\description{
This function implements the ADMM algorithm that optimize the following problem in C++

\eqn{\frac{1}{2}u^TAu-ub^T+\lambda ||u||_1}

where

\eqn{A = \frac{1}{n}X^TX + 0.01I}

and \eqn{||u||_1 = \sum_{j}|u_j|} is the LASSO penalty.

}
\usage{
cg_admm(X, b, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    Data matrix of n x p dimension
  }

  \item{b}{
    Parameter vector of length p
  }

  \item{lambda}{
    Tuning parameter
  }
}
\details{
This ADMM function will perform with a preset tolerence = 1e-4 and max number of iteration = 1000.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Estimated vector u of length p.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yu Bin Chen, Rachit Sabharwal, Caroline Schaefer, Pagna Sok
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#set parameters for simulation
rho <- 0.9
n <- 200
p <- 10
#lambda value
lambda <- 0.01


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

#calling ADMM
cg_admm(X, b, lambda, tol=1e-4, maxit = 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
