% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxmeg_m.R
\name{coxmeg_m}
\alias{coxmeg_m}
\title{Fit a Cox mixed-effects model for estimating HRs for a set of predictors}
\usage{
coxmeg_m(X, outcome, corr, type, FID = NULL, cov = NULL, tau = NULL,
  min_tau = 1e-04, max_tau = 5, eps = 1e-06, order = NULL,
  detap = NULL, opt = "bobyqa", score = FALSE, threshold = 0,
  solver = NULL, spd = TRUE, verbose = TRUE, mc = 100)
}
\arguments{
\item{X}{A matrix of the preidctors. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the predictors are columns.}

\item{outcome}{A matrix contains time (first column) and status (second column). The status is a binary variable (1 for failure / 0 for censored).}

\item{corr}{A relatedness matrix. Can be a matrix or a 'dgCMatrix' class in the Matrix package. Must be symmetric positive definite or symmetric positive semidefinite.}

\item{type}{A string indicating the sparsity structure of the relatedness matrix. Should be 'bd' (block diagonal), 'sparse', or 'dense'. See details.}

\item{FID}{An optional string vector of family ID. If provided, the data will be reordered according to the family ID.}

\item{cov}{An optional matrix of the covariates included in the null model for estimating the variance component. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the covariates are columns.}

\item{tau}{An optional positive value for the variance component. If \code{tau} is given, the function will skip estimating the variance component, and use the given \code{tau} to analyze the predictors.}

\item{min_tau}{An optional positive value indicating the lower bound in the optimization algorithm for the variance component \code{tau}. Default is 1e-4.}

\item{max_tau}{An optional positive value indicating the upper bound in the optimization algorithm for the variance component \code{tau}. Default is 5.}

\item{eps}{An optional positive value indicating the tolerance in the optimization algorithm. Default is 1e-6.}

\item{order}{An optional integer value starting from 0. Only valid when \code{dense=FALSE}. It specifies the order of approximation used in the inexact newton method. Default is NULL, which lets coxmeg choose an optimal order.}

\item{detap}{An optional string indicating whether to use approximation for log-determinant. Can be 'exact', 'diagonal' or 'slq'. Default is NULL, which lets the function select a method based on 'type' and other information. See details.}

\item{opt}{An optional logical value for the Optimization algorithm for tau. Can have the following values: 'bobyqa', 'Brent' or 'NM'. Default is 'bobyqa'.}

\item{score}{An optional logical value indicating whether to perform a score test. Default is FALSE.}

\item{threshold}{An optional non-negative value. If threshold>0, coxmeg_m will reestimate HRs for those SNPs with a p-value<threshold by first estimating a variant-specific variance component. Default is 0.}

\item{solver}{An optional bianry value that can be either 1 (Cholesky Decomposition using RcppEigen), 2 (PCG) or 3 (Cholesky Decomposition using Matrix). Default is NULL, which lets the function select a solver. See details.}

\item{spd}{An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. See details.}

\item{verbose}{An optional logical value indicating whether to print additional messages. Default is TRUE.}

\item{mc}{An optional integer value specifying the number of Monte Carlo samples used for approximating the log-determinant. Only valid when \code{dense=TRUE} and \code{detap='slq'}. Default is 100.}
}
\value{
beta: The estimated coefficient for each predictor in X.

HR: The estimated HR for each predictor in X.

sd_beta: The estimated standard error of beta.

p: The p-value of each SNP.

tau: The estimated variance component.

iter: The number of iterations until convergence.
}
\description{
\code{coxmeg_m} first estimates the variance component under a null model with only cov, and then analyzing each predictor in X one by one.
}
\section{About \code{type}}{

'bd' is used for a block-diagonal relatedness matrix, or a sparse matrix the inverse of which is also sparse. 'sparse' is used for a general sparse relatedness matrix the inverse of which is not sparse.
}

\section{About \code{spd}}{

When \code{spd=TRUE}, the relatedness matrix is treated as SPD. If the matrix is SPSD or not sure, use \code{spd=FALSE}.
}

\section{About \code{solver}}{

When \code{solver=1,3}/\code{solver=2}, Cholesky decompositon/PCG is used to solve the linear system. When \code{solver=3}, the solve function in the Matrix package is used, and when \code{solver=1}, it uses RcppEigen:LDLT to solve linear systems. When \code{type='dense'}, it is recommended to set \code{solver=2} to have better computational performance.
}

\section{About \code{detap}}{

When \code{detap='exact'}, the exact log-determinant is computed for estimating the variance component. Specifying \code{detap='diagonal'} uses diagonal approximation, and is only effective for a sparse relatedness matrix. Specifying \code{detap='slq'} uses stochastic lanczos quadrature approximation.
}

\examples{
library(Matrix)
library(MASS)
library(coxmeg)

## simulate a block-diagonal relatedness matrix
tau_var <- 0.2
n_f <- 100
mat_list <- list()
size <- rep(10,n_f)
offd <- 0.5
for(i in 1:n_f)
{
  mat_list[[i]] <- matrix(offd,size[i],size[i])
  diag(mat_list[[i]]) <- 1
}
sigma <- as.matrix(bdiag(mat_list))
n <- nrow(sigma)

## simulate random effects and outcomes
x <- mvrnorm(1, rep(0,n), tau_var*sigma)
myrates <- exp(x-1)
y <- rexp(n, rate = myrates)
cen <- rexp(n, rate = 0.02 )
ycen <- pmin(y, cen)
outcome <- cbind(ycen,as.numeric(y <= cen))

## simulate genotypes
g = matrix(rbinom(n*5,2,0.5),n,5)

## The following command will first estimate the variance component without g, 
## and then use it to estimate the HR for each preditor in g.
re = coxmeg_m(g,outcome,sigma,type='bd',tau=0.5,detap='diagonal',order=1)
re
}
\keyword{Cox}
\keyword{mixed-effects}
\keyword{model}
