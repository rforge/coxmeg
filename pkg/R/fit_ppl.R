
#' Estimate HRs using PPL given a known variance component (tau)
#'
#' \code{fit_ppl} returns estimates of HRs and their p-values given a known variance component (tau).
#' 
#' @section About \code{solver}:
#' When solver=1/solver=2, Cholesky decompositon/PCG is used to solve the linear system. However, when \code{dense=FALSE} and \code{eigen=FALSE}, the solve function in the Matrix package is used regardless of \code{solver}. 
#' @section About \code{invchol}:
#' Cholesky decomposition using \code{invchol=TRUE} is generally (but not always) much faster to invert a relatedness matrix (e.g., a block-diagonal matrix). But for some types of sparse matrices (e.g., a banded AR(1) matrix with rho=0.9), it sometimes can be very slow. In such cases, \code{invchol=FALSE} is can be used. 
#' 
#' @param tau A positive scalar. A variance component given by the user. Default is 0.5.
#' @param X A matrix of the preidctors. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the predictors are columns.
#' @param outcome A matrix contains time (first column) and status (second column). The status is a binary variable (1 for failure / 0 for censored).
#' @param corr A relatedness matrix. Can be a matrix or a 'dgCMatrix' class in the Matrix package. Must be symmetric positive definite or symmetric positive semidefinite.
#' @param FID An optional string vector of family ID. If provided, the data will be reordered according to the family ID.
#' @param eps An optional positive value indicating the tolerance in the optimization algorithm. Default is 1e-6.
#' @param dense An optional logical value indicating whether the relatedness matrix is dense. Default is FALSE.
#' @param spd An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. 
#' @param solver An optional bianry value taking either 1 or 2. Default is 1. See details.
#' @param verbose An optional logical value indicating whether to print additional messages. Default is TRUE.
#' @param order An optional integer value starting from 0. Only valid when dense=FALSE. It specifies the order of approximation used in the inexact newton method. Default is 1.
#' @param eigen An optional logical value Only effective when \code{dense=FALSE}. It indicates whether to use RcppEigen:LDLT to solve linear systems. Default is TRUE.
#' @param invchol An optional logical value. Only effective when \code{dense=FALSE}. If TRUE, sparse Cholesky decomposition is used to compute the inverse of the relatedness matrix. Otherwise, sparse LU is used.
#' @return beta: The estimated coefficient for each predictor in X.
#' @return HR: The estimated HR for each predictor in X.
#' @return sd_beta: The estimated standard error of beta.
#' @return iter: The number of iterations until convergence.
#' @return ppl: The PPL when the convergence is reached.
#' @keywords Cox mixed-effects model
#' @export fit_ppl
#' @examples
#' library(Matrix)
#' library(MASS)
#' library(coxmeg)
#' 
#' ## simulate a block-diagonal relatedness matrix
#' tau_var <- 0.2
#' n_f <- 100
#' mat_list <- list()
#' size <- rep(10,n_f)
#' offd <- 0.5
#' for(i in 1:n_f)
#' {
#'   mat_list[[i]] <- matrix(offd,size[i],size[i])
#'   diag(mat_list[[i]]) <- 1
#' }
#' sigma <- as.matrix(bdiag(mat_list))
#' n <- nrow(sigma)
#' 
#' ## simulate random effexts and outcomes
#' x <- mvrnorm(1, rep(0,n), tau_var*sigma)
#' myrates <- exp(x-1)
#' y <- rexp(n, rate = myrates)
#' cen <- rexp(n, rate = 0.02 )
#' ycen <- pmin(y, cen)
#' outcome <- cbind(ycen,as.numeric(y <= cen))
#' 
#' ## fit the ppl
#' re = fit_ppl(x,outcome,sigma,tau=0.5,order=1,eigen=TRUE,dense=FALSE)
#' re

fit_ppl <- function(X,outcome,corr,tau=0.5,FID=NULL,eps=1e-06,order=1,eigen=TRUE,dense=FALSE,solver=1,spd=TRUE,verbose=TRUE,invchol=TRUE){

  if(eps<0)
  {eps <- 1e-06}
  
  ## family structure
  if(is.null(FID)==FALSE)
  {
    ord <- order(FID)
    FID <- as.character(FID[ord])
    X <- as.matrix(X[ord,,drop = FALSE])
    outcome <- as.matrix(outcome[ord,,drop = FALSE])
    corr <- corr[ord,ord,drop = FALSE]
  }else{
    X <- as.matrix(X)
    outcome <- as.matrix(outcome)
  }
  
  min_d <- min(outcome[which(outcome[,2]==1),1])
  rem <- which((outcome[,2]==0)&(outcome[,1]<min_d))
  if(length(rem)>0)
  {
    outcome <- outcome[-rem, ,drop = FALSE]
    X <- as.matrix(X[-rem,,drop = FALSE])
    corr <- corr[-rem,-rem,drop = FALSE]
  }
  
  if(verbose==TRUE)
  {message(paste0('Remove ', length(rem), ' subjects censored before the first failure.'))}
  
  x_sd = which(as.vector(apply(X,2,sd))>0)
  x_ind = length(x_sd)
  if(x_ind==0)
  {stop("The predictors are all constants after the removal of subjects.")}else{
    k <- ncol(X)
    if(x_ind<k)
    {
      warning(paste0(k-x_ind," predictor(s) is/are removed because they are all constants after the removal of subjects."))
      X = X[,x_sd,drop=FALSE]
      k <- ncol(X)
    }
  }
  
  n <- nrow(outcome)
  if(min(outcome[,2] %in% c(0,1))<1)
  {stop("The status should be either 0 (censored) or 1 (failure).")}
  u <- rep(0,n)
  beta <- rep(0,k)
  
  d_v <- outcome[,2]
  
  ## risk set matrix
  ind <- order(outcome[,1])
  ind <- as.matrix(cbind(ind,order(ind)))
  rk <- rank(outcome[ind[,1],1],ties.method='min')
  n1 <- sum(d_v>0)
  
  rs <- rs_sum(rk-1,d_v[ind[,1]])
  if(spd==FALSE)
  {
    # minei = eigs_sym(corr, 1, which = "SM")
    # if(minei$values < -1e-10)
    # {
    #  stop(paste0("The relatedness matrix has negative eigenvalues (", minei$values,")."))
    # }
    if(max(colSums(corr))<1e-10)
    {
      stop("The relatedness matrix has a zero eigenvalue with an eigenvector of 1s.")
    }
    rk_cor = matrix.rank(as.matrix(corr),method='chol')
    spsd = FALSE
    if(rk_cor<n)
    {spsd = TRUE}
    if(verbose==TRUE)
    {message(paste0('The sample size included is ',n,'. The rank of the relatedness matrix is ', rk_cor))}
    
  }else{
    spsd = FALSE
    rk_cor = n
    if(verbose==TRUE)
    {message(paste0('The sample size included is ',n,'.'))}
  }
  
  nz <- nnzero(corr)
  sparsity = -1
  if(nz>(n*n/2))
  {dense <- TRUE}
  
  if(dense==TRUE)
  {
    if(verbose==TRUE)
    {message('The relatedness matrix is treated as dense.')}
    corr = as.matrix(corr)
    if(spsd==FALSE)
    {
      corr = chol(corr)
      corr = as.matrix(chol2inv(corr))
    }else{
      corr <- ginv(corr)
    }
    inv <- TRUE
    sigma_i_s = corr
    corr <- s_d <- NULL
    si_d <- as.vector(diag(sigma_i_s))
    res <- irls_ex(beta, u, tau, si_d, sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=FALSE,sigma_s=NULL,s_d=NULL,eigen=eigen,solver=solver)
  }else{
    
    if(verbose==TRUE)
    {message('The covariance matrix is treated as sparse.')}
    
    corr <- as(corr, 'dgCMatrix')
    # if(eigen==FALSE)
    # {
    #  corr <- as(corr, 'symmetricMatrix')
    #}
  
    if(spsd==FALSE)
    {
      if(invchol==TRUE)
      {
        sigma_i_s <- Matrix::chol2inv(Matrix::chol(corr))
      }else{
        sigma_i_s <- Matrix::solve(corr)
      }
    }else{
      sigma_i_s = eigen(corr)
      if(min(sigma_i_s$values) < -1e-10)
      {
        stop("The relatedness matrix has negative eigenvalues.")
      }
      sigma_i_s = sigma_i_s$vectors%*%(c(1/sigma_i_s$values[1:rk_cor],rep(0,n-rk_cor))*t(sigma_i_s$vectors))
    }
    
    nz_i <- nnzero(sigma_i_s)
    si_d <- NULL
    sparsity = nz_i/nz
    if((nz_i>nz)&&(spsd==FALSE))
    {
      inv <- FALSE
      s_d <- as.vector(Matrix::diag(corr))
    }else{
      inv <- TRUE
      si_d <- as.vector(Matrix::diag(sigma_i_s))
    }
    
    sigma_i_s <- as(sigma_i_s,'dgCMatrix')
    if(eigen==FALSE)
    {
      sigma_i_s = Matrix::forceSymmetric(sigma_i_s)
    }
    
    if(inv==TRUE)
    {corr <- s_d <- NULL}
    res <- irls_fast_ap(beta, u, tau, si_d, sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,order,det=FALSE,detap=TRUE,sigma_s=corr,s_d=s_d,eigen=eigen,solver=solver,sparsity=sparsity)
    
  }
  
  re = list(beta=res$beta,HR=exp(res$beta),sd_beta=sqrt(diag(as.matrix(res$v11))),iter=res$iter,ppl=res$ll)
  return(re)
}

