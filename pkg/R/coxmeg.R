
#' Fit a Cox mixed-effects model
#'
#' \code{coxmeg} returns estimates of the variance component, the HRs and p-values for the predictors.
#' 
#' @section About \code{spd}:
#' When \code{spd=TRUE}, the relatedness matrix is treated as SPD. If the matrix is SPSD or not sure, use \code{spd=FALSE}.
#' @section About \code{solver}:
#' When solver=1/solver=2, Cholesky decompositon/PCG is used to solve the linear system. However, when \code{dense=FALSE} and \code{eigen=FALSE}, the solve function in the Matrix package is used regardless of \code{solver}. When \code{dense=TRUE}, it is recommended to set \code{solver=2} to have better computational performance. 
#' @section About \code{invchol}:
#' Cholesky decomposition using \code{invchol=TRUE} is generally (but not always) much faster to invert a relatedness matrix (e.g., a block-diagonal matrix). But for some types of sparse matrices (e.g., a banded AR(1) matrix with rho=0.9), it sometimes can be very slow. In such cases, \code{invchol=FALSE} is can be used. 
#' @param outcome A matrix contains time (first column) and status (second column). The status is a binary variable (1 for events / 0 for censored).
#' @param corr A relatedness matrix. Can be a matrix or a 'dgCMatrix' class in the Matrix package. Must be symmetric positive definite or symmetric positive semidefinite.
#' @param X An optional matrix of the preidctors with fixed effects. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the predictors are columns. 
#' @param FID An optional string vector of family ID. If provided, the data will be reordered according to the family ID.
#' @param eps An optional positive scalar indicating the tolerance in the optimization algorithm. Default is 1e-6.
#' @param min_tau An optional positive scalar indicating the lower bound in the optimization algorithm for the variance component \code{tau}. Default is 1e-4.
#' @param max_tau An optional positive scalar indicating the upper bound in the optimization algorithm for the variance component \code{tau} Default is 5.
#' @param dense An optional logical scalar indicating whether the relatedness matrix is dense. Default is FALSE.
#' @param opt An optional logical scalar for the Optimization algorithm for tau. Can have the following values: 'bobyqa', 'Brent' or 'NM'. Default is 'bobyqa'.
#' @param spd An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. See details.
#' @param detap An optional logical scalar indicating whether to use approximation for log-determinant. Default is TRUE.
#' @param solver An optional bianry scalar taking either 1 or 2. Default is 1. See details.
#' @param order An optional integer scalar starting from 0. Only valid when \code{dense=FALSE}. It specifies the order of approximation used in the inexact newton method. Default is 1.
#' @param eigen An optional logical scalar. Only effective when \code{dense=FALSE}. It indicates whether to use RcppEigen:LDLT to solve linear systems. Default is TRUE.
#' @param verbose An optional logical scalar indicating whether to print additional messages. Default is TRUE.
#' @param mc An optional integer scalar specifying the number of Monte Carlo samples used for approximating the log-determinant. Only valid when dense=TRUE and detap=TRUE. Default is 100.
#' @param invchol An optional logical value. Only effective when \code{dense=FALSE}. If TRUE, sparse Cholesky decomposition is used to compute the inverse of the relatedness matrix. Otherwise, sparse LU is used.
#' @return beta: The estimated coefficient for each predictor in X.
#' @return HR: The estimated HR for each predictor in X.
#' @return sd_beta: The estimated standard error of beta.
#' @return iter: The number of iterations until convergence.
#' @return tau: The estimated variance component.
#' @return int_ll: The marginal likelihood (-2*log(lik)) of tau evaluated at the estimate of tau.
#' @return rank: The rank of the relatedness matrix.
#' @return nsam: Actual sample size.
#' @keywords Cox mixed-effects model
#' @export coxmeg
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
#' ## simulate random effects and outcomes
#' x <- mvrnorm(1, rep(0,n), tau_var*sigma)
#' myrates <- exp(x-1)
#' y <- rexp(n, rate = myrates)
#' cen <- rexp(n, rate = 0.02 )
#' ycen <- pmin(y, cen)
#' outcome <- cbind(ycen,as.numeric(y <= cen))
#' 
#' ## fit a Cox mixed-effects model
#' re = coxmeg(outcome,sigma,order=1,eigen=TRUE,dense=FALSE)
#' re

coxmeg <- function(outcome,corr,X=NULL,FID=NULL,eps=1e-6, min_tau=1e-04,max_tau=5,order=1,detap=TRUE,opt='bobyqa',eigen=TRUE,dense=FALSE,solver=1,spd=TRUE,verbose=TRUE, mc=100,invchol=TRUE){
  
  if(eps<0)
  {eps <- 1e-6}
  
  ## family structure
  if(is.null(FID)==FALSE)
  {
    ord <- order(FID)
    FID <- as.character(FID[ord])
    X <- as.matrix(X[ord,])
    outcome <- as.matrix(outcome[ord,])
    corr <- corr[ord,ord]
  }else{
    if(is.null(X)==FALSE)
    {
      X <- as.matrix(X)
    }
    outcome <- as.matrix(outcome)
  }
  
  min_d <- min(outcome[which(outcome[,2]==1),1])
  rem <- which((outcome[,2]==0)&(outcome[,1]<min_d))
  if(length(rem)>0)
  {
    outcome <- outcome[-rem,,drop = FALSE]
    if(is.null(X)==FALSE)
    {
      X <- as.matrix(X[-rem,,drop = FALSE])
    }
    corr <- corr[-rem,-rem,drop = FALSE]
  }
  
  if(verbose==TRUE)
  {message(paste0('Remove ', length(rem), ' subjects censored before the first failure.'))}
  
  n <- nrow(outcome)
  if(min(outcome[,2] %in% c(0,1))<1)
  {stop("The status should be either 0 (censored) or 1 (failure).")}
  
  u <- rep(0,n)
  if(is.null(X)==FALSE)
  {
    x_sd = which(as.vector(apply(X,2,sd))>0)
    x_ind = length(x_sd)
    if(x_ind==0)
    {
      warning("The predictors are all constants after the removal of subjects.")
      k <- 0
      beta <- numeric(0)
      X <- matrix(0,0,0)
    }else{
      k <- ncol(X)
      if(x_ind<k)
      {
        warning(paste0(k-x_ind," predictor(s) is/are removed because they are all constants after the removal of subjects."))
        X = X[,x_sd,drop=FALSE]
        k <- ncol(X)
      }
      beta <- rep(0,k)
    }
  }else{
    k <- 0
    beta <- numeric(0)
    X <- matrix(0,0,0)
  }
  
  tau <- 0.5
  
  d_v <- outcome[,2]
  
  ## risk set matrix
  ind <- order(outcome[,1])
  ind <- as.matrix(cbind(ind,order(ind)))
  rk <- rank(outcome[ind[,1],1],ties.method='min')
  n1 <- sum(d_v>0)
  
  rs <- rs_sum(rk-1,d_v[ind[,1]])
  if(spd==FALSE)
  {
    if(max(colSums(corr))<1e-10)
    {
      stop("The relatedness matrix has a zero eigenvalue with an eigenvector of 1s.")
    }
    rk_cor = matrix.rank(as.matrix(corr),method='chol')
    spsd = FALSE
    if(rk_cor<n)
    {spsd = TRUE}
    
    if(verbose==TRUE)
    {message(paste0('There is/are ',k,' predictors. The sample size included is ',n,'. The rank of the relatedness matrix is ', rk_cor))}
    
  }else{
    spsd = FALSE
    rk_cor = n
    if(verbose==TRUE)
    {message(paste0('There is/are ',k,' predictors. The sample size included is ',n,'.'))}
  }
  
  nz <- nnzero(corr)
  sparsity = -1
  if(nz>(n*n/2))
  {
    dense <- TRUE
  }
  
  rad = NULL
  
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
    if(verbose==TRUE)
    {message('The relatedness matrix is inverted.')}
    inv <- TRUE
    sigma_i_s = corr
    corr = s_d = NULL
    if(detap==TRUE)
    {
      rad = rbinom(n*mc,1,0.5)
      rad[rad==0] = -1
      rad = matrix(rad,n,mc)/sqrt(n)
    }
    si_d <- as.vector(diag(sigma_i_s))
  }else{
    if(verbose==TRUE)
    {message('The relatedness matrix is treated as sparse.')}
    corr <- as(corr, 'dgCMatrix')
    #if(eigen==FALSE)
    #{
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
      if(detap==FALSE)
      {
        si_d <- as.vector(Matrix::diag(sigma_i_s))
      }
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
  }
  
  marg_ll = 0
  new_t = switch(
    opt,
    'bobyqa' = bobyqa(tau, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,lower=min_tau,upper=max_tau,eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity),
    'Brent' = optim(tau, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,lower=min_tau,upper=max_tau,method='Brent',eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity),
    'NM' = optim(tau, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,method='Nelder-Mead',eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity),
    stop("The argument opt should be bobyqa, Brent or NM.")
  )
  marg_ll = new_t$value
  if(opt=='bobyqa')
  {iter <- new_t$iter}else{
    iter <- new_t$counts
  }
  
  tau_e <- new_t$par
  if(tau_e==min_tau)
  {warning(paste0("The estimated variance component equals the lower bound (", min_tau, "), probably suggesting no random effects."))}
  
  if(dense==TRUE)
  {
    re <- irls_ex(beta, u, tau_e, si_d, sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=FALSE,sigma_s=NULL,s_d=NULL,eigen=eigen,solver=solver)
  }else{
    re <- irls_fast_ap(beta, u, tau_e, si_d, sigma_i_s, X, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order,det=FALSE,detap=detap,sigma_s=corr,s_d=s_d,eigen=eigen,solver=solver,sparsity=sparsity)
  }

  if(k>0)
  {
    HR = exp(re$beta)
    sdb = sqrt(diag(as.matrix(re$v11)))
    p <- pchisq(re$beta^2/diag(re$v11),1,lower.tail=FALSE)
  }else{
    HR=sdb=p=NULL
  }
  
  res <- list(beta=re$beta,HR=HR,sd_beta=sdb,p=p,tau=tau_e,iter=iter,rank=rk_cor,nsam=n,int_ll=marg_ll)
  return(res)
}

