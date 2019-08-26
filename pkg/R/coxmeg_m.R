
#' Fit a Cox mixed-effects model for estimating HRs for a set of predictors
#'
#' \code{coxmeg_m} first estimates the variance component under a null model with only cov, and then analyzing each predictor in X one by one.
#' 
#' @section About \code{spd}:
#' When \code{spd=TRUE}, the relatedness matrix is treated as SPD. If the matrix is SPSD or not sure, use \code{spd=FALSE}.
#' @section About \code{solver}:
#' When solver=1/solver=2, Cholesky decompositon/PCG is used to solve the linear system. However, when \code{dense=FALSE} and \code{eigen=FALSE}, the solve function in the Matrix package is used regardless of \code{solver}. When \code{dense=TRUE}, it is recommended to set \code{solver=2} to have better computational performance.
#' @section About \code{invchol}:
#' Cholesky decomposition using \code{invchol=TRUE} is generally (but not always) much faster to invert a relatedness matrix (e.g., a block-diagonal matrix). But for some types of sparse matrices (e.g., a banded AR(1) matrix with rho=0.9), it sometimes can be very slow. In such cases, \code{invchol=FALSE} is can be used. 
#' 
#' @param outcome A matrix contains time (first column) and status (second column). The status is a binary variable (1 for failure / 0 for censored).
#' @param corr A relatedness matrix. Can be a matrix or a 'dgCMatrix' class in the Matrix package. Must be symmetric positive definite or symmetric positive semidefinite.
#' @param X A matrix of the preidctors. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the predictors are columns.
#' @param cov An optional matrix of the covariates included in the null model for estimating the variance component. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the covariates are columns. 
#' @param FID An optional string vector of family ID. If provided, the data will be reordered according to the family ID.
#' @param tau An optional positive value for the variance component. If \code{tau} is given, the function will skip estimating the variance component, and use the given \code{tau} to analyze the predictors.
#' @param eps An optional positive value indicating the tolerance in the optimization algorithm. Default is 1e-6.
#' @param min_tau An optional positive value indicating the lower bound in the optimization algorithm for the variance component \code{tau}. Default is 1e-4.
#' @param max_tau An optional positive value indicating the upper bound in the optimization algorithm for the variance component \code{tau}. Default is 5.
#' @param dense An optional logical value indicating whether the relatedness matrix is dense. Default is FALSE.
#' @param opt An optional logical value for the Optimization algorithm for tau. Can have the following values: 'bobyqa', 'Brent' or 'NM'. Default is 'bobyqa'.
#' @param spd An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. See details.
#' @param detap An optional logical value indicating whether to use approximation for log-determinant. Default is TRUE.
#' @param solver An optional bianry value taking either 1 or 2. Default is 1. See details.
#' @param score An optional logical value indicating whether to perform a score test. Default is FALSE.
#' @param order An optional integer value starting from 0. Only valid when \code{dense=FALSE}. It specifies the order of approximation used in the inexact newton method. Default is 1.
#' @param eigen An optional logical value. Only effective when \code{dense=FALSE}. It indicates whether to use RcppEigen:LDLT to solve linear systems. Default is TRUE.
#' @param verbose An optional logical value indicating whether to print additional messages. Default is TRUE.
#' @param threshold An optional non-negative value. If threshold>0, coxmeg_m will reestimate HRs for those SNPs with a p-value<threshold by first estimating a variant-specific variance component. Default is 0.
#' @param mc An optional integer value specifying the number of Monte Carlo samples used for approximating the log-determinant. Only valid when dense=TRUE and detap=TRUE. Default is 100.
#' @param invchol An optional logical value. Only effective when \code{dense=FALSE}. If TRUE, sparse Cholesky decomposition is used to compute the inverse of the relatedness matrix. Otherwise, sparse LU is used.
#' @return beta: The estimated coefficient for each predictor in X.
#' @return HR: The estimated HR for each predictor in X.
#' @return sd_beta: The estimated standard error of beta.
#' @return p: The p-value of each SNP.
#' @return tau: The estimated variance component.
#' @return iter: The number of iterations until convergence.
#' @keywords Cox mixed-effects model
#' @export coxmeg_m
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
#' ## simulate genotypes
#' g = matrix(rbinom(n*5,2,0.5),n,5)
#' 
#' ## The following command will first estimate the variance component without g, 
#' ## and then use it to estimate the HR for each preditor in g.
#' re = coxmeg_m(g,outcome,sigma,tau=0.5,order=1,eigen=TRUE,dense=FALSE)
#' re


coxmeg_m <- function(X,outcome,corr,FID=NULL,cov=NULL,tau=NULL,min_tau=1e-04,max_tau=5,eps=1e-6,order=1,detap=TRUE,opt='bobyqa',eigen=TRUE,score=FALSE,dense=FALSE,threshold=0,solver=1,spd=TRUE,verbose=TRUE,mc=100,invchol=TRUE){
  
  if(eps<0)
  {eps <- 1e-6}
  X = as.matrix(X)
  if(is.null(cov)==FALSE)
  {cov = as.matrix(cov)}
  
  ## family structure
  if(is.null(FID)==FALSE)
  {
    ord <- order(FID)
    FID <- as.character(FID[ord])
    X <- as.matrix(X[ord,])
    outcome <- as.matrix(outcome[ord,])
    corr <- corr[ord,ord]
    if(is.null(cov)==FALSE)
    {cov <- as.matrix(cov[ord,])}
  }else{
    X <- as.matrix(X)
    outcome <- as.matrix(outcome)
    if(is.null(cov)==FALSE)
    {cov <- as.matrix(cov)}
  }
  
  min_d <- min(outcome[which(outcome[,2]==1),1])
  rem <- which((outcome[,2]==0)&(outcome[,1]<min_d))
  if(length(rem)>0)
  {
    outcome <- outcome[-rem,,drop = FALSE]
    X <- as.matrix(X[-rem,,drop = FALSE])
    corr <- corr[-rem,-rem,drop = FALSE]
    if(is.null(cov)==FALSE)
    {cov <- as.matrix(cov[-rem,,drop = FALSE])}
  }
  
  if(verbose==TRUE)
  {message(paste0('Remove ', length(rem), ' subjects censored before the first failure.'))}
  
  n <- nrow(outcome)
  if(min(outcome[,2] %in% c(0,1))<1)
  {stop("The status should be either 0 (censored) or 1 (failure).")}
  
  p <- ncol(X)
  u <- rep(0,n)
  if(is.null(cov)==FALSE)
  {
    x_sd = which(as.vector(apply(cov,2,sd))>0)
    x_ind = length(x_sd)
    if(x_ind==0)
    {
      warning("The covariates are all constants after the removal of subjects.")
      k <- 0
      beta <- numeric(0)
      cov <- matrix(0,0,0)
    }else{
      k <- ncol(cov)
      if(x_ind<k)
      {
        warning(paste0(k-x_ind," covariate(s) is/are removed because they are all constants after the removal of subjects."))
        cov = cov[,x_sd,drop=FALSE]
        k <- ncol(cov)
      }
      beta <- rep(0,k)
    }
    
  }else{
    k <- 0
    beta <- numeric(0)
    cov <- matrix(0,0,0)
  }
  
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
    {message(paste0('There is/are ',k,' covariates. The sample size included is ',n,'. The rank of the relatedness matrix is ', rk_cor))}
  }else{
    spsd = FALSE
    rk_cor = n
    if(verbose==TRUE)
    {message(paste0('There is/are ',k,' covariates. The sample size included is ',n,'.'))}
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
      # corr = as.matrix(corr)
      rad = rbinom(n*mc,1,0.5)
      rad[rad==0] = -1
      rad = matrix(rad,n,mc)/sqrt(n)
    }
    si_d <- as.vector(diag(sigma_i_s))
  }else{
    if(verbose==TRUE)
    {message('The relatedness matrix is treated as sparse.')}
    corr <- as(corr, 'dgCMatrix')
    # if(eigen==FALSE)
    # {
    #  corr <- as(corr, 'symmetricMatrix')
    # }
    
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
    if(verbose==TRUE)
    {message('The relatedness matrix is inverted.')}
    
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
  
  if(is.null(tau))
  {
    tau_e = 0.5
    new_t = switch(
      opt,
      'bobyqa' = bobyqa(tau_e, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=cov, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,lower=min_tau,upper=max_tau,eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity),
      'Brent' = optim(tau_e, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=cov, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,lower=min_tau,upper=max_tau,method='Brent',eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity),
      'NM' = optim(tau_e, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=cov, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,method='Nelder-Mead',eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity),
      stop("The argument opt should be bobyqa, Brent or NM.")
    )
    
    if(opt=='bobyqa')
    {iter <- new_t$iter}else{
      iter <- new_t$counts
    }
    
    tau_e <- new_t$par
    if(tau_e==min_tau)
    {warning(paste0("The estimated variance component equals the lower bound (", min_tau, "), probably suggesting no random effects."))}
    
  }else{
    if(tau < 0)
    {stop("The variance component must be positive.")}
    tau_e = tau
  }
  
  snpval = which(apply(X,2,sd)>0)
  if(verbose==TRUE)
  {message(paste0('The variance component is estimated. Start analyzing SNPs...'))}
  
  if(score==FALSE)
  {
  
    c_ind <- c()
    if(k>0)
    {
      X <- cbind(X,cov)
      c_ind <- (p+1):(p+k)
    }

    beta <- rep(1e-16,k+1)
    u <- rep(0,n)
    sumstats <- data.frame(beta=rep(NA,p),HR=rep(NA,p),sd_beta=rep(NA,p),p=rep(NA,p))
    
    if(dense==TRUE)
    {
      cme_re <- sapply(snpval, function(i)
      {
        tryCatch({
        res <- irls_ex(beta, u, tau_e, si_d, sigma_i_s, as.matrix(X[,c(i,c_ind)]), eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=FALSE,sigma_s=NULL,s_d=NULL,eigen=eigen,solver=solver)
        c(res$beta[1],res$v11[1,1])},
        warning = function(war){
          message(paste0('The estimation may not converge for predictor ',i))
          c(NA,NA)
          },
        error = function(err){
          message(paste0('The estimation failed for predictor ',i))
          c(NA,NA)
          }
        )
      })
    }else{
      cme_re <- sapply(snpval, function(i)
      {
        tryCatch({
        res <- irls_fast_ap(beta, u, tau_e, si_d, sigma_i_s, as.matrix(X[,c(i,c_ind)]), eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order,det=FALSE,detap=detap,sigma_s=corr,s_d=s_d,eigen=eigen,solver=solver,sparsity=sparsity)
        c(res$beta[1],res$v11[1,1])
        },
        warning = function(war){
          message(paste0('The estimation may not converge for predictor ',i))
          c(NA,NA)
         },
        error = function(err){
          message(paste0('The estimation failed for predictor ',i))
          c(NA,NA)
          }
        )
      })
    }
    
    sumstats$beta[snpval] <- cme_re[1,]
    sumstats$HR[snpval] = exp(sumstats$beta[snpval])
    sumstats$sd_beta[snpval] <- sqrt(cme_re[2,])
    sumstats$p[snpval] <- pchisq(sumstats$beta[snpval]^2/cme_re[2,],1,lower.tail=FALSE)
    
    top = which(sumstats$p<threshold)
    if(length(top)>0)
    {
      tau = tau_e
      cme_re <- sapply(top, function(i)
      {
        if(opt=='bobyqa')
        {
          new_t <- bobyqa(tau, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=as.matrix(X[,c(i,c_ind)]), d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=TRUE,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,lower=min_tau,upper=max_tau,eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity)
        }else{
          if(opt=='Brent')
          {
            new_t <- optim(tau, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=as.matrix(X[,c(i,c_ind)]), d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=TRUE,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,lower=min_tau,upper=max_tau,method='Brent',eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity)
          }else{
            new_t <- optim(tau, mll, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=as.matrix(X[,c(i,c_ind)]), d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=TRUE,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,method='Nelder-Mead',eigen=eigen,dense=dense,solver=solver,rad=rad,sparsity=sparsity)
          }
        }
        tau_s <- new_t$par
        if(dense==TRUE)
        {
          res <- irls_ex(beta, u, tau_s, si_d, sigma_i_s, as.matrix(X[,c(i,c_ind)]), eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=FALSE,sigma_s=NULL,s_d=NULL,eigen=eigen,solver=solver)
        }else{
          res <- irls_fast_ap(beta, u, tau_s, si_d, sigma_i_s, as.matrix(X[,c(i,c_ind)]), eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order,det=FALSE,detap=detap,sigma_s=corr,s_d=s_d,eigen=eigen,solver=solver,sparsity=sparsity)
        }
        c(tau_s,res$beta[1],res$v11[1,1])
      })
      sumstats$tau = sumstats$beta_exact = sumstats$sd_beta_exact = sumstats$p_exact = NA
      sumstats$tau[top] = cme_re[1,]
      sumstats$beta_exact[top] = cme_re[2,]
      sumstats$sd_beta_exact[top] = sqrt(cme_re[3,])
      sumstats$p_exact[top] <- pchisq(sumstats$beta_exact[top]^2/cme_re[3,],1,lower.tail=FALSE)
    }
  }else{
    beta <- rep(1e-16,k)
    u <- rep(0,n)
    if(dense==FALSE)
    {
      model_n <- irls_fast_ap(beta, u, tau_e, si_d, sigma_i_s, cov, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order,det=FALSE,detap=detap,sigma_s=corr,s_d=s_d,eigen=eigen,solver=solver,sparsity=sparsity)
    }else{
      model_n <- irls_ex(beta, u, tau_e, si_d, sigma_i_s, cov, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,det=FALSE,detap=detap,sigma_s=NULL,s_d=NULL,eigen=eigen,solver=solver)
    }
    u <- model_n$u
    if(k>0)
    {
      gamma <- model_n$beta
      eta_v <- cov%*%gamma+u
    }else{
      eta_v <- u
    }
    
    w_v <- as.vector(exp(eta_v))
    s <- as.vector(cswei(w_v,rs$rs_rs-1,ind-1,1))
    a_v <- d_v/s
    a_v_p <- a_v[ind[,1]]
    a_v_2 <- as.vector(a_v_p*a_v_p)
    a_v_p <- a_v_p[a_v_p>0]
    b_v <- as.vector(cswei(a_v,rs$rs_cs-1,ind-1,0))
    bw_v <- as.vector(w_v)*as.vector(b_v)
    deriv <- d_v - bw_v
    
    dim_v = k + n
    brc = (k+1):dim_v
    v = matrix(NA,dim_v,dim_v)
    v[brc,brc] = -wma_cp(w_v,rs$rs_cs_p-1,ind-1,a_v_p)
    diag(v[brc,brc]) = diag(v[brc,brc]) + bw_v
    if(k>0)
    {
      temp = t(cov)%*%v[brc,brc]
      v[1:k,1:k] = temp%*%cov
      v[1:k,brc] = temp
      v[brc,1:k] = t(temp)
    }
    v[brc,brc] = v[brc,brc] + as.matrix(sigma_i_s/tau_e)
    v = chol(v)
    v = chol2inv(v)
    t_st <- score_test(deriv,bw_v,w_v,rs$rs_rs-1,rs$rs_cs-1,rs$rs_cs_p-1,ind-1,a_v_p,a_v_2,tau_e,v,cov,as.matrix(X[,snpval]))
    pv <- pchisq(t_st,1,lower.tail=FALSE)
    
    sumstats <- data.frame(score_test=rep(NA,p),p=rep(NA,p))
    sumstats$score_test[snpval]=t_st
    sumstats$p[snpval]=pv
  }
  res = list(summary=sumstats,tau=tau_e,rank=rk_cor,nsam=n)
  return(res)
}

