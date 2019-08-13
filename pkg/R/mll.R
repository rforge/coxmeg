mll <- function(tau, beta, u,si_d,sigma_i_s, X, d_v, ind, rs_rs, rs_cs,rs_cs_p,det,detap,inv,rk_cor,sigma_s=NULL,s_d=NULL,eps=1e-6,order=0,eigen=TRUE,dense=FALSE,solver=1,rad=NULL,sparsity=-1){
  
  if(dense==TRUE)
  {
    re <- irls_ex(beta, u, tau, si_d, as.matrix(sigma_i_s), X, eps, d_v, ind, rs_rs, rs_cs,rs_cs_p,det=det,detap=detap,sigma_s=sigma_s,s_d=NULL,eigen=eigen,solver=solver,rad=rad)
  }else{
    if(inv==TRUE){
      re <- irls_fast_ap(beta, u, tau, si_d, sigma_i_s, X, eps, d_v, ind, rs_rs, rs_cs,rs_cs_p,order,det=det,detap=detap,eigen=eigen,solver=solver,sparsity=sparsity)
    }else{
      re <- irls_fast_ap(beta, u, tau, si_d, sigma_i_s, X, eps, d_v, ind, rs_rs, rs_cs,rs_cs_p,order,det=det,detap=detap,sigma_s=sigma_s,s_d=s_d,eigen=eigen,solver=solver,sparsity=sparsity)
    } 
  }

  return(rk_cor*log(tau) + re$logdet[1] + (-2)*re$ll)
}
