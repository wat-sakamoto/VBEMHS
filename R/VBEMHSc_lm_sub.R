### Subroutines for VBEMHSc_lm (conjugate prior)

# Compute the lower bound for log marginal posterior
VBEMHSc_lm.lb <- function(n,p,sigma2,tau2,q_beta,omega,tr_S,s2_estim,
                          tau_log_prior,a=0.5,b=0.5,c_ig=0,d=0)
{
  lb <- switch(s2_estim, 
               mode = (n+2*c_ig+2)*log(sigma2)+(n-tr_S+2*c_ig+2), # marginal posterior mode
               mean = (n+2*c_ig)*log(sigma2)+(n-tr_S+2*c_ig)    # marginal posterior mean
  )      
  lb <- lb + p*log(tau2) 
  lb <- -lb/2 

  lb <- lb - q_beta$logdet 
  lb <- lb + (b+0.5)*(sum(log(omega))+p*(1-log(b+0.5))) + p*(lgamma(b+0.5)-lgamma(b))
  lb <- lb - (a+b)*sum(log(omega+1)) + p*(lgamma(a+b)-lgamma(a))
  lb <- lb + tau_log_prior(tau2)
  
  return(lb)
}  

# Initial search for the global shrinkage parameter (gs=1/tau2)
VBEMHSc_lm.gs <- function(X,y,lgs_range=c(-2,4),lgs_step=1,s2_estim="mean",s2_correct=FALSE,
                         tau2_scale=1,a=0.5,b=0.5,c_ig=0,d=0,tol=1e-6,verbose=FALSE)
{
  # Marginal posterior estimate for sigma2
  sigma2_mpe <- function(n,X,y,beta,tr_S,s2_estim,c_ig,d)
  {
    SS <- as.vector(crossprod(y, y-X %*% beta))
    sigma2 <- switch(s2_estim, 
                     mode = (SS+2*d)/(n-tr_S+2*c_ig+2), # marginal posterior mode
                     mean = (SS+2*d)/(n-tr_S+2*c_ig)    # marginal posterior mean
    )
    return(sigma2)
  }

  # Prior for tau2
  if (tau2_scale>0){
    tau_log_prior <- function(tau2){ -log(1+tau2/tau2_scale) } # half_Cauchy
  } else {
    tau_log_prior <- function(tau2){ -log(tau2) } # non-informative (improper)
  }
  
  n <- length(y)
  p <- ncol(X)
  XX <- crossprod(X)
  Xy <- crossprod(X,y)
  omega <- rep(1,p)
  
  lb_fun <- function(gs){
    prec <- gs*omega
    q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
    beta <- q_beta$mean
    var_b <- q_beta$var
    if (s2_correct) tr_S <- p-sum(var_b*prec) else tr_S <- 0
    sigma2 <- sigma2_mpe(n,X,y,beta,tr_S,s2_estim,c_ig,d)
    tau2 <- 1/gs
    Eb2 <- sigma2*var_b+beta^2
    lb <- VBEMHSc_lm.lb(n,p,sigma2,tau2,q_beta,omega,tr_S,s2_estim,
                       tau_log_prior,a,b,c_ig,d)
    return(lb)
  }
  
  lb_seq <- NULL
  lgs_seq <- seq(lgs_range[1],lgs_range[2],by=lgs_step)
  
  for (l in seq(along.with=lgs_seq)){
    gs <- 10^lgs_seq[l]
    lb <- lb_fun(gs)
    lb_seq <- c(lb_seq, lb)
    if(verbose) {
      print(paste("k=", 0, 
                  " sigma2/tau2=",sprintf("%f", gs),
                  " lower bound=",round(lb, 5)
      ))
    }
  }
  lgs_opt <- lgs_seq[which.max(lb_seq)]
  
  lb_opt <- optimize(
    f=function(x) lb_fun(10^x), 
    lower=max(lgs_opt-lgs_step, min(lgs_range)),
    upper=min(lgs_opt+lgs_step, max(lgs_range)),
    maximum=TRUE, tol=tol)
  
  gs_opt <- 10^lb_opt$maximum
  lb <- lb_opt$objective
  if(verbose) {
    print(paste("k=", 0, 
                " sigma2/tau2=",sprintf("%f", gs_opt),
                " lower bound=",round(lb, 5)
    ))
  }
  return(list(gs_opt=gs_opt, lb=lb))
}

# VBEM algorithm (to update the global and local shrinkage parameters)
VBEMHSc_lm.update <- function(X,y,gs_ini,gs_fix=FALSE,ls_fix=FALSE,s2_estim="mean",s2_correct=FALSE,
                              tau2_scale=1,a=0.5,b=0.5,c_ig=0,d=0,max_iter=1000,tol=1e-6,verbose=FALSE)
{
  
  # Marginal posterior estimate for sigma2
  sigma2_mpe <- function(n,X,y,beta,tr_S,s2_estim,c_ig,d)
  {
    SS <- as.vector(crossprod(y, y-X %*% beta))
    sigma2 <- switch(s2_estim, 
                     mode = (SS+2*d)/(n-tr_S+2*c_ig+2), # marginal posterior mode
                     mean = (SS+2*d)/(n-tr_S+2*c_ig)    # marginal posterior mean
    )
    return(sigma2)
  }
  
  # VBEM update for sigma2
  sigma2_mpe.update <- function(n,p,X,y,beta,sigma2_old,tr_S,s2_estim,c_ig,d)
  {
    SS <- as.vector(crossprod(y, y-X %*% beta))
    sigma2 <- switch(s2_estim, 
                     mode = (SS+sigma2_old*p+2*d)/(n+p-tr_S+2*c_ig+2), # marginal posterior mode
                     mean = (SS+sigma2_old*p+2*d)/(n+p-tr_S+2*c_ig)    # marginal posterior mean
    )
    return(sigma2)
  }
  
  # VBEM update for tau2
  tau2_update <- function(p,sigma2,Eob2,tau2_scale=1)
  {
    if (tau2_scale>0){
      A <- (p+2)/tau2_scale
      B <- p-Eob2/tau2_scale/sigma2
      D <- (p+Eob2/tau2_scale/sigma2)^2+8*Eob2/tau2_scale/sigma2
      tau2 <- (-B+sqrt(D))/A*0.5
    } else {
      tau2 <- Eob2/sigma2/(p+2)
    }
    return(tau2)
  }
  
  # Prior for tau2
  if (tau2_scale>0){
    tau_log_prior <- function(tau2){ -log(1+tau2/tau2_scale) } # half_Cauchy
  } else {
    tau_log_prior <- function(tau2){ -log(tau2) } # non-informative (improper)
  }
  
  conv <- function(par_old,par_new,tol) {
    abs(par_new-par_old)<tol 
  }
  
  if(verbose) {
    if (gs_fix & !ls_fix) print("Updating the local shrinkage parameters ...")
    else if (!gs_fix & ls_fix) print("Updating the global shrinkage parameters ...")
    else print("Updating the global and local shrinkage parameters ...")
  }
  
  n <- length(y)
  p <- ncol(X)
  XX <- crossprod(X)
  Xy <- crossprod(X,y)
  omega <- rep((b/a),p)

  # Initialize beta
  if (missing(gs_ini)) gs_ini <- n
  prec <- gs_ini*omega
  q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
  beta <- q_beta$mean
  var_b <- q_beta$var
  
  # Initialize sigma2 and tau2
  if (s2_correct) tr_S <- p-sum(var_b*prec) else tr_S <- 0
  sigma2 <- sigma2_mpe(n,X,y,beta,tr_S,s2_estim,c_ig,d)
  Eb2 <- as.vector(sigma2*var_b+beta^2)       # E[beta^2]
  tau2 <- 1/gs_ini
  
  apb <- a+b
  beta_old <- beta
  sigma2_old <- sigma2
  tau2_old <- tau2

  for (k in 1:max_iter){
    
    if (!ls_fix){
      # Update q(nu) and q(omega)
      nu <- apb/(omega+1)
      omega <- (b+0.5)/(Eb2/(2*sigma2_old*tau2_old)+nu)
    }
    
    # Update q(beta) 
    prec <- omega/tau2_old
    q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
    beta <- q_beta$mean
    var_b <- q_beta$var
    Eb2 <- as.vector(sigma2_old*var_b+beta^2) # E[beta^2]
    
    # Update sigma2
    if (s2_correct) tr_S <- p-sum(var_b*prec) else tr_S <- 0
    sigma2 <- sigma2_mpe.update(n,p,X,y,beta,sigma2_old,tr_S,s2_estim,c_ig,d)
    sigma2 <- max(sigma2,tol)

    # Update tau2
    if (!gs_fix) {
      Eob2 <- sum(omega*Eb2)         # sum of E[omega]E[beta^2]
      tau2 <- tau2_update(p,sigma2,Eob2,tau2_scale)
    }
    
    if(verbose) {
      lb <- VBEMHSc_lm.lb(n,p,sigma2,tau2,q_beta,omega,tr_S,s2_estim,tau_log_prior,a,b,c_ig,d)
      print(paste("k=", k, 
                  " tau2=",sprintf("%f", tau2), 
                  " sigma2=",sprintf("%f", sigma2), 
                  " lower bound=",round(lb, 5)
      ))
    }
    
    if (all(
      conv(beta_old,beta,tol),
      conv(sigma2_old,sigma2,tol),
      conv(log10(tau2_old),log10(tau2),tol)
    )) break
    
    beta_old <- beta
    sigma2_old <- sigma2
    tau2_old <- tau2
  }
  if (k==max_iter){
    print("Update stopped because the number of iterations exceeded max_iter.")
  }
  
  if (s2_correct) {
    tr_S <- p-sum(var_b*prec)
  } else {
    tr_S <- 0
  }
  lb <- VBEMHSc_lm.lb(n,p,sigma2,tau2,q_beta,omega,tr_S,s2_estim,tau_log_prior,a,b,c_ig,d)
  df_s2 <- 2*c_ig+n-tr_S
  
  list( 
    beta=beta, se_beta=sqrt(var_b*sigma2), V_beta=q_beta$V*sigma2, 
    omega=omega, sigma2=sigma2, tau2=tau2, df_s2=df_s2, lb=lb, iter=k
  )
}


# VBEM for several fixed values of the global shrinkage parameter
VBEMHSc_lm.trace <- function(X,y,lgs_range,lgs_step,ls_fix=FALSE,s2_estim="mean",s2_correct=FALSE,
                            tau2_scale=1,a=0.5,b=0.5,c_ig=0,d=0,max_iter=1000,tol=1e-6,verbose=FALSE)
{
  lgs_seq <- seq(lgs_range[1],lgs_range[2],by=lgs_step)

  beta_mat <- NULL
  se_beta_mat <- NULL
  omega_mat <- NULL
  sigma2_seq <- NULL
  tau2_seq <- NULL
  df_s2_seq <- NULL
  lb_seq <- NULL
  iter_seq <- NULL
  
  for (k in seq(along.with=lgs_seq)){
    gs <- 10^lgs_seq[k]
    res <- VBEMHSc_lm.update(X,y,gs,gs_fix=TRUE,ls_fix,s2_estim,s2_correct,tau2_scale,a,b,c_ig,d,max_iter,tol)
    beta <- as.vector(res$beta)
    se_beta <- as.vector(res$se_beta)
    omega <- as.vector(res$omega)
    beta_mat <- rbind(beta_mat,beta)
    se_beta_mat <- rbind(se_beta_mat,se_beta)
    omega_mat <- rbind(omega_mat,omega)
    sigma2_seq <- c(sigma2_seq,res$sigma2)
    tau2_seq <- c(tau2_seq,res$tau2)
    df_s2_seq <- c(df_s2_seq,res$df)
    lb_seq <- c(lb_seq,res$lb)
    iter_seq <- c(iter_seq,res$iter)
    if(verbose) {
      print(paste("tau2=",sprintf("%f", res$tau2), 
                  " sigma2=",sprintf("%f", res$sigma2), 
                  " lower bound=",round(res$lb, 5)
      ))
    }
  }
  row.names(beta_mat) <- NULL
  row.names(se_beta_mat) <- NULL
  row.names(omega_mat) <- NULL

  list(beta=beta_mat,se_beta=se_beta_mat,omega=omega_mat,sigma2=sigma2_seq,tau2=tau2_seq,
       df_s2=df_s2_seq,lb=lb_seq,iter=iter_seq)
}
