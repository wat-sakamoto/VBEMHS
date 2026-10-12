# VBEMHS_lm_grid: Variational Bayesian EM Approach to Variable Selection 
# in linear regression model with (independent/conjugate) Horseshoe Prior
# Notice: no intercept term is included, so X and y are supposed to be centered.

# Compute the lower bound for log marginal posterior : Independent prior
VBEMHS_lm_grid.lb <- function(n,p,X,y,sigma2,tau2,q_beta,omega,tau_log_prior,
                              a=0.5,b=0.5,c_ig=0,d=0)
{
  beta <- q_beta$mean
  var_b <- q_beta$var
  SS <- as.vector(crossprod(y-X %*% beta)+2*d)

  lb <- (n-p+2*c_ig+2)*log(sigma2) + p*log(tau2)
  lb <- lb + SS/sigma2 + sum(omega*beta^2)/tau2
  lb <- -lb/2 
  
  lb <- lb - q_beta$logdet 
  lb <- lb + (b+0.5)*(sum(log(omega))+p*(1-log(b+0.5))) + p*(lgamma(b+0.5)-lgamma(b))
  lb <- lb - (a+b)*sum(log(omega+1)) + p*(lgamma(a+b)-lgamma(a))
  lb <- lb + tau_log_prior(tau2)
  
  return(lb)
} 

# Compute the lower bound for log marginal posterior : Conjugate prior
VBEMHSc_lm_grid.lb <- function(n,p,X,y,sigma2,tau2,q_beta,omega,tau_log_prior,
                              a=0.5,b=0.5,c_ig=0,d=0)
{
  beta <- q_beta$mean
  var_b <- q_beta$var
  SS <- as.vector(crossprod(y, y-X %*% beta)+2*d)

  lb <- (n+2*c_ig+2)*log(sigma2) + p*log(tau2) + SS/sigma2
  lb <- -lb/2 
  
  lb <- lb - q_beta$logdet 
  lb <- lb + (b+0.5)*(sum(log(omega))+p*(1-log(b+0.5))) + p*(lgamma(b+0.5)-lgamma(b))
  lb <- lb - (a+b)*sum(log(omega+1)) + p*(lgamma(a+b)-lgamma(a))
  lb <- lb + tau_log_prior(tau2)
  
  return(lb)
} 

# VBEM algorithm (sigma2 and tau2 fixed) : Independent prior
VBEMHS_lm_grid.fixed <- function(X,y,sigma2,tau2,ls_fix=FALSE,tau2_scale=1,
                          a=0.5,b=0.5,c_ig=0,d=0,max_iter=1000,tol=1e-6)
{
  if (tau2_scale>0){
    tau_log_prior <- function(tau2){ -log(1+tau2/tau2_scale) } # half_Cauchy
  } else {
    tau_log_prior <- function(tau2){ -log(tau2) } # non-informative (improper)
  }
  
  conv <- function(par_old,par_new,tol) {
    abs(par_new-par_old)<tol 
  }
  
  n <- length(y)
  p <- ncol(X)
  XX <- crossprod(X)
  Xy <- crossprod(X,y)
  omega <- rep((b/a),p)

  # Initialize q(beta) 
  prec <- omega*(sigma2/tau2)
  q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
  beta <- q_beta$mean
  var_b <- q_beta$var
  Eb2 <- as.vector(sigma2*var_b+beta^2)  # E[beta^2]

  apb <- a+b
  beta_old <- beta

  for (k in 1:max_iter){

    if (!ls_fix){
      # Update q(nu) and q(omega)
      nu <- apb/(omega+1)
      omega <- (b+0.5)/(Eb2/(2*tau2)+nu)
    }

    # Update q(beta) 
    prec <- omega*(sigma2/tau2)
    q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
    beta <- q_beta$mean
    var_b <- q_beta$var
    Eb2 <- as.vector(sigma2*var_b+beta^2) # E[beta^2]
    
    if (all(
      conv(beta_old,beta,tol)
    )) {
      beta <- beta
      break
    }
    beta_old <- beta
  }
  if (k==max_iter){
    print("Update stopped because the number of loops exceeded max_iter.")
  }
  lb <- VBEMHS_lm_grid.lb(n,p,X,y,sigma2,tau2,q_beta,omega,tau_log_prior,a,b,c_ig,d)

  list(beta=beta, omega=omega, lb=lb, loop=k)
}

# VBEM algorithm (sigma2 and tau2 fixed) : Conjugate prior
VBEMHSc_lm_grid.fixed <- function(X,y,sigma2,tau2,ls_fix=FALSE,tau2_scale=1,
                                 a=0.5,b=0.5,c_ig=0,d=0,max_iter=1000,tol=1e-6)
{
  if (tau2_scale>0){
    tau_log_prior <- function(tau2){ -log(1+tau2/tau2_scale) } # half_Cauchy
  } else {
    tau_log_prior <- function(tau2){ -log(tau2) } # non-informative (improper)
  }
  
  conv <- function(par_old,par_new,tol) {
    abs(par_new-par_old)<tol 
  }
  
  n <- length(y)
  p <- ncol(X)
  XX <- crossprod(X)
  Xy <- crossprod(X,y)
  omega <- rep((b/a),p)
  
  # Initialize q(beta) 
  prec <- omega/tau2
  q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
  beta <- q_beta$mean
  var_b <- q_beta$var
  Eb2 <- as.vector(sigma2*var_b+beta^2)  # E[beta^2]
  
  apb <- a+b
  beta_old <- beta
  
  for (k in 1:max_iter){
    
    if (!ls_fix){
      # Update q(nu) and q(omega)
      nu <- apb/(omega+1)
      omega <- (b+0.5)/(Eb2/(2*sigma2*tau2)+nu)
    }
    
    # Update q(beta) 
    prec <- omega/tau2
    q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
    beta <- q_beta$mean
    var_b <- q_beta$var
    Eb2 <- as.vector(sigma2*var_b+beta^2) # E[beta^2]
    
    if (all(
      conv(beta_old,beta,tol)
    )) {
      beta <- beta
      break
    }
    beta_old <- beta
  }
  if (k==max_iter){
    print("Update stopped because the number of loops exceeded max_iter.")
  }
  lb <- VBEMHSc_lm_grid.lb(n,p,X,y,sigma2,tau2,q_beta,omega,tau_log_prior,a,b,c_ig,d)

  list(beta=beta, omega=omega, lb=lb, loop=k)
}

# VBEM for several fixed values of tau
VBEMHS_lm_grid.trace <- function(X,y,independent=TRUE,
                            ls2_range=c(-3,0),ls2_step=0.5,
                            lt2_range=c(-3,3),lt2_step=0.5,ls_fix=FALSE,tau2_scale=1,
                            a=0.5,b=0.5,c_ig=0,d=0,max_iter=1000,tol=1e-6,verbose=FALSE)
{
  lt2_seq <- seq(lt2_range[1],lt2_range[2],by=lt2_step)
  ls2_seq <- seq(ls2_range[1],ls2_range[2],by=ls2_step)
  
  lb_mat <- matrix(0,length(lt2_seq),length(ls2_seq))
  loop_seq <- NULL

  for (k in seq(along.with=lt2_seq)){
    for (ks in seq(along.with=ls2_seq)){
    tau2 <- 10^lt2_seq[k]
    sigma2 <- 10^ls2_seq[ks]
    if (independent) res <- VBEMHS_lm_grid.fixed(X,y,sigma2,tau2,ls_fix,tau2_scale,a,b,c_ig,d,max_iter,tol)
    else res <- VBEMHSc_lm_grid.fixed(X,y,sigma2,tau2,ls_fix,tau2_scale,a,b,c_ig,d,max_iter,tol)
    lb_mat[k,ks] <- res$lb

    if(verbose) {
      print(paste("tau2=",sprintf("%f", tau2), 
                  " sigma2=",sprintf("%f", sigma2), 
                  " lower bound=",round(res$lb, 5)
      ))
    }
  }
  }
  list(sigma2=10^ls2_seq, tau2=10^lt2_seq, lb=lb_mat)
}

### Wrapper functions ###
# "ls2" stands for log10(sigma^2) 
# "lt2" stands for log10(tau^2) 

VBEMHS_lm_grid <- function(formula, data=NULL, standardize=TRUE,
                          ls2_range=c(-3,0),ls2_step=0.5, # Range of log10(sigma^2)    
                          lt2_range=c(-3,3),lt2_step=0.5, # Range of log10(tau^2)
                          prior=VBEMHS_lm.prior(), 
                          control=VBEMHS_lm.control(), verbose=FALSE)
{
  
  # Make sure the outcome variable is the first column
  mf = model.frame(formula = formula, data = data)
  tf = terms.formula(formula, data=data)
  
  # get the outcome
  y = as.matrix(mf[,1])
  
  # get the predictors and remove the intercept column
  X = model.matrix(tf, data=mf)
  if (all(X[,1]==1)) X = X[,-1,drop=FALSE]
  colnames_X <- colnames(X)
  
  n <- nrow(X)
  p <- ncol(X)
  if (length(y)!=n) stop("The length of y is not equal to the number of rows.")
  
  mean_y <- mean(y)
  mean_X <- apply(X,2,mean)
  
  # standardization
  if (standardize){
    sd_y <- sd(y)
    y <- (y-mean_y)/sd_y
    sd_X <- apply(X,2,sd)
    X <- sweep(sweep(X,2,mean_X),2,sd_X,"/")
  } else {
    sd_y <- 1
    sd_X <- rep(1,p)
  }
  
  independent <- prior$independent
  tau2_scale <- prior$tau2_scale
  ls_fix <- prior$ls_fix
  a <- prior$a
  b <- prior$b
  c_ig <- prior$c_ig
  d <- prior$d
  
  s2_estim <- control$s2_estim
#  lgs_range <- control$lgs_range
#  lgs_step <- control$lgs_step
  max_iter <- control$max_iter
  tol <- control$tol

  res <- VBEMHS_lm_grid.trace(X,y,independent,
                          ls2_range,ls2_step,lt2_range,lt2_step,
                          ls_fix,tau2_scale,a,b,c_ig,d,max_iter,tol,verbose)

  return(res)
}
