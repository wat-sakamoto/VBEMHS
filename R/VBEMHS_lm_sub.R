# Penalized least square estimate
VBEM.reg <- function(XX,Xy,prec)
{
  XXA <- XX+diag(prec)
  R <- chol(XXA)
  beta <- backsolve(R,forwardsolve(t(R),Xy))

  return(list(beta=as.vector(beta), R=R))
}

# Penalized least square estimate (Sherman-Morrison-Woodbury formula)
VBEM.reg_nltp <- function(n,X,Xy,prec) # in case n is less than p
{
  AX <- (1/prec)*t(X) # A X^t
  AXy <- (1/prec)*Xy  # A X^t y
  XAX <- crossprod(AX,t(X)) # X A X^t
  XAXy <- crossprod(AX,Xy)  # X A X^t y

  IXAX <- diag(1,n)+XAX
  R <- chol(IXAX)
  beta <- AXy-tcrossprod(AX,t(backsolve(R,forwardsolve(t(R),XAXy))))

  return(list(beta=as.vector(beta), R=R))
}

# Update q(beta)
VBEMHS_lm.beta <- function(n,p,X,XX,Xy,prec)
{
  if (n < p){
    reg <- VBEM.reg_nltp(n,X,Xy,prec)
  } else {
    reg <- VBEM.reg(XX,Xy,prec)
  }
  beta <- reg$beta
  R <- reg$R
  logdet <- sum(log(diag(R)))

  if (n < p) {  # SMW formula
    logdet <- logdet + sum(log(prec))/2
    AX <- (1/prec)*t(X)
    V <- diag(1/prec)-crossprod(forwardsolve(t(R),t(AX)))
  } else {
    V <- tcrossprod(backsolve(R,diag(1,p)))
  }
  var.b <- diag(V)

  list(mean=beta,     # (XX+A)^{-1} Xy
       var=var.b,     # diag of (XX+A)^{-1}
       V=V,           # (XX+A)^{-1}
       logdet=logdet) # (1/2)*log|XX+A|
}

# Update q(omega), omega=1/lambda^2
VBEMHS_lm.omega<- function(Eb2,tau2,nu,b=0.5)
{
  omega <- (b+0.5)/(Eb2/(2*tau2)+nu)
  return(omega)   # omega.bar=E[ 1/lambda^2 ]
}

# Compute the lower bound for log marginal posterior
VBEMHS_lm.lb <- function(n,p,X,y,sigma2,tau2,q_beta,omega,
                      tau_log_prior,a=0.5,b=0.5,c_ig=0,d=0)
{
  beta <- q_beta$mean
  var_b <- q_beta$var
  SS <- crossprod(y, y-X %*% beta)

  lb <- (n-p+2*c_ig+2)*log(sigma2) + p*log(tau2)
  lb <- lb + (SS + 2*d)/sigma2
  lb <- -lb/2

  lb <- lb - q_beta$logdet
  lb <- lb + (b+0.5)*(sum(log(omega))+p*(1-log(b+0.5))) + lgamma(b+0.5)
  lb <- lb - (a+b)*sum(log(omega+1))  + lgamma(a+b)
  lb <- lb + tau_log_prior(tau2)

  return(as.numeric(lb))
}

# VBEM algorithm (tau2 updated in iterations)
# "lt2" stands for log10(tau2)
VBEMHS_lm.tau2_update <- function(X,y,
                               lt2_range=c(-6,0),lt2_step=0.5,tau_log_prior,tau2_scale=1,bc=0,
                               a=0.5,b=0.5,c_ig=0,d=0,max_iter=1000,tol=1e-6,verbose=FALSE)
{
  # bc=2: sigma2 = |y-X beta|^2/tr((I-S)^2)
  sigma2_correct2 <- function(n,p,X,y,beta,V,prec,c_ig,d,tol)
  {
    SS <- crossprod(y-X %*% beta)
    tr_IS2 <- n-p+crossprod(prec, (V^2) %*% prec)
    sigma2 <- (SS+2*d)/(tr_IS2+2*c_ig+2)
    return(max(as.numeric(sigma2),tol))
  }

  # bc=1: sigma2 = y(y-X beta)/(n-tr_S)
  sigma2_correct <- function(n,p,X,y,beta,var_b,prec,c_ig,d,tol)
  {
    SS <- crossprod(y, y-X %*% beta)
    tr_S <- p-sum(var_b*prec)
    sigma2 <- (SS+2*d)/(n-tr_S+2*c_ig+2)
    return(max(as.numeric(sigma2),tol))
  }

  # bc=0: sigma2 as marginal posterior mode
  sigma2_update <- function(n,p,X,y,beta,var_b,prec,sigma2_old,c_ig,d,tol)
  {
    SS <- crossprod(y-X %*% beta)
    tr_S <- p-sum(var_b*prec)
    sigma2 <- (SS+sigma2_old*tr_S+2*d)/(n+2*c_ig+2)
    return(max(as.numeric(sigma2),tol))
  }

  tau2_update <- function(p,EB2O,tau2_range,tau2_scale){
    if (tau2_scale>0){
      A <- (p+2)/tau2_scale
      B <- p-EB2O/tau2_scale
      D <- (p+EB2O/tau2_scale)^2+8*EB2O/tau2_scale
      tau2 <- (-B+sqrt(D))/A*0.5
    } else {
      tau2 <- EB2O/(p+2)
    }
    tau2 <- min(max(tau2,tau2_range[1]),tau2_range[2])
    return(tau2)
  }

  conv <- function(par_old,par_new,tol) {
    abs(par_new-par_old)<tol
  }

  n <- length(y)
  p <- ncol(X)
  XX <- crossprod(X)
  Xy <- crossprod(X,y)

  omega_ini <- rep(n,p)

  # Initialize q(beta)
  prec <- omega_ini
  q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
  beta <- q_beta$mean
  var_b <- q_beta$var

  # Initialize sigma2
  sigma2_old <- as.numeric(crossprod(y-X %*% beta)/n)
  Eb2 <- sigma2_old*var_b+beta^2 # E[beta^2]

  # Initial search for tau2
  apb <- a+b
  nu <- apb/(omega_ini+1)
  lb_seq <- NULL
  lt2_seq <- seq(lt2_range[1],lt2_range[2],by=lt2_step)

  for (l in seq(along.with=lt2_seq)){
    tau2 <- 10^lt2_seq[l]

    # Update q(omega)
    omega <- VBEMHS_lm.omega(Eb2,tau2,nu,b)

    # Update q(beta)
    prec <- omega*(sigma2_old/tau2)
    q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
    beta <- q_beta$mean
    var_b <- q_beta$var

    # Update sigma2
    sigma2 <- switch(bc+1,
                     sigma2_update(n,p,X,y,beta,var_b,prec,sigma2_old,c_ig,d,tol), # bc=0
                     sigma2_correct(n,p,X,y,beta,var_b,prec,c_ig,d,tol),           # bc=1
                     sigma2_correct2(n,p,X,y,beta,q_beta$V,prec,c_ig,d,tol)        # bc=2
    )

    lb <- VBEMHS_lm.lb(n,p,X,y,sigma2,tau2,q_beta,omega,tau_log_prior,a,b,c_ig)
    lb_seq <- c(lb_seq, lb)
    if(verbose) {
      print(paste("k=", 0,
                  " tau2=",sprintf("%f", tau2),
                  " sigma2=",sprintf("%f", sigma2),
                  " lower bound=",round(lb, 5)
      ))
    }
  }
  lt2_opt <- lt2_seq[which.max(lb_seq)]
  tau2 <- 10^lt2_opt

  omega <- VBEMHS_lm.omega(Eb2,tau2,nu,b)
  prec <- omega*(sigma2_old/tau2)
  q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
  beta <- q_beta$mean
  var_b <- q_beta$var

  sigma2 <- switch(bc+1,
                   sigma2_update(n,p,X,y,beta,var_b,prec,sigma2_old,c_ig,d,tol), # bc=0
                   sigma2_correct(n,p,X,y,beta,var_b,prec,c_ig,d,tol),           # bc=1
                   sigma2_correct2(n,p,X,y,beta,q_beta$V,prec,c_ig,d,tol)        # bc=2
  )
  Eb2 <- sigma2*var_b+beta^2       # E[beta^2]

  beta_old <- beta
  sigma2_old <- sigma2
  tau2_old <- tau2
  tau2_range <- 10^c(min(lt2_range),max(lt2_range))

  for (k in 1:max_iter){

    # Update q(nu) and q(omega)
    nu <- apb/(omega+1)
    omega <- VBEMHS_lm.omega(Eb2,tau2_old,nu,b)

    # Update q(beta)
    prec <- omega*(sigma2_old/tau2_old)
    q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
    beta_new <- q_beta$mean
    var_b <- q_beta$var
    Eb2 <- sigma2_old*var_b+beta_new^2 # E[beta^2]

    # Update sigma2
    sigma2_new <- switch(bc+1,
                         sigma2_update(n,p,X,y,beta_new,var_b,prec,sigma2_old,c_ig,d,tol), # bc=0
                         sigma2_correct(n,p,X,y,beta_new,var_b,prec,c_ig,d,tol),           # bc=1
                         sigma2_correct2(n,p,X,y,beta_new,q_beta$V,prec,c_ig,d,tol)        # bc=2
    )

    # Update tau2
    EB2O <- sum(omega*Eb2)         # sum of E[1/lambda^2]E[beta^2]
    tau2_new <- tau2_update(p,EB2O,tau2_range,tau2_scale)

    if(verbose) {
      lb <- VBEMHS_lm.lb(n,p,X,y,sigma2_new,tau2_new,q_beta,omega,tau_log_prior,a,b,c_ig)
      print(paste("k=", k,
                  " tau2=",sprintf("%f", tau2_new),
                  " sigma2=",sprintf("%f", sigma2_new),
                  " lower bound=",round(lb, 5)
      ))
    }
    if (all(
      conv(beta_old,beta_new,tol),
      conv(sigma2_old,sigma2_new,tol),
      conv(log10(tau2_old),log10(tau2_new),tol)
    )) {
      beta <- beta_new
      sigma2 <- sigma2_new
      tau2 <- tau2_new
      break
    }

    beta_old <- beta_new
    sigma2_old <- sigma2_new
    tau2_old <- tau2_new
  }
  if (k==max_iter){
    print("Update stopped because the number of iterations exceeded max_iter.")
    beta <- beta_new
    sigma2 <- sigma2_new
    tau2 <- tau2_new
  }
  lb <- VBEMHS_lm.lb(n,p,X,y,sigma2,tau2,q_beta,omega,tau_log_prior,a,b,c_ig)

  list(beta=beta, se_beta=sqrt(q_beta$var*sigma2), V_beta=q_beta$V*sigma2, sigma2=sigma2,
       tau2=tau2, omega=omega, kappa=omega/(omega+n*tau2/sigma2), lb=lb, iter=k)
}

# VBEM algorithm (tau2 fixed)
VBEMHS_lm.tau2_fixed <- function(X,y,tau2,tau_log_prior,bc=0,
                              a=0.5,b=0.5,c_ig=0,d=0,max_iter=1000,tol=1e-6,verbose=FALSE)
{

  # bc=2: sigma2 = |y-X beta|^2/tr((I-S)^2)
  sigma2_correct2 <- function(n,p,X,y,beta,V,prec,c_ig,d,tol)
  {
    SS <- crossprod(y-X %*% beta)
    tr_IS2 <- n-p+crossprod(prec, (V^2) %*% prec)
    sigma2 <- (SS+2*d)/(tr_IS2+2*c_ig+2)
    return(max(as.numeric(sigma2),tol))
  }

  # bc=1: sigma2 = y(y-X beta)/(n-tr_S)
  sigma2_correct <- function(n,p,X,y,beta,var_b,prec,c_ig,d,tol)
  {
    SS <- crossprod(y, y-X %*% beta)
    tr_S <- p-sum(var_b*prec)
    sigma2 <- (SS+2*d)/(n-tr_S+2*c_ig+2)
    return(max(as.numeric(sigma2),tol))
  }

  # bc=0: sigma2 as marginal posterior mode
  sigma2_update <- function(n,p,X,y,beta,var_b,prec,sigma2_old,c_ig,d,tol)
  {
    SS <- crossprod(y-X %*% beta)
    tr_S <- p-sum(var_b*prec)
    sigma2 <- (SS+sigma2_old*tr_S+2*d)/(n+2*c_ig+2)
    return(max(as.numeric(sigma2),tol))
  }

  conv <- function(par_old,par_new,tol) {
    abs(par_new-par_old)<tol
  }

  n <- length(y)
  p <- ncol(X)

  XX <- crossprod(X)
  Xy <- crossprod(X,y)

  omega <- rep(n,p)

  # Initialize q(beta)
  prec <- omega
  q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
  beta <- q_beta$mean
  var_b <- q_beta$var

  # Initialize sigma2
  sigma2 <- as.numeric(crossprod(y-X %*% beta)/n)
  Eb2 <- sigma2*var_b+beta^2 # E[beta^2]

  # Update q(nu) and q(omega)
  apb <- a+b
  nu <- apb/(omega+1)
  omega <- VBEMHS_lm.omega(Eb2,tau2,nu,b)

  beta_old <- beta
  sigma2_old <- sigma2

  for (k in 1:max_iter){

    # Update q(nu) and q(omega)
    nu <- apb/(omega+1)
    omega <- VBEMHS_lm.omega(Eb2,tau2,nu,b)

    # Update q(beta)
    prec <- omega*(sigma2_old/tau2)
    q_beta <- VBEMHS_lm.beta(n,p,X,XX,Xy,prec)
    beta_new <- q_beta$mean
    var_b <- q_beta$var
    Eb2 <- sigma2_old*var_b+beta_new^2 # E[beta^2]

    # Update sigma2
    sigma2_new <- switch(bc+1,
                         sigma2_update(n,p,X,y,beta_new,var_b,prec,sigma2_old,c_ig,d,tol), # bc=0
                         sigma2_correct(n,p,X,y,beta_new,var_b,prec,c_ig,d,tol),           # bc=1
                         sigma2_correct2(n,p,X,y,beta_new,q_beta$V,prec,c_ig,d,tol)        # bc=2
    )

    if(verbose) {
      lb <- VBEMHS_lm.lb(n,p,X,y,sigma2_new,tau2,q_beta,omega,tau_log_prior,a,b,c_ig)
      print(paste("k=", k,
                  " tau2=",sprintf("%f", tau2),
                  " sigma2=",sprintf("%f", sigma2_new),
                  " lower bound=",round(lb, 5)
      ))
    }

    if (all(
      conv(beta_old,beta_new,tol),
      conv(sigma2_old,sigma2_new,tol)
    )) {
      beta <- beta_new
      sigma2 <- sigma2_new
      break
    }
    beta_old <- beta_new
    sigma2_old <- sigma2_new
  }
  if (k==max_iter){
    print("Update stopped because the number of iterations exceeded max_iter.")
    beta <- beta_new
    sigma2 <- sigma2_new
  }
  lb <- VBEMHS_lm.lb(n,p,X,y,sigma2,tau2,q_beta,omega,tau_log_prior,a,b,c_ig)

  list(beta=beta, se_beta=sqrt(q_beta$var*sigma2), V_beta=q_beta$V*sigma2, sigma2=sigma2,
       omega=omega, kappa=omega/(omega+n*tau2/sigma2), lb=lb, iter=k)
}


# VBEM for several fixed values of tau
VBEMHS_lm.trace <- function(X,y,lt2_range,lt2_step,tau_log_prior,bc=0,
                         a=0.5,b=0.5,c_ig=0,d=0,max_iter=1000,tol=1e-6,verbose=FALSE)
{
  lt2_seq <- seq(lt2_range[1],lt2_range[2],by=lt2_step)

  beta_mat <- NULL
  se_beta_mat <- NULL
  sigma2_seq <- NULL
  kappa_mat <- NULL
  lb_seq <- NULL
  iter_seq <- NULL

  for (k in seq(along.with=lt2_seq)){
    tau2 <- 10^lt2_seq[k]
    res <- VBEMHS_lm.tau2_fixed(X,y,tau2,tau_log_prior,bc,
                             a,b,c_ig,d,max_iter,tol)
    beta_mat <- rbind(beta_mat,res$beta)
    se_beta_mat <- rbind(se_beta_mat,res$se_beta)
    sigma2_seq <- c(sigma2_seq,res$sigma2)
    kappa_mat <- rbind(kappa_mat,res$kappa)
    lb_seq <- c(lb_seq,res$lb)
    iter_seq <- c(iter_seq,res$iter)
    if(verbose) {
      print(paste("tau2=",sprintf("%f", tau2),
                  " sigma2=",sprintf("%f", res$sigma2),
                  " lower bound=",round(res$lb, 5)
      ))
    }
  }
  list(tau2=10^lt2_seq,beta=beta_mat,se_beta=se_beta_mat,kappa=kappa_mat,sigma2=sigma2_seq,
       lb=lb_seq,iter=iter_seq)
}
