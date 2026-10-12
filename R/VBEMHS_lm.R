# VBEMHS_lm: Variational Bayesian EM Approach to Variable Selection 
# in linear regression model with (independent/conjugate) Horseshoe Prior
# Notice: no intercept term is included, so X and y are supposed to be centered.

#### Setting

VBEMHS_lm.prior <- function(
    independent=TRUE, # prior for beta and sigma2: independent/conjugate horseshoe prior
    tau2_scale=1,     # prior for tau2: half-Cauchy (>0) / non-informative (=0)
    ls_fix=FALSE,     # option for lambda (local shrinkage parameter) (TRUE: ridge regression)
    a=0.5, b=0.5,     # hyper-parameters for lambda: beta-prime prior (default: half-Cauchy)
    c_ig=0, d=0       # hyper-parameters for sigma2: inverse gamma prior
){
  list(independent=independent, tau2_scale=tau2_scale, ls_fix=ls_fix, 
       a=a, b=b, c_ig=c_ig, d=d)
}
  
VBEMHS_lm.control <- function(
    s2_estim="mean",                  # posterior estimate of sigma2 (mean/mode)
    lgs_range=c(-2,4), lgs_step=0.5,  # range of log10 global shrinkage parameter (for gs_trace=TRUE)
    max_iter=1000,                    # maximum number of iterations
    tol_first=1e-3, tol=1e-6          # tolerance for convergence 
){
  list(s2_estim=s2_estim, lgs_range=lgs_range, lgs_step=lgs_step, 
      max_iter=max_iter, tol_first=tol_first, tol=tol)  
}

#### Main

VBEMHS_lm <- function(formula, data=NULL, standardize=TRUE, 
                      gs_trace=FALSE,     # If TRUE, trace along the global shrinkage parameter
                      s2_correct=FALSE,   # If TRUE, correct for underestimation of sigma2
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
  lgs_range <- control$lgs_range
  lgs_step <- control$lgs_step
  max_iter <- control$max_iter
  tol_first <- control$tol_first
  tol <- control$tol

  if (gs_trace){
    # Trace for several fixed values of the global shrinkage parameter
    if (independent) {
      # Independent prior
      fit <- VBEMHS_lm.trace(X,y,lgs_range,lgs_step,ls_fix,s2_estim,s2_correct,
                             tau2_scale,a,b,c_ig,d,max_iter,tol,verbose)
    } else {
      # Conjugate prior
      fit <- VBEMHSc_lm.trace(X,y,lgs_range,lgs_step,ls_fix,s2_estim,s2_correct,
                              tau2_scale,a,b,c_ig,d,max_iter,tol,verbose)
    }
  } else {

    if (independent) {
      # Independent prior
      if (lgs_step>0) {
        fit_first <- VBEMHS_lm.gs(X,y,lgs_range,lgs_step,s2_estim,s2_correct,
                                  tau2_scale,a,b,c_ig,d,tol=tol_first,verbose)
        gs_ini <- fit_first$gs_opt
      } else {
        gs_ini <- n
      }
      fit <- VBEMHS_lm.update(X,y,gs_ini,gs_fix=FALSE,ls_fix,s2_estim,s2_correct,
                                tau2_scale,a,b,c_ig,d,max_iter,tol,verbose)
    } else {
      # Conjugate prior
      if (lgs_step>0) {
        fit_first <- VBEMHSc_lm.gs(X,y,lgs_range,lgs_step,s2_estim,s2_correct,
                                  tau2_scale,a,b,c_ig,d,tol=tol_first,verbose)
        gs_ini <- fit_first$gs_opt
      } else {
        gs_ini <- n
      }
      fit <- VBEMHSc_lm.update(X,y,gs_ini,gs_fix=FALSE,ls_fix,s2_estim,s2_correct,
                                 tau2_scale,a,b,c_ig,d,max_iter,tol,verbose)
    }

  }
  
  res <- list(X=X, y=y, sd_X=sd_X, sd_y=sd_y, independent=independent)
  if (gs_trace){
    res$tau2 <- fit$tau2
    res$beta_std <- fit$beta
    res$beta <- fit$beta * (sd_y/sd_X)
    res$beta0 <- as.vector(mean_y-tcrossprod(mean_X, res$beta))
    res$se_beta <- fit$se_beta * (sd_y/sd_X)
    res$omega <- fit$omega
    res$sigma2_std <- fit$sigma2 
    res$sigma2 <- fit$sigma2 * sd_y^2
    res$df_s2 <- fit$df_s2
    res$lb <- fit$lb
    res$iter <- fit$iter
    colnames(res$beta_std) <- colnames_X
    colnames(res$beta) <- colnames_X
    colnames(res$se_beta) <- colnames_X
    colnames(res$omega) <- colnames_X
  } else {
    res$beta_std <- as.vector(fit$beta) 
    res$beta <- as.vector(fit$beta) * (sd_y/sd_X)
    res$se_beta <- fit$se_beta * (sd_y/sd_X)
    res$V_beta <- fit$V_beta * tcrossprod(sd_y/sd_X)
    res$omega <- as.vector(fit$omega)
    res$beta0 <- as.vector(mean_y-mean_X %*% res$beta)
    res$sigma2_std <- fit$sigma2
    res$sigma2 <- fit$sigma2 * sd_y^2
    res$tau2 <- fit$tau2
    res$df_s2 <- fit$df_s2
    res$lb <- fit$lb
    res$iter <- fit$iter
    names(res$beta_std) <- colnames_X
    names(res$beta) <- colnames_X
    names(res$se_beta) <- colnames_X
    names(res$omega) <- colnames_X
    rownames(res$V_beta) <- colnames_X
    colnames(res$V_beta) <- colnames_X
  }
  return(res)
}

##### summary 

# Summary
VBEMHS_lm.summary <- function(res, CI=FALSE, alpha=0.05)
{
  res.df <- data.frame(
    beta.mean=res$beta, 
    beta.se=res$se_beta
  )
    
  # Credible intervals
  if (CI) {
    res.df$beta.lower <- res$beta+qt(alpha/2, res$df_s2)*res$se_beta
    res.df$beta.upper <- res$beta+qt(1-alpha/2, res$df_s2)*res$se_beta
  }

  # Probability of direction (0.5<=pd<1)
  t_stat <- ifelse(res$se_beta>0, res$beta/res$se_beta, 0)
  pd <- pt(abs(t_stat), res$df_s2)
  star <- ifelse(pd>0.9995, "***", 
           ifelse(pd>0.995, "**",
            ifelse(pd>0.975, "*",
             ifelse(pd>0.95, ".", ""))))
  res.df$t_stat=t_stat
  res.df$pd=pd
  res.df$star=star
  
  rownames(res.df) <- names(res$beta) 
  return(res.df)
}
