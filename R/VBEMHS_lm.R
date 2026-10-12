#' Prior setting in the sparse linear model
#'
#' A list to use as an argument of VBEMHS_lm
#'
#' @param independent logical.
#' Prior for beta and sigma2.
#' If \code{independent=TRUE} (default), the priors for \code{beta} and \code{sigma2} are independent.
#' If \code{independent=FALSE}, the priors for  \code{beta} and \code{sigma2} are dependent and conjugate.
#' @param tau2_scale non-negative scalar.
#' Specify the type of the prior for the global shrinkage parameter \code{tau2}.
#' If \code{tau2_scale}\eqn{>0}, the half-Cauchy prior
#' \deqn{
#'   \pi(\tau) \propto 1/(1+\tau^2/\tau_0^2),
#' }
#' is used, where \eqn{\tau_0^2} is specified by \code{tau2_scale} (default: 1).
#' If \code{tau2_scale=0}, the prior is non-informative.
#' @param ls_fix logical.
#' Specify whether the local shrinkage parameters \code{lambda} are fixed or estimated.
#' If \code{ls_fix=TRUE}, all elements of \code{lambda} are fixed to 1.
#' If \code{ls_fix=FALSE} (default), all elements of \code{lambda} are estimated.
#' @param a,b hyper-parameters in beta-prime prior \eqn{\beta^{\prime}(a,b)} for \code{lambda}.
#' @param c_ig,d parameters in inverse gamma prior \eqn{\mathrm{IG}(c,d)} for \code{sigma2}.
#' @export
#'
VBEMHS_lm.prior <- function(
    independent=TRUE, tau2_scale=1, ls_fix=FALSE, a=0.5, b=0.5, c_ig=0, d=0
){
  list(independent=independent, tau2_scale=tau2_scale, ls_fix=ls_fix, a=a, b=b, c_ig=c_ig, d=d)
}

#' Control parameters for VBEM algorithm
#'
#' A list to use as an argument of VBEMHS_lm
#'
#' @param s2_estim character. Specify how to estimate the variance parameter \code{sigma2}.
#' If \code{s2_estim="mean"} (default), \code{sigma2} is estimated by the harmonic posterior mean.
#' If \code{s2_estim="mode"}, \code{sigma2} is estimated by the posterior mode.
#' @param tau2_scale non-negative scalar.
#' If \code{tau2_scale}\eqn{>0}, the prior for tau2 is half-Cauchy with specified scale (default: 1).
#' If \code{tau2_scale=0}, the prior for tau2 is non-informative.
#' @param lgs_range The range of the global shrinkage parameter in \eqn{\log_{10}} scale.
#' Default \code{lgs_range=c(-2,4)}.
#' @param lgs_step The step of the grid search for the global shrinkage parameter in \eqn{\log_{10}} scale.
#' Default \code{lgs_step=0.5}.
#' @param max_iter integer. Specify the maximum number of iterations for the VBEM update. Default 1000.
#' @param tol_first Positive convergence tolerance for the initial search. Default \code{1e-3}.
#' @param tol Positive convergence tolerance for the VBEM update. Default \code{1e-6}.
#' The program judges that the iterations converge if the updates for all parameters are smaller than \code{tol}.
#'
#' @export
#'
VBEMHS_lm.control <- function(
    s2_estim="mean", lgs_range=c(-2,4), lgs_step=0.5, max_iter=1000, tol_first=1e-3, tol=1e-6
){
  list(s2_estim=s2_estim, lgs_range=lgs_range, lgs_step=lgs_step,
      max_iter=max_iter, tol_first=tol_first, tol=tol)
}

#' VBEM Algorithm with Horseshoe Prior in Sparse Linear Model
#'
#' Variational Bayesian EM algorithm for variable selection in sparse linear model
#' with (conjugate/independent) horseshoe prior.
#' The variational distribution is factorized using gamma mixture representation of the horseshoe prior.
#' No intercept term is included, so \code{X} and \code{y} are supposed to be centered.
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data a data frame, list or environment containing the variables in the model.
#' @param standardize logical.
#' If \code{standardize=TRUE} (default), both response and explanatory variables are standardized prior to fitting the model.
#' The final results are returned on the original scale.
#' @param gs_trace logical.
#' If \code{gs_trace=TRUE}, trace the results along the global shrinkage parameter. Default \code{FALSE}.
#' @param s2_correct logical.
#' If \code{s2_correct=TRUE}, correct for under-estimation of \code{sigma2}.
#' @param prior list for setting the priors. See \code{VBEMHS_lm.prior} for the default setting.
#' @param control list for setting the control parameters for the VBEM algorithm.
#' See \code{VBEMHS_lm.control} for the default setting.
#' @param verbose logical. It specifies whether to display some values at each iteration.
#'
#' @returns If \code{gs_trace=TRUE}, a list object containing the following components:
#' \itemize{
#'  \item \code{tau2} : a vector of \eqn{\tau^2} values used.
#'  \item \code{beta} : a matrix of posterior mean of \eqn{\beta_j}'s. Each row corresponds to a value of \eqn{\tau^2},
#'  and each column corresponds to an explanatory variable.
#'  \item \code{se_beta} : a matrix of posterior standard errors of \eqn{\beta_j}'s.
#'  \item \code{kappa} : a matrix of estimated shrinkage coefficients.
#'  \item \code{sigma2} : a vector of estimated \eqn{\sigma^2}.
#'  Each element corresponds to a value of \eqn{\tau^2}.
#'  \item \code{lb} : a vector of the estimated lower bound of the log posterior density (up to constant).
#'  \item \code{iter} : a vector of the numbers of iterations.
#' }
#' If \code{gs_trace=FALSE}, a list object containing the following components:
#' \itemize{
#'  \item \code{beta} : posterior means of \eqn{\beta_j}'s.
#'  \item \code{se_beta} : posterior standard errors of \eqn{\beta_j}'s.
#'  \item \code{V_beta} : estimated variance-covariance matrix of \eqn{\beta}'s.
#'  \item \code{sigma2} : an estimated \eqn{\sigma^2}.
#'  \item \code{tau2} : an estimated \eqn{\tau^2}.
#'  \item \code{omega} : estimated precision parameters \eqn{\omega_j}.
#'  \item \code{kappa} : estimated shrinkage coefficients.
#'  \item \code{lb} : estimated lower bound of the log posterior density (up to constant).
#'  \item \code{iter} : the numbers of iterations.
#' }
#' @export
#'
#' @examples
#' # Diabetes data
#' library(lars); data(diabetes)
#' diabetes.lm <- VBEMHS_lm(y ~ x, data=diabetes)
#' diabetes.trace <- VBEMHS_lm(y ~ x, data=diabetes, tau_trace=TRUE)
#'
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
