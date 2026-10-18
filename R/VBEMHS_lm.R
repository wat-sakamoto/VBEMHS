#' Prior Settings in Bayesian Sparse Linear Regression Models
#'
#' A list of parameters for settings of the priors, used as an argument of VBEMHS_lm.
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
#' @returns A list with components named as the arguments.
#' @seealso [VBEMHS_lm]
#' @export
#'
VBEMHS_lm.prior <- function(
    independent=TRUE, tau2_scale=1, ls_fix=FALSE, a=0.5, b=0.5, c_ig=0, d=0
){
  list(independent=independent, tau2_scale=tau2_scale, ls_fix=ls_fix, a=a, b=b, c_ig=c_ig, d=d)
}

#' Control Parameters for VBEM Algorithm in Bayesian Sparse Linear Regression Models
#'
#' A list of parameters for controlling the VBEM algorithm, used as an argument of VBEMHS_lm.
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
#' @returns A list with components named as the arguments.
#' @seealso [VBEMHS_lm]
#' @export
#'
VBEMHS_lm.control <- function(
    s2_estim="mean", lgs_range=c(-2,4), lgs_step=0.5, max_iter=1000, tol_first=1e-3, tol=1e-6
){
  list(s2_estim=s2_estim, lgs_range=lgs_range, lgs_step=lgs_step,
      max_iter=max_iter, tol_first=tol_first, tol=tol)
}

#' Variational Bayesian EM Algorithm with Horseshoe Prior in Bayesian Sparse Linear Models
#'
#' Variational Bayesian EM algorithm for variable selection and inference in Bayesian sparse linear models
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
#' If \code{gs_trace=TRUE}, the results are traced along the global shrinkage parameter.
#' @param s2_correct logical.
#' If \code{s2_correct=TRUE}, make correction for underestimation of \code{sigma2}.
#' @param prior a list of parameters for settings on the priors.
#' See [VBEMHS_lm.prior] for the default setting.
#' @param control a list of parameters for controlling the VBEM algorithm.
#' See [VBEMHS_lm.control] for the default setting.
#' @param verbose logical.
#' If \code{verbose=TRUE}, some parameter values are displayed in the console at each iteration.
#'
#' @returns If \code{gs_trace=TRUE}, a list object containing the following components:
#' \itemize{
#'  \item \code{beta_std} : a matrix of posterior means of \eqn{\beta_j}'s in standardized scale.
#'  Each row corresponds to a value of \eqn{\tau^2}, and each column corresponds to an explanatory variable.
#'  \item \code{beta} : a matrix of posterior means of \eqn{\beta_j}'s in original scale.
#'  \item \code{beta0} : a vector of posterior means of \eqn{\beta_0} in original scale.
#'  Each element corresponds to a value of \eqn{\tau^2}.
#'  \item \code{se_beta} : a matrix of posterior standard errors of \eqn{\beta_j}'s in original scale.
#'  \item \code{sigma2_std} : a vector of estimated error variance values \eqn{\sigma^2} in standardized scale.
#'  \item \code{sigma2} : a vector of estimated error variance values \eqn{\sigma^2} in original scale.
#'  \item \code{tau2} : a vector of global shrinkage parameters \eqn{\tau^2} values used.
#'  \item \code{omega} : a matrix of estimated precision parameters \eqn{\omega_j}'s in standardized scale.
#'  \item \code{df_s2} : a vector of the degrees of freedom for estimated \eqn{\sigma^2}, also used as the df of t-distribution.
#'  \item \code{lb} : a vector of the estimated lower bounds of the log posterior density (up to constant).
#'  \item \code{iter} : a vector of the numbers of iterations during the VBEM algorithm.
#' }
#'
#' If \code{gs_trace=FALSE}, a list object containing the following components:
#' \itemize{
#'  \item \code{beta_std} : a vector of posterior means of \eqn{\beta_j}'s in standardized scale.
#'  \item \code{beta} : a vector of posterior means of \eqn{\beta_j}'s in original scale.
#'  \item \code{beta0} : a posterior mean of \eqn{\beta_0} in original scale.
#'  \item \code{se_beta} : a vector of posterior standard errors of \eqn{\beta_j}'s in original scale.
#'  \item \code{V_beta} : an estimated variance-covariance matrix of \eqn{\beta}'s in original scale.
#'  \item \code{sigma2_std} : an estimated error variance \eqn{\sigma^2} in standardized scale.
#'  \item \code{sigma2} : an estimated error variance \eqn{\sigma^2} in original scale.
#'  \item \code{tau2} : an estimated global shrinkage parameter \eqn{\tau^2}.
#'  \item \code{omega} : a vector of estimated precision parameters \eqn{\omega_j} in standardized scale.
#'  \item \code{df_s2} : the degree of freedom for estimated \eqn{\sigma^2}, also used as the df of t-distribution.
#'  \item \code{lb} : estimated lower bound of the log posterior density (up to constant).
#'  \item \code{iter} : the number of iterations during the VBEM algorithm.
#' }
#' @seealso [VBEMHS_lm.prior], [VBEMHS_lm.control], [VBEMHS_lm.summary].
#' @examples
#' # Diabetes data
#' library(lars); data(diabetes)
#'
#' # Using the independent prior
#' diabetes.VBEMHS_lm <- VBEMHS_lm(y ~ x, data=diabetes)
#' diabetes.VBEMHS_lm.trace <- VBEMHS_lm(y ~ x, data=diabetes, gs_trace=TRUE)
#'
#' # Using the independent prior, with corrected estimate of sigma2
#' diabetes.VBEMHS_lm <- VBEMHS_lm(y ~ x, data=diabetes, s2_correct=TRUE)
#' diabetes.VBEMHS_lm.trace <- VBEMHS_lm(y ~ x, data=diabetes, gs_trace=TRUE, s2_correct=TRUE)
#'
#' # Using the conjugate prior, with sigma2 estimated by posterior mode
#' diabetes.VBEMHSc_lm <- VBEMHS_lm(y ~ x, data=diabetes,
#'      prior=VBEMHS_lm.prior(independent=FALSE),
#'      control=VBEMHS_lm.control(s2_estim="mode"))
#' @export
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

  res <- list()
#  res <- list(X=X, y=y, sd_X=sd_X, sd_y=sd_y, independent=independent)
  if (gs_trace){
    res$beta_std <- fit$beta
    res$beta <- fit$beta * (sd_y/sd_X)
    res$beta0 <- as.vector(mean_y-tcrossprod(mean_X, res$beta))
    res$se_beta <- fit$se_beta * (sd_y/sd_X)
    res$sigma2_std <- fit$sigma2
    res$sigma2 <- fit$sigma2 * sd_y^2
    res$tau2 <- fit$tau2
    res$omega <- fit$omega
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
    res$beta0 <- as.vector(mean_y-mean_X %*% res$beta)
    res$se_beta <- fit$se_beta * (sd_y/sd_X)
    res$V_beta <- fit$V_beta * tcrossprod(sd_y/sd_X)
    res$sigma2_std <- fit$sigma2
    res$sigma2 <- fit$sigma2 * sd_y^2
    res$tau2 <- fit$tau2
    res$omega <- as.vector(fit$omega)
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

#' Variational Bayesian Variable Selection and Inference in Bayesian Sparse Linear Models
#'
#' Output useful for variable selection and inference in Bayesian sparse linear models.
#'
#' @param object a list object that is a result of a call to [VBEMHS_lm].
#' @param CI logical. if \code{CI=TRUE}, credible intervals for \eqn{\beta_j}'s are included in the output.
#' @param alpha a real value \eqn{\alpha} for which the level for credible intervals is \eqn{1-\alpha}.
#' Must be a single number between 0 and 1.
#'
#' @returns a data frame containing the following components:
#' \item{beta.mean}{posterior mean of \eqn{\beta_j}'s.}
#' \item{beta.se}{posterior standard error of \eqn{\beta_j}'s.}
#' \item{beta.lower, beta.upper}{lower and uppe bound of the credible intervals of \eqn{\beta_j}'s, included only if \code{CI=TRUE}.}
#' \item{t_stat}{t statistics, computed by beta.mean/beta.se.}
#' \item{pd}{probability of direction, defined by \eqn{ \max\{ P(\beta_j<0), P(\beta_j>0) \} },
#' which corresponds to the two-sided p-value \eqn{2(1-\mbox{pd})} for the frequentist methods.}
#' \item{star}{character symbols that shows whether each coefficient is significant. <br>
#' *** : pd>0.9995, &emsp; ** : pd>0.995, &emsp; * : pd>0.975, &emsp; . : pd>0.95.}
#'
#' @seealso [VBEMHS_lm].
#' @examples
#' # Using the independent prior
#' diabetes.VBEMHS_lm <- VBEMHS_lm(y ~ x, data=diabetes)
#' diabetes.VBEMHS_lm.trace <- VBEMHS_lm(y ~ x, data=diabetes, gs_trace=TRUE)
#' VBEMHS_lm.summary(diabetes.VBEMHS_lm, CI=TRUE, alpha=0.05)
#' @export
#'
VBEMHS_lm.summary <- function(object, CI=FALSE, alpha=0.05)
{
  res.df <- data.frame(
    beta.mean=object$beta,
    beta.se=object$se_beta
  )

  # Credible intervals
  if (CI) {
    res.df$beta.lower <- object$beta+qt(alpha/2, object$df_s2)*object$se_beta
    res.df$beta.upper <- object$beta+qt(1-alpha/2, object$df_s2)*object$se_beta
  }

  # Probability of direction (0.5<=pd<1)
  t_stat <- ifelse(object$se_beta>0, object$beta/object$se_beta, 0)
  pd <- pt(abs(t_stat), object$df_s2)
  star <- ifelse(pd>0.9995, "***",
           ifelse(pd>0.995, "**",
            ifelse(pd>0.975, "*",
             ifelse(pd>0.95, ".", ""))))
  res.df$t_stat=t_stat
  res.df$pd=pd
  res.df$star=star

  rownames(res.df) <- names(object$beta)
  return(res.df)
}
