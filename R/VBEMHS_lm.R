#' VBEM Algorithm with Horseshor Prior in Sparse Linear Model
#'
#' Variational Bayesian EM algorithm with (independent) horseshoe prior
#' for variable selection in sparse linear model.
#' The variational distribution is factorized using gamma mixture representation of horseshoe prior.
#' No intercept term is included, so X and y are supposed to be centered.
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data a data frame, list or environment containing the variables in the model.
#' @param standardize logical.
#' If \code{standardize=TRUE} (default), both response and explanatory variables are standardized prior to fitting the model.
#' The results are always returned on the original scale.
#' @param tau_trace logical.
#' If \code{tau_trace=TRUE}, the VBEM algorithm is applied for fixed \eqn{\tau}'s on the grid sequence of \eqn{\log_{10}\tau^2} specified by \code{lt2_range} and \code{lt2_step}.
#' If \code{tau_trace=FALSE} (default), after initial search of \eqn{\tau} on the grid sequence of \eqn{\log_{10}\tau^2} specified by \code{lt2_range} and \code{lt2_step},
#' and then the VBEM algorithm is applied to find the posterior mode of \eqn{(\tau,\sigma^2)},
#' @param lt2_range a vector of length 2 that specifies the lower and upper bounds of \eqn{\log_{10}\tau^2}.
#' @param lt2_step a positive value that specifies the step size of the sequence of \eqn{\log_{10}\tau^2}.
#' @param tau2_scale If \code{tau2_scale>0}, the half-Cauchy prior
#' \deqn{
#'   \pi(\tau) \propto 1/(1+\tau^2/\tau_0^2),
#' }
#' is used for the prior on \eqn{\tau}, where \eqn{\tau_0^2} is specified by \code{tau2_scale}.
#' If \code{tau2_scale=0}, the improper prior \eqn{\pi(\tau) \propto 1/\tau} is used.
#' The default is \code{tau2_scale=1}.
#' @param bc If \code{bc=0} (default), bias correction is not used for the estimate of \eqn{\sigma^2}.
#' If \code{bc=1} or \code{bc=2}, bias correction (of the type BC1 or BC2, respectively) is used.
#' In case \eqn{p} is as great as or greater than \eqn{n}, either bias correction is recommended.
#' See our paper for detailed explanation.
#' @param a,b parameters in beta-prime prior \eqn{\beta^{\prime}(a,b)} for local shrinkage parameters \eqn{\lambda_j (j=1,\ldots,p)}.
#' @param c_ig,d parameters in inverse gamma prior \eqn{\mathrm{IG}(c,d)} for the error variance \eqn{\sigma^2}.
#' If \code{c_ig=0} and \code{d=0} (default), the improper prior \eqn{\pi(\sigma^2) \propto 1/\sigma^2} is used.
#' @param max_iter Maximum number of iterations in the VBEM algorithm.
#' @param tol Positive convergence tolerance (say \eqn{\delta}).
#' The program judges that The iterations converge when
#' \eqn{|\beta_j^{\mathrm{new}}-\beta_j^{\mathrm{old}}|<\delta}, \eqn{j=1,\ldots,p},
#' \eqn{|\log_{10}(\tau_{\mathrm{new}}^2)-\log_{10}(\tau_{\mathrm{old}}^2)|<\delta} and
#' \eqn{|\sigma_{\mathrm{new}}^2-\sigma_{\mathrm{old}}^2|<\delta}.
#' @param alpha A value \eqn{\alpha} between 0 and 1 that specifies the coverage probability \eqn{1-\alpha} of the credible intervals for \eqn{\beta_j}'s.
#' @param verbose logical. It specifies whether to display some values at each iteration.
#'
#' @returns If \code{tau_search=TRUE}, a list object containing the following components:
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
#' If \code{tau_search=FALSE}, a list object containing the following components:
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
#' library(lars); data(diabates)
#' diabetes.trace <- VBEMHS_lm(y ~ x, data=diabetes, tau_trace=TRUE)
#' diabetes.lm <- VBEMHS_lm(y ~ x, data=diabetes)
#'
VBEMHS_lm <- function(formula, data=NULL, standardize=TRUE, tau_trace=FALSE,
                   lt2_range=c(-6,0), lt2_step=0.5, tau2_scale=1, bc=0,
                   a=0.5, b=0.5, c_ig=0, d=0,
                   max_iter=1000, tol=1e-6, alpha=0.05, verbose=FALSE)
{
  if (tau2_scale>0){
    tau_log_prior <- function(tau2){ -log(1+tau2/tau2_scale) } # half_cauchy
  } else {
    tau_log_prior <- function(tau2){ -log(tau2) } # improper
  }

  # Make sure the outcome variable is the first column
  mf = model.frame(formula = formula, data = data)
  t = terms.formula(formula, data=data)

  # get the outcome
  y = as.matrix(mf[,1])

  # get the predictors and remove the intercept column
  X = model.matrix(t, data=mf)
  if (all(X[,1]==1)) X = X[,-1,drop=FALSE]
  colnames_X <- colnames(X)

  n <- nrow(X)
  p <- ncol(X)
  if (length(y)!=n) stop("The length of y is not equal to the number of rows.")

  # standardization
  if (standardize){
    mean_y <- mean(y)
    sd_y <- sd(y)
    y <- (y-mean_y)/sd_y
    mean_X <- apply(X,2,mean)
    sd_X <- apply(X,2,sd)
    X <- sweep(sweep(X,2,mean_X),2,sd_X,"/")
  } else {
    sd_y <- 1
    sd_X <- rep(1,p)
  }

  if (tau_trace){
    fit <- VBEMHS_lm.trace(X,y,lt2_range,lt2_step,tau_log_prior,bc,
                        a,b,c_ig,d,max_iter,tol,verbose)
  } else {
    fit <- VBEMHS_lm.tau2_update(X,y,lt2_range,lt2_step,tau_log_prior,tau2_scale,bc,
                              a,b,c_ig,d,max_iter,tol,verbose)
  }

  if (tau_trace){
    res$beta <- fit$beta * (sd_y/sd_X)
    res$se_beta <- fit$se_beta * (sd_y/sd_X)
    res$sigma2 <- fit$sigma2 * sd_y^2
    colnames(res$beta) <- colnames_X
    colnames(res$se_beta) <- colnames_X
    colnames(res$kappa) <- colnames_X
#    res$eff <- p-apply(res$kappa,1,sum)
  } else {
    res$beta <- fit$beta * (sd_y/sd_X)
    res$beta0 <- as.numeric(mean_y-mean_X %*% res$beta)
    res$se_beta <- fit$se_beta * (sd_y/sd_X)
    res$V_beta <- fit$V_beta * tcrossprod(sd_y/sd_X)
    res$sigma2 <- fit$sigma2 * sd_y^2
    res$t_stat <- res$beta/res$se_beta
    res$CI <- data.frame(lower=res$beta+qnorm(alpha/2)*res$se_beta,
                         upper=res$beta+qnorm(1-alpha/2)*res$se_beta)
#    res$formula <- formula
    names(res$beta) <- colnames_X
    names(res$se_beta) <- colnames_X
    names(res$omega) <- colnames_X
    names(res$kappa) <- colnames_X
    rownames(res$V_beta) <- colnames_X
    colnames(res$V_beta) <- colnames_X
  }
  return(res)
}
