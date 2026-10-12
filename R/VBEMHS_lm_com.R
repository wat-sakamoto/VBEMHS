#' VBEMHS_lm.pls
#'
#' Compute penalized least square estimate of beta
#'
#' @noRd
#'
VBEMHS_lm.pls <- function(XX,Xy,prec)
{
  XXA <- XX+diag(prec)
  R <- chol(XXA)
  beta <- backsolve(R,forwardsolve(t(R),Xy))
  # beta <- (XX+A)^{-1} Xy

  return(list(beta=beta, R=R))
}

#' VBEMHS_lm.pls.nltp
#'
#' Compute penalized least square estimate of beta (Sherman-Morrison-Woodbury formula)
#'
#' @noRd
#'
VBEMHS_lm.pls.nltp <- function(n,X,Xy,prec) # in case n is less than p
{
  AX <- (1/prec)*t(X) # A X^t
  AXy <- (1/prec)*Xy  # A X^t y
  XAX <- crossprod(t(X),AX) # X A X^t
  XAXy <- crossprod(AX,Xy)  # X A X^t y

  IXAX <- XAX+diag(1,n)
  R <- chol(IXAX)
  beta <- AXy-tcrossprod(AX,t(backsolve(R,forwardsolve(t(R),XAXy))))
  # beta <- AXy-AX (I+XAX)^{-1} XAXy

  return(list(beta=beta, R=R))
}

#' VBEMHS_lm.beta
#'
#' Update q(beta) in VBEM algorithm
#'
#' @noRd
#'
VBEMHS_lm.beta <- function(n,p,X,XX,Xy,prec)
{
  if (n < p){
    reg <- VBEMHS_lm.pls.nltp(n,X,Xy,prec)
  } else {
    reg <- VBEMHS_lm.pls(XX,Xy,prec)
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
