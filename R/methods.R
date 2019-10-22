#functions that implement the 3 methods. not for export


#' Non-parametric method
#'
#' Calculates weights under I = 1, using the non-parametric method
#'
#' @param Y outcome
#'
#' @return vector of weights under I = 1

non_parm <- function(Y, X, controls, wgt = NULL, lambda = 100){
  if (is.null(wgt)){
    wgt <- rep(1,length(Y))
  }
  #basline weight (without controls)
  base_wgt1 <- wgt * Y^(1/lambda)
  #add a constant vector to controls
  controls <- as.matrix(cbind(rep(1,nrow(controls)), controls))
  #define function to solve, plug in variables
  f <- function(b) {dev_moments(b, base_wgt1, controls, wgt)}
  #get solution
  b0 <- rootSolve::multiroot(f, rep(0,ncol(controls)))$root
  #calculate probs with solution
  wgt1 <- base_wgt1 * exp(controls %*% as.matrix(b0))
  return(wgt1)
}

#' Moment deviation
#'
#' Find the moment deviation for a given lagrange multiplier
#'
#' @param beta a lagrange multiplier
#' @param base basline weights
#' @param controls controls matrix (with a constant)
#' @param wgt original weights
#'
#' @return vector of moment deviations

dev_moments <- function(beta, base, controls, wgt){
  #calculate probs for a given beta
  p <- base * exp(controls %*% as.matrix(beta))
  #calculate the difference in means between the two distributions
  dev <- apply(as.vector(p - wgt) * controls, 2, sum)
  return(dev)
}
