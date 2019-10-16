#' Non-parametric method
#'
#' Calculates weights under I = 1, using the non-parametric method
#'
#' @param Y outcome
#'
#' @return vector of weights under I = 1

non_parm <- function(Y, X, controls, wgt, lambda){
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
