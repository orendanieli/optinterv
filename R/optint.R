#' Optimal intervention
#'
#' identifies the factors with the greatest potential to increase
#' a pre-specified outcome, using varius methods.
#'
#' @param Y outcome vector (must be numeric without NA's).
#' @param X data frame or matrix of factors to be considered.
#' @param controls data frame or matrix of factors to control for. these are factors
#'                 that we can't consider while looking for the optimal intervention
#'                 (e.g. gender).
#' @param wgt an optional vector of weights.
#' @param method the method to be used. one of "correlation" / "non-parametric" /
#'               "nearest-neighbors".
#' @param lambda the lagrange multiplier. also known as the shadow price of an
#'               intervention.
#' @param sigma distance penalty for the nearest-neighbors method.
#' @param n.boot number of bootstrap replications to use for the standard errors
#'               calculation.
#' @param sign.factor
#'
#' @param seed
#'
#' @return a optint object.

optint <- function(Y, X,
                  controls = NULL,
                  wgt = NULL,
                  method = "non-parametric",
                  lambda = 100,
                  sigma = 1,
                  n.boot = 1000,
                  sign.factor = 2/3,
                  seed = runif(1, 0, .Machine$integer.max),
                  ...){
  #for nn & non-par
  if (min(Y) <= 0){
    Y = exp(Y)
  }
}
