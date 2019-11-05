#functions that implement the 3 methods. not for export


#' Non-parametric method
#'
#' Calculates weights under I = 1, using the non-parametric method
#'
#' @inheritParams optint
#'
#' @return vector of weights under I = 1

non_parm <- function(Y, X, controls = NULL, wgt = rep(1, length(Y)), lambda = 100, ...){
  #basline weights (without controls)
  base_wgt1 <- wgt * Y^(1/lambda)
  if (is.null(controls)){
    return(base_wgt1)
  } else {
    #add a constant vector to controls
    controls <- as.matrix(cbind(rep(1, nrow(controls)), controls))
    #define function to solve, plug in variables
    f <- function(b) {dev_moments(b, base_wgt1, controls, wgt)}
    #get solution
    b0 <- rootSolve::multiroot(f, rep(0,ncol(controls)))$root
    #calculate probs with solution
    wgt1 <- base_wgt1 * exp(controls %*% as.matrix(b0))
    return(wgt1)
  }
}

#' Moment deviation
#'
#' Finds the moment deviation for a given lagrange multiplier
#'
#' @param beta a lagrange multiplier
#' @param base basline weights
#' @param controls controls matrix (with a constant)
#' @param wgt original weights
#'
#' @return vector of moment deviations

dev_moments <- function(beta, base, controls, wgt, ...){
  #calculate probs for a given beta
  p <- base * exp(controls %*% as.matrix(beta))
  #calculate the difference in means between the two distributions
  dev <- apply(as.vector(p - wgt) * controls, 2, sum)
  return(dev)
}

#' Nearest-neighbors method
#'
#' Calculates weights under I = 1, using the nearest-neighbors method
#'
#' @inheritParams optint
#'
#' @return vector of weights under I = 1

nn <- function(Y, X, controls = NULL, wgt = rep(1, length(Y)), lambda = 100, sigma = 1, ...){
  char_matrix <- cbind(X, controls)
  vcov <- cov.wt(char_matrix, wgt)$cov
  #find mehalanobis distance for all pairs
  dist <- distances::distances(char_matrix, normalize = vcov)
  #log weights matrix (log is for precision)
  lwgt_mat <- -.5 * (dist^2) * (sigma^-2) + (1/lambda) * log(Y) + log(wgt)
  #make max 100 for precision
  lwgt_mat <- lwgt_mat - max(lwgt_mat) + 100
  wgt_mat <- exp(lwgt_mat)
  #now make wgt_mat_ij sum to wgt for every column (i) (j is the raw index)
  #first, make each column (i) sum to 1
  si <- apply(wgt_mat, 2, sum)
  wgt_mat <-  t(t(wgt_mat) / si)
  #second, make each column sum to the sampling weight of that column (i)
  wgt_mat <- t(t(wgt_mat) * wgt)
  #sum rows (over i) and that's the overall weight on each observation
  wgt1 = apply(wgt_mat, 1, sum)
  return(wgt1)
}
