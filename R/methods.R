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
    return(wgt_adjust(controls, base_wgt1, wgt))
  }
}

#' Weights adjustment
#'
#' Adjust new weights so that the distribution of the control variables doesn't change
#'
#' @param base.wgt1 baseline weights under I = 1
#' @inheritParams optint
#'
#' @return vector of adjusted weights under I = 1
wgt_adjust <- function(controls, base.wgt1, wgt, ...){
  #add a constant vector to controls
  controls <- as.matrix(cbind(rep(1, nrow(controls)), controls))
  #define function to solve, plug in variables
  f <- function(b) {dev_moments(b, base.wgt1, controls, wgt)}
  #get solution
  b0 <- rootSolve::multiroot(f, rep(0,ncol(controls)))$root
  #calculate probs with solution
  wgt1 <- base.wgt1 * exp(controls %*% as.matrix(b0))
  return(wgt1)
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


#' Nearest-neighbors weights
#'
#' Calculates unadjusted weights under I = 1, using the nearest-neighbors method
#'
#' @inheritParams optint
#'
#' @return vector of unadjusted weights under I = 1

nn_wgt <- function(Y, X, controls = NULL, wgt = rep(1, length(Y)), lambda = 100, sigma = 1, ...){
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

#' Nearest-neighbors method
#'
#' Calculates adjusted weights under I = 1, using the nearest-neighbors method
#'
#' @param grp.size if number of observations in each control group is smaller than grp.size,
#'                 perform weights adjustment the same way as with the non-parametric method.
#'                 else, calculate weights seperatly for each control group.
#' @inheritParams optint
#'
#' @return vector of adjusted weights under I = 1

nn <- function(Y, X, controls = NULL, wgt = rep(1, length(Y)),
               lambda = 100, sigma = 1, grp.size = 30, ...){
  if (is.null(controls)){
    return(nn_wgt(Y, X, wgt = wgt, lambda =  lambda, sigma = sigma))
  }
  #check average size of control groups
  n_controls <- ncol(controls)
  #first, transform controls to list:
  control_list <- lapply(seq_len(n_controls), function(i) controls[,i])
  control_val <- lapply(control_list,  unique)
  n_grps <- prod(unlist(lapply(control_val, length)))
  n_obs <- length(Y)
  if (n_obs / n_grps < grp.size){
    base_wgt1 <- nn_wgt(Y, X, controls, wgt, lambda, sigma)
    return(wgt_adjust(controls, base_wgt1, wgt))
  }
  #all possible combinations of control groups:
  control_grps <- expand.grid(control_val)
  #check that each control group of size grp.size
  comb <- 1
  too_small <- F
  grp_flag <- NULL
  while(!too_small & comb <= n_grps){
    grp_flag <- cbind(grp_flag, rep(T, n_obs))
    for(i in 1:n_controls){
      grp_flag[,comb] <- grp_flag[,comb] & (controls[,i] == control_grps[comb,i])
    }
    too_small <- sum(grp_flag[,comb]) < grp.size
    comb <- comb + 1
  }
  if(too_small){
    base_wgt1 <- nn_wgt(Y, X, controls, wgt, lambda, sigma)
    return(wgt_adjust(controls, base_wgt1, wgt))
  } else {
    #transform grp_flag to list
    grp_flag <- lapply(seq_len(n_grps), function(i) grp_flag[,i])
    wgt_list <- lapply(grp_flag,
                       function(inc){nn_wgt(Y[inc], X[inc,], wgt =  wgt[inc],
                                            lambda = lambda, sigma = sigma)})
    wgt1 <- rep(0, n_obs)
    for(i in 1:n_grps){
      wgt1[grp_flag[[i]]] <- wgt_list[[i]]
    }
    return(wgt1)
  }
}






