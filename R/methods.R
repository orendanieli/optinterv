#functions that implement the 3 methods. not for export


#' Non-parametric method
#'
#' Calculates weights under I = 1, using the non-parametric method
#'
#' @inheritParams optint
#'
#' @return vector of weights under I = 1

non_parm <- function(Y, X, control = NULL, wgt = rep(1, length(Y)), lambda = 100, ...){
  #basline weights (without control)
  base_wgt1 <- wgt * Y^(1/lambda)
  if (is.null(control)){
    return(base_wgt1)
  } else {
    return(wgt_adjust(control, base_wgt1, wgt))
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
wgt_adjust <- function(control, base.wgt1, wgt, ...){
  #add a constant vector to control
  control <- as.matrix(cbind(rep(1, nrow(control)), control))
  #define function to solve, plug in variables
  f <- function(b) {dev_moments(b, base.wgt1, control, wgt)}
  #get solution
  b0 <- rootSolve::multiroot(f, rep(0,ncol(control)))$root
  #calculate probs with solution
  wgt1 <- base.wgt1 * exp(control %*% as.matrix(b0))
  return(as.vector(wgt1))
}


#' Moment deviation
#'
#' Finds the moment deviation for a given lagrange multiplier
#'
#' @param beta a lagrange multiplier
#' @param base basline weights
#' @param control control matrix (with a constant)
#' @param wgt original weights
#'
#' @return vector of moment deviations

dev_moments <- function(beta, base, control, wgt, ...){
  #calculate probs for a given beta
  p <- base * exp(control %*% as.matrix(beta))
  #calculate the difference in means between the two distributions
  dev <- apply(as.vector(p - wgt) * control, 2, sum)
  return(dev)
}


#' Nearest-neighbors weights
#'
#' Calculates unadjusted weights under I = 1, using the nearest-neighbors method
#'
#' @inheritParams optint
#'
#' @return vector of unadjusted weights under I = 1

nn_wgt <- function(Y, X, control = NULL, wgt = rep(1, length(Y)), lambda = 100, sigma = 1, ...){
  char_matrix <- cbind(X, control)
  vcov <- cov.wt(char_matrix, wgt)$cov
  #find mehalanobis distance for all pairs
  dist <- distances::distances(char_matrix, normalize = vcov)
  dist <- as.matrix(dist)
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
  wgt1 <- apply(wgt_mat, 1, sum)
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

nn <- function(Y, X, control = NULL, wgt = rep(1, length(Y)),
               lambda = 100, sigma = 1, grp.size = 30, ...){
  if (is.null(control)){
    return(nn_wgt(Y, X, wgt = wgt, lambda =  lambda, sigma = sigma))
  }
  #check average size of control groups
  n_control <- ncol(control)
  #first, transform control to list:
  control_list <- lapply(seq_len(n_control), function(i) control[,i])
  control_val <- lapply(control_list,  unique)
  n_grps <- prod(unlist(lapply(control_val, length)))
  n_obs <- length(Y)
  if (n_obs / n_grps < grp.size){
    base_wgt1 <- nn_wgt(Y, X, control, wgt, lambda, sigma)
    return(wgt_adjust(control, base_wgt1, wgt))
  }
  #all possible combinations of control groups:
  control_grps <- expand.grid(control_val)
  #check that each control group of size grp.size
  comb <- 1
  too_small <- F
  grp_flag <- NULL
  while(!too_small & comb <= n_grps){
    grp_flag <- cbind(grp_flag, rep(T, n_obs))
    for(i in 1:n_control){
      grp_flag[,comb] <- grp_flag[,comb] & (control[,i] == control_grps[comb,i])
    }
    too_small <- sum(grp_flag[,comb]) < grp.size
    comb <- comb + 1
  }
  if(too_small){
    base_wgt1 <- nn_wgt(Y, X, control, wgt, lambda, sigma)
    return(wgt_adjust(control, base_wgt1, wgt))
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


#' Partial Correlation
#'
#' Calculates correlation / covariance between Y and X, holding control constant
#'
#' @inheritParams optint
#'
#' @return data frame with partial correlations & covariance

par_cor = function(Y, X, control = NULL, wgt = rep(1, length(Y)), ...){
  if(is.null(control)){
    cors <- as.vector(weights::wtd.cors(X, Y, wgt))
    Y_sd <- sd(Y)
    X_sd <- apply(X, 2, sd)
    return(data.frame(cors, covs = Y_sd * X_sd * cors))
  } else {
    #residualize control
    Y <- lm(Y ~ control, weights = wgt)$residuals
    X <- lm(X ~ control, weights = wgt)$residuals
    cors <- as.vector(weights::wtd.cors(X, Y, wgt))
    Y_sd <- sd(Y)
    X_sd <- apply(X, 2, sd)
    return(data.frame(cors, covs = Y_sd * X_sd * cors))
  }
}




