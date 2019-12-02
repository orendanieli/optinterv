#' Outcome difference
#'
#' Calculates the diffrence between E(Y|I=1) and E(Y|I=0)
#'
#' @param wgt1 weights under I = 1
#' @inheritParams optint
#'
#' @return outcome difference
#'
outcome_diff <- function(Y, wgt1, wgt = rep(1, length(Y))){
  Y0 <- weighted.mean(Y, wgt)
  Y1 <- weighted.mean(Y, wgt1)
  return(Y1 - Y0)
}

#' Distance Between Distribution
#'
#' Calculate distance in RMSE between quantiles of distribution
#'
#' @param x variable.
#' @param n.quant number of quantiles.
#' @param wgt original weights.
#' @param wgt1 weights under I = 1.
#' @param p proportion of quantiles that need to be increase (decrease) in order to
#'          return a positive (negative) sign.
#' @param plot.sign if F returns RMSE, if T returns the sign of effect.
#'
#' @return scalar for distance. If sign = TRUE, returns +1 (-1) for increasing (decreasing)
#'         p of quantiles. Else returns 0

per_distance <- function(x, n.quant,
                         wgt,
                         wgt1,
                         p = 2/3,
                         plot.sign = F,
                         ...) {
  #fix weights to sum to data length
  wgt  <- (wgt / sum(wgt)) * length(wgt)
  wgt1 <- (wgt1 / sum(wgt1)) * length(wgt1)
  #get quantiles
  q0 <- wtd_bin(x, n.quant, wgt)
  q1 <- wtd_bin(x, n.quant, wgt1)
  distance <- mean(abs(q1 - q0))
  if (!plot.sign){
    return(distance)
  } else {
    sgn <- sign(q1 - q0)
    tot <- sum(sgn) / n.quant
    cutoff <- p - (1 - p)
    res <- ifelse(tot > cutoff, 1, ifelse(tot < -cutoff, -1 ,0))
    return(res)
  }
}

#' Weighted Bin
#'
#' Divide the data into equal size bins and calculate mean
#'
#' @param wgt vector of weigts.
#' @inheritParams per_distance
#'
#' @return vector of means.

wtd_bin <- function(x, n.quant, wgt){
  #set weights to sum to number of quants
  wgt <- wgt / sum(wgt) * n.quant
  #create a data set of x and weights
  dat <- data.frame(x = x, w = wgt)
  #sort
  dat <- dat[order(x),]
  #get beginning and end of each observation by bins
  dat$end <- cumsum(dat$w)
  dat$start <- c(0,dat$end[-nrow(dat)])
  #weighted average of each bin
  bin_ave <- function(b){
    weighted.mean(dat$x, pmin(pmax(dat$end - b, 0), 1) - pmin(pmax(dat$start - b, 0), 1))
  }
  res <- sapply(0:(n.quant - 1), bin_ave)
  return(res)
}

#' Kullback-Leibler Divergence
#'
#' Calculates Kullback-Liebler Divergence for two weight vectors.
#'
#' @inheritParams per_distance
#'
#' @return scalar of kullback-liebler divergence.

kl_dist_def <- function(wgt, wgt1, ...){
  #make weights sum to 1
  wgt <- wgt / sum(wgt)
  wgt1 <- wgt1 / sum(wgt1)
  return(sum(wgt1 * log(wgt1 / wgt)))
}

#' Kullback-Leibler Divergence
#'
#' Calculates Kullback-Liebler Divergence for two multivariate normal distributions.
#'
#' @inheritParams optint
#' @param ni difference in means (mu1 - mu0)
#'
#' @return scalar of kullback-liebler divergence.

kl_dist_cor <- function(X, wgt, ni, ...){
  vcov <- cov.wt(X, wgt)$cov
  return(t(ni) %*% vcov %*% ni)
}

#' Bootstrap Confidence Intervals
#'
#' Calculates bootstrap confidence intervals for matrix of bootstrap replicates
#'
#' @param boot.res matrix of bootstrap replicates
#' @param alpha significance level
#'
#' @return matrix of confidence intervals

boot_ci <- function(boot.res, alpha = 0.05){
  quant <- alpha / 2
  apply(boot.res, 2, function(x){quantile(x, probs = c(quant, 1 - quant))})
}


#'Permutation test

#'Test the null hypothesis P(X|I=0) = P(X|I=1), using permutation test.

#' @param estimates point estimates of the percentile distance between P(X|I=0) & P(X|I=1).
#' @param n.perm number of permutations to permute from wgt1.
#' @inheritParams optint
#'
#' @return vector of p values.

perm_test <- function(estimates, wgt, wgt1, X, n.quant, n.perm = 1000, ...){
  n <- length(wgt1)
  p <- ncol(X)
  #permute 'n.perm' permutations from wgt1
  perm_w1 <- sapply(1:n.perm, function(x){sample(wgt1, size = n, replace = F)})
  #for each variable, calculate per_distance for each permutation.
  dist_sample <- matrix(nrow = n.perm, ncol = p)
  for (i in 1:p){
    dist_sample[,i] <- apply(perm_w1, 2, function(w1){per_distance(X[,i], n.quant, wgt, w1)})
  }
  #calculate P(estimate < dist_sample) (=p value)
  p_val <- sapply(1:p, function(j){mean(estimates[j] < dist_sample[,j])})
  return(p_val)
}





