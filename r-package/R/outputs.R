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
  if(min(Y) >= 0){
    #transform to log scale
    pos_ind <- Y > 0
    Y <- Y[pos_ind]
    wgt <- wgt[pos_ind]
    wgt1 <- wgt1[pos_ind]
    Y <- log(Y)
  }
  Y0 <- weighted.mean(Y, wgt)
  Y1 <- weighted.mean(Y, wgt1)
  return(Y1 - Y0)
}


#' Distance Between Distributions
#'
#' Calculate distance in RMSE between quantiles of distributions
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
                         plot.sign = F) {
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
    #dont consider quantiles that weren't changed
    sgn <- sgn[sgn != 0]
    tot <- mean(sgn)
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
#' @importFrom stats weighted.mean

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

kl_dist_def <- function(wgt, wgt1){
  #make weights sum to 1
  wgt <- wgt / sum(wgt)
  wgt1 <- wgt1 / sum(wgt1)
  #if wgt1 == 0, the contribution of the corresponding term is 0.
  lratio <- ifelse(wgt1 == 0, 0, log(wgt1 / wgt))
  return(sum(wgt1 * lratio))
}

#' Kullback-Leibler Divergence
#'
#' Calculates Kullback-Liebler Divergence for two multivariate normal distributions.
#'
#' @inheritParams optint
#' @param ni difference in means (mu1 - mu0)
#'
#' @return scalar of kullback-liebler divergence.
#' @importFrom stats cov.wt

kl_dist_cor <- function(X, wgt, ni){
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
#' @importFrom stats quantile

boot_ci <- function(boot.res, alpha = 0.05){
  quant <- alpha / 2
  ci <- apply(boot.res, 2, function(x){quantile(x, probs = c(quant, 1 - quant))})
  return(ci)
}


cor_ci <- function(estimates, n, alpha = 0.05){
  z <- 0.5 * log((1 + estimates) / (1 - estimates))
  s <- sqrt(1 / (n - 3))
  up_ci <- z + qnorm(1 - alpha/2)*s
  low_ci <- z - qnorm(1 - alpha/2)*s
  up_ci <- (exp(2*up_ci) - 1) / (exp(2*up_ci) + 1)
  low_ci <- (exp(2*low_ci) - 1) / (exp(2*low_ci) + 1)
  return(rbind(low_ci, up_ci))
}

#' Permutation test

#' Test the null hypothesis P(X|I=0) = P(X|I=1), using permutation test.

#' @param estimates point estimates of the percentile distance between P(X|I=0) & P(X|I=1).
#' @param n.perm number of permutations to permute from wgt1.
#' @param func either "non_parm" or "nn". for "nn", weights are recalculated for each
#'             permutation, and thus Y and control are needed. the default is "non_parm",
#'             and Y and control aren't needed.
#' @param wgt1 weights under I = 1.
#' @inheritParams optint
#'
#' @return vector of p values.
#' @importFrom pbapply setpb startpb closepb
#' @importFrom boot boot

perm_test <- function(estimates, wgt, wgt1, X, n.quant, n.perm = 1000,
                      Y = NULL, control = NULL, func = "non_parm",
                      lambda = 100, sigma = 1, grp.size = 30){
  p <- ncol(X)
  if(func == "nn"){
    perm_func <- function(d, i){
      #print progress:
      setpb(pb, rep_count)
      rep_count <<- rep_count + 1
      wgt1 <- nn(Y, d[i,,drop = F], control, wgt, lambda = lambda,
                 sigma = sigma, grp.size = grp.size)
      apply(d[i,,drop = F], 2,
            function(x) per_distance(x, n.quant, wgt, wgt1))}
    } else {
      perm_func <- function(d, i){
        #print progress:
        setpb(pb, rep_count)
        rep_count <<- rep_count + 1
        apply(d[i,,drop = F], 2,
              function(x) per_distance(x, n.quant, wgt, wgt1))}
  }
  #necessary variable for pbapply:
  rep_count <- 1
  message("Calculating p-value:", "\n")
  #start report:
  pb <- startpb(min = 0, max = n.perm)
  res <- boot(X, perm_func, sim = "permutation", n.perm, stype = "i")
  closepb(pb)
  #necessary for pbapply:
  rep_count <- 1
  p_val <- rep(NA, p)
  for(i in 1:p){
    p_val[i] <- mean(estimates[i] < res$t[,i])
  }
  return(p_val)
}

#' @importFrom stats weighted.mean
mean_diff <- function(x, wgt, wgt1){
  diff <- weighted.mean(x, wgt1) - weighted.mean(x, wgt)
  return(diff)
}

#' Bootstrap (default)

#' Bootstrap function for the non-parametric and the nearest neighbor methods
#'
#' @param func a function for weights calculation (nn / non_parm).
#' @param Y the original outcome.
#' @param Y_pos outcome after exponential transformation (if needed).
#' @param X the original X matrix.
#' @param X_std X matrix after standardization.
#' @inheritParams optint
#'
#' @return a list - the output from the function 'boot()'.
#' @importFrom pbapply setpb startpb closepb
#' @importFrom boot boot

boot_default <- function(func, Y, Y_pos, X, X_std, control, wgt, n.quant,
                         lambda, sigma, grp.size, n.boot, quick){
  #create bootstrap function
  boot_func <- function(d, i){
    #print progress:
    setpb(pb, rep_count)
    rep_count <<- rep_count + 1
    w <- do.call(func, list(Y_pos[i], X_std[i,,drop = F], control[i,,drop = F],
                            wgt =  wgt[i], lambda =  lambda, sigma = sigma,
                            grp.size = grp.size))
    if(quick){
      diff <- apply(X[i,,drop = F], 2, function(v) mean_diff(v, wgt[i], w))
      return(diff)
    }
    dists <- apply(X_std[i,,drop = F], 2,
                   function(v) per_distance(v, n.quant, wgt[i], w))
    diff <- outcome_diff(Y[i], w, wgt[i])
    return(c(dists, diff))
  }
  #necessary variable for pbapply:
  rep_count <-  1
  message("Calculating standard errors:", "\n")
  #start report:
  pb <- startpb(min = 0, max = n.boot)
  res <- boot(1:length(Y), boot_func, n.boot, stype = "i")
  closepb(pb)
  #necessary for pbapply:
  rep_count <- 1
  return(res)
}



