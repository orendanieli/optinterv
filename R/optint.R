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
#' @param n.quant
#'
#' @param seed
#'
#' @return a optint object.

optint <- function(Y, X,
                  controls = NULL,
                  wgt = rep(1, length(Y)),
                  method = "non-parametric",
                  lambda = 100,
                  sigma = 1,
                  n.boot = 1000,
                  sign.factor = 2/3,
                  n.quant = 20,
                  seed = runif(1, 0, .Machine$integer.max),
                  ...){
  #prepare data
  X_std <- apply(X, 2, function(x, w = wgt){x / sqrt(Hmisc::wtd.var(x, w))})
  if(min(Y) <= 0){
    Y_pos <- exp(Y)
  } else {
    Y_pos <- Y
  }
  if (method == "non-parametric"){
    non_parm_boot <- function(d, i){
      w <- non_parm(Y_pos[i], X_std[i,], controls[i], wgt =  wgt[i], lambda =  lambda)
      dists <- apply(X_std[i,], 2, function(v) per_distance(v, n.quant, wgt[i], w))
      diff <- outcome_diff(Y[i], w, wgt[i])
      c(dists, diff)
    }
    res <- boot::boot(1:length(Y), non_parm_boot, n.boot, stype = "i")
    wgt1 <- non_parm(Y_pos, X_std, controls, wgt =  wgt, lambda =  lambda)
  }
  signs <- apply(X_std, 2, function(v) per_distance(v, n.quant, wgt, wgt1, sign.factor, T))
  n <- length(res$t0)
  estimates <- res$t0[-n]
  estimates_sd <- apply(res$t[,-n], 2, sd)
  stand_factor <- sd(estimates)
  output <- list(estimates = estimates / stand_factor,
                 estimates_sd = estimates_sd / stand_factor,
                 details = list(Y_diff = res$t0[n],
                                Y_diff_sd = sd(res$t[n]),
                                method = method,
                                lambda = lambda,
                                new_sample = cbind(X,controls,wgt1),
                                signs = signs,
                                stand_factor = stand_factor,
                                kl_distance = NULL))
  class(output) <- "optint"
  return(output)
}


