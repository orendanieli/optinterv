#' Optimal intervention
#'
#' identifies the factors with the greatest potential to increase
#' a pre-specified outcome, using varius methods.
#'
#' @param Y outcome vector (must be numeric without NA's).
#' @param X data frame or matrix of factors to be considered.
#' @param controls data frame or matrix of factors to control for. these are factors
#'                 that we can't consider while looking for the optimal intervention
#'                 (e.g. race).
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
#' @return an optint object.
#' @export

optint <- function(Y, X,
                  controls = NULL,
                  wgt = rep(1, length(Y)),
                  method = "non-parametric",
                  lambda = 100,
                  sigma = 1,
                  n.boot = 1000,
                  sign.factor = 2/3,
                  n.quant = length(Y) / 10,
                  seed = runif(1, 0, .Machine$integer.max),
                  ...){
  #prepare data
  X_std <- apply(X, 2, function(x, w = wgt){x / sqrt(Hmisc::wtd.var(x, w))})
  if(min(Y) <= 0){
    Y_pos <- exp(Y)
  } else {
    Y_pos <- Y
  }
  set.seed(seed)
  if (method != "correlations"){
    func <- ifelse(method == "non-parametric", "non_parm", "nn")
    boot_func <- function(d, i){
      w <- do.call(func, list(Y_pos[i], X_std[i,], controls[i,], wgt =  wgt[i], lambda =  lambda))
      dists <- apply(X_std[i,], 2, function(v) per_distance(v, n.quant, wgt[i], w))
      diff <- outcome_diff(Y[i], w, wgt[i])
      c(dists, diff)
    }
    res <- boot::boot(1:length(Y), boot_func, n.boot, stype = "i")
    wgt1 <- do.call(func, list(Y_pos, X_std, controls, wgt =  wgt, lambda =  lambda))
  }
  signs <- apply(X_std, 2, function(v) per_distance(v, n.quant, wgt, wgt1, sign.factor, T))
  n <- length(res$t0)
  estimates <- res$t0[-n]
  estimates_sd <- apply(res$t[,-n], 2, sd)
  stand_factor <- sd(estimates)
  output <- list(estimates = estimates / stand_factor,
                 estimates_sd = estimates_sd / stand_factor,
                 details = list(Y_diff = res$t0[n],
                                Y_diff_sd = sd(res$t[,n]),
                                method = method,
                                lambda = lambda,
                                new_sample = cbind(X, controls, wgt1),
                                signs = signs,
                                stand_factor = stand_factor,
                                kl_distance = kl_dist(wgt, wgt1)))
  class(output) <- "optint"
  return(output)
}


#' Summary for optint object
#'
#' Report results from an optint object
#'
#' @param object an optint object
#' @param r number of decimal places to use
#'
#' @export

summary.optint <- function(object, r = 5){
  x <- object
  est <- round(x$estimates, r)
  se <- round(x$estimates_sd, r)
  kl <- x$details$kl_distance
  out <- round(x$details$Y_diff, r)
  out_sd <- round(x$details$Y_diff_sd, r)
  out_t <- out / out_sd
  out_p <- 2 * pnorm(out_t, lower.tail = F)
  n <- length(est)
  var_names <- colnames(x$details$new_sample[,1:n])
  if(is.null(var_names))
    for(i in 1:(n))
      var_names[i] <- paste0("X",i)
  coeffs <- matrix(c(est, se), ncol = 2,
                   dimnames = list(var_names, c("Estimate","Std. error")))
  out_mat <- matrix(c(out, out_sd, out_t, out_p), ncol = 4,
                      dimnames = list("E(Y|I=1) - E(Y|I=0)",
                                      c("Estimate","Std. error", "t value", "P(>|t|)")))
  method <- x$details$method
  lambda <- x$details$lambda
  cat("Method:", method, ", Lambda =", lambda, "\n")
  coef_title <- ifelse(method == "correlations", "Raw Correlations:", "CDF Distances:")
  cat("\n",coef_title, "\n", "\n")
  print(coeffs)
  cat("\n", "The Kullback–Leibler divergence of P(X|I=0) from P(X|I=1) is:", kl, "\n")
  cat("\n", "Outcome Difference:", "\n", "\n")
  print(out_mat)
}




















