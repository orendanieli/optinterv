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
#' @param grp.size
#' @param n.boot number of bootstrap replications to use for the standard errors
#'               calculation.
#' @param sign.factor
#' @param alpha
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
                  grp.size = 30,
                  n.boot = 1000,
                  sign.factor = 2/3,
                  alpha = 0.05,
                  n.quant = length(Y) / 10,
                  seed = runif(1, 0, .Machine$integer.max),
                  ...){
  n <- ncol(X)
  set.seed(seed)
  if (method != "correlations"){
    #prepare data
    X_std <- apply(X, 2, function(x, w = wgt){x / sqrt(Hmisc::wtd.var(x, w))})
    if(min(Y) <= 0){
      Y_pos <- exp(Y)
    } else {
      Y_pos <- Y
    }
    func <- ifelse(method == "non-parametric", "non_parm", "nn")
    boot_func <- function(d, i){
      w <- do.call(func, list(Y_pos[i], X[i,], controls[i,], wgt =  wgt[i],
                              lambda =  lambda, sigma = sigma, grp.size = grp.size))
      dists <- apply(X_std[i,], 2, function(v) per_distance(v, n.quant, wgt[i], w))
      diff <- outcome_diff(Y[i], w, wgt[i])
      c(dists, diff)
    }
    res <- boot::boot(1:length(Y), boot_func, n.boot, stype = "i")
    wgt1 <- do.call(func, list(Y_pos, X_std, controls, wgt =  wgt,
                               lambda =  lambda, sigma = sigma, grp.size = grp.size))
    signs <- apply(X_std, 2, function(v) per_distance(v, n.quant, wgt, wgt1, sign.factor, T))
    kl_distance <- kl_dist.def(wgt, wgt1)
  } else {
    boot_func <- function(d, i){
      cor_cov <- par_cor(Y[i], X[i,], controls[i,], wgt[i])
      covs <- cor_cov$covs
      cors <- cor_cov$cors
      betas <- lm(Y[i] ~ X[i,] + controls[i,], weights = wgt[i])$coefficients
      diff <- (1/lambda) * (betas[2:(n+1)] %*% covs)
      c(cors, diff)
    }
    res <- boot::boot(1:length(Y), boot_func, n.boot, stype = "i")
    ni <- (1/lambda) * par_cor(Y, X, controls, wgt)$covs
    kl_distance <- kl_dist.cor(X, wgt, ni)
    X <- t(t(X) + ni)
    wgt1 <- wgt
    signs <- sign(res$t0[-(n + 1)])
  }
  estimates <- res$t0[-(n + 1)]
  estimates_sd <- apply(res$t[,-(n + 1)], 2, sd)
  ci <- boot_ci(res$t[,-(n + 1)], alpha)
  stand_factor <- sd(estimates)
  new_sample <- cbind(X, controls, wgt1)
  output <- list(estimates = estimates / stand_factor,
                 estimates_sd = estimates_sd / stand_factor,
                 details = list(Y_diff = res$t0[n + 1],
                                Y_diff_sd = sd(res$t[,n + 1]),
                                method = method,
                                lambda = lambda,
                                signs = signs,
                                ci = ci / stand_factor,
                                stand_factor = stand_factor,
                                kl_distance = kl_distance))
  if (method == "nearest-neighbors")
    output[["details"]][["sigma"]] <- sigma
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
  out_p <- 2 * pnorm(abs(out_t), lower.tail = F)
  n <- length(est)
  var_names <- colnames(x$details$new_sample[,1:n])
  coeffs <- matrix(c(est, se), ncol = 2,
                   dimnames = list(var_names, c("Estimate","Std. error")))
  out_mat <- matrix(c(out, out_sd, out_t, out_p), ncol = 4,
                      dimnames = list("E(Y|I=1) - E(Y|I=0)",
                                     c("Estimate","Std. error", "t value", "P(>|t|)")))
  method <- x$details$method
  lambda <- x$details$lambda
  cat("Method:", method, ", Lambda =", lambda)
  if (method == "nearest-neighbors")
    cat(", Sigma =", x$details$sigma)
  coef_title <- ifelse(method == "correlations", "Raw Correlations:", "CDF Distances:")
  cat("\n", "\n",coef_title, "\n", "\n")
  print(coeffs)
  cat("\n", "The Kullback–Leibler divergence of P(X|I=0) from P(X|I=1) is:", kl, "\n")
  cat("\n", "Outcome Difference:", "\n", "\n")
  print(out_mat)
}




















