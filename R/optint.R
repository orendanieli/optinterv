#' Optimal intervention
#'
#' identifies the factors with the greatest potential to increase
#' a pre-specified outcome, using varius methods.
#'
#' @param Y outcome vector (must be numeric without NA's).
#' @param X numeric data frame or matrix of factors to be considered.
#' @param control numeric data frame or matrix of factors to control for. these are factors
#'                 that we can't consider while looking for the optimal intervention
#'                 (e.g. race).
#' @param wgt an optional vector of weights.
#' @param method the method to be used. either "non-parametric" (default), "correlation" or
#'               "nearest-neighbors".
#' @param lambda the lagrange multiplier. also known as the shadow price of an
#'               intervention.
#' @param sigma distance penalty for the nearest-neighbors method.
#' @param grp.size for the nearest-neighbors method; if the number of examples in each
#'                 control group is smaller than grp.size, performs weight adjustment
#'                 using \code{\link[optinterv]{wgt_adjust}}. else,
#'                 calculate weights seperatly for each control group.
#' @param n.boot number of bootstrap replications to use for the standard errors /
#'               confidence intervals calculation.
#' @param sign.factor what proportion of quantiles should to be increased (decreased)
#'                    in order to return a positive (negative) sign? not relevant for
#'                    the correlation method (there the correlation sign is returned).
#' @param alpha significance level for the confidence intervals.
#' @param n.quant number of quantiles to use when calculating CDF distance.
#' @param n.perm number of permutations for the permutation test.
#' @param quick logical. if TRUE, returns only \eqn{E(X | I=1) - E(X | I=0)} as an estimate.
#'              this estimate is used by \code{\link{optint_by_group}}.
#' @param plot logical. if TRUE (default), the results are plotted by \code{\link{plot.optint}}.
#' @param seed the seed of the random number generator.
#'
#' @return an object of class "optint". This object is a list containing
#'         the folowing components:
#'  \item{estimates}{standardized point estimates (correlations for the
#'  correlation method and cdf distances otherwise).}
#'  \item{estimates_sd}{estimates standard deviation.}
#'  \item{details}{a list containing further details, such as:}
#'  \itemize{
#'   \item Y_diff - \eqn{E(Y | I=1) - E(Y | I=0)}.
#'   \item Y_diff_sd - standard deviation for Y_diff.
#'   \item method - the method used.
#'   \item lambda - the lagrange multiplier used.
#'   \item signs - signs (i.e. directions) for the estimates.
#'   \item p_value - p-values for the estimates.
#'   \item ci - a matrix of confidence intervals for the estimates.
#'   \item stand_factor - the standardization factor used to standardize the results.
#'   \item kl_distance - the Kullback–Leibler divergence of \eqn{P(X | I=0) from P(X | I=1)}.
#'   \item new_sample - a data frame containing X, control (if provided),
#'         wgt (the original weights) and wgt1 (the new weights under \eqn{I = 1}.)
#'  }
#'
#' In addition, the function \code{\link[Matrix]{summary}} can be used to
#' print a summary of the results.
#'
#' @examples
#'
#' @export

optint <- function(Y, X,
                  control = NULL,
                  wgt = rep(1, length(Y)),
                  method = "non-parametric",
                  lambda = 100,
                  sigma = 1,
                  grp.size = 30,
                  n.boot = 1000,
                  sign.factor = 2/3,
                  alpha = 0.05,
                  n.quant = length(Y) / 10,
                  n.perm = 1000,
                  quick = F,
                  plot = T,
                  seed = runif(1, 0, .Machine$integer.max),
                  ...){
  validate_data(Y, X, control, wgt = wgt)
  #create var names if missing
  n <- ncol(X)
  var_names <- colnames(X)
  if(is.null(var_names)){
    var_names <- as.character(rep(1:n))
    colnames(X) <- var_names
  }
  set.seed(seed)
  if (method != "correlations"){
    #prepare data
    X_std <- apply(X, 2, function(x, w = wgt){x / sqrt(Hmisc::wtd.var(x, w))})
    Y_pos <- prepare_Y(Y)
    func <- ifelse(method == "non-parametric", "non_parm", "nn")
    boot_func <- function(d, i){
      w <- do.call(func, list(Y_pos[i], X_std[i,], control[i,], wgt =  wgt[i],
                              lambda =  lambda, sigma = sigma, grp.size = grp.size))
      if(quick){
        diff <- apply(X[i,], 2, function(v) mean_diff(v, wgt[i], w))
        return(diff)
      }
      dists <- apply(X_std[i,], 2, function(v) per_distance(v, n.quant, wgt[i], w))
      diff <- outcome_diff(Y[i], w, wgt[i])
      return(c(dists, diff))
    }
    res <- boot::boot(1:length(Y), boot_func, n.boot, stype = "i")
    estimates <- res$t0[-(n + 1)]
    if(quick){
      estimates_sd <- apply(res$t[,-(n + 1)], 2, sd)
      return(data.frame(estimates = estimates,
                        estimates_sd = estimates_sd))
    }
    wgt1 <- do.call(func, list(Y_pos, X_std, control, wgt =  wgt,
                               lambda =  lambda, sigma = sigma, grp.size = grp.size))
    signs <- apply(X_std, 2,
                   function(v) per_distance(v, n.quant, wgt, wgt1, sign.factor, T))
    kl_distance <- kl_dist_def(wgt, wgt1)
    #estimates = apply(X_std, 2, function(v) per_distance(v, n.quant, wgt, wgt1))#
    p_val <- perm_test(estimates, wgt, wgt1, X, n.quant, n.perm)
    #return(p_val)#
  } else {
    if(!is.matrix(X))
      X <- as.matrix(X)
    if(quick){
      corrs <- par_cor(Y, X, control, wgt)
      return(data.frame(estimates = corrs$correlation, estimated_sd = corrs$std.err))
    }
    boot_func <- function(d, i){
      cor_cov <- par_cor(Y[i], X[i,], control[i,], wgt[i])
      covs <- cor_cov$covariance
      cors <- cor_cov$correlation
      betas <- lm(Y[i] ~ cbind(X[i,], control[i,]), weights = wgt[i])$coefficients
      diff <- (1/lambda) * (betas[2:(n+1)] %*% covs)
      c(cors, diff)
    }
    res <- boot::boot(1:length(Y), boot_func, n.boot, stype = "i")
    point_est <- par_cor(Y, X, control, wgt)
    ni <- (1/lambda) * point_est$covariance
    kl_distance <- kl_dist_cor(X, wgt, ni)
    X <- t(t(X) + ni)
    wgt1 <- wgt
    estimates <- point_est$correlation
    signs <- sign(estimates)
    p_val <- point_est$p.value
  }
  estimates_sd <- apply(res$t[,-(n + 1)], 2, sd)
  ci <- boot_ci(res$t[,-(n + 1)], alpha)
  stand_factor <- sd(estimates)
  new_sample <- cbind(X, control, wgt, wgt1) #return also wgt for plot_change()
  output <- list(estimates = estimates / stand_factor,
                 estimates_sd = estimates_sd / stand_factor,
                 details = list(Y_diff = res$t0[n + 1],
                                Y_diff_sd = sd(res$t[,n + 1]),
                                method = method,
                                lambda = lambda,
                                signs = signs,
                                p_value = p_val,
                                ci = ci / stand_factor,
                                stand_factor = stand_factor,
                                kl_distance = kl_distance,
                                new_sample = new_sample))
  if (method == "nearest-neighbors")
    output[["details"]][["sigma"]] <- sigma
  class(output) <- "optint"
  if(plot){
    plot(output, alpha = alpha)
  }
  return(output)
}

#' Optimal intervention, by group
#'
#' Similar to \code{\link{optint}}, identifies the factors with the greatest
#' potential to increase a pre-specified outcome for each group separately, and thus allowing
#' to detect heterogeneity between groups.
#'
#' @param group vector with group labels (i.e. grouping variable). the function
#'              \code{\link{optint}} implemented for each group separately.
#' @inheritParams optint
#'
#' @return an object of class "optint_by_group". This object is a list containing
#'         two components:
#'  \item{est}{a matrix of estimates (in their original units), for each group.
#'             here estimates are \eqn{E(X | I=1) - E(X | I=0)}, and they are
#'             used by \code{\link{plot.optint_by_group}}.}
#'  \item{sd}{estimates standard deviation.}
#'
#' @export

optint_by_group <- function(Y, X, group,
                            control = NULL,
                            wgt = rep(1, length(Y)),
                            method = "non-parametric",
                            lambda = 100,
                            sigma = 1,
                            grp.size = 30,
                            n.boot = 1000,
                            alpha = 0.05){
  validate_data(Y, X, control, wgt = wgt)
  validate_group(Y, group)
  group_names <- unique(group)
  var_names <- colnames(X)
  n_vars <- ncol(X)
  if(is.null(var_names)){
    var_names <- as.character(rep(1:n_vars))
  }
  #initialize empty matrixs for the output:
  estimates <- matrix(NA, nrow = n_vars, ncol = length(group_names),
                      dimnames = list(var_names, as.character(group_names)))
  sd <- estimates
  #implement optint() for each group separately.
  for(g in group_names){
    gr_inc <- which(group == g)
    res <- do.call("optint", list(Y[gr_inc], X[gr_inc,], control[gr_inc,],
                                     wgt[gr_inc], method, lambda, sigma, grp.size,
                                     n.boot, quick = T))
    estimates[,paste(g)] <- res$estimates
    sd[,paste(g)] <- res$estimates_sd
  }
  output <- list(est = estimates, sd = sd)
  class(output) <- "optint_by_group"
  plot(output, alpha = alpha)
  return(output)
}






#' Summary for optint object
#'
#' Report results from an optint object.
#'
#' @param object an optint object.
#' @param r number of decimal places to use.
#'
#' @export

summary.optint <- function(object, r = 5){
  x <- object
  est <- round(x$estimates, r)
  se <- round(x$estimates_sd, r)
  p_val <- round(x$details$p_value, r)
  kl <- x$details$kl_distance
  out <- round(x$details$Y_diff, r)
  out_sd <- round(x$details$Y_diff_sd, r)
  out_t <- out / out_sd
  out_p <- 2 * pnorm(abs(out_t), lower.tail = F)
  n <- length(est)
  var_names <- colnames(x$details$new_sample[,1:n])
  #add signs to var names:
  var_names <- add_sign(var_names, x$details$signs)
  coeffs <- matrix(c(est, se, p_val), ncol = 3,
                   dimnames = list(var_names, c("Estimate","Std. error","P-Value")))
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




















