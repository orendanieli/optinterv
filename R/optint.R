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
#' @param perm.test logical. if TRUE (default) performs permutation test and calculates p-values.
#' @param n.perm number of permutations for the permutation test.
#' @param quick logical. if TRUE, returns only \eqn{E(X | I=1) - E(X | I=0)} as an estimate.
#'              this estimate is used by \code{\link{optint_by_group}}.
#' @param plot logical. if TRUE (default), the results are plotted by \code{\link{plot.optint}}.
#' @param seed the seed of the random number generator.
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
#' In addition, the function \code{\link[base]{summary}} can be used to
#' print a summary of the results.
#'
#' @examples
#' # generate data
#' n <- 50
#' p <- 3
#' features <- matrix(rnorm(n*p), ncol = p)
#' men <- matrix(rbinom(n, 1, 0.5), nrow = n)
#' outcome <- 2*(features[,1] > 1) + men*pmax(features[,2], 0) + rnorm(n)
#' outcome <- as.vector(outcome)
#'
#' #find the optimal intervention using the non-parametric method:
#' imp_feat <- optint(Y = outcome, X = features, control = men,
#'                    method = "non-parametric", lambda = 10, plot = TRUE,
#'                    n.boot = 100, n.perm = 100)
#'
#' #by default, only the significant features are displayed
#' #(see ?plot.optint for further details).
#' #for customized variable importance plot, use plot():
#' plot(imp_feat, plot.vars = 3)
#'
#' #show summary of the results using summary():
#' summary(imp_feat)
#' @export
#' @importFrom stats sd lm predict runif


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
                  n.quant = floor(length(Y) / 10),
                  perm.test = T,
                  n.perm = 1000,
                  quick = F,
                  plot = T,
                  seed = runif(1, 0, .Machine$integer.max)){
  validate_data(Y, X, control, wgt = wgt)
  #create var names if missing
  n <- ncol(X)
  var_names <- colnames(X)
  if(is.null(var_names)){
    var_names <- as.character(rep(1:n))
    colnames(X) <- var_names
  }
  X <- as.matrix(X)
  #replace "method" with its abbreviation (crr / nn- / nr-)
  method <- as.character(abbreviate(method, 3))
  set.seed(seed)
  if (method != "crr"){
    #prepare data
    X_std <- apply(X, 2, function(x, w = wgt){x / sqrt(Hmisc::wtd.var(x, w))})
    Y_pos <- prepare_Y(Y)
    func <- ifelse(method == "nn-", "non_parm", "nn")
    res <- boot_default(func, Y, Y_pos, X, X_std, control, wgt, n.quant,
                         lambda, sigma, grp.size, n.boot, quick)
    estimates <- res$t0[-(n + 1)]
    if(quick){
      estimates_sd <- apply(res$t[,-(n + 1), drop = F], 2, sd)
      return(data.frame(estimates = estimates,
                        estimates_sd = estimates_sd))
    }
    wgt1 <- do.call(func, list(Y_pos, X_std, control, wgt =  wgt,
                               lambda =  lambda, sigma = sigma, grp.size = grp.size))
    signs <- apply(X_std, 2,
                   function(v) per_distance(v, n.quant, wgt, wgt1, sign.factor, T))
    kl_distance <- kl_dist_def(wgt, wgt1)
    estimates_sd <- apply(res$t[,-(n + 1), drop = F], 2, sd)
    ci <- boot_ci(res$t[,-(n + 1), drop = F], alpha)
    #estimates = apply(X_std, 2, function(v) per_distance(v, n.quant, wgt, wgt1))#
    if(perm.test){
      p_val <- perm_test(estimates, wgt, wgt1, X_std, n.quant, n.perm,
                         Y_pos, control, func)
    } else {
      p_val <- rep(NA, n)
    }
    Y_diff = res$t0[n + 1]
    Y_diff_sd = sd(res$t[,n + 1])
    #return(p_val)#
  } else {
    #correlations method
    #transform to matrix for lm:
    if(!is.null(control)){
      control <- as.matrix(control)
    }
    if(min(Y) >= 0 & length(unique(Y)) > 2){
      #transform to log
      pos_ind <- Y > 0
      if(!all(pos_ind)){
        Y <- Y[pos_ind]
        X <- X[pos_ind,,drop = F]
        control <- control[pos_ind,,drop = F]
        wgt <- wgt[pos_ind]
        warning("excluding observations with Y = 0")
      }
      Y <- log(Y)
    }
    res <- par_cor(Y, X, control, wgt)
    estimates <- res$correlation
    estimates_sd <- res$std.err
    if(quick){
      return(data.frame(estimates = estimates, estimated_sd = estimates_sd))
    }
    ni <- (1/lambda) * res$covariance
    kl_distance <- kl_dist_cor(X, wgt, ni)
    #outcome difference:
    ols <- lm(Y ~ cbind(X, control), weights = wgt)
    y_hat <- predict(ols)
    #make E(y_hat) = 0
    y_hat <- y_hat - weighted.mean(y_hat, wgt)
    n_obs <- length(Y)
    Y_diff <- (1 / lambda) * weighted.mean(y_hat^2, wgt)
    Y_diff_sd <- sqrt(Hmisc::wtd.var(y_hat^2, wgt) / (n_obs * lambda^2))
    #update X
    X <- t(t(X) + ni)
    wgt1 <- wgt
    signs <- sign(estimates)
    p_val <- res$p.value
    ci <- cor_ci(estimates, length(Y), alpha)
  }
  #we dont need standardization factor if we have only one variable
  stand_factor <- ifelse(n == 1, 1, sd(estimates))
  new_sample <- cbind(X, control, wgt, wgt1) #return also wgt for plot_change()
  method <- ifelse(method == "crr", "correlations",
                   ifelse(method == "nn-", "non-parametric", "nearest-neighbors"))
  output <- list(estimates = estimates / stand_factor,
                 estimates_sd = estimates_sd / stand_factor,
                 details = list(Y_diff = Y_diff,
                                Y_diff_sd = Y_diff_sd,
                                method = method,
                                lambda = lambda,
                                signs = signs,
                                p_value = p_val,
                                ci = ci / stand_factor,
                                stand_factor = stand_factor,
                                kl_distance = kl_distance,
                                new_sample = new_sample))
  if (method == "nr-")
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
#' @param plot logical. if TRUE (default), the results are plotted by
#'             \code{\link{plot.optint_by_group}}.
#' @inheritParams optint
#'
#' @return an object of class "optint_by_group". This object is a list containing
#'         two components:
#'  \item{est}{a matrix of estimates (in their original units), for each group.
#'             here estimates are \eqn{E(X | I=1) - E(X | I=0)}, and they are
#'             used by \code{\link{plot.optint_by_group}}.}
#'  \item{sd}{estimates standard deviation.}
#'
#' @examples
#' # generate data
#' n <- 50
#' p <- 3
#' features <- matrix(rnorm(n*p), ncol = p)
#' men <- matrix(rbinom(n, 1, 0.5), nrow = n)
#' outcome <- 2*(features[,1] > 1) + men*pmax(features[,2], 0) + rnorm(n)
#' outcome <- as.vector(outcome)
#'
#' #find the optimal intervention using the non-parametric method:
#' imp_feat <- optint(Y = outcome, X = features, control = men,
#'                    method = "non-parametric", lambda = 10, plot = TRUE,
#'                    n.boot = 100, n.perm = 100)
#'
#' #we can explore how the optimal intervention varies between genders using optint_by_group():
#' men <- as.vector(men)
#' imp_feat_by_gender <- optint_by_group(Y = outcome, X = features,
#'                                       group = men,
#'                                       method = "non-parametric",
#'                                       lambda = 10)
#'
#' #by default, only the significant features are displayed
#' #(see ?plot.optint_by_group for further details).
#' #for customized variable importance plot, use plot():
#' plot(imp_feat_by_gender, plot.vars = 3)
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
                            alpha = 0.05,
                            plot = TRUE){
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
  if(plot){
    plot(output, alpha = alpha)
  }
  return(output)
}


























