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
  #add stars as additional column:
  star <- gen_star(p_val)
  coeffs <- matrix(c(est, se, p_val, star), ncol = 4,
                   dimnames = list(var_names, c("Estimate","Std. error","P-Value"," ")))
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
  print(as.data.frame(coeffs))
  cat("\n", "The Kullback–Leibler divergence of P(X|I=0) from P(X|I=1) is:", kl, "\n")
  cat("\n", "Outcome Difference:", "\n", "\n")
  print(out_mat)
}


gen_star <- function(p_val){
  stars <- rep("", length(p_val))
  stars[p_val >= 0.05 & p_val < 0.1] <- "*"
  stars[p_val >= 0.01 & p_val < 0.05] <- "**"
  stars[p_val < 0.01] <- "***"
  return(stars)
}





