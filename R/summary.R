#' Summary for optint object
#'
#' Report results from an optint object.
#'
#' @param object an optint object.
#' @param r number of decimal places to use.
#'
#' @export
#' @importFrom stats pnorm

summary.optint <- function(object, r = 4){
  x <- object
  est <- round(x$estimates, r)
  se <- round(x$estimates_sd, r)
  p_val <- round(x$details$p_value, r)
  #add stars as additional column:
  star <- gen_star(p_val)
  #truncate p values:
  p_val[p_val == 0] <- paste0("<", 10^-r, sep = "")
  kl <- x$details$kl_distance
  out <- round(x$details$Y_diff, r)
  out_sd <- round(x$details$Y_diff_sd, r)
  out_t <- round(out / out_sd , r)
  out_p <- round(2 * pnorm(abs(out_t), lower.tail = F), r)
  #truncate p value:
  out_p <- ifelse(out_p == 0, paste0("<", 10^-r, sep = "") , out_p)
  n <- length(est)
  var_names <- colnames(x$details$new_sample[,1:n])
  #add signs to var names:
  var_names <- add_sign(var_names, x$details$signs)
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
  cat("---", "\n", "(1) Signif. codes: 0.1 '\u002A' 0.05 '\u002A\u002A' 0.01 '\u002A\u002A\u002A'")
  cat("\n", "(2) The Kullback–Leibler divergence of P(X|I=0) from P(X|I=1) is:", kl, "\n")
  cat("\n", "Outcome Difference (excluding zeros):", "\n", "\n")
  print(as.data.frame(out_mat))
  cat("---", "\n", "Note: Y is in log units")
}


gen_star <- function(p_val){
  stars <- rep("", length(p_val))
  stars[p_val >= 0.05 & p_val < 0.1] <- "\u002A"
  stars[p_val >= 0.01 & p_val < 0.05] <- "\u002A\u002A"
  stars[p_val < 0.01] <- "\u002A\u002A\u002A"
  return(stars)
}





