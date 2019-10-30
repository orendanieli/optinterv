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
#' @param n_quant number of quantiles.
#' @param wgt original weights.
#' @param wgt1 weights under I = 1.
#' @param p proportion of quantiles that need to be increase (decrease) in order to
#'          return a positive (negative) sign.
#'
#' @return dataframe with percentile distance & sign: +1 (-1) for increasing (decreasing)
#'         p of quantiles. Else returns 0.

per_distance <- function(x, n_quant = length(x),
                         wgt = rep(1,length(x)),
                         wgt1,
                         p = 2/3){
  #fix weights to sum to data length
  wgt  <- (wgt / sum(wgt)) * length(wgt)
  wgt1 <- (wgt1 / sum(wgt1)) * length(wgt1)
  #get quantiles
  q0 <- wtd_bin(x, n_quant, wgt)
  q1 <- wtd_bin(x, n_quant, wgt1)
  distance <- mean(abs(q1 - q0))
  sgn <- sign(q1 - q0)
  tot <- sum(sgn) / n_quant
  cutoff <- p - (1 - p)
  res <- data.frame(distance = distance,
                    sign = ifelse(tot > cutoff, 1, ifelse(tot < -cutoff, -1 ,0)))
  return(res)
}






