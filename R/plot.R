#' Variable Position
#'
#' Find which variables to plot.
#'
#' @inheritParams plot.optint
#' @return vector of variables incidents

var_pos <- function(object, plot.vars = "all"){
  x <- object
  n <- length(x$estimates)
  #estimates positions
  if(length(plot.vars) == 1){
    if(plot.vars == "all"){
      inc <- 1:n
    }
  } else {
    var_names <- colnames(x$details$new_sample[,1:n])
    inc <- which(var_names %in% plot.vars)
  }
  return(inc)
}

#' Plot optint object
#'
#' Produce variables important plot from an optint object.
#'
#' @param object an optint object.
#' @param plot.vars which variables to plot? "all" (default) plot all variables and
#'                  "sig" plot only siginficance variables. Another option is vector
#'                   with variables names to plot.
#' @param plot.ci logical. if TRUE (default) plot condifendece intervals. Otherwise
#'                plot only point estimates.
#' @export

#order variables & deal with negative ci
plot.optint <- function(object, plot.vars = "all", plot.ci = T, graph.col = 1, ...){
  x <- object
  inc <- var_pos(x, plot.vars)
  estimates <- x$estimates[inc]
  var_names <- colnames(x$details$new_sample[,inc])
  #absolute value of point estimates
  estimates <- abs(estimates)
  #order by magnitude
  inc <- inc[order(estimates)]
  estimates <- estimates[inc]
  var_names <- var_names[inc]
  #decide font size
  fsize <- ifelse(length(inc) > 20, 0.5, 0.9)
  #add sign to name
  sgn <- x$details$signs[inc]
  sgn_symbol <- ifelse(sgn > 0, " (+)", ifelse(sgn < 0, " (-)", " (±)"))
  var_names <- paste0(var_names, sgn_symbol)
  if(plot.ci){
    #find values for confidence intervals
    low_ci <- x$details$ci[1,inc]
    up_ci <- x$details$ci[2,inc]
    if(x$details$method == "correlations"){
      #flip ci for negative estimates
      neg_est <- x$details$signs[inc] < 0
      low_ci[neg_est] <- -1 * low_ci[neg_est]
      up_ci[neg_est] <- -1 * up_ci[neg_est]
    }
    #plot
    graph_bor <- c(min(min(low_ci), 0), max(up_ci))
    dotchart(estimates, var_names, xlim = graph_bor, pch = 19, cex = fsize,
             col = graph.col, col.axis = 2)
    #add CIs manually as lines
    j <- 1
    n <- length(inc)
    for (i in 1:n){
      lines(x = c(low_ci[i], up_ci[i]), y = c(j, j))
      j <- j + 1
    }
  } else {
    dotchart(estimates, var_names, xlim = c(0, max(estimates)), pch = 19, cex = fsize,
             col = .col, col.axis = 2)
  }
  #add 0
  abline(v = 0, lty = 2, col = 2)
  #add sample size
  legend("bottomright", paste0("N = ", nrow(x$details$new_sample)), cex = 0.8 * fsize)
}

plot_change <- function(object, plot.vars = "all", plot.ci = T, graph.col = 1, ...){
  x <- object
  inc <- var_pos(x, plot.vars)
  wgt <- x$details$new_sample[,"wgt"]
  wgt1 <- x$details$new_sample[,"wgt1"]
  #var_names <- colnames(x$details$new_sample[,inc])
  par(mfrow = c(length(inc), 2))
  for(i in inc){
    plot(density(x$details$new_sample[,i], weights = wgt))
    plot(density(x$details$new_sample[,i], weights = wgt1))
  }
}


