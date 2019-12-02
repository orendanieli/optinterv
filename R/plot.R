#' Variable Position
#'
#' Find which variables to plot.
#'
#' @inheritParams plot.optint
#' @return vector of variables incidents

var_pos <- function(object, plot.vars = "all", alpha, ...){
  x <- object
  n <- length(x$estimates)
  #estimates positions
  if(length(plot.vars) == 1){
    if(plot.vars == "all"){
      inc <- 1:n
    } else {
      inc <- which(x$details$p_value < alpha)
      if(length(inc) == 0){
        inc <- 1
        warning("There are no significant variables. displays the first variable")
      }
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
#' @param alpha significance level. only used if plot.vars = "sig.
#' @export

#order variables & deal with negative ci
plot.optint <- function(object, plot.vars = "sig", plot.ci = T,
                        graph.col = 1, alpha = 0.05, ...){
  x <- object
  inc <- var_pos(x, plot.vars, alpha)
  estimates <- x$estimates[inc]
  var_names <- colnames(x$details$new_sample)[inc]
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

plot_change <- function(object, plot.vars = "sig",
                        graph.col = 1, alpha = 0.05, ...){
  x <- object
  if(x$details$method == "correlations")
    stop("plot_change() isn't available for the correlations method")
  inc <- var_pos(x, plot.vars, alpha)
  wgt <- x$details$new_sample[,"wgt"]
  wgt1 <- x$details$new_sample[,"wgt1"]
  #make weights sum to 1
  wgt <- wgt / sum(wgt)
  wgt1 <- wgt1 / sum(wgt1)
  wgt_all <- c(wgt, wgt1)
  var_names <- colnames(x$details$new_sample)[inc]
  #prepare graph data
  n <- nrow(x$details$new_sample)
  gdata <- data.frame(w = wgt_all, inter = c(rep("Before", n), rep("After",n)))
  n_graphs <- length(inc)
  count <- 1
  print(var_names)
  for(i in inc){
    gdata <- cbind(var = rep(x$details$new_sample[,i], 2), gdata)
    graph <- densityplot(~ var, weights = w, groups = inter, data = gdata,
                         auto_key = T, ylab = "", xlab = var_names[count],
                         col = c("red", "blue"), plot.points = F)
    final <- count == n_graphs
    print(graph, split = c(count, 1, n_graphs, 1), more = !final)
    count <- count + 1
  }
}

#legend + histogram

