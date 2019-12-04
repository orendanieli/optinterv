#' Variable Position
#'
#' Find which variables to plot.
#'
#' @inheritParams plot.optint
#' @return vector of variables incidents

var_pos <- function(object, plot.vars = "sig", alpha, ...){
  x <- object
  n <- length(x$estimates)
  #estimates positions
  if(plot.vars == "all"){
    inc <- 1:n
  } else {
    if(plot.vars == "sig") {
      inc <- which(x$details$p_value < alpha)
      if(length(inc) == 0){
        inc <- 1
        warning("There are no significant variables. displays the first variable")
      }
    } else {
      var_names <- colnames(x$details$new_sample[,1:n])
      inc <- which(var_names %in% plot.vars)
    }
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
#' @param graph.col graph color/s.
#' @param alpha significance level. only used if plot.vars = "sig.
#' @export

plot.optint <- function(object, plot.vars = "sig", plot.ci = T,
                        graph.col = 1, alpha = 0.05, ...){
  x <- object
  inc <- var_pos(x, plot.vars, alpha)
  estimates <- x$estimates[inc]
  #absolute value of point estimates
  estimates <- abs(estimates)
  #order by magnitude
  inc <- inc[order(estimates)]
  estimates <- x$estimates[inc]
  var_names <- colnames(x$details$new_sample)[inc]
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
  legend("bottomright", paste0("N = ", nrow(x$details$new_sample)))#, cex = 0.8 * fsize)
}

#' Plot change in X distribution
#'
#' Plot denisty or histogram of X, before and after the intervention.
#'
#' @param print.sep logical. If TRUE (default) plot each graph seperatly.
#'                  This option is highly recommended for more then 6 variables.
#' @param line.type
#' @inheritParams plot.optint
#' @export

plot_change <- function(object, plot.vars = "sig",
                        graph.col = c("lightcoral", "cadetblue3"),
                        alpha = 0.05, print.sep = T, line.type = c(1,2), ...){
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
  #prepare graph data (for densityplot)
  n <- nrow(x$details$new_sample)
  interv <- c(rep("Before", n), rep("After",n))
  n_graphs <- length(inc)
  count <- 1
  #edit legend (for densityplot)
  leg <- list(space="top",
           lines=list(col=graph.col,  lty=line.type),
           text=list(c("Before","After")))
  for(i in inc){
    var <- x$details$new_sample[,i]
    #for binary variables, plot barchart
    if(all(unique(var) %in% c(0,1))){
      freq_bef <- weighted.mean(var, wgt)
      freq_aft <- weighted.mean(var, wgt1)
      gdata <- matrix(c(freq_bef, freq_aft, 1 - freq_bef, 1 - freq_aft),
                      ncol = 2, dimnames = list(c("Before", "After")))
      graph <- lattice::barchart(gdata, stack = T, horizontal = F,
                                 col = graph.col, xlab = var_names[count])
    } else {
      #for continouous variables, plot densityplot
      var <- rep(var, 2)
      graph <- lattice::densityplot(~ var, weights = wgt_all, groups = interv,
                                    ylab = "", xlab = var_names[count],
                                    col = graph.col, plot.points = T,
                                    key = leg, lwd = 3, lty = line.type)
    }
    if(print.sep | count == n_graphs){
      print(graph)
    } else {
      if(count > 4){
        stop("can't print more than 6 graphs. please set 'print.sep = T'")
      }
      final <- count == n_graphs | count == 4
      tot_col <- min(n_graphs, 4)
      print(graph, split = c(count, 1, tot_col, 1), more = !final)
    }
    count <- count + 1
  }
}

plot.optint_by_group <- function(object,
                                 plot.vars = "sig",
                                 graph.col = 1,
                                 alpha = 0.05, ...){
  est <- object$est
  sd <- object$sd
  #max absolute value for each group
  ext_est <- apply(est, 1, function(x) max(abs(x)))
  #which group is the max for each var
  raw.max = apply(vs,1,function(r) r[which.max(abs(r))])
  #t-state for group differences
  if (se){
    #difference divided by SE of the difference
    tstat = (vs[,1]-vs[,2])/(vse[,1]^2+vse[,2]^2)^0.5
  }
}



