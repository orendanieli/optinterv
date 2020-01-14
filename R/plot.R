#' Add signs to variable names
#'
#' @param names vector of variable names.
#' @param signs vector of signs (the same length as names).
add_sign <- function(names, signs){
  sgn_symbol <- ifelse(signs > 0, " (+)", ifelse(signs < 0, " (-)", " (\u00B1)"))
  var_names <- paste0(names, sgn_symbol)
  return(var_names)
}


#' Variable Position
#'
#' Find which variables to plot.
#'
#' @inheritParams plot.optint
#' @return vector of variables incidents

var_pos <- function(x, plot.vars = "sig", alpha, ...){
  n <- length(x$estimates)
  #estimates positions
  if(is.numeric(plot.vars)){
    inc <- 1:plot.vars
  } else {
    if(plot.vars == "sig") {
      inc <- which(x$details$p_value < alpha)
      if(length(inc) == 0){
        inc <- 1
        warning("There are no significant variables. display only the first one")
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
#' Produce variable importance plot from an optint object.
#'
#' @param x an optint object.
#' @param plot.vars which variables to plot? either a number (n) -
#'                  indicating to plot the first n variables,
#'                  "sig" (default) - plot only significant  variables, or a vector
#'                   with names of variables to plot.
#' @param plot.ci logical. if TRUE (default) plot confidence intervals. Otherwise
#'                plot only point estimates.
#' @param graph.col graph color/s.
#' @param alpha significance level for the confidence intervals. also
#'              used in order to determine which variables are significant.
#' @export

plot.optint <- function(x, plot.vars = "sig", plot.ci = T,
                        graph.col = 1, alpha = 0.05, ...){
  inc <- var_pos(x, plot.vars, alpha)
  estimates <- x$estimates[inc]
  #absolute value of point estimates
  estimates <- abs(estimates)
  #order by magnitude
  inc <- inc[order(estimates)]
  estimates <- estimates[order(estimates)]
  var_names <- colnames(x$details$new_sample)[inc]
  #decide font size
  fsize <- ifelse(length(inc) > 20, 0.5, 0.9)
  #add sign to name
  var_names <- add_sign(var_names, x$details$signs[inc])
  if(plot.ci){
    #find values for confidence intervals
    low_ci <- x$details$ci[1,inc]
    up_ci <- x$details$ci[2,inc]
    if(x$details$method == "correlations"){
      #flip ci for negative estimates
      neg_est <- x$details$signs[inc] < 0
      temp_ci <- -1 * low_ci[neg_est]
      low_ci[neg_est] <- -1 * up_ci[neg_est]
      up_ci[neg_est] <- temp_ci
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
    dotchart(estimates, var_names, xlim = c(0, max(estimates)), pch = 19,
             cex = fsize,
             col = .col, col.axis = 2)
  }
  #add 0
  abline(v = 0, lty = 2, col = 2)
  #add sample size
  title(xlab = paste0("N = ", nrow(x$details$new_sample)),
        adj = 0.99, line = -1, cex.lab = 1.5)
  #legend("bottomright", paste0("N = ", nrow(x$details$new_sample)))#, cex = 0.8 * fsize)
}

#' Plot the change in the distribution of X
#'
#' Illustrates how the intervention changes the distribution of X by plotting
#' barchart (for categorical variables) or denisty plot (for continouous variables),
#' before and after the intervention.
#'
#' @param n.val variable with more values than 'n.val' will be displayed by
#'              density plot, while variable with fewer values will be
#'              displayed by barchart.
#' @param line.type line type for \code{\link[lattice]{densityplot}}
#' @inheritParams plot.optint
#' @export

plot_change <- function(x, plot.vars = "sig",
                        graph.col = c("red", "blue"),
                        alpha = 0.05, line.type = c(1,2), n.val = 10, ...){
  if(x$details$method == "correlations")
    stop(paste("plot_change() isn't available for the correlations method.",
               "the distribution of x just shifts by 1/lambda * cov(x,y)"))
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
           lines=list(col=graph.col,  lty=line.type, lwd = 4),
           text=list(c("After","Before")))
  for(i in inc){
    var <- x$details$new_sample[,i]
    values <- unique(var)
    num_val <- length(values)
    #for categorical variables, plot barchart
    if(num_val <= n.val){
      values <- values[order(values)]
      #create data:
      gdata <- matrix(nrow = num_val, ncol = 2,
                      dimnames = list(as.character(values), c("Before", "After")))
      r_ind <- 1
      #calculate frequency for each value seperately
      for(v in values){
        gdata[r_ind,1] <- weighted.mean(var == v, wgt)
        gdata[r_ind,2] <- weighted.mean(var == v, wgt1)
        r_ind <- r_ind + 1
      }
      graph <- lattice::barchart(gdata, stack = F, horizontal = F, ylab = "",
                                 auto.key = list(space="top", columns= 2),
                                 par.settings = list(superpose.polygon = list(col = graph.col[c(2,1)])),
                                 xlab = var_names[count])

      print(graph)
    } else {
      #for continouous variables, plot densityplot
      #graph boarders:
      x_lim <- c(min(var), max(var))
      var <- rep(var, 2)
      graph <- lattice::densityplot(~ var, weights = wgt_all,
                                    groups = interv, xlim = x_lim,
                                    ylab = "", xlab = var_names[count],
                                    col = graph.col, plot.points = T,
                                    key = leg, lwd = 3, lty = line.type)
      print(graph)
    }
    count <- count + 1
  }
}

#' Plot optint object, by group
#'
#' Produce variables importance plot from an optint_by_group object.
#' This plot has several features:
#' \describe{
#'   \item{1.}{Estimates here are \eqn{E(X | I=1) - E(X | I=0)} and not cdf distances.}
#'   \item{2.}{Star is added to variable name if there is a
#'             significant difference between at least two groups.}
#'   \item{3.}{Estimates are standardized before they plotted (so different set of variables
#'             will have different standardization factor.)}}
#' @param x an optint_by_group object.
#' @param plot.vars which variables to plot? either a number (n) -
#'                  indicating to plot the first n variables,
#'                  "sig" (default) - plot only significant variables
#'                  (here significant means that variabe is signifcant for
#'                  at least one group, or that there is significant heterogeneity),
#'                  or a vector with names of variables to plot.
#' @inheritParams plot.optint
#' @export

plot.optint_by_group <- function(x,
                                 plot.vars = "sig",
                                 graph.col = NULL,
                                 alpha = 0.05, ...){
  est <- x$est
  sd <- x$sd
  n_group <- ncol(est)
  if(is.null(graph.col)){
    graph.col <- seq(2, n_group*2, 2)
  }
  z <- qnorm(1 - (alpha / 2))
  #t-state for group differences
  tstat <- (est[,-1] - est[,-n_group]) / sqrt(sd[,-1]^2 + sd[,-n_group]^2)
  tstat <- as.matrix(abs(tstat))
  #find which group is the max for each var (will be useful later)
  var_max <- apply(est, 1, function(x){ x[which.max(abs(x))] })
  #absolut estimates
  est <- abs(est)
  #confidence intervals:
  lower_ci <- est - sd*z
  upper_ci <- est + sd*z
  #min tstat by variable
  tstat_min <- apply(tstat, 1, min)
  var_names <- row.names(est)
  #which variables to plot?
  if(is.numeric(plot.vars)){
    inc <- 1:plot.vars
  } else {
    if(plot.vars == "sig"){
      #lower low_ci by variable
      lower_ci_min <- apply(lower_ci, 1, min)
      #take variables with significance difference between groups or with at
      #least one significance group
      inc <- which((lower_ci_min > 0 | tstat_min > z))
      if(length(inc) == 0){
        inc <- 1
        warning(paste("There are no variables with significance difference",
                      "between groups or with at least one significance group.",
                      "displays the first variable"))
      }
    } else {
      inc <- which(var_names %in% plot.vars)
    }
  }
  var_names <- var_names[inc]
  var_max <- var_max[inc]
  tstat_min <- tstat_min[inc]
  #sign is based on the higher point estimate
  sgn <- sign(var_max)
  #add sign to name
  var_names <- add_sign(var_names, sgn)
  #add star if at least one difference is significant
  var_names <- paste0(ifelse(tstat_min > z, "*", ""), var_names)
  #sort to put the highest values first
  inc <- inc[order(abs(var_max))]
  #go back to original estimates (not absolut):
  est <- x$est[inc,]
  sd <- sd[inc,]
  #if est & sd are vectors, transform them to matrix
  n_vars <- length(inc)
  if(n_vars == 1){
    est <- t(as.matrix(est))
    sd <- t(as.matrix(sd))
  }
  #normalize to SD units
  stand_factor <- sd(est)
  lower_ci <- (est - sd*z) / stand_factor
  upper_ci <- (est + sd*z) / stand_factor
  est <- est / stand_factor
  #borders of graph to display, to include all CI
  x_lim <- c(min(min(lower_ci), 0), max(upper_ci))
  #decide font size
  fsize <- ifelse(n_vars  > 20, 0.5, 0.9)
  #empty plot with correct borders
  dotchart(rep(80, n_vars), var_names,
           pch = 19, cex = fsize, col = 1,
           col.axis = 2,
           xlim = x_lim)
  #add points  and CIs
  for (j in 1:n_group){
    #soem ofset to get the groups not on top of each other
    ofst <- (j - 1) / (2 * n_group)
    points(est[,j], (1:n_vars) - ofst, col = graph.col[j],
                       pch = 19, cex = fsize)
    #add CIs
    for (i in 1:n_vars){
      lines(x = c(lower_ci[i,j], upper_ci[i,j]), y = c(i,i) - ofst,
            col = graph.col[j])
    }
  }
  #add legend
  #legend isnt stable. change it.
  group_names <- colnames(est)
  legend("bottomright", group_names, col = graph.col, pch = 19, cex = 1)
  #add vertical line at 0
  abline(v = 0, lty = 2)
}









