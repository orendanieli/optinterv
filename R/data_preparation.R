validate_data <- function(Y, X, control = NULL, wgt, ...){
  n <- length(Y)
  if(!inherits(X, c("matrix","data.frame")) || n != nrow(X)){
    stop(paste("X must be either matrix or data.frame,",
               "with the same number of examples as Y."))
  }
  if(!is.null(control)){
    if(!inherits(control, c("matrix","data.frame")) || n != nrow(control)){
      stop(paste("control must be either matrix or data.frame,",
                 "with the same number of examples as Y."))
    }
  }
  if(!inherits(wgt, "numeric") | n != length(wgt)){
    stop(paste("wgt must be numeric,",
               "with the same number of examples as Y."))
  }
  if(is.null(control)){
    dat <- cbind(Y, X, wgt)
  } else {
    dat <- cbind(Y, X, control, wgt)
  }
  if(any(is.na(dat))){
    stop("Data contains NA's.")
  }
  if(!is.numeric(as.matrix(dat))){
    stop(paste("Only numerical variables are allowed.",
               "For categorical variables, please use 1/0 encoding",
               "(see ?model.matrix for further details)"))
  }
  if(min(Y) < 0){
    message("Assume Y is in log units (since it includes negative values)")
  }
}

prepare_Y <- function(Y){
  if(min(Y) < 0 | length(unique(Y)) < 3){
    Y <- exp(Y)
  }
  return(Y)
}

validate_group <- function(Y, group){
  if(!is.vector(group)){
    stop("'group' must be a vector")
  }
  n <- length(Y)
  p <- length(group)
  if(!n==p){
    stop("Different length for Y and group")
  }
  num_group <- length(unique(group))
  if((n / num_group) < 30){
    warning(paste("number of groups =", num_group,
                  " ,but only ", n, "observations"))
  }
}




