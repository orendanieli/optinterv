validate_data <- function(Y, X, control, wgt){
  n <- length(Y)
  if(!inherits(X, c("matrix","data.frame")) | n != nrow(X)){
    stop(paste("X must be either matrix or data.frame,",
               "with the same number of examples as Y."))
  }
  if(!inherits(control, c("matrix","data.frame")) | n != nrow(control)){
    stop(paste("control must be either matrix or data.frame,",
               "with the same number of examples as Y."))
  }
  if(!inherits(wgt, "numeric") | n != length(wgt)){
    stop(paste("wgt must be numeric,",
               "with the same number of examples as Y."))
  }
  dat <- cbind(Y, X, control, wgt)
  if(any(is.na(dat))){
    stop("Data contains NA's.")
  }
  if(!is.numeric(as.matrix(dat))){
    stop(paste("Only numerical variables are allowed.",
               "For categorical variables, please use 1/0 encoding",
               "(see ?model.matrix for further details)"))
  }
}

prepare_Y <- function(Y){
  if(min(Y) < 0 | length(unique(Y)) < 3){
    Y <- exp(Y)
  }
  return(Y)
}




