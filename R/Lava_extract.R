`coefVar` <-
  function(x,...) UseMethod("coefVar")

`coefCov` <-
  function(x,...) UseMethod("coefCov")

`loadings` <- function(object, ...) UseMethod("loadings")

#' @title Extract the name or the position of the variance coefficients
coefVar.lvm <- function(x, index = FALSE){
  names.var <- paste(x$index$vars, x$index$vars, sep = ",")
  return(grep(paste(names.var, collapse = "|"), coef(x), value = index))
}

#' @title Extract the name or the position of the covariance coefficients
coefCov.lvm <- function(x, index = FALSE){
  names.cov <- setdiff(coef(x)[x$index$parBelongsTo$cov],
                       paste(x$index$vars, x$index$vars, sep = ",")
  )
  return(grep(paste(names.cov, collapse = "|"), coef(x), value = index))
}

#' @title Extract the summary table for the loadings
loadings.lvmfit <- function(x, col = NULL){
  expr <- paste("^",endogenous(x),"\\~", collapse = "|", sep = "")
  index.loadings <- grep(expr, names(coef(x)), value = TRUE)
  loadings <- summary(x)$coef[index.loadings,]
  if(is.null(col)){
    return(loadings)
  }else{
    return(loadings[,col])
  }
}

loadings.lvm.missing <- loadings.lvmfit