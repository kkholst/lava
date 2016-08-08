`getPath` <- function(x, ...) UseMethod("getPath")

`getPath.plvmfit` <- function(x, names = NULL, getCoef, getLambda, rm.duplicated = FALSE, ascending = TRUE, row = NULL) {
  
  if(isPath(x)){
    stop("plotPath.plvmfit: no penalization path in the plvmfit object \n",
         "set argument \'regularizationPath\' to 1 or 2 when calling estimate \n")
  }
  
  regPath <- x$regularizationPath
  if(rm.duplicated){
    indexChange <- regPath$indexChange
    indexChange[is.na(indexChange)] <- -1
    test.change <- which(diff(na.omit(indexChange))!=0)+1
    index <- sort(union(c(1,NROW(regPath)),intersect(which(!is.na(regPath$indexChange)), test.change)))
    regPath <- regPath[index, ,drop = FALSE]
  }
  if(ascending == TRUE){
    regPath <- regPath[order(regPath$lambda1.abs, decreasing = FALSE),,drop = FALSE]
  }else{
    regPath <- regPath[order(regPath$lambda1.abs, decreasing = TRUE),,drop = FALSE]
  }
  validNames <- names(regPath)
  
  if(!is.null(names)){
    
    if(all(names %in% validNames)){
      res <- regPath[,names,drop = FALSE]
    }else{
      stop("getPath.plvmfit: ",names," contains invalid names \n",
           "invalid names: ",paste(names[names %in% validNames == FALSE], collapse = "\" \""),"\n",
           "valid names: ",paste(validNames, collapse = "\" \""),"\n")
    }
    
  } else {
    
    if(missing(getLambda)){
      names.lambda <- c("lambda1.abs", "lambda1", "lambda2.abs", "lambda2")
    }else if(is.null(getLambda)){
      names.lambda <- NULL
    }else{
      validValues <- c("abs", "nabs", "lambda1", "lambda2", "lambda1.abs", "lambda2.abs")
      if(any(getLambda %in% validValues == FALSE)){
        stop("getPath.plvmfit: invalid value for \'getLambda\' \n",
             "getLambda: ",getLambda,"\n",
             "valid values: ",paste(validValues, collapse = "\" \""),"\n")
      }
      
      names.lambda <- sapply(getLambda, function(x){
        switch(x,
               "abs" = c("lambda1.abs", "lambda2.abs"),
               "nabs" = c("lambda1", "lambda2"),
               "lambda1" = c("lambda1"),
               "lambda2" = c("lambda2"),
               "lambda1.abs" = c("lambda1.abs"),
               "lambda2.abs" = c("lambda2.abs")
        )
      })
    }
    
    if(missing(getCoef)){
      names.coef <- setdiff(validNames, c("lambda1.abs", "lambda1", "lambda2.abs", "lambda2", "indexChange"))
    }else if(is.null(getCoef)){
      names.coef <- NULL
    }else{
      validValues <- c("penalized", "npenalized", "coef0", "coefn0")
      if(getCoef %in% validValues == FALSE){
        stop("getPath.plvmfit: invalid value for \'getCoef\' \n",
             "getCoef: ",getCoef,"\n",
             "valid values: ",paste(validValues, collapse = "\" \""),"\n")
      }
      
      names.penalized <- intersect(validNames, x$penalty$names.penaltyCoef)
      names.coef <- switch(getCoef,
                           "penalized" = names.penalized,
                           "npenalized" = setdiff(names(coef(x)),names.penalized),
                           "coef0" = {regPath$coef0 <- rowSums(abs(regPath[,names.penalized, drop = FALSE] == 0)) ; "coef0"},
                           "coefn0" = {regPath$coefn0 <- rowSums(abs(regPath[,names.penalized, drop = FALSE] != 0)) ; "coefn0"}
      )
    }
    
    res <- regPath[,c(names.lambda, names.coef),drop = FALSE]
  }
    
  if(is.null(row)){
    return(res)
  }else{
    return(res[row,,drop = FALSE])
  }
}

`setPath<-` <- function (x, ..., value) {
  UseMethod("setPath<-", x)
}

`setPath<-.plvmfit` <- function(x, row = NULL, names = NULL, value) {
  
  if(is.null(names)){
    names.value <- names(value)
    
    if(length(names.value)>0){
      names <- names(value)  
    }else{
      stop("setPath<-: argument names must be specified \n")
    }
  }
  
  if(is.null(row)){
    x$regularizationPath[, names] <- value
  }else{
    x$regularizationPath[row, names] <- value
  }
  
  return(x)
}

`getLambda` <- function(x, ...) UseMethod("getLambda")

`getLambda.plvmfit` <- function(x, lambda1 = TRUE, lambda2 = FALSE, abs = TRUE) {
  
  name <- NULL
  if(lambda1){name <- c(name, "lambda1")}
  if(lambda2){name <- c(name, "lambda2")}
  if(abs){name <- paste(name, "abs", sep = ".")}
  return(getPath(x, names = name))
 
}

`isPath` <- function(x, ...) UseMethod("isPath")

`isPath.plvmfit` <- function(x) {
  
  return( is.null(x$regularizationPath) )
  
}







# 
# ### smooth path
# regPath <- x$regularizationPath
# n.coef0 <- rowSums(abs(regPath[,x$penalty$names.penaltyCoef, drop = FALSE] == 0))
# nmax.coef0 <- max(nCoef)
# 
# nCoefChange <- sapply(seq(0, n.coef0), function(x){which(nCoef==x)[1]})
# diffCoef <- diff(nCoef)
# lambda <- regPath$lambda1.abs
# diffLambda <- diff(lambda)
# 
# if(any(diffLambda == 0)){
#   
#   indexdiffLambda0 <- which(diffLambda == 0)
#   iterIndex <- 1
#   
#   while(iterIndex <= length(indexdiffLambda0)){
#     
#     test.smallLambda <- diffLambda[indexdiffLambda0[iterIndex]+1]
#     test.reverse <- diffCoef[indexdiffLambda0[iterIndex]+1] == - diffCoef[indexdiffLambda0[iterIndex]]
#     if(test.smallLambda && test.reverse){
#       lambda[indexdiffLambda0[iterIndex]+1]
#       regPath[-(indexdiffLambda0[iterIndex]+2),]
#     }
#     
#     
#     
#   }
#   
#   
#   
#   for(iterIndex in nCoefChange){
#     
#   }
#   
#   cond2 <- diffCoef
# }
