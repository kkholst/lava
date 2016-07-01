`getPath` <- function(x, ...) UseMethod("getPath")

`isPath` <- function(x, ...) UseMethod("isPath")

`getPath<-` <- function (x, ..., value) {
  UseMethod("getPath<-", x)
}

`isPath.plvmfit` <- function(x) {
  
  return( is.null(x$regularizationPath) )
  
}

`getPath.plvmfit` <- function(x, type = NULL, row = NULL) {
  
  validNames <- names(x$regularizationPath)
  
  if(is.null(type)){
    res <- x$regularizationPath
  }else if(type == "coef"){
    res <- x$regularizationPath[,names(coef(x)),drop = FALSE]
  }else if(type %in% validNames){
    res <- x$regularizationPath[,type,drop = FALSE]
  }else{
    stop("getPath.plvmfit: ",type," is an invalid name \n",
         "valid names: \"coef\" ",paste(validNames, collapse = "\" \""),"\n")
  }
  
  
  if(is.null(row)){
    return(res)
  }else{
    return(res[row,,drop = FALSE])
  }
}



`getPath<-.plvmfit` <- function(x, row = NULL, names = NULL, value) {
  
  if(is.null(names)){
    names.value <- names(value)
    
    if(length(names.value)>0){
      names <- names(value)  
    }else{
      stop("getPath<-: argument names must be specified \n")
    }
  }
  
  if(is.null(row)){
    x$regularizationPath[, names] <- value
  }else{
    x$regularizationPath[row, names] <- value
  }
  
  return(x)
}