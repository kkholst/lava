`penPath` <- function(x, ...) UseMethod("penPath")

`penPath<-` <- function (x, ..., value) {
  UseMethod("penPath<-", x)
}

`penPath.plvmfit` <- function(x, type = NULL, row = NULL) {
  
  validNames <- names(x$opt$message)
  
  if(is.null(type)){
    res <- x$opt$message
  }else if(type == "coef"){
    res <- x$opt$message[,names(coef(x)),drop = FALSE]
  }else if(type %in% validNames){
    res <- x$opt$message[,type,drop = FALSE]
  }else{
    stop("penPath.plvmfit: ",type," is an invalid name \n",
         "valid names: \"coef\" ",paste(validNames, collapse = "\" \""),"\n")
  }
  
  
  if(is.null(row)){
    return(res)
  }else{
    return(res[row,,drop = FALSE])
  }
}



`penPath<-.plvmfit` <- function(x, row = NULL, names = NULL, value) {
  
  if(is.null(names)){
    names.value <- names(value)
    
    if(length(names.value)>0){
      names <- names(value)  
    }else{
      stop("penPath<-: argument names must be specified \n")
    }
  }
  
  if(is.null(row)){
    x$opt$message[, names] <- value
  }else{
    x$opt$message[row, names] <- value
  }
  
  return(x)
}