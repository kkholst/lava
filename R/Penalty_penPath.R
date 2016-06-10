`penPath` <- function(x, ...) UseMethod("penPath")

`penPath.plvmfit` <- function(x, ...) {
  
  return(x$opt$message)
  
}

