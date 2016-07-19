`plotPath` <- function(x, ...) UseMethod("plotPath")

`plotPath.plvmfit` <- function(x, lambda = "lambda1.abs", 
                               add.line = TRUE, line.size = 2,
                               add.point = FALSE, point.size = 2) {
  
  if(isPath(x)){
    stop("plotPath.plvmfit: no penalization path in the plvmfit object \n",
         "set argument \'regularizationPath\' to 1 or 2 when calling estimate \n")
  }
  
  dt.Path <- data.table::melt(getPath(x, getCoef = "penalized", getLambda = lambda), 
                              measure=coef0(x, tol = NULL, penalized = TRUE, value = FALSE), 
                              value.name = "value", variable.name = "coefficient")
  
  ggPath <- ggplot(dt.Path, aes_string(y = "value", x = lambda, group = "coefficient", col = "coefficient"))
  if(add.line){ggPath <- ggPath + geom_line(size = line.size)}
  if(add.line){ggPath <- ggPath + geom_point(size = point.size)}
  return(ggPath)
}
