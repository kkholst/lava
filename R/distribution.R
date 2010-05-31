
"distribution<-" <- function(x,...,value) UseMethod("distribution<-")
"distribution" <- function(x,...,value) UseMethod("distribution")

"distribution<-.lvm" <- function(x,variable,...,value) {
  if (class(variable)[1]=="formula") {
    myvars <- all.vars(variable)
    if (length(value)!=length(myvars) & length(value)!=1) stop("Wrong number of values")
    for (i in 1:length(myvars))
      if (length(value)==1) {
        distribution(x,myvars[i],...) <- value
      } else {
        distribution(x,myvars[i],...) <- value[[i]]
      }
    return(x)
  }  
  var. <- paste("\"", variable, "\"", sep="")
  mytext <- paste("c(", paste(paste(var.,"=",expression(value),sep=""),collapse=","),")")
  nodeRenderInfo(Graph(x))$"distribution" <- eval(parse(text=mytext))
  x
}
"distribution.lvm" <- function(x,var,...) {
  nodeRenderInfo(Graph(x))$distribution[var]
}


