##' Generic method for adding attributes to object or components of a class
##'
##' @title Add attribute to class
##' @param x Object 
##' @param ... Additional arguments
##' @author Klaus K. Holst
##' @export
`addattr` <-
function(x,...) UseMethod("addattr")

##' @S3method addattr lvm
`addattr.lvm` <-
function(x, attr, var=NULL, val=TRUE, fun="nodeRenderInfo",debug=FALSE,...) {
  if (!is.null(var)) {
    Graph(x) <- addattr(Graph(x), attr=attr, var=var, val=val, fun=fun, debug=debug)
    return(x)
  } else {
    addattr(Graph(x), attr=attr, var=var, val=val, fun=fun)
  }
}

##' @S3method addattr graphNEL
`addattr.graphNEL` <-
function(x, attr, var=NULL, val=TRUE,fun="nodeRenderInfo",debug=FALSE,...) {
  if (is.null(var)) {
    f <- eval(call(fun,x))
    if (is.null(val))
      attrvar <- names(f[[attr]])
    else
      attrvar <- names(f[[attr]])[which(val==f[[attr]])]
    return(attrvar)
  }
  if (is.character(val)) 
    myexpr <- paste("list(",attr,"=c(", paste("\"",var,"\"=\"",val,"\"" , sep="", collapse=", "), "))", sep="")
  else
    myexpr <- paste("list(",attr,"=c(", paste("\"",var,"\"=",val, sep="", collapse=", "), "))", sep="")
  Debug(list("str=",myexpr),debug)
  eval(parse(text=paste(fun,"(x) <- ",myexpr,sep="")))
  return(x)
}
