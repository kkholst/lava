`addattr` <-
function(x,...) UseMethod("addattr")

`addattr.lvm` <-
function(x, attr, var=NULL, val=TRUE, fun="nodeRenderInfo",debug=FALSE,...) {
  if (!is.null(var)) {
    Graph(x) <- addattr(Graph(x), attr=attr, var=var, val=val, fun=fun, debug=debug)
    return(x)
  } else {
    addattr(Graph(x), attr=attr, var=var, val=val, fun=fun)
  }
}

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
